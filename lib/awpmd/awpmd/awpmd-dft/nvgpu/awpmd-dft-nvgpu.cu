#include "../DataTypes.hpp"
#include "../utils/pls.hpp"
#include "../xcfunctionals/LDA.hpp"
#include "../xcfunctionals/LSDA.hpp"
#include "../xcfunctionals/ModLDA.hpp"
#include "CudaKernels.hpp"
#include "awpmd-dft-nvgpu.hpp"
#include "cuApproxTableObject.hpp"
#include "gpu_api.h"
#include "gpu_functors.hpp"
#include "gpu_integration.hpp"
#include "numeric_gpu.cuh"
#include <chrono>
#include <iostream>

bool XCEnergy_nvgpu::check_texture_object_support(){
  cudaDeviceProp cuda_dev_prop{};
  int current_device;
  cudaGetDevice(&current_device);
  cudaGetDeviceProperties(&cuda_dev_prop, current_device);

  return m_config.use_xc_tables && cuda_dev_prop.major >= 3;
}

class XCApproximation {
public:
  static std::unique_ptr<gpuApproxTablesObject> &gpu_approximation_table() {
    static std::unique_ptr<gpuApproxTablesObject> gpu_approximation_table_ =
        nullptr;
    return gpu_approximation_table_;
  }
};

XCEnergy::XCResult
XCEnergy_nvgpu::energy(std::vector<class WavePacket> const &wp_spin_up,
                       std::vector<class WavePacket> const &wp_spin_down,
                       bool calc_force) {
  if (m_config.approximation->Type == ApproxType::T_VOID) {
    return XCResult{0.0f, 0.0f, {}};
  }

  auto total_packet_count = wp_spin_down.size() + wp_spin_up.size();
  m_packet_count = total_packet_count;

  host_allocation(m_packet_count);
  std::vector<WavePacketInfo> wp;
  std::transform(wp_spin_up.begin(), wp_spin_up.end(), std::back_inserter(wp),
                 [calc_force](WavePacket const &p) {
                   return WavePacketInfo{p.get_r().get_ptr(), p.get_width(), ElectronSpin::E_UP, calc_force};
                 });
  std::transform(wp_spin_down.begin(), wp_spin_down.end(),
                 std::back_inserter(wp), [calc_force](WavePacket const &p) {
        return WavePacketInfo{p.get_r().get_ptr(), p.get_width(), ElectronSpin::E_DOWN, calc_force};
      });
  copy_packet_to_host(wp);

  XCEnergy::Energy xcenergy{};
  std::vector<std::vector<float>> derivatives;
  if (!m_config.use_adaptive_mesh) {
    std::tie(xcenergy.potential, xcenergy.kinetic) =
        regular_mesh_integration(derivatives, calc_force);
  } else {
    std::tie(xcenergy.potential, xcenergy.kinetic) =
        adaptive_mesh_integration(derivatives, calc_force, wp);
  }

  xcenergy.potential *= units().Hartree2Energy;
  xcenergy.kinetic *= units().Hartree2Energy;

  for (auto &va : derivatives)
    for (auto &a : va) {
      a *= units().Hartree2Energy * units().Distance2Bohr;
    }

  return XCResult(xcenergy.potential, xcenergy.kinetic, std::move(derivatives));
}

void XCEnergy_nvgpu::init_xc_approximation() {
  Logger::Info("Init XCApprox on GPU");
  IApproximation **ptr_ref;

  if (m_config.approximation_device != nullptr) {
    ::__cudft_remove_xc_approximation<<<1, 1>>>(m_config.approximation_device);
    m_config.approximation_device = nullptr;
  }

  SAFECALL(cudaMalloc(&ptr_ref, sizeof(IApproximation *)),
           "Error in cudaMalloc for xc_approx");

  switch (m_config.approximation->Type) {
    case ApproxType::T_LDA:
      ::__cudft_init_xc_approximation<LDA>
      <<<1, 1>>>(ptr_ref, 0.738558766f, -0.01554534543482745f, 20.4562557f);
      break;
    case ApproxType::T_LDA_2:
      ::__cudft_init_xc_approximation<ModLDA><<<1, 1>>>(ptr_ref);
      break;
    case ApproxType::T_LSDA:
      ::__cudft_init_xc_approximation<LSDA><<<1, 1>>>(ptr_ref);
      break;
    case ApproxType::T_VOID:
      ::__cudft_init_xc_approximation<VoidApproximation><<<1, 1>>>(ptr_ref);
    case ApproxType::T_DUMMY:
      ::__cudft_init_xc_approximation<DummyApproximation><<<1, 1>>>(ptr_ref);
      break;
  }
  SAFECALL(cudaThreadSynchronize(), "Error in thread sync for xc_approx");

  IApproximation **ptr = new IApproximation *[1];
  SAFECALL(cudaMemcpy(ptr, ptr_ref, sizeof(IApproximation *),
                      cudaMemcpyDeviceToHost),
           "Error in cudaMemcpy for xc_approx");

  Logger::Info("XCApprox GPU Ptr:", *ptr);
  m_config.approximation_device = *ptr;

//#ifndef CUDA_NO_TEXTURE_OBJ
  if (check_texture_object_support()) {
    if (m_config.approximation->Type == ApproxType::T_LSDA) {
      Logger::Info("Texture object is available. Switch to table xc-functional.");

      auto approx_table_props = ApproximationTableProps();
      approx_table_props.density_range[0] = static_cast<float>(m_config.packet_number) * 0.0001f;
      approx_table_props.density_range[1] = static_cast<float>(m_config.packet_number);

      XCApproximation::gpu_approximation_table() =
          std::unique_ptr<gpuApproxTablesObject>(new gpuApproxTablesObject(
              m_config.approximation, m_config.approximation_device,
              approx_table_props));

      ::__cudft_init_xc_approximation<TableApproximation<gpuApproxTablesObject>>
      <<<1, 1>>>(ptr_ref,
                 *(XCApproximation::gpu_approximation_table().get()));
      SAFECALL(cudaThreadSynchronize(), "Error in thread sync for xc_approx");
      SAFECALL(cudaMemcpy(ptr, ptr_ref, sizeof(IApproximation *), cudaMemcpyDeviceToHost),
               "Error in cudaMemcpy for xc_approx");
      Logger::Info("New XCApprox GPU Ptr:", *ptr);
      m_config.approximation_device = *ptr;
    }
  }
//#endif

  Logger::Info("Check...");
  ::__cudft_xc_check<<<1, 8>>>(m_config.approximation_device);
  SAFECALL(cudaThreadSynchronize(), "Error in thread sync for check xc_approx");

  Logger::Info("Done.");
}

XCEnergy_nvgpu::XCEnergy_nvgpu(unsigned int packet_number,
                               DFTConfig const &meshConfig, long gpu_id) {
  Logger::ModuleName("AWPMD-DFT GPU");

  Logger::Warning("AWPMD-DFT VERSION FROM", __DATE__);

  cudaSetDevice(gpu_id);

  m_packet_count = packet_number;
  Logger::Info("Init dft module on GPU =", gpu_id,
               " Packet number =", m_packet_count);

  setConfig(meshConfig);

  Logger::Info("Allocate",
               m_packet_count * sizeof(GaussPacket<float>) / (1024.0 * 1024.0),
               "Mb in host memory for packet.");
  host_allocation(m_packet_count);

  Logger::Info("Allocate",
               m_packet_count * sizeof(GaussPacket<float>) / (1024.0 * 1024.0),
               "Mb in device memory for packet.");
  device_allocation(m_packet_count);
}

XCEnergy_nvgpu::~XCEnergy_nvgpu() {
  if (XCApproximation::gpu_approximation_table() != nullptr)
    XCApproximation::gpu_approximation_table()->free();
}

void XCEnergy_nvgpu::device_allocation(size_t new_size) {

  if (!m_config.use_adaptive_mesh) {
    if (raw_energy.count() == 0) {
      raw_energy = Mesh_nvgpu(m_config.mesh_size.size.as_struct.x,
                              m_config.mesh_size.size.as_struct.y,
                              m_config.mesh_size.size.as_struct.z);

      raw_kinetic = Mesh_nvgpu(m_config.mesh_size.size.as_struct.x,
                               m_config.mesh_size.size.as_struct.y,
                               m_config.mesh_size.size.as_struct.z);

      rho = ElectronDensity<Mesh_nvgpu>(
          (m_config.mesh_size.size.as_struct.x + 1) * 2,
          m_config.mesh_size.size.as_struct.y * 2 + 2,
          m_config.mesh_size.size.as_struct.z * 2 + 2);
    }

    if (m_config.calc_force && raw_derivatives.m_depth < new_size) {
      raw_derivatives = Mesh_nvgpu(PacketType::derivatives_count(),
                                   m_config.mesh_size.size.as_struct.x *
                                   m_config.mesh_size.size.as_struct.y *
                                   m_config.mesh_size.size.as_struct.z,
                                   new_size);
    }
  }

  if (dev_mem_size < new_size) {
    SAFECALL(cudaFree(device_mem), "Error in memory reallocation.");
    SAFECALL(cudaMalloc(&device_mem, new_size * sizeof(GaussPacket<float>)),
             "Error in allocation memory for packet");
    dev_mem_size = new_size;
  }
}

std::tuple<double, double> XCEnergy_nvgpu::adaptive_mesh_integration(std::vector<std::vector<float>> &derivatives,
    bool calc_force, const vector<WavePacketInfo> &wavepackets) {

  Integrator_nvgpu<MeshCell> energy_integration_engine_;

  integrator_.refine_mesh(
      MeshCell::RangeType(m_config.mesh_start.size.as_struct), MeshCell::RangeType(m_config.mesh_fin.size.as_struct),
      [this](MeshCell::RangeType a, MeshCell::RangeType b) {
        return this->refine_mesh_condition(host_mem.begin(), host_mem.end(), a, b);
      },
      MeshhWorldTopo{m_config.nodes, m_config.node_rank});

  energy_functor_sic<MeshCell, PacketType> energy_eval{};
  energy_eval.approx = m_config.approximation_device;

  auto xck_energy = integrator_.integrate<NumericType<float, float>>(
      energy_integration_engine_, host_mem, energy_eval);

  if (calc_force) {
    Integrator_nvgpu<ForceMeshCell> force_integration_engine_;

    force_integrator_.clear();
    auto nbins = m_config.force_cell_bins;

    for (auto packet_index = 0; packet_index < host_mem.size();
         ++packet_index) {
      if (wavepackets[packet_index].calc_force) {
        auto const &packet = host_mem[packet_index];

        float volume_dh =
            packet.width * static_cast<PacketType::fields_type>(m_config.max_distance);

        ForceMeshCell::RangeType pbegin{packet.r[0] - volume_dh,
                                        packet.r[1] - volume_dh,
                                        packet.r[2] - volume_dh};
        ForceMeshCell::RangeType pend{packet.r[0] + volume_dh,
                                      packet.r[1] + volume_dh,
                                      packet.r[2] + volume_dh};

        force_integrator_.refine_linked_mesh(pbegin, pend, nbins, packet_index);
      }
    }

    derivatives = std::move(
        std::vector<std::vector<float>>(m_packet_count, std::vector<float>(PacketType::derivatives_count(), 0)));

    auto derivative_eval = derivatives_functor_sic<ForceMeshCell, PacketType>::from_energy_functor(energy_eval);

    auto batch_size = nbins * nbins * nbins;
    auto forces =
        force_integrator_.integrate<std::vector<NumericType<float, float, float, float>>>(
                force_integration_engine_, host_mem, batch_size, derivative_eval);

    auto force_it = forces.begin();
    for (auto packet_index = 0u; packet_index < m_packet_count;
         ++packet_index) {
      float dw = 0.0f, dx = 0.0f, dy = 0.0f, dz = 0.0f;
      if (wavepackets[packet_index].calc_force) {
        thrust::tie(dx, dy, dz, dw) = *force_it;
        ++force_it;
      }

      auto round_float = [](float v) {
        return std::abs(static_cast<double>(v)) <= std::numeric_limits<float>::epsilon() ? 0.0 : static_cast<double>(v);
      };

      derivatives[packet_index][0] += round_float(dx);
      derivatives[packet_index][1] += round_float(dy);
      derivatives[packet_index][2] += round_float(dz);
      derivatives[packet_index][3] += round_float(dw);
    };
  }

  return std::make_tuple<double, double>(thrust::get<0>(xck_energy),
                                         thrust::get<1>(xck_energy));
}

XCEnergy::XCResult
XCEnergy_nvgpu::energy(const vector<WavePacketInfo> &wavepackets,
                       bool calc_force) {
  if (m_config.approximation->Type == ApproxType::T_VOID) {
    return XCResult{0.0f, 0.0f, {}};
  }

  m_packet_count = wavepackets.size();

  host_allocation(m_packet_count);
  copy_packet_to_host(wavepackets);

  XCEnergy::Energy xcenergy{};
  std::vector<std::vector<float>> derivatives;
  if (!m_config.use_adaptive_mesh) {
    std::tie(xcenergy.potential, xcenergy.kinetic) =
        regular_mesh_integration(derivatives, calc_force);
  } else {
    std::tie(xcenergy.potential, xcenergy.kinetic) =
        adaptive_mesh_integration(derivatives, calc_force, wavepackets);
  }

  xcenergy.potential *= units().Hartree2Energy;
  xcenergy.kinetic *= units().Hartree2Energy;

  for (auto &va : derivatives)
    for (auto &a : va) {
      a *= units().Hartree2Energy * units().Distance2Bohr;
    }

  return XCResult(xcenergy.potential, xcenergy.kinetic, std::move(derivatives));
}

// NOT IMPL YET!
std::tuple<double, double> XCEnergy_nvgpu::regular_mesh_integration(
    std::vector<std::vector<float>> &derivatives, bool ) {
  throw std::string("Regular mesh not supported!");
  /*device_allocation(m_packet_count);

  thrust::copy(host_mem.begin(), host_mem.end(),
  thrust::device_ptr<PacketType>(device_mem));

  compute_density(device_mem);
  auto xcenergy = this->compute_energy(device_mem);

  derivatives = std::move(((calc_force) ? this->compute_derivatives(device_mem)
  : std::vector<std::vector<float>>{})); return std::tuple<double,
  double>(xcenergy.potential, xcenergy.kinetic);*/
}

XCEnergy::Energy XCEnergy_nvgpu::compute_energy(PacketType *) {
  throw std::string("Regular mesh not supported!");
  /*KernelExecutionConfiguration config;
  config.positive = {m_spin_up_index.start, m_spin_up_index.end};

  config.space_delta = {m_config.mesh_step.size.as_struct.x,
  m_config.mesh_step.size.as_struct.y, m_config.mesh_step.size.as_struct.z};
  config.space_shift = {(float) m_config.mesh_start[0], (float)
  m_config.mesh_start[1], (float) m_config.mesh_start[2]};

  std::tie(config.width, config.height, config.depth) =
  std::tie(m_config.mesh_size.size.as_array[0],
                                                                 m_config.mesh_size.size.as_array[1],
                                                                 m_config.mesh_size.size.as_array[2]);
  config.packets_count = m_packet_count;

  dim3 block(8, 8, 4);
  dim3 grid(m_config.mesh_size.size.as_struct.x / block.x,
  m_config.mesh_size.size.as_struct.y / block.y,
            m_config.mesh_size.size.as_struct.z / block.z);

  unsigned long shared_mem_size = block.x * block.y * block.z *
  sizeof(PacketType);

  raw_energy.clear();
  raw_kinetic.clear();

  auto CuError = launch_kernel_eng(KernelConfig{grid, block, shared_mem_size},
  m_config.approximation_device, rho, raw_energy, raw_kinetic, config);

  if (CuError != cudaSuccess) {
    Logger::Error("Error in start cuda kernels. [", __FUNCTION__, __LINE__, "]",
  cudaGetErrorString(CuError)); std::terminate();
  }

  SAFECALL(cudaThreadSynchronize(), "Error in thread sync.");

  auto make_reduction = [this](float *ptr) -> float {
    thrust::device_ptr<float> raw_ptr(ptr);
    return thrust::reduce(raw_ptr, raw_ptr + raw_energy.count(), 0.0f,
  thrust::plus<float>());
  };

  auto potential_energy = make_reduction(raw_energy.dev_ptr());
  auto kinetic_energy = make_reduction(raw_kinetic.dev_ptr());

  return {potential_energy, kinetic_energy};*/
}

std::vector<std::vector<float>>
XCEnergy_nvgpu::compute_derivatives(PacketType const *) {
  throw std::string("Regular mesh not supported!");
  /*KernelExecutionConfigurationDerivs config{};

  config.space_delta = {m_config.mesh_step.size.as_struct.x,
                        m_config.mesh_step.size.as_struct.y,
                        m_config.mesh_step.size.as_struct.z};

  config.space_shift = {(float) m_config.mesh_start[0],
                        (float) m_config.mesh_start[1],
                        (float) m_config.mesh_start[2]};

  std::tie(config.width, config.height, config.depth) =
  std::tie(m_config.mesh_size.size.as_array[0],
                                                                 m_config.mesh_size.size.as_array[1],
                                                                 m_config.mesh_size.size.as_array[2]);
  config.packets_count = m_packet_count;

  dim3 block(8, 8, 4);
  dim3 grid(m_config.mesh_size.size.as_struct.x / block.x,
  m_config.mesh_size.size.as_struct.y / block.y,
            m_config.mesh_size.size.as_struct.z / block.z);

  unsigned long shared_mem_size = block.x * block.y * block.z *
  sizeof(PacketType);

  raw_derivatives.clear();
  auto cuError = launch_kernel_derivs({grid, block, shared_mem_size},
  m_config.approximation_device, rho, packets, raw_derivatives, config);

  if (cuError != cudaSuccess) {
    Logger::Error("Error in start cuda kernels. [", __FUNCTION__, __LINE__, "]",
  cudaGetErrorString(cuError)); std::terminate();
  }

  SAFECALL(cudaThreadSynchronize(), "Error in thread sync.");

  auto make_reduction = [this](float *begin, size_t count) -> float {
    thrust::device_ptr<float> raw_begin(begin);
    thrust::device_ptr<float> raw_end(begin + count);
    return thrust::reduce(raw_begin, raw_end, 0.0f, thrust::plus<float>());
  };

  std::vector<std::vector<float>> derivatives(m_packet_count,
  std::vector<float>(PacketType::derivatives_count(), 0));

  for (auto i = 0u; i < m_packet_count; ++i) {
    auto space_base = raw_derivatives.m_height;
    auto particle_shift = space_base * raw_derivatives.m_width;
    auto ptr = raw_derivatives.dev_ptr() + i * particle_shift;
    for (auto k = 0u; k < raw_derivatives.m_width; ++k) {
      derivatives[i][k] = make_reduction(ptr + k * space_base, space_base);
    }
  }

  return derivatives;*/
}

void XCEnergy_nvgpu::compute_density(XCEnergy::PacketType *) {
  throw std::string("Regular mesh not supported!");
  /*KernelExecutionConfiguration config{};
  config.positive = {m_spin_up_index.start, m_spin_up_index.end};

  config.space_delta = {m_config.mesh_step.size.as_struct.x / 2.0f,
  m_config.mesh_step.size.as_struct.y / 2.0f,
                        m_config.mesh_step.size.as_struct.z / 2.0f};
  config.space_shift = {
      (float) m_config.mesh_start[0],
      (float) m_config.mesh_start[1],
      (float) m_config.mesh_start[2]
  };

  std::tie(config.width, config.height, config.depth) =
  std::tie(rho.total.m_width, rho.total.m_height, rho.total.m_depth);
  config.packets_count = m_packet_count;

  dim3 block(8, 8, 4);
  dim3 grid(config.width / block.x + 1, config.height / block.y + 1,
            config.depth / block.z + 1);

  unsigned long shared_mem_size = block.x * block.y * block.z *
  sizeof(PacketType);
  __CuDftDensityCalculation << < grid, block, shared_mem_size >> > (packets,
  rho, config); auto CuError = cudaGetLastError();

  if (CuError != cudaSuccess) {
    Logger::Error("Error in start cuda kernels. [", __FUNCTION__, __LINE__, "]",
  cudaGetErrorString(CuError)); std::terminate();
  }

  SAFECALL(cudaThreadSynchronize(), "Error in thread sync.");*/
}