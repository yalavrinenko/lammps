#include "numeric_gpu.cuh"
#include "cuApproxTableObject.hpp"
#include "gpu_exceptions.h"
#include "gpu_functors.hpp"
#include "gpu_integration.hpp"


#include "../DataTypes.hpp"
#include "../xcfunctionals/LDA.hpp"
#include "../xcfunctionals/LSDA.hpp"
#include "../xcfunctionals/ModLDA.hpp"
#include "../xcfunctionals/TableApproximation.hpp"
#include "xc_init_kernels.hpp"
#include "awpmd-dft-nvgpu.hpp"
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

void XCEnergy_nvgpu::init_xc_approximation() {
  Logger::Info("Init XCApprox on GPU");
  IApproximation **ptr_ref;

  if (m_config.approximation_device != nullptr) {
    ::nvgpu_dft_remove_xc_approximation<<<1, 1>>>(m_config.approximation_device);
    m_config.approximation_device = nullptr;
  }

  SAFECALL(cudaMalloc(&ptr_ref, sizeof(IApproximation *)),
           "Error in cudaMalloc for xc_approx");

  switch (m_config.approximation->Type) {
    case ApproxType::T_LDA:
      ::nvgpu_dft_init_xc_approximation<LDA>
      <<<1, 1>>>(ptr_ref, 0.738558766f, -0.01554534543482745f, 20.4562557f);
      break;
    case ApproxType::T_LDA_2:
      ::nvgpu_dft_init_xc_approximation<ModLDA><<<1, 1>>>(ptr_ref);
      break;
    case ApproxType::T_LSDA:
      ::nvgpu_dft_init_xc_approximation<LSDA><<<1, 1>>>(ptr_ref);
      break;
    case ApproxType::T_VOID:
      ::nvgpu_dft_init_xc_approximation<VoidApproximation><<<1, 1>>>(ptr_ref);
    case ApproxType::T_DUMMY:
      ::nvgpu_dft_init_xc_approximation<DummyApproximation><<<1, 1>>>(ptr_ref);
      break;
  }
  SAFECALL(cudaDeviceSynchronize(), "Error in thread sync for xc_approx");

  auto **ptr = new IApproximation *[1];
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

      ::nvgpu_dft_init_xc_approximation<TableApproximation<gpuApproxTablesObject>>
      <<<1, 1>>>(ptr_ref,
                 *(XCApproximation::gpu_approximation_table()));
      SAFECALL(cudaDeviceSynchronize(), "Error in thread sync for xc_approx");
      SAFECALL(cudaMemcpy(ptr, ptr_ref, sizeof(IApproximation *), cudaMemcpyDeviceToHost),
               "Error in cudaMemcpy for xc_approx");
      Logger::Info("New XCApprox GPU Ptr:", *ptr);
      m_config.approximation_device = *ptr;
    }
  }
//#endif

  Logger::Info("Check...");
  ::nvgpu_dft_xc_check<<<1, 8>>>(m_config.approximation_device);
  SAFECALL(cudaDeviceSynchronize(), "Error in thread sync for check xc_approx");

  Logger::Info("Done.");
}

XCEnergy_nvgpu::XCEnergy_nvgpu(unsigned int packet_number,
                               DFTConfig const &meshConfig, long gpu_id) {
  Logger::ModuleName("AWPMD-DFT GPU");

  Logger::Warning("AWPMD-DFT VERSION FROM", __DATE__);

  cudaSetDevice(static_cast<int>(gpu_id));

  m_packet_count = packet_number;
  Logger::Info("Init dft module on GPU =", gpu_id,
               " Packet number =", m_packet_count);

  setConfig(meshConfig);

  Logger::Info("Allocate",
               static_cast<double>(m_packet_count) * sizeof(GaussPacket<float>) / (1024.0 * 1024.0),
               "Mb in host memory for packet.");
  host_allocation(m_packet_count);

  Logger::Info("Allocate",
               static_cast<double>(m_packet_count) * sizeof(GaussPacket<float>) / (1024.0 * 1024.0),
               "Mb in device memory for packet.");
  device_allocation(m_packet_count);
}

XCEnergy_nvgpu::~XCEnergy_nvgpu() {
  if (XCApproximation::gpu_approximation_table() != nullptr)
    XCApproximation::gpu_approximation_table()->free();
}

void XCEnergy_nvgpu::device_allocation(size_t new_size) {

  if (dev_mem_size < new_size) {
    SAFECALL(cudaFree(device_mem), "Error in memory reallocation.");
    SAFECALL(cudaMalloc(&device_mem, new_size * sizeof(GaussPacket<float>)),
             "Error in allocation memory for packet");
    dev_mem_size = new_size;
  }
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
        regular_mesh_integration(derivatives, wavepackets, calc_force);
  } else {
    std::tie(xcenergy.potential, xcenergy.kinetic) =
        adaptive_mesh_integration(derivatives, wavepackets, calc_force);
  }

  xcenergy.potential *= units().Hartree2Energy;
  xcenergy.kinetic *= units().Hartree2Energy;

  for (auto &va : derivatives)
    for (auto &a : va) {
      a *= units().Hartree2Energy * units().Distance2Bohr;
    }

  return {xcenergy.potential, xcenergy.kinetic, std::move(derivatives)};
}


std::tuple<double, double> XCEnergy_nvgpu::adaptive_mesh_integration(std::vector<std::vector<float>> &derivatives,
                                                                     const vector<WavePacketInfo> &wavepackets,
                                                                     bool calc_force) {

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
    compute_forces(derivatives, wavepackets, energy_eval);
  }

  return std::make_tuple<double, double>((double)thrust::get<0>(xck_energy),
                                         (double)thrust::get<1>(xck_energy));
}

std::tuple<double, double> XCEnergy_nvgpu::regular_mesh_integration(
    std::vector<std::vector<float>> &derivatives, const vector<WavePacketInfo> &wavepackets, bool calc_force) {

  Integrator_nvgpu<MeshCell> energy_integration_engine_;
  integrator_.refine_mesh<float>(MeshCell::RangeType(m_config.mesh_start.size.as_struct),
                          MeshCell::RangeType(m_config.mesh_fin.size.as_struct),
                          m_config.mesh_step.size.as_struct);

  energy_functor_sic<MeshCell, PacketType> energy_eval{};
  energy_eval.approx = m_config.approximation_device;

  auto xck_energy = integrator_.integrate<NumericType<float, float>>(
      energy_integration_engine_, host_mem, energy_eval);

  if (calc_force) {
    compute_forces(derivatives, wavepackets, energy_eval);
  }

  return std::make_tuple<double, double>((double)thrust::get<0>(xck_energy),
                                         (double)thrust::get<1>(xck_energy));
}

void
XCEnergy_nvgpu::compute_forces(std::vector<std::vector<float>> &derivatives, const vector<WavePacketInfo> &wavepackets,
                               energy_functor_sic<MeshCell, PacketType> const &energy_eval) {
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
    }; //???

    derivatives[packet_index][0] += static_cast<float>(round_float(dx));
    derivatives[packet_index][1] += static_cast<float>(round_float(dy));
    derivatives[packet_index][2] += static_cast<float>(round_float(dz));
    derivatives[packet_index][3] += static_cast<float>(round_float(dw));
  };
}
