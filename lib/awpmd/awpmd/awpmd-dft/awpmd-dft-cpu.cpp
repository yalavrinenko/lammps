//
// Created by cheshire on 10.03.17.
//
#include "awpmd-dft-cpu.hpp"
#include "DataTypes.hpp"
#include "numeric/numeric.h"
#include <cstring>
#include <omp.h>

using std::cout;
using std::endl;

XCEnergy_cpu::XCEnergy_cpu(unsigned int packet_number, DFTConfig const &meshConfig) {
  Logger::ModuleName("AWPMD-DFT CPU");

  m_packet_count = packet_number;
  Logger::Info("Init dft module. Packet number =", m_packet_count);

  setConfig(meshConfig);

  host_allocation(packet_number);
  device_allocation(packet_number);

//  throw std::string("CPU Code is not tested!");
}

void XCEnergy_cpu::device_allocation(size_t ) {
  device_mem = &host_mem[0];
  dev_mem_size = host_mem.size();
}

std::tuple<double, double>
XCEnergy_cpu::regular_mesh_integration(IApproximation const *, std::vector<std::vector<double>> &derivatives,
                                       bool calc_force) {
  double3 r{};
  double3 dh = {m_config.mesh_step.size.as_struct.x, m_config.mesh_step.size.as_struct.y,
                m_config.mesh_step.size.as_struct.z};
  auto &sh = m_config.mesh_start.size.as_struct;

  double xcenergy = 0;
  double kinetic = 0;

  for (auto i = 0u; i < m_config.mesh_size.size.as_struct.x; ++i)
    for (auto j = 0u; j < m_config.mesh_size.size.as_struct.y; ++j)
      for (auto k = 0u; k < m_config.mesh_size.size.as_struct.z; ++k) {

        r = {double(sh.x + i * dh.x), double(sh.y + j * dh.y), double(sh.z + k * dh.z)};
        double xc, ke;
        std::tie(xc, ke) = IntegrationMethods::simpson_integration<NumericType<double, double>>(
            [this](double x, double y, double z) {
              return xcenergy_function(x, y, z);
            }, r, dh);

        xcenergy += xc;
        kinetic += ke;

        //XC-Energy derivs by args_i
        //dE/dAi = dF/d(rho) * dOii/dAi
        if (calc_force) {
          for (auto packet_index = 0u; packet_index < m_packet_count; ++packet_index) {
            using output_type = NumericType<double, double, double, double>;
            double dw, dx, dy, dz;
            std::tie(dx, dy, dz, dw) = IntegrationMethods::simpson_integration<output_type>(
                [this, packet_index](double x, double y, double z) {
                  return this->force_function(x, y, z, packet_index);
                }, r, dh);
            derivatives[packet_index][0] += dx;
            derivatives[packet_index][1] += dy;
            derivatives[packet_index][2] += dz;
            derivatives[packet_index][3] += dw;
          }
        }
      }
  return std::tuple<double, double>{xcenergy, kinetic};
}

std::pair<double, double> XCEnergy_cpu::rho_function(double3 const &r) const {
  double rho_v_up = 0, rho_v_down = 0;
  for (const auto & packet : host_mem) {
    if (packet.rdot(r) < (packet.width * 5.0) * (packet.width * 5.0)) {
      if (packet.spin == ElectronSpin::E_UP)
        rho_v_up += packet.overlap_at(r);
      else
        rho_v_down += packet.overlap_at(r);
    }
  }
  double fraction_v;
  if (rho_v_up + rho_v_down < std::numeric_limits<double>::epsilon())
    fraction_v = 0.0;
  else
    fraction_v = (rho_v_up - rho_v_down) / (rho_v_up + rho_v_down);

  return std::pair<double, double>{rho_v_up + rho_v_down, fraction_v};
}

NumericType<double, double> XCEnergy_cpu::xcenergy_function(double x, double y, double z) const {
  IApproximation *approx = m_config.approximation;
  auto r = double3{x, y, z};

  auto sic_xc = 0.0;
  auto sic_kin = 0.0;

  double rho_v_up = 0, rho_v_down = 0;

  for (const auto & packet : host_mem) {
    if (packet.rdot(r) < (packet.width * 5.0) * (packet.width * 5.0)) {
      auto overlap = packet.overlap_at(r);
      sic_xc += approx->energy(overlap, 1);
      sic_kin += approx->kinetic(overlap, 1);
      if (packet.spin == ElectronSpin::E_UP)
        rho_v_up += overlap;
      else
        rho_v_down += overlap;
    }
  }
  auto fraction_v = (rho_v_up - rho_v_down) / (rho_v_up + rho_v_down);
  auto rho_total = rho_v_up + rho_v_down;

  auto xc = approx->energy(rho_total, fraction_v) - sic_xc;
  auto kin = approx->kinetic(rho_total, fraction_v) - sic_kin;

  return NumericType<double, double>{xc, kin};
}

std::tuple<double, double>
XCEnergy_cpu::adaptive_mesh_integration(IApproximation const *, std::vector<std::vector<double>> &derivatives,
                                        bool calc_force) {
  integrator_.refine_mesh(MeshCell::MeshPoint(m_config.mesh_start.size.as_struct),
                          MeshCell::MeshPoint(m_config.mesh_fin.size.as_struct), [this](MeshCell::MeshPoint const &a, MeshCell::MeshPoint const &b) {
        return this->refine_mesh_condition(host_mem.begin(), host_mem.end(), a, b);
      }, {m_config.nodes, m_config.node_rank});

  auto xck_energy = integrator_.integrate<NumericType<double, double>>(integration_engine_, [this](double x, double y, double z) {
    return this->xcenergy_function(x, y, z);
  });

  if (calc_force) {
    for (auto packet_index = 0u; packet_index < m_packet_count; ++packet_index) {
      double dw, dx, dy, dz;
      std::tie(dx, dy, dz, dw) = integrator_.integrate<NumericType<double, double, double, double>>(
          integration_engine_,
          [this, packet_index](double x, double y, double z) {
            return this->force_function(x, y, z, packet_index);
          });
      derivatives[packet_index][0] += dx;
      derivatives[packet_index][1] += dy;
      derivatives[packet_index][2] += dz;
      derivatives[packet_index][3] += dw;
    }
  }

  return xck_energy;
}

NumericType<double, double, double, double>
XCEnergy_cpu::force_function(double x, double y, double z, size_t packet_index) const {
  auto rho_ = rho_function({x, y, z});
  auto r = double3{x, y, z};
  auto overlap = host_mem[packet_index].overlap_at(r);

  auto spin = host_mem[packet_index].spin;

  auto dw = m_config.approximation->derives(rho_.first, rho_.second, host_mem[packet_index].dw(r), spin)
      - m_config.approximation->derives(overlap, static_cast<float>(spin), host_mem[packet_index].dw(r), spin);
  auto dx = m_config.approximation->derives(rho_.first, rho_.second, host_mem[packet_index].dr<0>(r), spin);
  auto dy = m_config.approximation->derives(rho_.first, rho_.second, host_mem[packet_index].dr<1>(r), spin);
  auto dz = m_config.approximation->derives(rho_.first, rho_.second, host_mem[packet_index].dr<2>(r), spin);
  return {dx, dy, dz, dw};
}

std::pair<std::vector<XCEnergy_cpu::cell_density>, double>
XCEnergy_cpu::build_density_map(std::vector<WavePacketInfo> const &wavepackets) {
  host_allocation(wavepackets.size());
  device_allocation(wavepackets.size());

  copy_packet_to_host(wavepackets);

  integrator_.refine_mesh(MeshCell::MeshPoint(m_config.mesh_start.size.as_struct),
                          MeshCell::MeshPoint(m_config.mesh_fin.size.as_struct), [this](MeshCell::MeshPoint const &a, MeshCell::MeshPoint const &b) {
        return this->refine_mesh_condition(host_mem.begin(), host_mem.end(), a, b);
      }, {m_config.nodes, m_config.node_rank});


  vector<XCEnergy_cpu::cell_density> density_array(integrator_.mesh_cells().size());
  std::transform(integrator_.mesh_cells().begin(), integrator_.mesh_cells().end(), density_array.begin(),
                 [](AdaptiveMeshCell<double> const &cell) {
                   return XCEnergy_cpu::cell_density{cell, 0.0, 0.0};
                 });

  auto overall_density = build_density_map(wavepackets, density_array);
  return {density_array, overall_density};
}

double XCEnergy_cpu::build_density_map(std::vector<WavePacketInfo> const &wavepackets,
                                           vector<cell_density> &cells, bool center_only) {
  host_allocation(wavepackets.size());
  device_allocation(wavepackets.size());

  copy_packet_to_host(wavepackets);

  double overall_rho = 0.0;

  for (auto &cell : cells) {
    auto f = [this, center_only](double3 a, double3 b) {
      double rho_v_up = 0, rho_v_down = 0;
      for (auto & packet : host_mem) {
        if (!center_only) {
          if (packet.spin == ElectronSpin::E_UP)
            rho_v_up += packet.I(a, b);
          else
            rho_v_down += packet.I(a, b);
        }
        else {
          rho_v_up += packet.in_range(a, b);
        }
      }
      double fraction_v;
      if (rho_v_up + rho_v_down < std::numeric_limits<double>::epsilon())
        fraction_v = 0.0;
      else
        fraction_v = (rho_v_up - rho_v_down) / (rho_v_up + rho_v_down);

      return NumericType<double, double>{rho_v_up + rho_v_down, fraction_v};
    };

    double3 begin = {cell.cell.begin().x, cell.cell.begin().y, cell.cell.begin().z};
    double3 end = {cell.cell.end().x, cell.cell.end().y, cell.cell.end().z};

    auto rho_ = f(begin, end);

    auto dh = cell.cell.dh();
    double dv = dh.x * dh.y * dh.z;
    cell.rho = std::get<0>(rho_) / (dv);
    cell.fraction = std::get<1>(rho_);

    overall_rho += cell.rho;
  }

  return overall_rho;
}

XCEnergy::XCResult XCEnergy_cpu::energy(std::vector<WavePacketInfo> const &wavepackets, bool calc_force) {
  if (m_config.approximation->Type == ApproxType::T_VOID)
    return XCResult(0, 0, {});

  m_packet_count = wavepackets.size();
  std::vector<std::vector<double>> derivatives((calc_force) ? m_packet_count : 0,
                                               std::vector<double>(PacketType::derivatives_count(), 0));

  host_allocation(m_packet_count);
  device_allocation(m_packet_count);

  //Initialization. Copy wavepacket to host memory
  copy_packet_to_host(wavepackets);

  double xcenergy{}, kinetic{};

  if (!m_config.use_adaptive_mesh) {
    std::tie(xcenergy, kinetic) = this->regular_mesh_integration(m_config.approximation, derivatives, calc_force);
  } else {
    std::tie(xcenergy, kinetic) = this->adaptive_mesh_integration(m_config.approximation, derivatives, calc_force);
  }

  xcenergy *= units().Hartree2Energy;
  kinetic *= units().Hartree2Energy;

  for (auto &va: derivatives)
    for (auto &a:va) {
      //if (std::abs(a) < 1e-7) a = 0; //bad fix
      a *= units().Hartree2Energy * units().Distance2Bohr;
    }

  std::vector<std::vector<float>> float_derivatives(derivatives.size(),
                                                    std::vector<float>(PacketType::derivatives_count(), 0));

  for (auto i = 0ul; i < derivatives.size(); ++i) {
    std::transform(derivatives[i].begin(), derivatives[i].end(), float_derivatives[i].begin(),
                   [](double &v) { return v; });
  }

  last_result = XCResult{static_cast<float>(xcenergy),
                         static_cast<float>(kinetic),
                         std::move(float_derivatives)};
  return last_result;
}
