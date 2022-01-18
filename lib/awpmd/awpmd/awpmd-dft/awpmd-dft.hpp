//
// Created by cheshire on 09.03.17.
//

#ifndef AWPMD_DFT_AWPMD_DFT_HPP
#define AWPMD_DFT_AWPMD_DFT_HPP

#include "DataTypes.hpp"
#include "utils/DerivativesFunction.hpp"
#include <math.h>

#include "utils/CCP4_DMap.h"
#include "utils/SpaceMesh.hpp"
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

class XCEnergy {
public:
  struct WavePacketInfo{
    //class WavePacket packet{};
    double const *coord;
    double width;
    ElectronSpin spin{};
    bool calc_force = false;

//    WavePacketInfo(WavePacket packet, ElectronSpin spin, bool calc_force):
//        packet{std::move(packet)}, spin{spin}, calc_force{calc_force}{
//    }

    WavePacketInfo(double const *coord, double width, ElectronSpin spin, bool calc_force):
        coord{coord}, width{width},
        spin{spin}, calc_force{calc_force}{
    }
  };

  struct Energy {
    double potential;
    double kinetic;
  };

  class XCResult {
  public:
    Energy eng{};
    double &energy = eng.potential;
    double &kinetic_energy = eng.kinetic;

    std::vector<std::vector<float>> derivatives;

    XCResult() = default;

    XCResult(double _potential_energy, double _kinetic_energy, std::vector<std::vector<float>> _deriv) :
        eng{_potential_energy, _kinetic_energy}, derivatives{std::move(_deriv)} {
    }

    XCResult &operator=(XCResult const &r) {
      this->eng = r.eng;
      this->derivatives = r.derivatives;
      return *this;
    }
  };

  virtual XCEnergy::XCResult energy(std::vector<WavePacketInfo> const &wavepackets, bool calc_force) = 0;

  void setConfig(DFTConfig const &config) {
    m_config = config;

    Logger::Info("DFT config.");
    for (unsigned int &dim_space_size : m_config.mesh_size.size.as_array) {
      if (dim_space_size % 16 != 0) {
        dim_space_size += 16;
        dim_space_size -= (dim_space_size % 16);
      }
    }

    Logger::Info("Mesh topo (X, Y, Z):");

    for (auto i : {0, 1, 2}){
      m_config.mesh_step.size.as_array[i] = (m_config.mesh_fin.size.as_array[i] - m_config.mesh_start.size.as_array[i]) / m_config.mesh_size.size.as_array[i];
      Logger::Info("\t", m_config.mesh_start.size.as_array[i], " --- ", m_config.mesh_fin.size.as_array[i]);
    }

    if (m_config.use_adaptive_mesh) {
      Logger::Info("Use adaptive mesh.");
    } else {
      Logger::Info("Mesh size:", m_config.mesh_size.size.as_struct.x, 'x', m_config.mesh_size.size.as_struct.y, 'x',
                   m_config.mesh_size.size.as_struct.z);
      Logger::Info("Mesh step:", m_config.mesh_step.size.as_struct.x, 'x', m_config.mesh_step.size.as_struct.y, 'x',
                   m_config.mesh_step.size.as_struct.z);
    }
    Logger::Info("Approximation:", (int) m_config.approximation->Type);

    this->init_xc_approximation();
  }

  DFTConfig const &get_dft_config() const {
    return m_config;
  }

  XCResult const &get_last() const {
    return last_result;
  }

  virtual ~XCEnergy() = default;

  UnitsSystem<float> &units() { return m_config.units; }

protected:

  using PacketType = GaussPacket<float>;

  virtual void init_xc_approximation() = 0;

  void host_allocation(size_t new_size) {
    host_mem.resize(new_size);
  }

  void copy_packet_to_host(std::vector<struct WavePacketInfo> const &wp) {
    std::transform(wp.begin(), wp.end(), host_mem.begin(), [this](WavePacketInfo const &i){
      return PacketType(i.coord, i.width, units(), i.spin);
    });
  }

  template<typename range_t, typename iter_t>
  bool refine_mesh_condition(iter_t const& packets_begin, iter_t const& packets_end, range_t a, range_t b) const {
    auto check_packet = [this](PacketType const &p, range_t a, range_t b) {
      auto const min_dh = p.width * m_config.min_cell;

      auto const max_shift = (p.width * m_config.max_distance) * (p.width * m_config.max_distance);

      using real_t = decltype(a.x);
      auto is_in_range = [](real_t a, real_t x, real_t b) { return a <= x && x <= b; };
      bool in_range = is_in_range(a.x, p.r[0], b.x) &&
                      is_in_range(a.y, p.r[1], b.y) &&
                      is_in_range(a.z, p.r[2], b.z);

      bool not_min = b.x - a.x > min_dh;

      range_t c = {static_cast<real_t >((a.x + b.x) / 2.0),
                   static_cast<real_t >((a.y + b.y) / 2.0),
                   static_cast<real_t >((a.z + b.z) / 2.0)};
      bool near_packet = (p.rdot(c) - static_cast<real_t >(2.0 * (b.x - a.x) * (b.x - a.x)) < max_shift);

      return (in_range || near_packet) && not_min;
    };

    for (auto packet = packets_begin; packet != packets_end; ++packet) {
      if (check_packet(*packet, a, b))
        return true;
    }
    return false;

  }

  std::vector<GaussPacket<float>> host_mem;

  GaussPacket<float> *device_mem = nullptr;
  size_t dev_mem_size{0};

  unsigned int m_packet_count{};
  DFTConfig m_config;

  XCResult last_result;
};

#endif //AWPMD_DFT_AWPMD_DFT_HPP
