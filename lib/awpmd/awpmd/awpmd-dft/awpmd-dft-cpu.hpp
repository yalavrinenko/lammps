//
// Created by cheshire on 09.03.17.
//

#ifndef AWPMD_DFT_AWPMD_DFT_CPU_HPP
#define AWPMD_DFT_AWPMD_DFT_CPU_HPP

#include "awpmd-dft.hpp"
#include "numeric/adaptive_mesh_integrator.hpp"
#include "numeric/cpu_integrator.hpp"

class XCEnergy_cpu : public XCEnergy {
public:
  XCEnergy_cpu(unsigned int packet_number, DFTConfig const &meshConfig);

  XCResult energy(std::vector<WavePacketInfo> const &wavepackets, bool calc_force) override;

  ~XCEnergy_cpu() override = default;

  struct cell_density{
    AdaptiveMeshCell<double> cell;
    double rho;
    double fraction;

    template <class T>
    cell_density(T const &cell, double rho, double fraction) : cell(cell), rho{rho}, fraction{fraction}{}

    cell_density() = default;
  };

  std::pair<std::vector<cell_density>, double>
  build_density_map(std::vector<WavePacketInfo> const &wavepackets);

  double
  build_density_map(std::vector<WavePacketInfo> const &wavepackets, std::vector<cell_density> &cells, bool center_only = false);

protected:

  void init_xc_approximation() override {}

  void device_allocation(size_t new_size);

  std::tuple<double, double>
  regular_mesh_integration(IApproximation const *approx, std::vector<std::vector<double>> &derivatives, bool calc_force);

  std::tuple<double, double>
  adaptive_mesh_integration(IApproximation const *approx, std::vector<std::vector<double>> &derivatives, bool calc_force);

  using MeshCell = AdaptiveMeshCell<double>;
  AdaptiveMeshIntegrator<MeshCell, Integrator_cpu<MeshCell>> integrator_;

  Integrator_cpu<MeshCell> integration_engine_;

private:
  std::pair<double, double> rho_function(double3 const &r) const;

  NumericType<double, double> xcenergy_function(double x, double y, double z) const;

  NumericType<double, double, double, double> force_function(double x, double y, double z, size_t packet_index) const;
};

#endif //AWPMD_DFT_AWPMD_DFT_CPU_HPP
