//
// Created by yalavrinenko on 29.11.2019.
//
#ifdef COMPUTE_CLASS

ComputeStyle(denprof,ComputeDensityAwpmd)

#else
#ifndef LAMMPS_COMPUTE_ELECTRON_DENSITY_H
#define LAMMPS_COMPUTE_ELECTRON_DENSITY_H
#include <array>
#include "compute.h"
#include <awpmd-dft-cpu.hpp>
#include <region.h>
#include <memory.h>

class ComputeDensityAwpmd: public LAMMPS_NS::Compute{
public:
  ComputeDensityAwpmd(LAMMPS_NS::LAMMPS* lmp, int argc, char** argv);

  double compute_scalar() override;

  void compute_array() override;

  void init() override;

  ~ComputeDensityAwpmd() override;

protected:
  struct plane_info {
    double x, y, z;
    double dx, dy, dz;
    int Nx, Ny, Nz;
  };

  void create_cell_list(plane_info const &info);
  void make_packets();

  std::vector<XCEnergy_cpu::cell_density> cells_;

  LAMMPS_NS::Region *region_ = nullptr;
  std::array<int, 3> nbins_{{1, 1, 1}};
  std::array<double, 3> axis_{};
  std::array<bool, 3> vary_axis_{{false, false, false}};
  std::vector<WavePacket> packets_;

  double scalef_ = 1.0;
  bool use_center_ = false;

  DFTConfig config_;

  std::unique_ptr<XCEnergy_cpu> xc_energy_ = nullptr;

};

#endif //LAMMPS_COMPUTE_ELECTRON_DENSITY_H
#endif
