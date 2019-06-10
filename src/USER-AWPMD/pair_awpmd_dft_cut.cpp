//
// Created by yalavrinenko on 10.06.19.
//

#include "pair_awpmd_dft_cut.h"
#include <atom.h>
#include <force.h>
#include <DataTypes.hpp>
#include <LSDA.hpp>
#include <domain.h>
#include <comm.h>
#include <wpmd_split.h>
#include <style_pair.h>

LAMMPS_NS::PairAWPMD_DFTCut::PairAWPMD_DFTCut(LAMMPS_NS::LAMMPS *lammps) : PairAWPMDCut(lammps) {
  delete []pvector;
  nextra = 7;
  pvector = new double[nextra];

  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });
  const double SPACE_MESH_SCALE = 1.5;
  DFTConfig mesh_config;
  mesh_config.packet_number = electron_count;
  mesh_config.calc_derivs = false;
  mesh_config.approximation = new LSDA();
  mesh_config.mesh_size.size.as_struct = {100, 100, 100};

  double3 domain_size {(domain->boxhi[0] - domain->boxlo[0]) * SPACE_MESH_SCALE,
                       (domain->boxhi[1] - domain->boxlo[1]) * SPACE_MESH_SCALE,
                       (domain->boxhi[2] - domain->boxlo[2]) * SPACE_MESH_SCALE};

  uint3 grid_size {static_cast<unsigned int>(comm->procgrid[0]),
                   static_cast<unsigned int>(comm->procgrid[1]),
                   static_cast<unsigned int>(comm->procgrid[2])};

  uint3 my_grid_pos {static_cast<unsigned int>(comm->myloc[0]),
                     static_cast<unsigned int>(comm->myloc[1]),
                     static_cast<unsigned int>(comm->myloc[2])};

  mesh_config.space_size = {domain_size.x / grid_size.x,
                            domain_size.y / grid_size.y,
                            domain_size.z / grid_size.z};

  mesh_config.space_shift = {domain->boxlo[0] * SPACE_MESH_SCALE + my_grid_pos.x * mesh_config.space_size.x,
                             domain->boxlo[1] * SPACE_MESH_SCALE + my_grid_pos.y * mesh_config.space_size.y,
                             domain->boxlo[2] * SPACE_MESH_SCALE + my_grid_pos.z * mesh_config.space_size.z};

  xc_energy_ = new XCEnergy(electron_count, mesh_config);

  UnitsScale.distance_to_bohr = 1.0 / (0.52917721092 * force->angstrom);
  UnitsScale.hartree_to_energy = 627.509474; //only for real
}

void LAMMPS_NS::PairAWPMD_DFTCut::compute(int i, int i1) {
  PairAWPMDCut::compute(i, i1);

  auto energy = xc_energy_->energy(wpmd->wp[0], wpmd->wp[1], {});
  output.like_vars.xc_energy = UnitsScale.hartree_to_energy * energy.eng.potential;
  output.like_vars.kinetic_energy = UnitsScale.hartree_to_energy * energy.eng.kinetic;

  force->pair->eng_coul += output.like_vars.xc_energy + output.like_vars.kinetic_energy;
  pvector[5] = output.like_vars.xc_energy;
  pvector[6] = output.like_vars.kinetic_energy;
}
