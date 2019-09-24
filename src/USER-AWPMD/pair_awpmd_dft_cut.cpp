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
#include <awpmd-dft-cpu.hpp>
#include "neigh_list.h"
#include <DerivativesFunction.hpp>

LAMMPS_NS::PairAWPMD_DFTCut::PairAWPMD_DFTCut(LAMMPS_NS::LAMMPS *lammps) : PairAWPMDCut(lammps) {
  delete []pvector;
  nextra = 6;
  pvector = new double[nextra];

  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });

  set_units();

  xc_energy_ = new XCEnergy_cpu(electron_count, make_dft_config());
}

void LAMMPS_NS::PairAWPMD_DFTCut::compute(int _i, int _i1) {
  auto one_h = force->mvh2r;
  PairAWPMDCut::compute(_i, _i1);

  auto electrons_count = atom->nlocal + atom->nghost;
  e_sup.clear();
  e_sup.reserve(electrons_count);
  e_sdown.clear();
  e_sdown.reserve(electrons_count);

  for (auto i = 0u; i < electrons_count; ++i){
    if (std::abs(atom->spin[i]) == 1){
      if (atom->spin[i] == 1)
        e_sup.emplace_back(packets[i]);
      else
        e_sdown.emplace_back(packets[i]);
    }
  }

  auto energy = xc_energy_->energy(e_sup, e_sdown, {});
  output.like_vars.xc_energy = energy.eng.potential;
  output.like_vars.kinetic_energy = energy.eng.kinetic;
  force->pair->eng_coul += output.like_vars.xc_energy + output.like_vars.kinetic_energy;

    packet.int2phys_der<eq_second>(dx, dx, dp, dw, pw, 1. / force->mvh2r);
    auto scale = UnitsScale.hartree_to_energy * UnitsScale.distance_to_bohr;
    for (auto k = 0u; k < 3; ++k)
      atom->f[lmp_index][k] += -dx[k] * scale;
    atom->erforce[lmp_index] += *dw * scale;
    atom->ervelforce[lmp_index] += *pw * scale;
  };

  auto packet_up = e_sup.begin();
  auto packet_down = e_sdown.begin();
  auto derive = energy.derivatives.begin();
  for (auto i = 0u; i < electrons_count; ++i) {
    if (std::abs(atom->spin[i]) == 1) {
      derive_to_phys_unit(*packet_up, *derive, i);
      ++packet_up;
      ++derive;
    } else {
      derive_to_phys_unit(*packet_down, *derive, i);
      ++packet_down;
      ++derive;
    }
  }
//  output.like_vars.xc_energy = energy.eng.potential;
//  output.like_vars.kinetic_energy = energy.eng.kinetic;
//  force->pair->eng_coul += output.like_vars.xc_energy + output.like_vars.kinetic_energy;

  //std::cout << "ENG:" << comm->me << "\t" << energy.energy << "\t" << energy.kinetic_energy << "\t" << electrons_count << std::endl;
  pvector[4] = output.like_vars.xc_energy;
  pvector[5] = output.like_vars.kinetic_energy;
}

DFTConfig LAMMPS_NS::PairAWPMD_DFTCut::make_dft_config() {
  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });
  const double SPACE_MESH_SCALE = 1.5 * UnitsScale.distance_to_bohr;
  DFTConfig mesh_config;
  mesh_config.packet_number = electron_count;
  mesh_config.calc_derivs = false;
  mesh_config.approximation = new LSDA();

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

  mesh_config.mesh_size.size.as_struct = {100 / grid_size.x, 100 / grid_size.y, 100 / grid_size.z};

  return mesh_config;
}

LAMMPS_NS::PairAWPMD_DFTCut::PairAWPMD_DFTCut(LAMMPS_NS::LAMMPS *lammps, XCEnergy* xc_energy_ptr) : PairAWPMDCut(lammps) {
  delete []pvector;
  nextra = 7;
  pvector = new double[nextra];

  xc_energy_ = xc_energy_ptr;

  set_units();
}

void LAMMPS_NS::PairAWPMD_DFTCut::settings(int i, char **pString) {
  PairAWPMDCut::settings(i, pString);
  overlap_derivs = DerivsFunction::GetFunctions();
}
