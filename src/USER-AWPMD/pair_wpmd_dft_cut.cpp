//
// Created by yalavrinenko on 10.06.19.
//
#ifdef AWPMD_ENABLE_DFT
#include "pair_wpmd_dft_cut.h"
#include <atom.h>
#include <force.h>
#include <DataTypes.hpp>
#include <LSDA.hpp>
#include <domain.h>
#include <comm.h>
#include <wpmd_split.h>
#include <awpmd-dft-cpu.hpp>
#include "neigh_list.h"
#include <DerivativesFunction.hpp>
#include <awpmd-dft.hpp>
#include <cstring>

LAMMPS_NS::PairAWPMD_DFTCut::PairAWPMD_DFTCut(LAMMPS_NS::LAMMPS *lammps) : PairAWPMD_DFTCut(lammps, nullptr) {
}

void LAMMPS_NS::PairAWPMD_DFTCut::compute(int _i, int _i1) {
  auto one_h = force->mvh2r;
  PairWPMD::compute(_i, _i1);

  auto electrons_count = atom->nlocal + atom->nghost;
  electrons.clear();

  double self_ee = 0;
  for (auto i = 0u; i < electrons_count; ++i){
    if (std::abs(atom->spin[i]) == 1){
      electrons.emplace_back(packets[i], ElectronSpin(atom->spin[i]), (calc_force_ && i < atom->nlocal) );
    }
  }

  auto energy = xc_energy_->energy(electrons, calc_force_);
  output.like_vars.xc_energy = energy.eng.potential ;
  output.like_vars.kinetic_energy = energy.eng.kinetic;
  force->pair->eng_coul += output.like_vars.xc_energy + output.like_vars.kinetic_energy;

  pvector[4] = output.like_vars.xc_energy;
  pvector[5] = output.like_vars.kinetic_energy;

  if (calc_force_) {
    auto force_it = energy.derivatives.begin();
    for (auto i = 0; i < electrons_count; ++i) {
      if (std::abs(atom->spin[i]) == 1 && i < atom->nlocal) {
        tally_electron_force(i, *(force_it++));
      }
    }
  }
}

DFTConfig LAMMPS_NS::PairAWPMD_DFTCut::make_dft_config(int nargs, char **pString) {
  bool is_daptive_mesh = false;
  unsigned int MeshSize = 50;

  DFTConfig mesh_config;
  auto get_next_float = [pString, this](size_t i) {
    return force->numeric(FLERR, pString[i+1]);
  };
  for (int i = 1; i < nargs; i++){
    if (std::strcmp(pString[i], "adaptive") == 0) {
      is_daptive_mesh = true;
    }

    if (std::strcmp(pString[i], "min_cell_size") == 0)
      mesh_config.min_cell = get_next_float(i);

    if (std::strcmp(pString[i], "max_distance") == 0)
      mesh_config.max_distance = get_next_float(i);

    if (std::strcmp(pString[i], "regular") == 0) {
      is_daptive_mesh = false;
      MeshSize = get_next_float(i);
    }

    if (std::strcmp(pString[i], "dynamic") == 0)
      calc_force_ = true;

    if (std::strcmp(pString[i], "force_mesh_bins") == 0){
      mesh_config.force_cell_bins = get_next_float(i);
    }
  }

  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });

  mesh_config.units.Hartree2Energy =  627.509474;
  mesh_config.units.Distance2Bohr = 1.0 / (0.52917721092 * force->angstrom);

  const double SPACE_MESH_SCALE = mesh_config.units.Distance2Bohr;

  mesh_config.packet_number = electron_count;
  mesh_config.approximation = new LSDA();
  mesh_config.use_adaptive_mesh = is_daptive_mesh;
  mesh_config.calc_force = calc_force_;
  mesh_config.node_rank = comm->me;
  mesh_config.nodes = comm->nprocs;

  double3 domain_size {domain->boxhi[0] - domain->boxlo[0],
                       domain->boxhi[1] - domain->boxlo[1],
                       domain->boxhi[2] - domain->boxlo[2]};

  uint grid_size[] = {static_cast<unsigned int>(comm->procgrid[0]),
                   static_cast<unsigned int>(comm->procgrid[1]),
                   static_cast<unsigned int>(comm->procgrid[2])};

  uint my_grid_pos[] = {static_cast<unsigned int>(comm->myloc[0]),
                     static_cast<unsigned int>(comm->myloc[1]),
                     static_cast<unsigned int>(comm->myloc[2])};

  double block_size[] = {domain_size.x / grid_size[0],
                            domain_size.y / grid_size[1],
                            domain_size.z / grid_size[2]};

  double const EDGE_MULT[] = {SPACE_MESH_SCALE, SPACE_MESH_SCALE, SPACE_MESH_SCALE};

  double const EDGE_BORDER_MULT = 3.0;

  for (auto i = 0; i < 3; ++i){
    mesh_config.mesh_start.size.as_array[i] = domain->boxlo[i] + my_grid_pos[i] * block_size[i];
    mesh_config.mesh_fin.size.as_array[i] = mesh_config.mesh_start.size.as_array[i] + block_size[i];

    if (my_grid_pos[i] == 0)
      mesh_config.mesh_start.size.as_array[i] *= EDGE_BORDER_MULT;

    if (my_grid_pos[i] == grid_size[i] - 1)
      mesh_config.mesh_fin.size.as_array[i] *= EDGE_BORDER_MULT;

    mesh_config.mesh_start[i] *= EDGE_MULT[i];
    mesh_config.mesh_fin[i] *= EDGE_MULT[i];
  }

  mesh_config.mesh_size.size.as_struct = {MeshSize / grid_size[0], MeshSize / grid_size[1], MeshSize / grid_size[2]};

  return mesh_config;
}

LAMMPS_NS::PairAWPMD_DFTCut::PairAWPMD_DFTCut(LAMMPS_NS::LAMMPS *lammps, XCEnergy* xc_energy_ptr) : PairWPMD(lammps) {
  delete []pvector;
  nextra = 6;
  pvector = new double[nextra];

  xc_energy_ = xc_energy_ptr;
}

void LAMMPS_NS::PairAWPMD_DFTCut::settings(int i, char **pString) {
  PairWPMD::settings(i, pString);

  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });

  xc_energy_ = new XCEnergy_cpu(electron_count, make_dft_config(i, pString));
}

void LAMMPS_NS::PairAWPMD_DFTCut::tally_electron_force(unsigned electron_id, std::vector<float> const& force_array) {
  atom->f[electron_id][0] += force_array[0];
  atom->f[electron_id][1] += force_array[1];
  atom->f[electron_id][2] += force_array[2];
  atom->erforce[electron_id] += force_array[3];
}

double LAMMPS_NS::PairAWPMD_DFTCut::wpmd_kinetic() const {
  double ke = 0.0;
  auto dot_v = [this](size_t i){
    return atom->v[i][0] * atom->v[i][0] +
           atom->v[i][1] * atom->v[i][1] +
           atom->v[i][2] * atom->v[i][2];
  };
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->spin[i] != 0)
      ke += atom->mass[atom->type[i]] * (atom->ervel[i] * atom->ervel[i] + dot_v(i));

  return ke * 0.5 * force->mvv2e;
}
#endif
