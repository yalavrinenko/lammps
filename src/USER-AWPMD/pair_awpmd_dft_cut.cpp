//
// Created by yalavrinenko on 10.06.19.
//
#ifdef AWPMD_ENABLE_DFT
#include "pair_awpmd_dft_cut.h"
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

LAMMPS_NS::PairAWPMD_DFTCut::PairAWPMD_DFTCut(LAMMPS_NS::LAMMPS *lammps) : PairAWPMDCut(lammps) {
  delete []pvector;
  nextra = 6;
  pvector = new double[nextra];
}

void LAMMPS_NS::PairAWPMD_DFTCut::compute(int _i, int _i1) {
  auto one_h = force->mvh2r;
  PairAWPMDCut::compute(_i, _i1);

  auto electrons_count = atom->nlocal + atom->nghost;
  e_sup.clear();
  e_sdown.clear();

  double self_ee = 0;
  for (auto i = 0u; i < electrons_count; ++i){
    if (std::abs(atom->spin[i]) == 1){
      if (atom->spin[i] == 1)
        e_sup.emplace_back(packets[i]);
      else
        e_sdown.emplace_back(packets[i]);

      //self_ee += this->wpmd->interaction_ee_single(packets[i], packets[i], nullptr, nullptr);
    }
  }

  auto energy = xc_energy_->energy(e_sup, e_sdown, calc_force_);
  output.like_vars.xc_energy = energy.eng.potential + self_ee ;
  output.like_vars.kinetic_energy = energy.eng.kinetic;
  force->pair->eng_coul += output.like_vars.xc_energy + output.like_vars.kinetic_energy;

  pvector[4] = output.like_vars.xc_energy;
  pvector[5] = output.like_vars.kinetic_energy;

  if (calc_force_) {
    auto force_sup_it = energy.derivatives.up();
    auto force_sdown_it = energy.derivatives.down();
    for (auto i = 0; i < electrons_count; ++i) {
      if (std::abs(atom->spin[i]) == 1) {
        if (atom->spin[i] == 1) {
          tally_electron_force(i, *(force_sup_it++));
        } else {
          tally_electron_force(i, *(force_sdown_it++));
        }
      }
    }
  }
}

DFTConfig LAMMPS_NS::PairAWPMD_DFTCut::make_dft_config(int nargs, char **pString) {
  bool is_daptive_mesh = false;
  unsigned int MeshSize = 50;

  DFTConfig mesh_config;
  for (int i = 1; i < nargs; i++){
    if (std::strcmp(pString[i], "adaptive") == 0) {
      is_daptive_mesh = true;
    }

    if (std::strcmp(pString[i], "min_cell_size") == 0)
      mesh_config.min_cell = force->numeric(FLERR, pString[i+1]);

    if (std::strcmp(pString[i], "max_distance") == 0)
      mesh_config.max_distance = force->numeric(FLERR, pString[i+1]);

    if (std::strcmp(pString[i], "regular") == 0) {
      is_daptive_mesh = false;
      MeshSize = force->numeric(FLERR, pString[i + 1]);
    }

    if (std::strcmp(pString[i], "dynamic") == 0)
      calc_force_ = true;
  }

  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });
  const double SPACE_MESH_SCALE = UnitsScale.distance_to_bohr * ((is_daptive_mesh) ? 1.5 : 1.5);

  mesh_config.packet_number = electron_count;
  mesh_config.approximation = new LSDA();
  mesh_config.use_adaptive_mesh = is_daptive_mesh;
  mesh_config.calc_force = calc_force_;
  mesh_config.node_rank = comm->me;

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

  mesh_config.mesh_size.size.as_struct = {MeshSize / grid_size.x, MeshSize / grid_size.y, MeshSize / grid_size.z};

  mesh_config.units.Hartree2Energy =  627.509474;
  mesh_config.units.Distance2Bohr = 1.0 / (0.52917721092 * force->angstrom);

  return mesh_config;
}

LAMMPS_NS::PairAWPMD_DFTCut::PairAWPMD_DFTCut(LAMMPS_NS::LAMMPS *lammps, XCEnergy* xc_energy_ptr) : PairAWPMDCut(lammps) {
  delete []pvector;
  nextra = 6;
  pvector = new double[nextra];

  xc_energy_ = xc_energy_ptr;
}

void LAMMPS_NS::PairAWPMD_DFTCut::settings(int i, char **pString) {
  PairAWPMDCut::settings(i, pString);

  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });

  xc_energy_ = new XCEnergy_cpu(electron_count, make_dft_config(i, pString));
}

void LAMMPS_NS::PairAWPMD_DFTCut::tally_electron_force(unsigned electron_id, std::vector<float> const& force_array) {
  atom->f[electron_id][0] += force_array[0];
  atom->f[electron_id][1] += force_array[1];
  atom->f[electron_id][2] += force_array[2];
  atom->erforce[electron_id] += force_array[3];

//  for (auto &f: force_array)
//    std::cout << f << " ";
//  std::cout << '\n';
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
