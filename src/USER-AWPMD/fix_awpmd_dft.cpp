//
// Created by yalavrinenko on 29.05.19.
//

#include "fix_awpmd_dft.h"
#include "atom.h"
#include <LSDA.hpp>
#include <DataTypes.hpp>
#include <comm.h>
#include <domain.h>
#include <force.h>
#include <style_fix.h>

int LAMMPS_NS::FixDftAwpmd::setmask() {
  return LAMMPS_NS::FixConst::PRE_REVERSE;
}

LAMMPS_NS::FixDftAwpmd::FixDftAwpmd(LAMMPS_NS::LAMMPS *lammps, int i, char **pString) : Fix(lammps, i, pString) {
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

void LAMMPS_NS::FixDftAwpmd::pre_reverse(int, int) {
  auto one_h=force->mvh2r;

  auto electrons_count = atom->nlocal + atom->nghost;
  std::vector<WavePacket> e_sup, e_sdown;
  e_sup.reserve(electrons_count), e_sdown.reserve(electrons_count);

  for (auto i = 0u; i < electrons_count; ++i){
    if (std::abs(atom->spin[i]) == 1){
      WavePacket packet;

      double width = atom->eradius[i] * UnitsScale.distance_to_bohr;
      Vector_3 r{atom->x[i][0], atom->x[i][1], atom->x[i][2]}, p{atom->v[i][0], atom->v[i][1], atom->v[i][2]};
      r *= UnitsScale.distance_to_bohr;
      p *= UnitsScale.distance_to_bohr * one_h * atom->mass[i];

      double pw = atom->ervel[i];
      pw *= UnitsScale.distance_to_bohr * one_h * atom->mass[i];

      packet.init(width, r, p, pw);

      if (atom->spin[i] == 1)
        e_sup.emplace_back(packet);
      else
        e_sdown.emplace_back(packet);
    }
  }

}
