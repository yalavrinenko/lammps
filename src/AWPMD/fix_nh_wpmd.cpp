//
// Created by yalavrinenko on 09.02.2022.
//

#include "fix_nh_wpmd.h"

#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

LAMMPS_NS::FixNHWpmd::FixNHWpmd(LAMMPS_NS::LAMMPS *lammps, int argc, char **argv) :
    FixNH(lammps, argc, argv){
  if (!atom->wavepacket_flag)
    error->all(FLERR,"Fix nve/awpmd requires atom style wavepacket");
}

void LAMMPS_NS::FixNHWpmd::nve_x() {
  FixNH::nve_x();

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      if (atom->spin[i] != 0) {
        atom->eradius[i] += dtv * atom->ervel[i];
      }
    }
  }
}

void LAMMPS_NS::FixNHWpmd::nve_v() {
  FixNH::nve_v();

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      if (atom->spin[i] != 0) {
        double dtfm = dtf / atom->mass[atom->type[i]];
        if (abs(atom->spin[i]) != 0)
          atom->ervel[i] += -dtfm * atom->erforce[i];
      }
    }
  }
}

void LAMMPS_NS::FixNHWpmd::nh_v_temp() {
  FixNH::nh_v_temp();

  for (int i = 0; i < atom->nlocal; i++)
    if (atom->mask[i] & groupbit)
      if (atom->spin[i] != 0) atom->ervel[i] *= factor_eta;
}
