//
// Created by yalavrinenko on 19.09.2019.
//

#include "fix_nvt_awpmd.h"
#include "fix_nve_awpmd.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

LAMMPS_NS::FixNVTAwpmd::FixNVTAwpmd(LAMMPS_NS::LAMMPS *lammps, int argc, char **argv) :
    FixNVT(lammps, argc, argv){
  if (!atom->wavepacket_flag)
    error->all(FLERR,"Fix nve/awpmd requires atom style wavepacket");
}

void LAMMPS_NS::FixNVTAwpmd::initial_integrate(int i) {
  if (atom->mass || atom->rmass) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        double dtfm = dtf / ((atom->mass) ? atom->mass[atom->type[i]] : atom->rmass[atom->type[i]]);
        if (atom->spin[i] != 0) {
          atom->ervel[i] += -dtfm * atom->erforce[i];
          atom->eradius[i] += dtv * atom->ervel[i];
        }
      }
    }
  }
  FixNVT::initial_integrate(i);
}

void LAMMPS_NS::FixNVTAwpmd::final_integrate() {
  if (atom->mass || atom->rmass) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        double dtfm = dtf / ((atom->mass) ? atom->mass[atom->type[i]] : atom->rmass[atom->type[i]]);
        if (abs(atom->spin[i]) != 0) {
          atom->ervel[i] += -dtfm * atom->erforce[i];
        }
      }
    }
  }
  FixNVT::final_integrate();
}
