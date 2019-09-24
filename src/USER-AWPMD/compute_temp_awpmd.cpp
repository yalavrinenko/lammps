//
// Created by yalavrinenko on 24.09.2019.
//

#include "compute_temp_awpmd.h"
#include <force.h>
#include <atom.h>
#include <update.h>
#include <error.h>

double LAMMPS_NS::ComputeTempAwpmd::compute_scalar() {
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;
  unsigned nelectrons = 0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && atom->spin[i] != 0) {
        t += atom->ervel[i] * atom->ervel[i] * rmass[i];
        ++nelectrons;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && atom->spin[i] != 0) {
        t += atom->ervel[i] * atom->ervel[i] * mass[type[i]];
        ++nelectrons;
      }
  }

  MPI_Allreduce(&t, &scalar,1,MPI_DOUBLE,MPI_SUM,world);

  //Ekpw = 1 / (2 * ne * k * T)
  auto tscale_factor = force->mvv2e / (force->boltz * nelectrons);
  auto pw_temp_scalar = scalar * tscale_factor;
  auto temp = ComputeTemp::compute_scalar();
  scalar = temp + pw_temp_scalar;
  return scalar;
}
