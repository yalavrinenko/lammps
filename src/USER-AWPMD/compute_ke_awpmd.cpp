//
// Created by yalavrinenko on 24.09.2019.
//

#include "compute_ke_awpmd.h"
#include <atom.h>
#include <update.h>
#include <force.h>
double LAMMPS_NS::ComputeKEAwpmd::compute_scalar() {
  invoked_scalar = update->ntimestep;

  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double ke = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && atom->spin[i] != 0)
        ke += rmass[atom->type[i]] * (atom->ervel[i] * atom->ervel[i]);
  } else {
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && atom->spin[i] != 0)
        ke += atom->mass[atom->type[i]] * (atom->ervel[i] * atom->ervel[i]);
  }
  double tmp_scalar = 0;
  MPI_Allreduce(&ke,&tmp_scalar,1,MPI_DOUBLE,MPI_SUM,world);
  tmp_scalar *= 0.5 * force->mvv2e;
  scalar = tmp_scalar + ComputeKE::compute_scalar();
  return scalar;
}
