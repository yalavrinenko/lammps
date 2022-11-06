//
// Created by Ilya Valuev on 05.11.2022.
//

#include "compute_awpmd_ke.h"
#include "wpmd_split.h"
#include <cstring>
#include <domain.h>
#include <error.h>
#include <atom.h>
#include <force.h>
#include <update.h>
#include <DataTypes.hpp>


double ComputeAwpmdKE::compute_scalar() {
  
  if (invoked_scalar != update->ntimestep) {

    scalar = 0;
    double result = wppair->awpmd()->get_kin_energy();
    MPI_Allreduce(&result, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

    invoked_scalar = update->ntimestep;
  }

  return scalar;
}

ComputeAwpmdKE::ComputeAwpmdKE(LAMMPS_NS::LAMMPS *lmp, int argc, char **argv) : Compute(lmp, argc, argv) {
  //array_flag = true;
  scalar_flag = true;
  wppair = dynamic_cast<LAMMPS_NS::PairAWPMD *>(force->pair);
  if(!wppair)
    error->all(FLERR, "awpmd_ke compute requires pair_awpmd");
}





ComputeAwpmdKE::~ComputeAwpmdKE() {
  //memory->destroy(array);
}
