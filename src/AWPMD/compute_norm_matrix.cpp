//
// Created by Ilya Valuev on 20.10.2022.
//

#include "compute_norm_matrix.h"
#include "wpmd_split.h"
#include <cstring>
#include <domain.h>
#include <error.h>
#include <atom.h>
#include <force.h>
#include <update.h>
#include <DataTypes.hpp>


double ComputeNormAwpmd::compute_scalar() {
  
  if (invoked_scalar != update->ntimestep) {

    scalar = 0;
    double result = wppair->awpmd()->norm_matrix_detl(0) + wppair->awpmd()->norm_matrix_detl(1) ;
    MPI_Allreduce(&result, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

    invoked_scalar = update->ntimestep;
  }

  return scalar;
}

ComputeNormAwpmd::ComputeNormAwpmd(LAMMPS_NS::LAMMPS *lmp, int argc, char **argv) : Compute(lmp, argc, argv) {
  //array_flag = true;
  scalar_flag = true;
  wppair = dynamic_cast<LAMMPS_NS::PairAWPMD *>(force->pair);
  if(!wppair)
    error->all(FLERR, "Norm matrix compute requires pair_awpmd");
}





ComputeNormAwpmd::~ComputeNormAwpmd() {
  //memory->destroy(array);
}
