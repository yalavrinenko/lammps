//
// Created by yalavrinenko on 24.09.2019.
//
#ifdef COMPUTE_CLASS

ComputeStyle(ke/awpmd,ComputeKEAwpmd)

#else
#ifndef LAMMPS_COMPUTE_KE_AWPMD_H
#define LAMMPS_COMPUTE_KE_AWPMD_H
#include "../compute_ke.h"

namespace LAMMPS_NS{
  class ComputeKEAwpmd: public ComputeKE{
  public:
    ComputeKEAwpmd(LAMMPS* lmp, int argc, char** argv) : ComputeKE(lmp, argc, argv){}

    double compute_scalar() override;
  };
}

#endif //LAMMPS_COMPUTE_KE_AWPMD_H
#endif
