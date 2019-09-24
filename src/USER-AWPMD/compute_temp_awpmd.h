//
// Created by yalavrinenko on 24.09.2019.
//
#ifdef COMPUTE_CLASS

ComputeStyle(temp/awpmd,ComputeTempAwpmd)

#else
#ifndef LAMMPS_COMPUTE_TEMP_AWPMD_H
#define LAMMPS_COMPUTE_TEMP_AWPMD_H
#include "../compute_temp.h"

namespace LAMMPS_NS{
  class ComputeTempAwpmd: public ComputeTemp{
  public:
    ComputeTempAwpmd(LAMMPS* lmp, int argc, char** argv): ComputeTemp(lmp, argc, argv)  {}

    double compute_scalar() override;
  };
}

#endif //LAMMPS_COMPUTE_TEMP_AWPMD_H
#endif
