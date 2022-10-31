//
// Created by Ilya Valuev on 20.10.2022.
//
#ifdef COMPUTE_CLASS

ComputeStyle(normmatr,ComputeNormAwpmd)

#else
#ifndef LAMMPS_COMPUTE_NORM_MATRIX_H
#define LAMMPS_COMPUTE_NORM_MATRIX_H
#include "compute.h"
#include <memory.h>
#include <array>
#include "pair_awpmd_cut.h"

class ComputeNormAwpmd: public LAMMPS_NS::Compute{
public:
  ComputeNormAwpmd(LAMMPS_NS::LAMMPS* lmp, int argc, char** argv);

  double compute_scalar() override;

  //void compute_array() override;

  void init() override {};

  ~ComputeNormAwpmd() override;

  AWPMD_split *awpmd() { return wppair->awpmd(); }

protected:
  LAMMPS_NS::PairAWPMD *wppair = nullptr;

};

#endif //LAMMPS_COMPUTE_NORM_MATRIX_H
#endif
