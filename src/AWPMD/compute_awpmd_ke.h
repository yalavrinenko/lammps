//
// Created by Ilya Valuev on 05.11.2022.
//
#ifdef COMPUTE_CLASS

ComputeStyle(awpmd_ke,ComputeAwpmdKE)

#else
#ifndef LAMMPS_COMPUTE_AWPMD_KE_H
#define LAMMPS_COMPUTE_AWPMD_KE_H
#include "compute.h"
#include <memory.h>
#include <array>
#include "pair_awpmd_cut.h"

/// Class to compute kinetic energy for AWPMD
class ComputeAwpmdKE: public LAMMPS_NS::Compute{
public:
  ComputeAwpmdKE(LAMMPS_NS::LAMMPS* lmp, int argc, char** argv);

  double compute_scalar() override;

  //void compute_array() override;

  void init() override {};

  ~ComputeAwpmdKE() override;

  AWPMD_split *awpmd() { return wppair->awpmd(); }

protected:
  LAMMPS_NS::PairAWPMD *wppair = nullptr;

};

#endif //LAMMPS_COMPUTE_AWPMD_KE_H
#endif
