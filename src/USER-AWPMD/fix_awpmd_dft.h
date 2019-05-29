//
// Created by yalavrinenko on 29.05.19.
//
#ifdef FIX_CLASS
FixStyle(dft/awpmd,FixDftAwpmd)
#else
#ifndef LAMMPS_FIX_AWPMD_DFT_H
#define LAMMPS_FIX_AWPMD_DFT_H

#include <limits>
#include "fix.h"
#include "pair.h"
#include "random_park.h"
#include "compute.h"
#include "variable.h"
#include "mc_utils.h"

namespace LAMMPS_NS{
  class FixDftAwpmd: public Fix{
  public:
    FixDftAwpmd(class LAMMPS *lammps, int i, char **pString);

  public:
    int setmask() override;
  };
}

#endif //LAMMPS_FIX_AWPMD_DFT_H
#endif
