//
// Created by yalavrinenko on 09.02.2022.
//
#ifdef FIX_CLASS

FixStyle(npt/wpmd,FixNPTWpmd)

#else
#ifndef LAMMPS_FIX_NPT_WPMD_H
#define LAMMPS_FIX_NPT_WPMD_H
#include "fix_nh_wpmd.h"

namespace LAMMPS_NS{
  class FixNPTWpmd: public FixNHWpmd{
  public:
    FixNPTWpmd(struct LAMMPS *lammps, int i, char **pString);
  };
}

#endif //LAMMPS_FIX_NPT_WPMD_H
#endif