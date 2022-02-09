//
// Created by yalavrinenko on 09.02.2022.
//
#ifdef FIX_CLASS

FixStyle(nph/wpmd,FixNPHWpmd)

#else
#ifndef LAMMPS_FIX_NPH_WPMD_H
#define LAMMPS_FIX_NPH_WPMD_H

#include "fix_nh_wpmd.h"

namespace LAMMPS_NS{
  class FixNPHWpmd: public FixNHWpmd{
  public:
    FixNPHWpmd(struct LAMMPS *lammps, int i, char **pString);
  };
}

#endif //LAMMPS_FIX_NPH_WPMD_H
#endif
