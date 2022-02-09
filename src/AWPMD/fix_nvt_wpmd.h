//
// Created by yalavrinenko on 19.09.2019.
//
#ifdef FIX_CLASS

FixStyle(nvt/wpmd,FixNVTWpmd)

#else
#ifndef LAMMPS_FIX_NVT_AWPMD_H
#define LAMMPS_FIX_NVT_AWPMD_H

#include "fix_nh_wpmd.h"
namespace LAMMPS_NS {
  class FixNVTWpmd : public FixNHWpmd {
  public:
    FixNVTWpmd(class LAMMPS *, int, char **);

    ~FixNVTWpmd() override = default;
  };
}

#endif //LAMMPS_FIX_NVT_AWPMD_H
#endif
