//
// Created by yalavrinenko on 19.09.2019.
//
#ifdef FIX_CLASS

FixStyle(nvt/awpmd,FixNVTAwpmd)

#else
#ifndef LAMMPS_FIX_NVT_AWPMD_H
#define LAMMPS_FIX_NVT_AWPMD_H

#include "../fix_nvt.h"
namespace LAMMPS_NS {
  class FixNVTAwpmd : public FixNVT {
  public:
    FixNVTAwpmd(class LAMMPS *, int, char **);

    void final_integrate() override;

    void initial_integrate(int i) override;
  };
}

#endif //LAMMPS_FIX_NVT_AWPMD_H
#endif
