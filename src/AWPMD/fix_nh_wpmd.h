//
// Created by yalavrinenko on 09.02.2022.
//

#ifndef LAMMPS_FIX_NH_WPMD_H
#define LAMMPS_FIX_NH_WPMD_H

#include <fix_nh.h>

namespace LAMMPS_NS{
  class FixNHWpmd: public FixNH{
  public:
    FixNHWpmd(class LAMMPS *, int, char **);
  protected:
    void nve_x() override;

    void nve_v() override;

    void nh_v_temp() override;
  };
}

#endif //LAMMPS_FIX_NH_WPMD_H
