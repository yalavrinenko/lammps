//
// Created by yalavrinenko on 29.05.19.
//

#include "fix_awpmd_dft.h"

int LAMMPS_NS::FixDftAwpmd::setmask() {
  return LAMMPS_NS::FixConst::PRE_REVERSE;
}

LAMMPS_NS::FixDftAwpmd::FixDftAwpmd(LAMMPS_NS::LAMMPS *lammps, int i, char **pString) : Fix(lammps, i, pString) {

}
