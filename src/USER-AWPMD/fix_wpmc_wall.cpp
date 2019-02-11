//
// Created by yalavrinenko on 11.02.19.
//

#include "fix_wpmc_wall.h"

LAMMPS_NS::FixWallAwpmd::FixWallAwpmd(LAMMPS_NS::LAMMPS *lammps, int i, char **pString) : Fix(lammps, i, pString) {
  std::cout << "FIX WALL" << std::endl;
}

LAMMPS_NS::FixWallAwpmd::~FixWallAwpmd() {
  std::cout << "UNFIX WALL" << std::endl;
}

int LAMMPS_NS::FixWallAwpmd::setmask() {
  return 0;
}
