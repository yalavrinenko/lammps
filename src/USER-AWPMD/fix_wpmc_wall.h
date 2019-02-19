//
// Created by yalavrinenko on 11.02.19.
//
#ifdef FIX_CLASS

FixStyle(wall/awpmd,FixWallAwpmd)

#else
#ifndef LAMMPS_FIX_WPMC_WALL_H
#define LAMMPS_FIX_WPMC_WALL_H

#include "fix.h"
#include <box_hamiltonian.h>

namespace LAMMPS_NS {
  class PairAWPMDCut;
  class FixWallAwpmd : public Fix {
  public:
    FixWallAwpmd(class LAMMPS *lammps, int i, char **pString);

    ~FixWallAwpmd() override;

    int setmask() override;

  private:
    BoxHamiltonian construct_box(char** pString);

    class PairAWPMDCut* m_pair;
  };
}

#endif //LAMMPS_FIX_WPMC_WALL_H
#endif
