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
#include <memory>

namespace LAMMPS_NS {
  class PairAWPMDCut;
  class FixWallAwpmd : public Fix {
  public:
    FixWallAwpmd(class LAMMPS *lammps, int i, char **pString);

    ~FixWallAwpmd() override;

    int setmask() override;

  private:
  public:
    void pre_reverse(int i, int i1) override;

  private:
    std::unique_ptr<BoxHamiltonian> construct_box(char **pString, double half_box_size);

    class PairAWPMDCut* m_pair;

    std::unique_ptr<BoxHamiltonian> box = nullptr;
    double wall_energy = 0;
  public:
    double compute_scalar() override;
  };
}

#endif //LAMMPS_FIX_WPMC_WALL_H
#endif
