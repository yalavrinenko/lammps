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
#include <awpmd-dft.hpp>

namespace LAMMPS_NS{
  class FixDftAwpmd: public Fix{
  public:
    FixDftAwpmd(class LAMMPS *lammps, int i, char **pString);

    int setmask() override;

    ~FixDftAwpmd() override {
      delete xc_energy_;
    }

    void pre_reverse(int i, int i1) override;

  protected:
    XCEnergy* xc_energy_;

    struct {
      double distance_to_bohr = 1.0;
      double hartree_to_energy = 1.0;
    } UnitsScale;
  public:
    double compute_vector(int i) override;

  protected:
    union {
      struct {
        double xc_energy;
        double kinetic_energy;
      } like_vars;

      double like_vector[sizeof(like_vars) / sizeof(double)];
    } output;
  };
}

#endif //LAMMPS_FIX_AWPMD_DFT_H
#endif
