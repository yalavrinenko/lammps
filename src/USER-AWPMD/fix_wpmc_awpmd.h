//
// Created by yalavrinenko on 28.12.18.
//
#ifdef FIX_CLASS

FixStyle(wpmc/awpmd,FixWPMCAwpmd)

#else

#ifndef LAMMPS_FIX_WMPC_AWPMD_H
#define LAMMPS_FIX_WMPC_AWPMD_H

#include <limits>
#include "fix.h"
#include "pair.h"
#include "random_park.h"
#include "compute.h"
#include "variable.h"
#include "mc_utils.h"

namespace LAMMPS_NS {

  class FixWPMCAwpmd : public Fix {
  public:
    FixWPMCAwpmd(class LAMMPS *, int, char **);

    void final_integrate() override;

    int setmask() override {
      int mask = 0;
      mask |= LAMMPS_NS::FixConst::PRE_FORCE;
      mask |= LAMMPS_NS::FixConst::FINAL_INTEGRATE;
      return mask;
    }

    void init() override;

    double memory_usage() override;

    ~FixWPMCAwpmd() override;

    double compute_vector(int i) override;

    void pre_force(int i) override;

  protected:

    void init_mc_steppers(int argc, char** argv);

    union {
      struct {
        double accept_flag;
        double accepted_energy;
        double step_energy;
        double accepted_count;
        double rejected_count;
      } like_vars;

      double like_vector[sizeof(like_vars) / sizeof(double)];
    } output;

    RanPark* random;

    MCStepperSet steppers;

    double energy_old = std::numeric_limits<double>::max();
    int v_id = -1;

    double target_temperature = 1.0;
  };

}
#endif //LAMMPS_FIX_WMPC_AWPMD_H

#endif