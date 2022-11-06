//
// Created by yalavrinenko on 28.12.18.
//
#ifdef FIX_CLASS

FixStyle(mc/wpmd,FixMCAwpmd)

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
#include "pair_awpmd_cut.h"

namespace LAMMPS_NS {

  class FixMCAwpmd : public Fix {
  public:
    FixMCAwpmd(class LAMMPS *, int, char **);

    void final_integrate() override;

    int setmask() override {
      int mask = 0;
      mask |= LAMMPS_NS::FixConst::PRE_FORCE;
      mask |= LAMMPS_NS::FixConst::FINAL_INTEGRATE;
      //mask |= LAMMPS_NS::FixConst::INITIAL_INTEGRATE;
      return mask;
    }
    void initial_integrate(int i) override;

    void init() override;

    double memory_usage() override;

    ~FixMCAwpmd() override;

    double compute_vector(int i) override;

    void pre_force(int i) override;

  protected:

    void update_ghosts();

    void init_mc_steppers(int argc, char** argv);

    union {
      struct {
        double accept_flag;
        double accepted_energy;
        double step_energy;
        double accepted_count;
        double rejected_count;
        double stepper_id;
      } like_vars;

      double like_vector[sizeof(like_vars) / sizeof(double)];
    } output;

    MCStepperSet steppers;

    Compute *temp, *pe, *norm, *ke;

    double energy_old = std::numeric_limits<double>::max();

    double target_temperature = 1.0;

    PairAWPMD *awpmd = nullptr;
    bool is_first = true;
    bool use_norm = false;
    bool use_awpmd_ke = false;  // use awpmd-calculated kinetic energies for electrons instead of \sum p^2/(2m)
    
  };

}
#endif //LAMMPS_FIX_WMPC_AWPMD_H

#endif