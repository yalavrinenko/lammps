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
      mask |= LAMMPS_NS::FixConst::INITIAL_INTEGRATE;
      mask |= LAMMPS_NS::FixConst::FINAL_INTEGRATE;
      return mask;
    }

    void init() override;

    double memory_usage() override;

    void initial_integrate(int i) override;

    ~FixWPMCAwpmd() override;

    double compute_vector(int i) override;

  protected:
    void reject();

    int index_by_tag(int tag) const;

    bool test_step(double dE, double prefactor) const;

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

    union{
      struct {
        mc_step coord;
        mc_step vel;
        mc_step width;
        mc_step pwidth;
      } vars;

      mc_step v[sizeof(vars) / sizeof(mc_step)];
    } steppers;

    size_t const vars_count = sizeof(steppers) / sizeof(mc_step);
    size_t current_stepper = -1;

    double energy_old = std::numeric_limits<double>::max();
    int v_id = -1;

    double target_temperature = 1.0;
  };

}
#endif //LAMMPS_FIX_WMPC_AWPMD_H

#endif