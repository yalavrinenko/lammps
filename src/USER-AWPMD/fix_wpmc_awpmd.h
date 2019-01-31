//
// Created by yalavrinenko on 28.12.18.
//
#ifdef FIX_CLASS

FixStyle(wpmc/awpmd,FixWPMCAwpmd)

#else

#ifndef LAMMPS_FIX_WMPC_AWPMD_H
#define LAMMPS_FIX_WMPC_AWPMD_H
#include "fix.h"
#include "pair.h"
#include "random_park.h"

namespace LAMMPS_NS {

  class FixWPMCAwpmd : public Fix {
  public:
   FixWPMCAwpmd(class LAMMPS *, int, char **);

   int setmask() override{
     int mask = 0;
     mask |= LAMMPS_NS::FixConst::INITIAL_INTEGRATE;
     mask |= LAMMPS_NS::FixConst::POST_FORCE;
     return mask;
   }

    void init() override;

    double memory_usage() override;

    void initial_integrate(int i) override;

    void post_force(int i) override;

    ~FixWPMCAwpmd() override;

    double compute_vector(int i) override;

  protected:
   double dtv, dtf;
   double *step_respa;
   int mass_require;

   Pair* awpmd_pair;
  };

}
#endif //LAMMPS_FIX_WMPC_AWPMD_H

#endif