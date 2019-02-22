//
// Created by yalavrinenko on 28.12.18.
//
#ifdef FIX_CLASS

FixStyle(wpmc/awpmd,FixNVEAwpmd)

#else

#ifndef LAMMPS_FIX_WMPC_AWPMD_H
#define LAMMPS_FIX_WMPC_AWPMD_H
#include "fix.h"
#include "pair_awpmd_cut.h"

namespace LAMMPS_NS {

  class FixWPMCAwpmd : public Fix {
  public:
   FixWPMCAwpmd(class LAMMPS *, int, char **);

   int setmask() override{
     return 0;
   }

   void init() override {}

   void initial_integrate(int) override {}

   void final_integrate() override {}

   void initial_integrate_respa(int, int, int) override {}

   void final_integrate_respa(int, int) override {}

   void reset_dt() override {}

  protected:
   double dtv, dtf;
   double *step_respa;
   int mass_require;

   PairAWPMDCut* awpmd_pair;
  };

}
#endif //LAMMPS_FIX_WMPC_AWPMD_H

#endif