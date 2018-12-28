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

   int setmask();

   virtual void init();

   virtual void initial_integrate(int);

   virtual void final_integrate();

   void initial_integrate_respa(int, int, int);

   void final_integrate_respa(int, int);

   void reset_dt();

  protected:
   double dtv, dtf;
   double *step_respa;
   int mass_require;

   PairAWPMDCut *awpmd_pair;
  };

}
#endif //LAMMPS_FIX_WMPC_AWPMD_H

#endif