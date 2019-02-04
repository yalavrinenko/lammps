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

namespace LAMMPS_NS {

  class FixWPMCAwpmd : public Fix {
  public:
   FixWPMCAwpmd(class LAMMPS *, int, char **);

   void final_integrate() override;

   int setmask() override{
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
   union {
     struct {
       double accept_flag;
       double accepted_energy;
       double step_energy;
     } like_vars;

     double like_vector[sizeof(like_vars) / sizeof(double)];
   } output;

   struct {
     RanPark *coord = nullptr;
     RanPark *particle_index = nullptr;
     RanPark *step_approve_prob = nullptr;
   } uniform;

   struct v_single_particle {
     double coord[3];
     double vel[3];
   };
   v_single_particle state_old;

   void save_state(int p_index);

   void set_state(int p_index, int tag, v_single_particle const &state);

   int index_by_tag(int tag) const;

   double energy_old = std::numeric_limits<double>::max();

   double max_displacement = 0.01;

   int current_particles_index = -1;
   int current_particles_tag = -1;

   int v_id = -1;
  };

}
#endif //LAMMPS_FIX_WMPC_AWPMD_H

#endif