//
// Created by yalavrinenko on 28.12.18.
//

#include <style_fix.h>
#include "fix_wpmc_awpmd.h"
#include "atom.h"
#include "force.h"
#include "error.h"
#include "memory.h"
#include <random>
#include <cmath>
#include "modify.h"

namespace LAMMPS_NS {
  FixWPMCAwpmd::FixWPMCAwpmd(LAMMPS_NS::LAMMPS *lmp, int narg, char **args) :
      Fix(lmp, narg, args){
    //if (!atom->wavepacket_flag)
      //error->all(FLERR, "Fix wpmc/awpmd requires atom style wavepacket");

    vector_flag = 1;
    size_vector = sizeof(output) / sizeof(double);

    global_freq = 1;
    extvector = 0;
    time_depend = 1;

    nevery = 1;

    random = new RanPark(lmp, 42);

    v_id = input->variable->find(args[3]);
    if (v_id == -1)
      error->all(FLERR, "Fix wpmc/awpmd requires a valid variable style");

    target_temperature = force->numeric(FLERR, args[4]) * force->boltz;
    output.like_vars.accepted_count = output.like_vars.rejected_count = 0.0;

    auto debug_stepper_count = 8u;
    for (auto i = 0u; i < debug_stepper_count; ++i){
      steppers.add(lmp, (stepper_type)i);
      steppers.get(i).max_shift = 0.1;
      steppers.get(i).engine.setT(target_temperature);
    }

  }

  void FixWPMCAwpmd::init() {
  }

  double FixWPMCAwpmd::memory_usage() {
    return sizeof(energy_old) + sizeof(random) + sizeof(output);
  }


  FixWPMCAwpmd::~FixWPMCAwpmd() = default;

  double FixWPMCAwpmd::compute_vector(int i) {
    return output.like_vector[i];
  }

  void FixWPMCAwpmd::final_integrate() {
    auto energy_new = input->variable->compute_equal(v_id);
    this->output.like_vars.accept_flag = steppers.current().engine.test(energy_new - energy_old, 1.);

    if (output.like_vars.accept_flag == 1){
      //accepted
      energy_old = energy_new;
    } else {
      //rejected
      steppers.current().system->restore((size_t)atom->nlocal);
    }

    output.like_vars.step_energy = energy_new;
    output.like_vars.accepted_count += output.like_vars.accept_flag;
    output.like_vars.accepted_energy += energy_old;
    output.like_vars.rejected_count += (output.like_vars.accept_flag == 0);

  }

  void FixWPMCAwpmd::pre_force(int i) {
    steppers.current().system->save((size_t)atom->nlocal);
    steppers.current().system->make((size_t)atom->nlocal, steppers.current());
  }


}
