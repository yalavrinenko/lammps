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

    for (auto i = 0u; i < vars_count; ++i){
      steppers.v[i].active = true;
      steppers.v[i].shift = 0.1;
      steppers.v[i].type = (mc_step::mc_type)i;
    }

    target_temperature = force->numeric(FLERR, args[4]) * force->boltz;
    output.like_vars.accepted_count = output.like_vars.rejected_count = 0.0;
  }

  void FixWPMCAwpmd::init() {
  }

  double FixWPMCAwpmd::memory_usage() {
    return sizeof(energy_old) + sizeof(random) + sizeof(output) + sizeof(steppers);
  }

  void FixWPMCAwpmd::initial_integrate(int i) {
    auto npart = atom->nlocal;
    auto particle_index = static_cast<int>(random->uniform() * (npart - 1));

    for (auto &stepper : steppers.v){
      if (stepper.active){
        stepper.index = particle_index;
        stepper.tag = atom->tag[particle_index];

        double **src = nullptr;

        switch (stepper.type){
          case mc_step::mc_type ::coord:
          case mc_step::mc_type ::vel:
            src = (stepper.type == mc_step::mc_type::coord) ? atom->x : atom->v;
            std::copy(src[particle_index], src[particle_index] + 3, stepper.old);

            for (auto k = 0u; k < 3; ++k)
              src[particle_index][k] += (2.0 * random->uniform() - 1.0) * stepper.shift;

            break;

          case mc_step::mc_type ::width:
            if (atom->spin[stepper.index] != 0) {
              stepper.old[0] = atom->eradius[particle_index];
              atom->eradius[particle_index] += (2.0 * random->uniform() - 1.0) * stepper.shift;
            }
            break;

          case mc_step::mc_type ::pwidth:
            if (atom->spin[stepper.index] != 0) {
              stepper.old[0] = atom->ervel[particle_index];
              atom->ervel[particle_index] += (2.0 * random->uniform() - 1.0) * stepper.shift;
            }
            break;
        };
      }
    }

  }

  FixWPMCAwpmd::~FixWPMCAwpmd() = default;

  double FixWPMCAwpmd::compute_vector(int i) {
    return output.like_vector[i];
  }

  void FixWPMCAwpmd::final_integrate() {
    auto energy_new = input->variable->compute_equal(v_id);
    auto accept_local = test_step(energy_new - energy_old, 1.0);

    auto accept_all = accept_local;
    if (accept_all){
      output.like_vars.accept_flag = 1;
      ++output.like_vars.accepted_count;
      energy_old = energy_new;
    } else {
      reject();
      output.like_vars.accept_flag = 0;
      ++output.like_vars.rejected_count;
    }

    output.like_vars.step_energy = energy_new;
    output.like_vars.accepted_energy = energy_old;
  }

  void FixWPMCAwpmd::reject() {
    for (auto &stepper : steppers.v) {
      if (stepper.active) {
        if (stepper.tag != atom->tag[stepper.index]) {
          stepper.index = index_by_tag(stepper.tag);
        }

        double **src = nullptr;
        switch (stepper.type) {
          case mc_step::mc_type::coord:
          case mc_step::mc_type::vel:
            src = (stepper.type == mc_step::mc_type::coord) ? atom->x : atom->v;
            std::copy(std::begin(stepper.old), std::end(stepper.old), src[stepper.index]);
            break;

          case mc_step::mc_type::width:
            if (atom->spin[stepper.index] != 0)
              atom->eradius[stepper.index] = stepper.old[0];
            break;

          case mc_step::mc_type::pwidth:
            if (atom->spin[stepper.index] != 0)
              atom->ervel[stepper.index] = stepper.old[0];
            break;
        };
      }
    }

  }

  int FixWPMCAwpmd::index_by_tag(int tag) const {
    auto npart = atom->nlocal + atom->nghost;
    for (auto i = 0u; i < npart; ++i)
      if (atom->tag[i] == tag)
        return i;
    return -1;
  }

  bool FixWPMCAwpmd::test_step(double dE, double prefactor) const {
    if (dE < 0)
      return true;
    else {
      double r1 = prefactor * exp(-dE / target_temperature);
      double r2 = random->uniform();

      if (r1 < 1) {
        if (r1 >= r2)
          return true;
      } else
        return true;
    }
    return false;
  }

}
