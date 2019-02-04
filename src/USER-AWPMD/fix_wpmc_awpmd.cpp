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

    uniform.coord = new RanPark(lmp, 42);
    uniform.particle_index = new RanPark(lmp, 42);
    uniform.step_approve_prob = new RanPark(lmp, 42);

    v_id = input->variable->find(args[3]);
    if (v_id == -1)
      error->all(FLERR, "Fix wpmc/awpmd requires a valid variable style");
  }

  void FixWPMCAwpmd::init() {
  }

  double FixWPMCAwpmd::memory_usage() {
    return sizeof(energy_old) + sizeof(uniform) + sizeof(output) + sizeof(state_old);
  }

  void FixWPMCAwpmd::initial_integrate(int i) {
    auto npart = atom->nlocal;

    current_particles_index = static_cast<int>(uniform.particle_index->uniform() * (npart - 1));
    current_particles_tag = atom->tag[current_particles_index];

    save_state(current_particles_index);

    v_single_particle s_new{};
    for (auto k = 0; k < 3; ++k){
      s_new.coord[k] = atom->x[current_particles_index][k] + (2.0 * uniform.coord->uniform()  - 1.0) * max_displacement;
      s_new.vel[k] = atom->v[current_particles_index][k] + (2.0 * uniform.coord->uniform()  - 1.0) * max_displacement;
    }

    set_state(current_particles_index, current_particles_tag, s_new);
  }

  FixWPMCAwpmd::~FixWPMCAwpmd() = default;

  double FixWPMCAwpmd::compute_vector(int i) {
    return output.like_vector[i];
  }

  void FixWPMCAwpmd::final_integrate() {
    auto reservoir_temperature = 1.0;
    auto beta = 1.0/(force->boltz*reservoir_temperature);

    auto energy_new = input->variable->compute_equal(v_id);

    if (uniform.step_approve_prob->uniform() < exp(beta*(energy_old - energy_new))){
      output.like_vars.accept_flag = 1;
      energy_old = energy_new;
    } else {
      set_state(current_particles_index, current_particles_tag, state_old);
      output.like_vars.accept_flag = 0;
    }

    output.like_vars.step_energy = energy_new;
    output.like_vars.accepted_energy = energy_old;
  }

  void FixWPMCAwpmd::save_state(int p_index) {
    for (auto k = 0; k < 3; ++k){
      state_old.coord[k] = atom->x[p_index][k];
      state_old.vel[k] = atom->v[p_index][k];
    }
  }

  void FixWPMCAwpmd::set_state(int p_index, int tag, const FixWPMCAwpmd::v_single_particle &state) {
    if (tag != atom->tag[p_index]){
      p_index = index_by_tag(tag);
    }

    if (p_index == -1)
      error->all(FLERR, "Invalid particle index");

    for (auto k = 0; k < 3; ++k){
      atom->x[p_index][k] = state.coord[k];
      atom->v[p_index][k] = state.vel[k];
    }
  }

  int FixWPMCAwpmd::index_by_tag(int tag) const {
    auto npart = atom->nlocal + atom->nghost;
    for (auto i = 0u; i < npart; ++i)
      if (atom->tag[i] == tag)
        return i;
    return -1;
  }

}
