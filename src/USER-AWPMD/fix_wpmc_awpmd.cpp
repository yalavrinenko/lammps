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
#include <cstring>
#include "comm.h"
#include "atom_vec.h"
#include <future>
using namespace std::string_literals;

#include "timemetrics.h"

namespace LAMMPS_NS {
  FixWPMCAwpmd::FixWPMCAwpmd(LAMMPS_NS::LAMMPS *lmp, int narg, char **args) :
      Fix(lmp, narg, args) {
    //if (!atom->wavepacket_flag)
    //error->all(FLERR, "Fix wpmc/awpmd requires atom style wavepacket");

    vector_flag = 1;
    size_vector = sizeof(output) / sizeof(double);

    global_freq = 1;
    extvector = 0;
    time_depend = 1;

    nevery = 1;

    v_id = input->variable->find(args[3]);
    if (v_id == -1)
      error->all(FLERR, "Fix wpmc/awpmd requires a valid variable style");

    target_temperature = force->numeric(FLERR, args[4]) * force->boltz;
    output.like_vars.accepted_count = output.like_vars.rejected_count = 0.0;

    init_mc_steppers(narg, args);
  }

  void FixWPMCAwpmd::init() {
  }

  double FixWPMCAwpmd::memory_usage() {
    return sizeof(energy_old) + sizeof(output);
  }

  FixWPMCAwpmd::~FixWPMCAwpmd() = default;

  double FixWPMCAwpmd::compute_vector(int i) {
    return output.like_vector[i];
  }

  void FixWPMCAwpmd::final_integrate() {
    throw std::runtime_error("It's a wrong step!");
    auto energy_new = input->variable->compute_equal(v_id);
    printf("OUTRY:%lf\n", energy_new);
    this->output.like_vars.accept_flag = steppers.current().engine.test(energy_new - energy_old, 1.);

    if (output.like_vars.accept_flag == 1) {
      energy_old = energy_new;
    } else {
      steppers.current().restore((size_t) atom->nlocal);
    }

    output.like_vars.step_energy = energy_new;
    output.like_vars.accepted_energy = energy_old;

    output.like_vars.accepted_count += output.like_vars.accept_flag;
    output.like_vars.rejected_count += (output.like_vars.accept_flag == 0);

    steppers.current().adjust();
    steppers.next();
  }

  void FixWPMCAwpmd::pre_force(int i) {
    TimeMetrics tm;
    auto energy_new = input->variable->compute_equal(v_id);
    tm.make_tick("1");
    steppers.current().save((size_t) atom->nlocal);
    tm.make_tick("2");
    steppers.current().make((size_t) atom->nlocal);
    tm.make_tick("3");
    if (comm->nprocs > 1)
      update_ghosts();
  }

  void FixWPMCAwpmd::init_mc_steppers(int argc, char **argv) {
    unsigned const ARG_SHIFT = 5u;

    auto electron_filter = [this](int index) { return atom->spin[index] != 0; };
    auto ion_filter = [this](int index) { return atom->spin[index] == 0; };
    unsigned long engine_seed =  std::random_device{}();
    if (comm->nprocs > 1)
      MPI_Bcast(&engine_seed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    for (auto i = ARG_SHIFT; i < argc; ++i) {
      auto random_seed = std::abs((int) std::random_device{}());
      if (!std::strcmp(argv[i], "ix")) {
        steppers.add(lmp, stepper_type::ion_r, random_seed, engine_seed).assign_subsystem(
            std::make_unique<MC3DVectorSystem>(atom->x, ion_filter));
      } else if (!std::strcmp(argv[i], "ex")) {
        steppers.add(lmp, stepper_type::electron_r, random_seed, engine_seed).assign_subsystem(
            std::make_unique<MC3DVectorSystem>(atom->x, electron_filter));
      } else if (!std::strcmp(argv[i], "ev")) {
        steppers.add(lmp, stepper_type::electron_p, random_seed, engine_seed).assign_subsystem(
            std::make_unique<MC3DVectorSystem>(atom->v, electron_filter));
      } else if (!std::strcmp(argv[i], "ew")) {
        steppers.add(lmp, stepper_type::electron_w, random_seed, engine_seed).assign_subsystem(
            std::make_unique<MCScalarSystem>(atom->eradius, electron_filter));
      } else if (!std::strcmp(argv[i], "ewp")) {
        steppers.add(lmp, stepper_type::electron_pw, random_seed, engine_seed).assign_subsystem(
            std::make_unique<MCScalarSystem>(atom->ervel, electron_filter));
      } else if (!std::strcmp(argv[i], "ec0")) {
        throw std::logic_error("Not impl yet.");
      } else if (!std::strcmp(argv[i], "ec1")) {
        throw std::logic_error("Not impl yet.");
      } else if (!std::strcmp(argv[i], "iv")) {
        steppers.add(lmp, stepper_type::ion_p, random_seed, engine_seed).assign_subsystem(
            std::make_unique<MC3DVectorSystem>(atom->v, ion_filter));
      } else {
        error->all(FLERR, ("Invalid stepper name"s + argv[i]).c_str());
      }
      steppers.get(i - ARG_SHIFT).max_shift = 0.1;
      steppers.get(i - ARG_SHIFT).engine.setT(target_temperature);
    }
  }

  void FixWPMCAwpmd::update_ghosts() {
    std::unordered_map<int, int> tag_to_index;
    auto ghost_map = std::async(std::launch::async, [&tag_to_index, this]() {
      for (auto i = atom->nlocal; i < atom->nghost; ++i)
        tag_to_index[atom->tag[i]] = i;
    });

    auto particle_data = std::move(steppers.current().pack(atom->nlocal, atom->tag));
    auto data_size = particle_data.size();

    int recv_size[comm->nprocs];
    MPI_Allgather(&data_size, 1, MPI_INT, recv_size, 1, MPI_INT, MPI_COMM_WORLD);

    int displace[comm->nprocs];
    displace[0] = 0;
    auto total_size = recv_size[comm->nprocs - 1];
    for (auto i = 1; i < comm->nprocs; ++i) {
      displace[i] = displace[i - 1] + recv_size[i - 1];
      total_size += recv_size[i - 1];
    }

    double recv_buf[total_size];
    MPI_Allgatherv(particle_data.data(), data_size, MPI_DOUBLE, recv_buf, recv_size, displace, MPI_DOUBLE, MPI_COMM_WORLD);

    ghost_map.wait();
    auto unpacked = steppers.current().unpack(recv_buf, total_size, tag_to_index);
  }

  void FixWPMCAwpmd::initial_integrate(int i) {
    if (is_first) {
      is_first = false;
      return;
    }
    auto energy_new = input->variable->compute_equal(v_id);
    this->output.like_vars.accept_flag = steppers.current().engine.test(energy_new - energy_old, 1.);

    if (output.like_vars.accept_flag == 1) {
      energy_old = energy_new;
    } else {
      steppers.current().restore((size_t) atom->nlocal);
    }

    output.like_vars.step_energy = energy_new;
    output.like_vars.accepted_energy = energy_old;

    output.like_vars.accepted_count += output.like_vars.accept_flag;
    output.like_vars.rejected_count += (output.like_vars.accept_flag == 0);
    output.like_vars.stepper_id = steppers.current_stepped_id();

    steppers.current().adjust();
    steppers.next();
  }
}
