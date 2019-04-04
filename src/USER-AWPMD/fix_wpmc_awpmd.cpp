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
using namespace std::string_literals;

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
    auto energy_new = input->variable->compute_equal(v_id);
    this->output.like_vars.accept_flag = steppers.current().engine.test(energy_new - energy_old, 1.);

    double collective_decision;
    MPI::COMM_WORLD.Allreduce(&(output.like_vars.accept_flag), &collective_decision, 1, MPI_DOUBLE, MPI_SUM);

    if (collective_decision > 0)
      this->output.like_vars.accept_flag = true;

    if (output.like_vars.accept_flag == 1){
      energy_old = energy_new;
    } else {
      steppers.current().restore((size_t)atom->nlocal);
    }

    output.like_vars.step_energy = energy_new;
    output.like_vars.accepted_energy = energy_old;

    output.like_vars.accepted_count += output.like_vars.accept_flag;
    output.like_vars.rejected_count += (output.like_vars.accept_flag == 0);

    steppers.current().adjust();
    steppers.next();
  }

  void FixWPMCAwpmd::pre_force(int i) {
    steppers.current().save((size_t)atom->nlocal);
    steppers.current().make((size_t)atom->nlocal);
    //update_ghosts();
  }

  void FixWPMCAwpmd::init_mc_steppers(int argc, char **argv) {
    unsigned const ARG_SHIFT = 5u;

    auto electron_filter = [this](int index) { return atom->spin[index] != 0; };
    auto ion_filter =[this](int index) { return atom->spin[index] == 0; };
    unsigned engine_seed = 1239553;

    for (auto i = ARG_SHIFT; i < argc; ++i) {
      auto random_seed = std::abs((int)std::random_device{}());
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
    auto* avec = atom->avec;
    double send_buf[(avec->size_border + avec->size_velocity + 2) * atom->nlocal];
    int send_shift = 0;
    for (auto i = 0; i < atom->nlocal; ++i)
      send_shift += avec->pack_exchange(i, &send_buf[send_shift]);

    int recv_size[comm->nprocs];
    MPI::COMM_WORLD.Allgather(&send_shift, 1, MPI_INT, recv_size, 1, MPI_INT);

    int displace[comm->nprocs];
    displace[0] = 0;
    for (auto i = 1; i < comm->nprocs; ++i)
      displace[i] = displace[i-1] + recv_size[i - 1];

    int total_recv = (avec->size_border + avec->size_velocity + 2) * (atom->nlocal + atom->nghost);
    double recv_buf[total_recv];
    MPI::COMM_WORLD.Allgatherv(send_buf, send_shift, MPI_DOUBLE, recv_buf, recv_size, displace, MPI_DOUBLE);

    auto tmp_nlocal = atom->nlocal;
    auto tmp_nghost = atom->nghost;

    std::unordered_map<int, int> tag_to_index;
    for (auto i = tmp_nlocal; i < tmp_nlocal + tmp_nghost; ++i)
      tag_to_index[atom->tag[i]] = i;

    atom->nlocal += atom->nghost;

    int unpuck_shift = 0;
    while (unpuck_shift < total_recv) {
      unpuck_shift += avec->unpack_exchange(&recv_buf[unpuck_shift]);
      if (tag_to_index.count(atom->tag[atom->nlocal - 1])){
        avec->copy(atom->nlocal - 1, tag_to_index[atom->tag[atom->nlocal - 1]], 0);
      }
      --atom->nlocal;
    }

    atom->nlocal = tmp_nlocal;
    atom->nghost = tmp_nghost;
  }

}
