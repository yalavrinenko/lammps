//
// Created by yalavrinenko on 28.12.18.
//

#include <style_fix.h>
#include "fix_mc_wpmd.h"
#include "pair_awpmd_cut.h"
#include "wpmd_split.h"
#include "group.h"
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
#include <memory>
#include <compute.h>
#include <input.h>

namespace LAMMPS_NS {

  namespace {
    template<typename T, typename ... TArgs>
    std::unique_ptr<T> make_unique(TArgs &&... args) {
      return std::unique_ptr<T>{new T{std::forward<TArgs>(args)...}};
    }
  }

  FixMCAwpmd::FixMCAwpmd(LAMMPS_NS::LAMMPS *lmp, int narg, char **args) :
      Fix(lmp, narg, args) {
    //if (!atom->wavepacket_flag)
    //error->all(FLERR, "Fix wpmc/awpmd requires atom style wavepacket");

    vector_flag = 1;
    size_vector = sizeof(output) / sizeof(double);

    global_freq = 1;
    extvector = 0;
    time_depend = 1;

    nevery = 1;
    
    
    target_temperature = utils::numeric(FLERR, args[3], true, lmp) * force->boltz;
    output.like_vars.accepted_count = output.like_vars.rejected_count = 0.0;

    init_mc_steppers(narg, args);
    
    awpmd = dynamic_cast<PairAWPMD *>(force->pair);
    if (awpmd) {
      ke= modify->get_compute_by_id("thermo_awpmd_ke");
      if (!ke) {
        modify->add_compute("thermo_awpmd_ke all awpmd_ke");
        ke = modify->get_compute_by_id("thermo_awpmd_ke");
        if (!ke)
          error->all(FLERR, "Can't create compute for kinetic energy in awpmd");
      }
      use_awpmd_ke = true; // use awpmd-calculated kinetic energies instead of \sum p^2/(2m)
      awpmd->need_force = false; // don't need force for MC
    }
    if (use_norm) { 
      if (awpmd)
        awpmd->need_norm = true;
      norm = modify->get_compute_by_id("thermo_norm");
      if (!norm) {
        modify->add_compute("thermo_norm all normmatr");
        norm = modify->get_compute_by_id("thermo_norm");
        if (!norm)
          error->all(FLERR, "Can't create compute for normmatr, use pair/awpmd for norm");
      }
    }
    if (use_awpmd_ke) { // need a compute for ion kinetic energy
      temp = modify->get_compute_by_id("thermo_temp_ion");
      if (!temp) { // trying to define compute
        if(group->find("ion")<0)
          error->message(FLERR, "Can't create compute for ion kinetic energy, assumed constant. Add 'ion' group with all ions to fix this.");
        else {
          modify->add_compute("thermo_temp_ion ion temp");
          temp = modify->get_compute_by_id("thermo_temp_ion");
        } 
      }
    }
    else
      temp = modify->compute[modify->find_compute("thermo_temp")];

    pe = modify->compute[modify->find_compute("thermo_pe")];
  }

  void FixMCAwpmd::init() {
  }

  double FixMCAwpmd::memory_usage() {
    return sizeof(energy_old) + sizeof(output);
  }

  FixMCAwpmd::~FixMCAwpmd() = default;

  double FixMCAwpmd::compute_vector(int i) {
    return output.like_vector[i];
  }

  void FixMCAwpmd::final_integrate() {
    double energy_new;
    double nmlog = 0;
    if (use_norm)
      nmlog = norm->compute_scalar();
  
    double T = target_temperature;
    double classic_dof = 0.; // number of ion degrees of freedom, to be used when we don't have a compute for ion ke
    if (temp) { // ion part for awpmd_ke
      T = temp->compute_scalar();
      classic_dof = temp->dof;
    }
    else {
      if (awpmd)
        classic_dof = (double)(3 * awpmd->awpmd()->ni - 3);
      if (classic_dof < 0)
        classic_dof = 0.;
    }
    double quantum_ke = 0.;
    if (use_awpmd_ke)
      quantum_ke = ke->compute_scalar();
    
    
    if (awpmd)
      energy_new = T * 0.5 * classic_dof * force->boltz - target_temperature * 0.5 * nmlog +
       awpmd->awpmd()->get_energy();
    else
      energy_new = T * 0.5 * classic_dof * force->boltz - target_temperature * 0.5 * nmlog   +
      quantum_ke + pe->compute_scalar(); //input->variable->compute_equal(v_id);

    

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

  void FixMCAwpmd::pre_force(int i) {
    steppers.current().save((size_t) atom->nlocal);
    steppers.current().make((size_t) atom->nlocal);
    if (comm->nprocs > 1)
      update_ghosts();
  }

  void FixMCAwpmd::init_mc_steppers(int argc, char **argv) {
    unsigned const ARG_SHIFT = 4u;

    auto electron_filter = [this](int index) { return this->atom->mask[index] && atom->spin[index] != 0; };
    auto ion_filter = [this](int index) { return this->atom->mask[index] && atom->spin[index] == 0; };
    unsigned long engine_seed = std::random_device{}();
    if (comm->nprocs > 1)
      MPI_Bcast(&engine_seed, 1, MPI_UNSIGNED_LONG, 0, world);

    for (int i = ARG_SHIFT; i < argc; ++i) {
      auto random_seed = std::abs((int) std::random_device{}() / 2) + 1;
      if (!std::strcmp(argv[i], "ix")) {
        steppers.add(lmp, stepper_type::ion_r, random_seed, engine_seed).assign_subsystem(
            make_unique<MCVectorSystem<3>>(atom->x, ion_filter));
      } else if (!std::strcmp(argv[i], "ex")) {
        steppers.add(lmp, stepper_type::electron_r, random_seed, engine_seed).assign_subsystem(
            make_unique<MCVectorSystem<3>>(atom->x, electron_filter));
      } else if (!std::strcmp(argv[i], "ev")) {
        steppers.add(lmp, stepper_type::electron_p, random_seed, engine_seed).assign_subsystem(
            make_unique<MCVectorSystem<3>>(atom->v, electron_filter));
      } else if (!std::strcmp(argv[i], "ew")) {
        steppers.add(lmp, stepper_type::electron_w, random_seed, engine_seed).assign_subsystem(
            make_unique<MCScalarSystem>(atom->eradius, electron_filter));
      } else if (!std::strcmp(argv[i], "ewp")) {
        steppers.add(lmp, stepper_type::electron_pw, random_seed, engine_seed).assign_subsystem(
            make_unique<MCScalarSystem>(atom->ervel, electron_filter));
      } else if (!std::strcmp(argv[i], "ec_re")) {
        auto c_re_proj = [](double **src, unsigned i, unsigned j) -> double & {
          return src[i][0];
        };
        steppers.add(lmp, stepper_type::electron_c, random_seed, engine_seed).assign_subsystem(
            make_unique<MCVectorSystem<1, decltype(c_re_proj)>>(atom->cs, electron_filter, c_re_proj));
      } else if (!std::strcmp(argv[i], "ec_im")) {
        auto c_im_proj = [](double **src, unsigned i, unsigned j) -> double & {
          return src[i][1];
        };
        steppers.add(lmp, stepper_type::electron_c, random_seed, engine_seed).assign_subsystem(
            make_unique<MCVectorSystem<1, decltype(c_im_proj)>>(atom->cs, electron_filter, c_im_proj));
      } else if (!std::strcmp(argv[i], "iv")) {
        steppers.add(lmp, stepper_type::ion_p, random_seed, engine_seed).assign_subsystem(
            make_unique<MCVectorSystem<3>>(atom->v, ion_filter));
      }
      else if (!std::strcmp(argv[i], "norm")) {
        use_norm = true;
        continue;
      } else {
        error->all(FLERR, (std::string{"Invalid stepper name "} + argv[i]).c_str());
      }
      steppers.get(i - ARG_SHIFT).max_shift = 0.1;
      steppers.get(i - ARG_SHIFT).engine.setT(target_temperature);
    }
  }

  void FixMCAwpmd::update_ghosts() {
    std::unordered_map<int, int> tag_to_index;

    for (auto i = atom->nlocal; i < atom->nlocal + atom->nghost; ++i)
      tag_to_index[atom->tag[i]] = i;

    auto particle_data = steppers.current().pack(atom->nlocal, atom->tag);
    auto data_size = static_cast<int>(particle_data.size());

    std::vector<int> recv_size(comm->nprocs);
    MPI_Allgather(&data_size, 1, MPI_INT, &recv_size[0], 1, MPI_INT, world);

    std::vector<int> displace(comm->nprocs);
    displace[0] = 0;
    auto total_size = recv_size[comm->nprocs - 1];
    for (auto i = 1; i < comm->nprocs; ++i) {
      displace[i] = displace[i - 1] + recv_size[i - 1];
      total_size += recv_size[i - 1];
    }

    vector<double> recv_buf(total_size);
    MPI_Allgatherv(particle_data.data(), data_size, MPI_DOUBLE, &recv_buf[0], &recv_size[0], &displace[0], MPI_DOUBLE,
                   world);

    //ghost_map.wait();
    steppers.current().unpack(&recv_buf[0], total_size, tag_to_index);
  }

  void FixMCAwpmd::initial_integrate(int i) {
  }
}
