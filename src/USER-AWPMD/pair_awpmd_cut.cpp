//
// Created by yalavrinenko on 16.03.2020.
//

#include "pair_awpmd_cut.h"
#include <wpmd_split.h>
#include "pair_wpmd_cut.h"
#include "atom.h"
#include "update.h"
#include "min.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include <cstring>

LAMMPS_NS::WavepacketPairCommon::awpmd_energies LAMMPS_NS::PairAWPMD::compute_energy_force() {
  awpmd_ions ions;
  awpmd_electrons electrons;
  std::vector<Vector_3> fi;

  std::tie(ions, electrons) = this->make_packets();
  this->init_wpmd(ions, electrons);

  if (wpmd->ni)
    fi.resize(static_cast<unsigned long>(wpmd->ni));

  wpmd->interaction(0x1u | 0x4u | 0x10u, fi.data());
  wpmd->forces2phys();

  auto coul_energy = wpmd->get_energy() - electron_ke_;

  double **f = atom->f;
  //tally ion force
  for (auto const &ion : ions) {
    auto &i_lmp = ion.lmp_index;
    auto &i_wpmd = ion.wpmd_index;
    for (auto k : {0, 1, 2})
      f[i_lmp][k] = fi[i_wpmd][k];
  }

  //tally electron force
  for (auto const &electron : electrons) {
    for (auto const &packets : electron.second) {
      auto i_lmp = packets.lmp_index;
      auto i_wpmd = packets.wpmd_index;

      int s = atom->spin[i_lmp] > 0 ? 0 : 1;
      Vector_3 fv;
      Vector_3 vforce;

      double erforce;
      double ervelforce;

      Vector_2 csforce;
      wpmd->get_wp_force(s, i_wpmd, &fv, &vforce,
                         &erforce,
                         &ervelforce, &csforce, 0);

      for (auto k : {0, 1, 2})
        f[i_lmp][k] = fv[k];
      atom->erforce[i_lmp] = erforce;
      atom->ervelforce[i_lmp] = ervelforce;
    }
  }
  awpmd_energies output;
  output.ee = coul_energy; //ee -energy. Coul only
  output.ei = 0.0; //ei - energy. Coul only
  output.ii = 0.0; //ii - ii-energy. Coul only
  output.ke = 0.0; //ps^2/(2.0 * me)
  output.ee_w = 0.0; //1/s^2
  return output;
}

LAMMPS_NS::PairAWPMD::awpmd_packets LAMMPS_NS::PairAWPMD::make_packets() const {
  int *spin = atom->spin;
  int *etag = atom->tag;

  awpmd_ions ions;
  awpmd_electrons electrons{};

  auto insert_particle = [&ions, &electrons, spin, etag, this](unsigned index) {
    if (spin[index] == 0) {
      ions.emplace_back(index, 0, etag[index]);
    } else if (spin[index] == 1 || spin[index] == -1) {
      electrons[etag[index]].emplace_back(index, 0, etag[index]);
    } else {
      error->all(FLERR, fmt("Invalid spin value (%d) for particle %d !", spin[index], index));
    }
  };

  for (int i = 0; i < atom->nlocal + atom->nghost; ++i)
    insert_particle(static_cast<unsigned int>(i));

  return LAMMPS_NS::PairAWPMD::awpmd_packets{std::move(ions), std::move(electrons)};
}

void LAMMPS_NS::PairAWPMD::init_wpmd(awpmd_ions &ions, awpmd_electrons &electrons) {
  int newton_pair = force->newton_pair;

  wpmd->newton_pair = newton_pair;
  wpmd->reset();

  double **x = atom->x;
  double *q = atom->q;
  int *spin = atom->spin;
  int *type = atom->type;
  double **v = atom->v;

  int nlocal = atom->nlocal;

  std::sort(ions.begin(), ions.end());
  for (auto &ion_index : ions) {
    auto &insert_index = ion_index.lmp_index;
    ion_index.wpmd_index = (unsigned) wpmd->add_ion(q[insert_index], Vector_3(x[insert_index][0], x[insert_index][1],
                                                                              x[insert_index][2]),
                                                    (insert_index < nlocal ? atom->tag[insert_index]
                                                                           : -atom->tag[insert_index]));
  }

  electron_ke_ = 0.0;
  for (auto &electron : electrons) {
    std::sort(electron.second.begin(), electron.second.end());
    auto &main_packet_index = electron.second.begin()->lmp_index;
    int s = spin[main_packet_index] > 0 ? 0 : 1;
    wpmd->add_electron(s);
    for (auto &e_split_index : electron.second) {
      auto &insert_index = e_split_index.lmp_index;
      if (spin[insert_index] != spin[main_packet_index])
        error->all(FLERR,
                   fmt("WP splits for one electron should have the same spin (at particles %d, %d)!", insert_index,
                       main_packet_index));

      double m = atom->mass ? atom->mass[type[insert_index]] : force->e_mass;
      Vector_3 xx = Vector_3(x[insert_index][0], x[insert_index][1], x[insert_index][2]);
      Vector_3 rv = Vector_3(v[insert_index][0], v[insert_index][1], v[insert_index][2]);

      double pv = m * atom->ervel[insert_index];
      Vector_2 cc = Vector_2(atom->cs[2 * insert_index], atom->cs[2 * insert_index + 1]);

      e_split_index.wpmd_index = (unsigned) wpmd->add_split(xx, rv, atom->eradius[insert_index], pv, cc, m,
                                                            atom->q[insert_index],
                                                            (insert_index < nlocal ? atom->tag[insert_index]
                                                                                   : -atom->tag[insert_index]));
      electron_ke_ += (insert_index < nlocal) ? wpmd->wp[s][e_split_index.wpmd_index].get_p().norm2() * (wpmd->h2_me / 2.0) : 0.0;
    }
  }
}

void LAMMPS_NS::PairAWPMD::settings(int narg, char **arg) {
  WavepacketPairCommon::settings(narg, arg);
  if (narg < 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = force->numeric(FLERR, arg[0]);

  ermscale = 1.;
  width_pbc = 0.;

  wpmd->calc_ei = wpmd->calc_ii = wpmd->calc_ee = true;

  for (int i = 1; i < narg; i++) {
    // reading commands
    if (!strcmp(arg[i], "hartree"))
      wpmd->approx = AWPMD::HARTREE;
    else if (!strcmp(arg[i], "dproduct"))
      wpmd->approx = AWPMD::DPRODUCT;
    else if (!strcmp(arg[i], "uhf"))
      wpmd->approx = AWPMD::UHF;
    else if (!strcmp(arg[i], "free"))
      wpmd->constraint = AWPMD::NONE;
    else if (!strcmp(arg[i], "fix")) {
      wpmd->constraint = AWPMD::FIX;
      i++;
      if (i >= narg)
        error->all(FLERR, "Setting 'fix' should be followed by a number in awpmd/cut");
      wpmd->w0 = force->numeric(FLERR, arg[i]);
    } else if (!strcmp(arg[i], "harm")) {
      wpmd->constraint = AWPMD::HARM;
      i++;
      if (i >= narg)
        error->all(FLERR, "Setting 'harm' should be followed by a number in awpmd/cut");
      wpmd->w0 = force->numeric(FLERR, arg[i]);
      wpmd->set_harm_constr(wpmd->w0);
    } else if (!strcmp(arg[i], "pbc")) {
      i++;
      if (i >= narg)
        error->all(FLERR, "Setting 'pbc' should be followed by a number in awpmd/cut");
      width_pbc = force->numeric(FLERR, arg[i]);
    } else if (!strcmp(arg[i], "relax"))
      wpmd->constraint = AWPMD::RELAX;
    else if (!strcmp(arg[i], "ermscale")) {
      i++;
      if (i >= narg)
        error->all(FLERR, "Setting 'ermscale' should be followed by a number in awpmd/cut");
      ermscale = force->numeric(FLERR, arg[i]);
    }
  }
}
