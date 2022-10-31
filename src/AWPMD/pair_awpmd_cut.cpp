// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ilya Valuev (JIHT, Moscow, Russia)
------------------------------------------------------------------------- */

#include "pair_awpmd_cut.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "min.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_wpmd_cut.h"
#include "update.h"
#include <cstring>
#include <wpmd_split.h>

LAMMPS_NS::WavepacketPairCommon::awpmd_energies LAMMPS_NS::PairAWPMD::compute_energy_force()
{
  awpmd_ions ions;
  awpmd_electrons electrons;
  std::vector<Vector_3> fi;

  wpmd->norm_mode = AWPMD_split::NORMALIZE;

  std::tie(ions, electrons) = this->make_packets();
  this->init_wpmd(ions, electrons);

  if (wpmd->ni) fi.resize(static_cast<unsigned long>(wpmd->ni));

  wpmd->interaction(0x1u | 0x4u | 0x10u | 0x20u, fi.data());
  wpmd->forces2phys();

  auto coul_energy = wpmd->get_energy() - electron_ke_;

  double **f = atom->f;
  //tally ion force
  for (auto const &ion : ions) {
    auto &i_lmp = ion.lmp_index;
    auto &i_wpmd = ion.wpmd_index;
    for (auto k : {0, 1, 2}) f[i_lmp][k] = fi[i_wpmd][k];
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
      wpmd->get_wp_force(s, i_wpmd, &fv, &vforce, &erforce, &ervelforce, &csforce, 0);

      for (auto k : {0, 1, 2}) f[i_lmp][k] = fv[k];
      atom->erforce[i_lmp] = erforce;
      atom->ervelforce[i_lmp] = ervelforce;
    }
  }
  awpmd_energies output;
  output.ee = coul_energy;    //ee -energy. Coul only
  output.ei = 0.0;            //ei - energy. Coul only
  output.ii = 0.0;            //ii - ii-energy. Coul only
  output.ke = 0.0;            //ps^2/(2.0 * me)
  output.ee_w = 0.0;          //1/s^2
  return output;
}

LAMMPS_NS::PairAWPMD::awpmd_packets LAMMPS_NS::PairAWPMD::make_packets() const
{
  int *spin = atom->spin;
  int *tag = atom->tag;
  int *etag = atom->etag;

  // check electrons with auto-assigned tags
  std::vector<int> newetag;
  if (wp_per_electron > 1) {
    std::vector<std::pair<int, int>> tagv;
    for (int i = 0; i < atom->nlocal + atom->nghost; ++i)
      if (spin[i] && !etag[i])    // searching for electrons with non-assigned tags
        tagv.push_back(std::make_pair(tag[i], i));

    if (tagv.size()) {
      std::sort(tagv.begin(), tagv.end());

      newetag.resize(atom->nlocal + atom->nghost);
      for (int i = 0; i < atom->nlocal + atom->nghost; ++i) newetag[i] = etag[i];
      for (
          size_t i = 0; i < tagv.size();
          i++)    // assuming wavepackets were created using groups with size of multiple of wp_per_electron
        newetag[tagv[i].second] = tag[tagv[i].second] / wp_per_electron;
      etag = &newetag[0];    // the new vector has all electron tags assigned
    }
  }
  awpmd_ions ions;
  awpmd_electrons electrons{};

  auto insert_particle = [&ions, &electrons, spin, etag, tag, this](unsigned index) {
    if (spin[index] == 0) {
      ions.emplace_back(index, 0, tag[index]);
    } else if (spin[index] == 1 || spin[index] == -1) {
      if (!etag[index])    // efficient solution for wp_per_electron =1, no preprocessing required
        etag[index] = tag[index];
      electrons[etag[index]].emplace_back(index, 0, etag[index]);
    } else {
      error->all(FLERR,
                 fmt::format("Invalid spin value ({}) for particle {} !", spin[index], index));
    }
  };

  for (int i = 0; i < atom->nlocal + atom->nghost; ++i)
    insert_particle(static_cast<unsigned int>(i));

  return LAMMPS_NS::PairAWPMD::awpmd_packets{std::move(ions), std::move(electrons)};
}

void LAMMPS_NS::PairAWPMD::init_wpmd(awpmd_ions &ions, awpmd_electrons &electrons)
{
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
                                                    (insert_index < (unsigned) nlocal ? atom->tag[insert_index]
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
        error->all(
            FLERR,
            fmt::format(
                "WP splits for one electron should have the same spin (at particles {}, {})!",
                insert_index, main_packet_index));

      double m = atom->mass ? atom->mass[type[insert_index]] : force->e_mass;
      Vector_3 xx = Vector_3(x[insert_index][0], x[insert_index][1], x[insert_index][2]);
      Vector_3 rv = Vector_3(v[insert_index][0], v[insert_index][1], v[insert_index][2]);

      double pv = m * atom->ervel[insert_index];
      Vector_2 cc = Vector_2(
          atom->cs[insert_index][0],
          atom->cs[insert_index][1]);    //atom->cs[2*insert_index], atom->cs[2*insert_index + 1]

      e_split_index.wpmd_index = (unsigned) wpmd->add_split(xx, rv, atom->eradius[insert_index], pv, cc, m,
                                                            atom->q[insert_index],
                                                            (insert_index < (unsigned) nlocal ? atom->tag[insert_index]
                                                                                   : -atom->tag[insert_index]));
      electron_ke_ += (insert_index < (unsigned) nlocal)
          ? wpmd->wp[s][e_split_index.wpmd_index].get_p().norm2() * (wpmd->h2_me / 2.0)
          : 0.0;
    }
  }
}

void LAMMPS_NS::PairAWPMD::settings(int narg, char **arg)
{
  WavepacketPairCommon::settings(narg, arg);
  if (narg < 1) error->all(FLERR, "Illegal pair_style command");

  if (!comm->ghost_velocity)
    error->all(
        FLERR,
        "pair_style requires ghost_velocity flag. Add comm_modify vel yes to your input script.");

  auto numeric = [this](char const *str) {
    return utils::numeric(FLERR, str, false, lmp);
  };

  cut_global = numeric(arg[0]);

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
      if (i >= narg) error->all(FLERR, "Setting 'fix' should be followed by a number in awpmd/cut");
      wpmd->w0 = numeric(arg[i]);
    } else if (!strcmp(arg[i], "harm")) {
      wpmd->constraint = AWPMD::HARM;
      i++;
      if (i >= narg)
        error->all(FLERR, "Setting 'harm' should be followed by a number in awpmd/cut");
      wpmd->w0 = numeric(arg[i]);
      wpmd->set_harm_constr(wpmd->w0);
    } else if (!strcmp(arg[i], "pbc")) {
      i++;
      if (i >= narg) error->all(FLERR, "Setting 'pbc' should be followed by a number in awpmd/cut");
      width_pbc = numeric(arg[i]);
    } else if (!strcmp(arg[i], "relax"))
      wpmd->constraint = AWPMD::RELAX;
    else if (!strcmp(arg[i], "ermscale")) {
      i++;
      if (i >= narg)
        error->all(FLERR, "Setting 'ermscale' should be followed by a number in awpmd/cut");
      ermscale = numeric(arg[i]);
    } else if (!strcmp(arg[i], "wp_per_electron")){
      ++i;
      wp_per_electron = utils::inumeric(FLERR, arg[i], false, lmp);
    }
  }
}
