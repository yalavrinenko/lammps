/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_awpmd_cut.h"
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

#include <wpmd_split.h>


#include <chrono>
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairAWPMDCut::PairAWPMDCut(LAMMPS *lmp) : Pair(lmp) {
  single_enable = 0;

  nmax = 0;
  min_var = nullptr;
  min_varforce = nullptr;
  nextra = 5;
  pvector = new double[nextra];

  ermscale = 1.;
  width_pbc = 0.;
  wpmd = new AWPMD_split();

  half_box_length = 0;
}

/* ---------------------------------------------------------------------- */

PairAWPMDCut::~PairAWPMDCut() {
  delete[] pvector;
  memory->destroy(min_var);
  memory->destroy(min_varforce);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }

  delete wpmd;
}


struct cmp_x {
  double **xx;
  double tol;

  cmp_x(double **xx_ = nullptr, double tol_ = 1e-12) : xx(xx_), tol(tol_) {}

  bool operator()(const pair<int, int> &left, const pair<int, int> &right) const {
    if (left.first == right.first) {
      double d = xx[left.second][0] - xx[right.second][0];
      if (d < -tol)
        return true;
      else if (d > tol)
        return false;
      d = xx[left.second][1] - xx[right.second][1];
      if (d < -tol)
        return true;
      else if (d > tol)
        return false;
      d = xx[left.second][2] - xx[right.second][2];
      return d < -tol;
    } else
      return left.first < right.first;
  }
};

/* ---------------------------------------------------------------------- */

PairAWPMDCut::awpmd_packets PairAWPMDCut::make_packets() const {
  int *spin = atom->spin;
  //int *etag = atom->etag;
  int *etag = atom->tag;

  awpmd_ions ions;
  awpmd_electrons electrons{};

  auto insert_particle = [&ions, &electrons, spin, etag, this](unsigned index) {
    if (spin[index] == 0) {
      ions.emplace_back(awpmd_pair_index{index, 0});
    } else if (spin[index] == 1 || spin[index] == -1) {
      electrons[etag[index]].emplace_back(awpmd_pair_index{index, 0});
    } else {
      error->all(FLERR, fmt("Invalid spin value (%d) for particle %d !", spin[index], index));
    }
  };

  for (int i = 0; i < atom->nlocal + atom->nghost; ++i)
    insert_particle(static_cast<unsigned int>(i));

  return LAMMPS_NS::PairAWPMDCut::awpmd_packets{std::move(ions), std::move(electrons)};
}


void PairAWPMDCut::init_wpmd(awpmd_ions &ions, awpmd_electrons &electrons) {
  int newton_pair = force->newton_pair;

  if (width_pbc < 0)
    wpmd->Lextra = 2 * half_box_length;
  else
    wpmd->Lextra = width_pbc;

  wpmd->newton_pair = newton_pair;

  // prepare the solver object
  wpmd->reset();
  wpmd->set_pbc(nullptr); // not required for LAMMPS

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *spin = atom->spin;
  int *type = atom->type;
  double **v = atom->v;

  int nlocal = atom->nlocal;

  for (auto &ion_index : ions) {
    auto &insert_index = ion_index.lmp_index;
    ion_index.wpmd_index = (unsigned) wpmd->add_ion(q[insert_index], Vector_3(x[insert_index][0], x[insert_index][1],
                                                                              x[insert_index][2]),
                                                    (insert_index < nlocal ? atom->tag[insert_index]
                                                                           : -atom->tag[insert_index]));
  }

  electron_ke_ = 0.0;
  for (auto &electron : electrons) {
    auto &main_packet_index = electron.second.front().lmp_index;
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
      electron_ke_ += (insert_index < nlocal) ? m * (rv * rv) / 2.0 : 0.0;

      double pv = ermscale * m * atom->ervel[insert_index];
      Vector_2 cc = Vector_2(atom->cs[2 * insert_index], atom->cs[2 * insert_index + 1]);

      e_split_index.wpmd_index = (unsigned) wpmd->add_split(xx, rv, atom->eradius[insert_index], pv, cc, m,
                                                            atom->q[insert_index],
                                                            (insert_index < nlocal ? atom->tag[insert_index]
                                                                                   : -atom->tag[insert_index]));
    }
  }

}

void PairAWPMDCut::compute(int eflag, int vflag) {
  // pvector = [KE, Pauli, ecoul, radial_restraint]
  for (int i = 0; i < 5; i++) pvector[i] = 0.0;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0; //??

  awpmd_ions ions;
  awpmd_electrons electrons;

  std::tie(ions, electrons) = this->make_packets();
  this->init_wpmd(ions, electrons);

  std::vector<Vector_3> fi;
  if (wpmd->ni)
    fi.resize(static_cast<unsigned long>(wpmd->ni));

  auto begin = std::chrono::high_resolution_clock::now();
  wpmd->interaction(0x1 | 0x4 | 0x10, fi.data());
  auto eng = std::chrono::high_resolution_clock::now();
  std::cerr << std::chrono::duration_cast<std::chrono::milliseconds>(eng - begin).count() << std::endl;

  auto full_coul_energy = wpmd->get_energy() - electron_ke_ * force->mvv2e;

  double **f = atom->f;

  for (auto const &ion : ions) {
    auto &i_lmp = ion.lmp_index;
    auto &i_wpmd = ion.wpmd_index;
    f[i_lmp][0] = fi[i_wpmd][0];
    f[i_lmp][0] = fi[i_wpmd][1];
    f[i_lmp][0] = fi[i_wpmd][2];
  }

  for (auto const &electron : electrons) {
    for (auto const &packets : electron.second) {
      auto &i_lmp = packets.lmp_index;
      auto &i_wpmd = packets.wpmd_index;

      int s = atom->spin[i_lmp] > 0 ? 0 : 1;
      wpmd->get_wp_force(s, i_wpmd, (Vector_3 *) f[i_lmp], (Vector_3 *) (atom->vforce + 3 * i_lmp),
                         atom->erforce + i_lmp,
                         atom->ervelforce + i_lmp, (Vector_2 *) (atom->csforce + 2 * i_lmp));
    }
  }

  // update LAMMPS energy
  if (eflag_either) {
    if (eflag_global) {
      eng_coul += full_coul_energy;

      // pvector = [KE, Pauli, ecoul, radial_restraint]
      pvector[0] = wpmd->Ee[0] + wpmd->Ee[1];
      pvector[2] = wpmd->Eii + wpmd->Eei[0] + wpmd->Eei[1] + wpmd->Eee;
      pvector[1] = pvector[0] + pvector[2] - wpmd->Edk - wpmd->Edc - wpmd->Eii;  // All except diagonal terms
      pvector[3] = wpmd->Ew;
      pvector[4] = wpmd->Ebord + wpmd->Ebord_ion;
    }

    if (eflag_atom) {
      // transfer per-atom energies here
      for (auto const &ion : ions) {
        auto &i_lmp = ion.lmp_index;
        auto &i_wpmd = ion.wpmd_index;
        eatom[i_lmp] = wpmd->Eiep[i_wpmd] + wpmd->Eiip[i_wpmd];
      }

      for (auto const &electron : electrons) {
        for (auto const &packets : electron.second) {
          auto &i_lmp = packets.lmp_index;
          auto &i_wpmd = packets.wpmd_index;

          int s = atom->spin[i_lmp] > 0 ? 0 : 1;
          eatom[i_lmp] = wpmd->Eep[s][i_wpmd] + wpmd->Eeip[s][i_wpmd] + wpmd->Eeep[s][i_wpmd] + wpmd->Ewp[s][i_wpmd];
        }
      }
    }
  }

  if (vflag_fdotr) {
    virial_fdotr_compute();
    if (flexible_pressure_flag)
      virial_eradius_compute();
  }
}

/* ----------------------------------------------------------------------
   electron width-specific contribution to global virial
------------------------------------------------------------------------- */

void PairAWPMDCut::virial_eradius_compute() {
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  double e_virial;
  int *spin = atom->spin;

  // sum over force on all particles including ghosts

  if (neighbor->includegroup == 0) {
    int nall = atom->nlocal + atom->nghost;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i] * eradius[i] / 3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }

    // neighbor includegroup flag is set
    // sum over force on initial nfirst particles and ghosts

  } else {
    int nall = atom->nfirst;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i] * eradius[i] / 3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }

    nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i] * eradius[i] / 3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }
  }
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairAWPMDCut::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
}

/* ---------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
// the format is: pair_style awpmd/cut [<global_cutoff|-1> [command1] [command2] ...]
// commands:
// [hartree|dproduct|uhf]  -- quantum approximation level (default is hartree)
// [free|pbc <length|-1>|fix <w0|-1>|relax|harm <w0>] -- width restriction (default is free)
// [ermscale <number>]  -- scaling factor between electron mass and effective width mass (used for equations of motion only) (default is 1)
// [flex_press]  -- set flexible pressure flag
// -1 for length means default setting (L/2 for cutoff and L for width PBC)

void PairAWPMDCut::settings(int narg, char **arg) {
  if (narg < 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = force->numeric(FLERR, arg[0]);

  ermscale = 1.;
  width_pbc = 0.;

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
    } else if (!strcmp(arg[i], "flex_press"))
      flexible_pressure_flag = 1;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
// pair settings are as usual
void PairAWPMDCut::coeff(int narg, char **arg) {
  if (narg < 2 || narg > 3) error->all(FLERR, "Incorrect args for pair coefficients");

  /*if(domain->xperiodic == 1 || domain->yperiodic == 1 ||
    domain->zperiodic == 1) {*/
  double delx = domain->boxhi[0] - domain->boxlo[0];
  double dely = domain->boxhi[1] - domain->boxlo[1];
  double delz = domain->boxhi[2] - domain->boxlo[2];
  half_box_length = 0.5 * MIN(delx, MIN(dely, delz));
  //}
  if (cut_global < 0)
    cut_global = half_box_length;

  if (!allocated)
    allocate();
  else {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  double cut_one = cut_global;
  if (narg == 3) cut_one = force->numeric(FLERR, arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAWPMDCut::init_style() {
  // error and warning checks

  if (!atom->q_flag || !atom->spin_flag ||
      !atom->eradius_flag || !atom->erforce_flag)  // TO DO: adjust this to match approximation used
    error->all(FLERR, "Pair awpmd/cut requires atom attributes "
                      "q, spin, eradius, erforce");

  /*
  if(vflag_atom){ // can't compute virial per atom
    //warning->
    error->all(FLERR,"Pair style awpmd can't compute per atom virials");
  }*/

  // add hook to minimizer for eradius and erforce

  if (update->whichflag == 2)
    int ignore = update->minimize->request(this, 1, 0.01);

  // make sure to use the appropriate timestep when using real units

  /*if (update->whichflag == 1) {
    if (force->qqr2e == 332.06371 && update->dt == 1.0)
      error->all(FLERR,"You must lower the default real units timestep for pEFF ");
  }*/

  // need a half neigh list and optionally a granular history neigh list

  //int irequest = neighbor->request(this,instance_me);

  //if (atom->tag_enable == 0)
  //  error->all(FLERR,"Pair style reax requires atom IDs");

  //if (force->newton_pair == 0)
  //error->all(FLERR,"Pair style awpmd requires newton pair on");

  //if (strcmp(update->unit_style,"real") != 0 && comm->me == 0)
  //error->warning(FLERR,"Not using real units with pair reax");

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->newton = 2;

  if (force->e_mass == 0. || force->hhmrr2e == 0. || force->mvh2r == 0.)
    error->all(FLERR,
               "Pair style awpmd requires e_mass and conversions hhmrr2e, mvh2r to be properly set for unit system");

  wpmd->me = force->e_mass;
  wpmd->h2_me = force->hhmrr2e / force->e_mass;
  wpmd->one_h = force->mvh2r;
  wpmd->coul_pref = force->qqrd2e;

  wpmd->calc_ii = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAWPMDCut::init_one(int i, int j) {
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairAWPMDCut::write_restart(FILE *fp) {
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) fwrite(&cut[i][j], sizeof(double), 1, fp);
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairAWPMDCut::read_restart(FILE *fp) {
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j], sizeof(int), 1, fp);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) fread(&cut[i][j], sizeof(double), 1, fp);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairAWPMDCut::write_restart_settings(FILE *fp) {
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairAWPMDCut::read_restart_settings(FILE *fp) {
  if (comm->me == 0) {
    fread(&cut_global, sizeof(double), 1, fp);
    fread(&offset_flag, sizeof(int), 1, fp);
    fread(&mix_flag, sizeof(int), 1, fp);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   returns pointers to the log() of electron radius and corresponding force
   minimizer operates on log(radius) so radius never goes negative
   these arrays are stored locally by pair style
------------------------------------------------------------------------- */

void PairAWPMDCut::min_xf_pointers(int ignore, double **xextra, double **fextra) {
  // grow arrays if necessary
  // need to be atom->nmax in length
  int nvar = atom->nmax * (3 + 1 + 1 + 2);  // w(1), vel(3),  pw(1), cs(2)

  if (nvar > nmax) {
    memory->destroy(min_var);
    memory->destroy(min_varforce);
    nmax = nvar;
    memory->create(min_var, nmax, "pair:min_var");
    memory->create(min_varforce, nmax, "pair:min_varforce");
  }

  *xextra = min_var;
  *fextra = min_varforce;
}

/* ----------------------------------------------------------------------
   minimizer requests the log() of electron radius and corresponding force
   calculate and store in min_eradius and min_erforce
------------------------------------------------------------------------- */

void PairAWPMDCut::min_xf_get(int ignore) {
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  double **v = atom->v;
  double *vforce = atom->vforce;
  double *ervel = atom->ervel;
  double *ervelforce = atom->ervelforce;
  double *cs = atom->cs;
  double *csforce = atom->csforce;

  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (spin[i]) {
      min_var[7 * i] = log(eradius[i]);
      min_varforce[7 * i] = eradius[i] * erforce[i];
      for (int j = 0; j < 3; j++) {
        min_var[7 * i + 1 + 3 * j] = v[i][j];
        min_varforce[7 * i + 1 + 3 * j] = vforce[3 * i + j];
      }
      min_var[7 * i + 4] = ervel[i];
      min_varforce[7 * i + 4] = ervelforce[i];
      min_var[7 * i + 5] = cs[2 * i];
      min_varforce[7 * i + 5] = csforce[2 * i];
      min_var[7 * i + 6] = cs[2 * i + 1];
      min_varforce[7 * i + 6] = csforce[2 * i + 1];

    } else {
      for (int j = 0; j < 7; j++)
        min_var[7 * i + j] = min_varforce[7 * i + j] = 0.0;
    }
}

/* ----------------------------------------------------------------------
   propagate the minimizer values to the atom values
------------------------------------------------------------------------- */

void PairAWPMDCut::min_x_set(int ignore) {
  double *eradius = atom->eradius;
  double **v = atom->v;
  double *ervel = atom->ervel;
  double *cs = atom->cs;

  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (spin[i]) {
      eradius[i] = exp(min_var[7 * i]);
      for (int j = 0; j < 3; j++)
        v[i][j] = min_var[7 * i + 1 + 3 * j];
      ervel[i] = min_var[7 * i + 4];
      cs[2 * i] = min_var[7 * i + 5];
      cs[2 * i + 1] = min_var[7 * i + 6];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairAWPMDCut::memory_usage() {
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom * 6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}

double PairAWPMDCut::ghost_energy() {
  awpmd_ions ions;
  awpmd_electrons electrons;

  std::tie(ions, electrons) = this->make_packets();

  auto &local = atom->nlocal;

  ions.erase(std::remove_if(ions.begin(), ions.end(), [&local](auto const &v) { return v.lmp_index < local; }),
             ions.end());
  for (auto it = electrons.begin(); it != electrons.end();) {
    auto &v = it->second;
    v.erase(std::remove_if(v.begin(), v.end(), [&local](auto const &i) { return i.lmp_index < local; }), v.end());
    if (v.empty()) {
      it = electrons.erase(it);
    } else
      ++it;
  }

  this->init_wpmd(ions, electrons);

  std::vector<Vector_3> fi;
  if (wpmd->ni)
    fi.resize(static_cast<unsigned long>(wpmd->ni));

  wpmd->interaction(0x1 | 0x4 | 0x10, fi.data());

  return wpmd->get_energy();
}

AWPMD_split *PairAWPMDCut::awpmd() {
  return wpmd;
}