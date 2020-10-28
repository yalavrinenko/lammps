//
// Created by yalavrinenko on 16.03.2020.
//
#include "WavepacketPairCommon.h"
#include <cstdio>
#include <cstring>
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

void LAMMPS_NS::WavepacketPairCommon::compute(int eflag, int vflag) {
  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0; //??

  energy_components_ = this->compute_energy_force();
  interaction_energy_ += energy_components_.sum();

  if (eflag_global) {
    eng_coul += interaction_energy_;

    pvector[0] = energy_components_.ii;
    pvector[1] = energy_components_.ee;
    pvector[2] = energy_components_.ei;
    pvector[3] = energy_components_.ke + energy_components_.ee_w;
  }

  interaction_energy_ = 0.0;
  if (vflag_fdotr) {
    virial_fdotr_compute();
    if (flexible_pressure_flag)
      virial_eradius_compute();
  }
}

LAMMPS_NS::WavepacketPairCommon::WavepacketPairCommon(LAMMPS_NS::LAMMPS *lmp) : Pair(lmp){
  single_enable = 0;
  nextra = 4;
  pvector = new double[nextra];

  wpmd = new AWPMD_split();
}

LAMMPS_NS::WavepacketPairCommon::~WavepacketPairCommon() {
  delete[] pvector;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
  delete wpmd;
}

void LAMMPS_NS::WavepacketPairCommon::settings(int narg, char **arg) {
  if (narg < 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = force->numeric(FLERR, arg[0]);

  wpmd->calc_ei = wpmd->calc_ii = wpmd->calc_ee = true;

  for (int i = 1; i < narg; i++) {
    if (!strcmp(arg[i], "flex_press"))
      flexible_pressure_flag = 1;
    else if (!strcmp(arg[i], "disable_ii")) {
      wpmd->calc_ii = false;
      error->warning(FLERR, "Ion-ion interaction disabled.");
    } else if (!strcmp(arg[i], "disable_ei")) {
      wpmd->calc_ei = false;
      error->warning(FLERR, "Electron-ion interaction disabled.");
    } else if (!strcmp(arg[i], "disable_ee")) {
      wpmd->calc_ee = false;
      error->warning(FLERR, "Electron-electron interaction disabled.");
    }
  }
}

void LAMMPS_NS::WavepacketPairCommon::coeff(int narg, char **arg) {
  if (narg < 2 || narg > 3) error->all(FLERR, "Incorrect args for pair coefficients");

  double delx = domain->boxhi[0] - domain->boxlo[0];
  double dely = domain->boxhi[1] - domain->boxlo[1];
  double delz = domain->boxhi[2] - domain->boxlo[2];
  auto half_box_length = 0.5 * MIN(delx, MIN(dely, delz));

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

void LAMMPS_NS::WavepacketPairCommon::init_style() {
  if (!atom->q_flag || !atom->spin_flag ||
      !atom->eradius_flag || !atom->erforce_flag)  // TO DO: adjust this to match approximation used
    error->all(FLERR, "Pair awpmd/cut requires atom attributes "
                      "q, spin, eradius, erforce");

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->newton = 0;

  if (force->e_mass == 0. || force->hhmrr2e == 0. || force->mvh2r == 0.)
    error->all(FLERR,
               "Pair style awpmd requires e_mass and conversions hhmrr2e, mvh2r to be properly set for unit system");

  wpmd->me = force->e_mass;
  wpmd->h2_me = force->hhmrr2e / force->e_mass;
  wpmd->one_h = force->mvh2r;
  wpmd->coul_pref = force->qqrd2e;
  wpmd->mvv2e = force->mvv2e;
}

double LAMMPS_NS::WavepacketPairCommon::init_one(int i, int j) {
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);

  return cut[i][j];
}

void LAMMPS_NS::WavepacketPairCommon::write_restart(FILE *fp) {
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) fwrite(&cut[i][j], sizeof(double), 1, fp);
    }
}

void LAMMPS_NS::WavepacketPairCommon::read_restart(FILE *fp) {
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

void LAMMPS_NS::WavepacketPairCommon::write_restart_settings(FILE *fp) {
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

void LAMMPS_NS::WavepacketPairCommon::read_restart_settings(FILE *fp) {
  if (comm->me == 0) {
    fread(&cut_global, sizeof(double), 1, fp);
    fread(&offset_flag, sizeof(int), 1, fp);
    fread(&mix_flag, sizeof(int), 1, fp);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

double LAMMPS_NS::WavepacketPairCommon::memory_usage() {
  double bytes = static_cast<double>(maxeatom) * sizeof(double);
  bytes += static_cast<double>(maxvatom) * 6 * sizeof(double);
  return bytes;
}

void LAMMPS_NS::WavepacketPairCommon::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
}

void LAMMPS_NS::WavepacketPairCommon::virial_eradius_compute() {
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
std::vector<WavePacket> const &
LAMMPS_NS::WavepacketPairCommon::electrons_packets() const {
  return packets;
}
