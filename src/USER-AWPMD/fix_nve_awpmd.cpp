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
#include <cstring>
#include "fix_nve_awpmd.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

#include <wpmd_split.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEAwpmd::FixNVEAwpmd(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (!atom->wavepacket_flag)
    error->all(FLERR,"Fix nve/awpmd requires atom style wavepacket");
}

/* ----------------------------------------------------------------------
   allow for only per-type  mass
------------------------------------------------------------------------- */

void FixNVEAwpmd::initial_integrate(int vflag)
{
  FixNVE::initial_integrate(vflag);
  if (atom->mass || atom->rmass) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        double dtfm = dtf / ((atom->mass) ? atom->mass[atom->type[i]] : atom->rmass[atom->type[i]]);
        if (atom->spin[i] != 0) {
          atom->ervel[i] += -dtfm * atom->erforce[i];
          atom->eradius[i] += dtv * atom->ervel[i];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEAwpmd::final_integrate(){
  FixNVE::final_integrate();
  if (atom->mass || atom->rmass) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit) {
        double dtfm = dtf / ((atom->mass) ? atom->mass[atom->type[i]] : atom->rmass[atom->type[i]]);
        if (abs(atom->spin[i]) != 0) {
          atom->ervel[i] += -dtfm * atom->erforce[i];
        }
      }
    }
  }
}