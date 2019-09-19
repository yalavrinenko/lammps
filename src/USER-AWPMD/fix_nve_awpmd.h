/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing author: Ilya Valuev (JIHT RAS)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nve/awpmd,FixNVEAwpmd)

#else

#ifndef LMP_FIX_NVE_awpmd_H
#define LMP_FIX_NVE_awpmd_H

#include "../fix_nve.h"
#include "pair_awpmd_cut.h"


namespace LAMMPS_NS {

class FixNVEAwpmd : public FixNVE {
 public:
  FixNVEAwpmd(class LAMMPS *, int, char **);
  void initial_integrate(int) override;
  void final_integrate() override;
};

}

#endif
#endif
