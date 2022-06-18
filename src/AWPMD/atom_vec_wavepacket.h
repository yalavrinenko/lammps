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


#ifdef ATOM_CLASS

AtomStyle(wavepacket,AtomVecWavepacket)

#else

#ifndef LMP_ATOM_VEC_WAVEPACKET_H
#define LMP_ATOM_VEC_WAVEPACKET_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecWavepacket : public AtomVec {
public:
  AtomVecWavepacket(class LAMMPS *);

  void grow_pointers() override;
  void force_clear(int, size_t) override;
  void create_atom_post(int) override;
  void data_atom_post(int) override;
  int property_atom(const std::string &) override;
  void pack_property_atom(int, double *, int, int) override;

private:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;

  ///\en spin: -1 or 1 for electron, 0 for ion (compatible with eff)
  int *spin;
  ///\en charge: must be specified in the corresponding units (-1 for electron in real units, eff compatible)
  double *q;
  ///\en width of the wavepacket (compatible with eff)
  double *eradius;
  ///\en width velocity for the wavepacket (compatible with eff)
  double *ervel;
  ///\en (generalized) force on width  (compatible with eff)
  double *erforce;

  // AWPMD- specific:
  ///\en electron tag: must be the same for the WPs belonging to the same electron
  int *etag;
  ///\en wavepacket split coefficients: cre, cim, size is 2*N
  double **cs;
  ///\en force on wavepacket split coefficients: re, im, size is 2*N
  double **csforce;
  ///\en (generalized) force on velocity, size is 3*N
  double **vforce;
   ///\en (generalized) force on radius velocity, size is N
  double *ervelforce;
};

}

#endif
#endif
