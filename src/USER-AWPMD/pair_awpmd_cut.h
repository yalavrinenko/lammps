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


#ifdef PAIR_CLASS

PairStyle(awpmd/cut,PairAWPMDCut)

#else

#ifndef LMP_PAIR_AWPMD_CUT_H
#define LMP_PAIR_AWPMD_CUT_H

#include "pair.h"
#include <vector>
#include <unordered_map>

class AWPMD_split;


namespace LAMMPS_NS {

  class PairAWPMDCut : public Pair {
   friend class FixNVEAwpmd;

  public:
   explicit PairAWPMDCut(class LAMMPS *);

   ~PairAWPMDCut() override;

   void compute(int, int) override;

   void settings(int, char **) override;

   void coeff(int, char **) override;

   void init_style() override;

   void min_pointers(double **, double **);

   double init_one(int, int) override;

   void write_restart(FILE *) override;

   void read_restart(FILE *) override;

   void write_restart_settings(FILE *) override;

   void read_restart_settings(FILE *) override;

   void min_xf_pointers(int, double **, double **) override;

   void min_xf_get(int) override;

   void min_x_set(int) override;

   double memory_usage() override;

  private:

   struct awpmd_pair_index{
     unsigned lmp_index{};
     unsigned wpmd_index{};
   };

   using awpmd_ions = std::vector<awpmd_pair_index >;
   using awpmd_electrons = std::unordered_map<unsigned, std::vector<awpmd_pair_index>>;
   using awpmd_packets = std::tuple<awpmd_ions, awpmd_electrons>;

   awpmd_packets make_packets() const;

   //SOME PROP FUNCTION
   double ghost_energy();

   void init_wpmd(awpmd_ions &ions, awpmd_electrons &electrons);

   int flexible_pressure_flag;
   double cut_global;
   double **cut;


   int nmax; // number of additional variables for minimizer
   double *min_var, *min_varforce; // additional variables for minimizer

   void allocate();

   void virial_eradius_compute();


   AWPMD_split *wpmd; // solver oybject
   double ermscale; // scale of width mass for motion
   double width_pbc; // setting for width pbc
   double half_box_length; // calculated by coeff function
  };

}

#endif
#endif
