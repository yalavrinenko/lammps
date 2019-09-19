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
#include <map>
#include <vector_3.h>

class AWPMD_split;
class WavePacket;

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

    double init_one(int, int) override;

    void write_restart(FILE *) override;

    void read_restart(FILE *) override;

    void write_restart_settings(FILE *) override;

    void read_restart_settings(FILE *) override;

    void min_xf_pointers(int, double **, double **) override;

    void min_xf_get(int) override;

    void min_x_set(int) override;

    double memory_usage() override;

    AWPMD_split *awpmd();

    std::vector<WavePacket> const& electrons_packets() const {
      return packets;
    }
  protected:

    struct awpmd_energies{
      double ke{};
      double ee{};
      double ei{};
      double ii{};

      double sum() const {
        return ke + ee + ei + ii;
      }
    };

    struct awpmd_pair_index {
      int tag{};
      unsigned lmp_index{};
      unsigned wpmd_index{};

      awpmd_pair_index(unsigned l_index, unsigned w_index, int tag) : lmp_index(l_index), wpmd_index(w_index),
                                                                           tag(tag) {}

      bool operator < (awpmd_pair_index const &v) const{
        return this->tag < v.tag;
      }
    };

    using awpmd_ions = std::vector<awpmd_pair_index>;
    using awpmd_electrons = std::map<unsigned, std::vector<awpmd_pair_index>>;
    using awpmd_packets = std::tuple<awpmd_ions, awpmd_electrons>;

    awpmd_packets make_packets() const;

    void init_wpmd(awpmd_ions &ions, awpmd_electrons &electrons);

    awpmd_energies compute_pair();

    int flexible_pressure_flag;
    double cut_global;
    double **cut;


    int nmax; // number of additional variables for minimizer
    double *min_var, *min_varforce; // additional variables for minimizer

    void update_energy(double full_coul_energy, awpmd_ions const &ions, awpmd_electrons const &electrons);

    void update_force(awpmd_ions const &ions, awpmd_electrons const &electrons, std::vector<Vector_3> const &fi);

    void allocate();

    void virial_eradius_compute();

    AWPMD_split *wpmd; // solver oybject
    double ermscale; // scale of width mass for motion
    double width_pbc; // setting for width pbc
    double half_box_length; // calculated by coeff function
    double electron_ke_;

    double const time_unit = 1; //10.12;
    double const ev_to_energy = 1; //1.0 / 27.211386 * 627.509474;

    std::vector<WavePacket> packets;

  private:
    void check_with_native_wpmd(double coul_energy);
  };

}

#endif
#endif
