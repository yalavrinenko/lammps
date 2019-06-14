//
// Created by yalavrinenko on 10.06.19.
//
#ifdef PAIR_CLASS

PairStyle(awpmd/dft/cut,PairAWPMD_DFTCut)

#else
#ifndef LAMMPS_PAIR_AWPMD_DFT_CUT_H
#define LAMMPS_PAIR_AWPMD_DFT_CUT_H

#include "pair_awpmd_cut.h"
#include <awpmd-dft.hpp>
#include <force.h>

namespace LAMMPS_NS {
  class PairAWPMD_DFTCut : public PairAWPMDCut{
  public:
    explicit PairAWPMD_DFTCut(class LAMMPS *lammps);

    PairAWPMD_DFTCut(class LAMMPS *lammps, XCEnergy* xc_energy_ptr);

    void compute(int i, int i1) override;

  protected:
    DFTConfig make_dft_config();

    void set_units(){
      UnitsScale.distance_to_bohr = 1.0 / (0.52917721092 * force->angstrom);
      UnitsScale.hartree_to_energy = 627.509474; //only for real
    }

    XCEnergy* xc_energy_;

    struct {
      double distance_to_bohr = 1.0;
      double hartree_to_energy = 1.0;
    } UnitsScale;

    union {
      struct {
        double xc_energy;
        double kinetic_energy;
      } like_vars;

      double like_vector[sizeof(like_vars) / sizeof(double)];
    } output{};
  };
}

#endif //LAMMPS_PAIR_AWPMD_DFT_CUT_H
#endif