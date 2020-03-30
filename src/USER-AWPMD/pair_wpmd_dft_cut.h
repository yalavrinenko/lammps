//
// Created by yalavrinenko on 10.06.19.
//
#ifdef AWPMD_ENABLE_DFT
#ifdef PAIR_CLASS

PairStyle(wpmd/dft/cut,PairAWPMD_DFTCut)

#else
#ifndef LAMMPS_PAIR_AWPMD_DFT_CUT_H
#define LAMMPS_PAIR_AWPMD_DFT_CUT_H

#include "pair_wpmd_cut.h"
#include <awpmd-dft.hpp>
#include <force.h>
#include <vector>

namespace LAMMPS_NS {
  class PairAWPMD_DFTCut : public PairWPMD{
  public:
    explicit PairAWPMD_DFTCut(class LAMMPS *lammps);

    PairAWPMD_DFTCut(class LAMMPS *lammps, XCEnergy* xc_energy_ptr);

    void compute(int i, int i1) override;

    void settings(int i, char **pString) override;

    ~PairAWPMD_DFTCut() override{
      delete  xc_energy_;
    }
  protected:
    DFTConfig make_dft_config(int i, char **pString);

    bool calc_force_ = false;

    void tally_electron_force(unsigned electron_id, std::vector<float> const& force_array);

    double wpmd_kinetic() const;

    XCEnergy* xc_energy_ = nullptr;

    union {
      struct {
        double xc_energy;
        double kinetic_energy;
      } like_vars;

      double like_vector[sizeof(like_vars) / sizeof(double)];
    } output{};

    std::vector<XCEnergy::WavePacketInfo> electrons;
  };
}

#endif //LAMMPS_PAIR_AWPMD_DFT_CUT_H
#endif
#endif