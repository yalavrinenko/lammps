//
// Created by yalavrinenko on 16.03.2020.
//

#ifdef PAIR_CLASS

PairStyle(wpmd/cut,PairWPMD)

#else

#ifndef LAMMPS_PAIRWPMD_H
#define LAMMPS_PAIRWPMD_H
#include "WavepacketPairCommon.h"
#include <vector>

namespace LAMMPS_NS {
  class PairWPMD: public WavepacketPairCommon {
  public:
    explicit PairWPMD(class LAMMPS *lmp): WavepacketPairCommon(lmp) {}
    void settings(int i, char **pString) override;

  protected:
    awpmd_energies compute_energy_force() override;

    bool is_pbc_{false};
  };
}

#endif //LAMMPS_PAIRWPMD_H
#endif
