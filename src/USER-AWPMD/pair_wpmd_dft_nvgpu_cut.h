//
// Created by yalavrinenko on 10.06.19.
//
#ifdef AWPMD_ENABLE_GPU
#ifdef PAIR_CLASS

PairStyle(wpmd/dft-nvgpu/cut,PairAWPMD_DFT_NVGPUCut)

#else
#ifndef LAMMPS_PAIR_AWPMD_DFT_NVGPU_CUT_H
#define LAMMPS_PAIR_AWPMD_DFT_NVGPU_CUT_H

#include "pair_wpmd_dft_cut.h"

namespace LAMMPS_NS {
  class PairAWPMD_DFT_NVGPUCut : public PairAWPMD_DFTCut {
  public:
    explicit PairAWPMD_DFT_NVGPUCut(class LAMMPS *lammps);

    void settings(int i, char **pString) override;
  };

}

#endif //LAMMPS_PAIR_AWPMD_DFT_NVGPU_CUT_H
#endif
#endif