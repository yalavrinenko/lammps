//
// Created by yalavrinenko on 10.06.19.
//

#include "pair_awpmd_dft_nvgpu_cut.h"
#include <GPU-CUDA/awpmd-dft-nvgpu.hpp>
#include <atom.h>

LAMMPS_NS::PairAWPMD_DFT_NVGPUCut::PairAWPMD_DFT_NVGPUCut(LAMMPS_NS::LAMMPS *lammps) :
  PairAWPMD_DFTCut(lammps, nullptr) {
  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });

  xc_energy_ = new XCEnergy_nvgpu(electron_count, make_dft_config());
}
