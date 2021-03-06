//
// Created by yalavrinenko on 10.06.19.
//
#ifdef AWPMD_ENABLE_GPU
#include "pair_wpmd_dft_nvgpu_cut.h"
#include <nvgpu/awpmd-dft-nvgpu.hpp>
#include <atom.h>
#include <comm.h>
#include <cstring>
#include <utils/DerivativesFunction.hpp>

LAMMPS_NS::PairAWPMD_DFT_NVGPUCut::PairAWPMD_DFT_NVGPUCut(LAMMPS_NS::LAMMPS *lammps) :
  PairAWPMD_DFTCut(lammps, nullptr) {
}

void LAMMPS_NS::PairAWPMD_DFT_NVGPUCut::settings(int argc, char **pString) {
  LAMMPS_NS::PairWPMD::settings(argc, pString);

  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });
  int gpu_per_node = 1;
  for (auto i = 0; i < argc; ++i){
    if (!std::strcmp(pString[i], "gppn")){
      gpu_per_node = static_cast<int>(force->numeric(FLERR, pString[i+1]));
    }
  }

  auto gpu_num = comm->me % gpu_per_node;
  xc_energy_ = new XCEnergy_nvgpu(electron_count, make_dft_config(argc, pString), gpu_num);
}
#endif