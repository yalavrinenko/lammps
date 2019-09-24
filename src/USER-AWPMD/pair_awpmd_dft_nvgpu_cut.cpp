//
// Created by yalavrinenko on 10.06.19.
//
#ifdef AWPMD_ENABLE_GPU
#include "pair_awpmd_dft_nvgpu_cut.h"
#include <GPU-CUDA/awpmd-dft-nvgpu.hpp>
#include <DerivativesFunction.hpp>
#include <atom.h>
#include <comm.h>
#include <cstring>

LAMMPS_NS::PairAWPMD_DFT_NVGPUCut::PairAWPMD_DFT_NVGPUCut(LAMMPS_NS::LAMMPS *lammps) :
  PairAWPMD_DFTCut(lammps, nullptr) {
}

void LAMMPS_NS::PairAWPMD_DFT_NVGPUCut::settings(int argc, char **pString) {
  PairAWPMDCut::settings(argc, pString);

  auto electron_count = std::count_if(atom->spin, atom->spin + atom->nlocal + atom->nghost,
                                      [](auto &spin) { return std::abs(spin) == 1; });

  int gpu_per_node = 1;

  for (auto i = 0; i < argc; ++i){
    if (!std::strcmp(pString[i], "gppn")){
      gpu_per_node = force->numeric(FLERR, pString[i+1]);
    }
  }

  overlap_derivs = DerivsFunction_NVGPU::GetFunctions();
  xc_energy_ = new XCEnergy_nvgpu(electron_count, make_dft_config(), comm->me % gpu_per_node);
}
#endif