//
// Created by cheshire on 21.05.17.
//

#ifndef AWPMD_DFT_CUDAKERNELS_HPP
#define AWPMD_DFT_CUDAKERNELS_HPP

template<class AType, class ... ArgsType>
__global__ void nvgpu_dft_init_xc_approximation(IApproximation **ptr, ArgsType ... args) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    (*ptr) = new AType(args...);
  }
};

__global__ void nvgpu_dft_remove_xc_approximation(IApproximation *ptr) {
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    delete ptr;
  }
};

__global__ void nvgpu_dft_xc_check(IApproximation *ptr) {
  auto index = threadIdx.x;

  auto rho = 0.26f * static_cast<float>(index), spin_up = 0.5f, derivs = 1.0f;
  auto eng_p = ptr->energy(rho, spin_up);
  auto eng_k = ptr->kinetic(rho, spin_up);
  auto deriv_eup = ptr->derives(rho, spin_up, derivs, ElectronSpin::E_UP);
  auto deriv_down = ptr->derives(rho, spin_up, derivs, ElectronSpin::E_DOWN);

  printf(
      "Thread %d: Energy (pot, kin) for rho = %0.6f, spin_up_f = %0.6f and deriv_coef = %0.6f: (%0.6f, %0.6f), %0.6f %0.6lf\n",
      index, rho, spin_up, derivs, eng_p, eng_k, deriv_eup, deriv_down);
};
#endif //AWPMD_DFT_CUDAKERNELS_HPP
