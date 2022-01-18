//
// Created by cheshire on 11.05.17.
//

#ifndef AWPMD_DFT_DERIVATIVESFUNCTION_HPP
#define AWPMD_DFT_DERIVATIVESFUNCTION_HPP

#include "../DataTypes.hpp"
#include "vector"

class DerivsFunction {
public:
  template<int deriv_index>
  __host__ __device__ static float dr(GaussPacket<float> const &b, float3 &r) {
    return b.dr<deriv_index>(r);
  }

  __host__ __device__ static float dw(GaussPacket<float> const &b, float3 &r) {
    return b.dw(r); // * (-3.0f / b.width + 3.0f / (b.width * b.width * b.width) * b.rdot(r));
  }

  static std::vector<DerivFunction> GetFunctions();
};

class DerivsFunction_NVGPU {
public:
  __device__ static float EnergyByA_Re(GaussPacket<float> &b, float3 &r) {
    return 1.0; //2.0f * ((float)b.coeff.lambda.real() - b.dot(ufloat3(r), ufloat3(r))) * b.at(r);
  }

  __device__ static float EnergyByA_Im(GaussPacket<float> &b, float3 &r) {
    return 0.0f;//-2.0f * (float)b.coeff.nu.imag() * b.at(r);
  }

  template<unsigned int index>
  __device__ static float EnergyByB_Re(GaussPacket<float> &b, float3 &r) {
    return 1.0; //2.0f * b.at(r) * ( b.coeff.u.Re[index]  + ( (index == 0) ? r.x: (index == 1) ? r.y : r.z) );
  }

  template<unsigned int index>
  __device__ static float EnergyByB_Im(GaussPacket<float> &b, float3 &r) {
    return 0.0f; //2.0f * b.at(r) * b.coeff.v.Re[index];
  }

  static std::vector<DerivFunction> GetFunctions();
};

#endif //AWPMD_DFT_DERIVATIVESFUNCTION_HPP
