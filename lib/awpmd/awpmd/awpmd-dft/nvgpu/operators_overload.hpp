//
// Created by cheshire on 18.05.17.
//

#ifndef AWPMD_DFT_CUDALAUNCHCONFIG_HPP
#define AWPMD_DFT_CUDALAUNCHCONFIG_HPP

#include "../xcfunctionals/ApproxList.hpp"
#include <cstdlib>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#define SHUTUP_CLION //Remove in release
#ifdef SHUTUP_CLION

#include "cuda_runtime_api.h"

#endif

__host__ __device__ uint3 operator+(const uint3 &a, const uint3 &b);

__host__ __device__ uint3 operator*(const uint3 &a, const uint3 &b);

__host__ __device__ uint3 operator*(const uint3 &a, const dim3 &b);

__host__ __device__ uint3 operator*(const dim3 &a, const uint3 &b);

#endif //AWPMD_DFT_CUDALAUNCHCONFIG_HPP
