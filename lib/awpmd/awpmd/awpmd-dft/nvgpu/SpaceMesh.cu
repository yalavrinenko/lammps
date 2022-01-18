//
// Created by cheshire on 20.04.17.
//
#include "../DataTypes.hpp"
#include "../utils/SpaceMesh.hpp"

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <iostream>

#include <curand_mtgp32_kernel.h>
#include <device_launch_parameters.h>

#include <iostream>
#include <cstring>
#include "gpu_exceptions.h"

__host__ __device__ unsigned int inline Mesh_nvgpu::address(unsigned int i, unsigned int j, unsigned int k) const {
  return k * (m_width * m_height) + i * m_height + j;
}

Mesh_nvgpu::Mesh_nvgpu(unsigned int width, unsigned int height, unsigned int depth) : m_width(width), m_height(height),
                                                                                      m_depth(depth) {
  grid_size = (width * height * depth);
  Logger::Info("Create GPU Mesh. Size:", width, 'x', height, 'x', depth);
  Logger::Info("Allocate", grid_size * sizeof(float) / (1024.0 * 1024.0), "Mb of host memory for mesh.");
  m_host_storage = new float[grid_size];

  Logger::Info("Allocate", grid_size * sizeof(float) / (1024.0 * 1024.0), "Mb of device memory for mesh.");
  SAFECALL(
      cudaMalloc(&m_device_storage, grid_size * sizeof(float)),
      "Error in allocation device memory."
  );
}

Mesh_nvgpu::Mesh_nvgpu(const Mesh_nvgpu &m) {
  m_width = m.m_width;
  m_height = m.m_height;
  m_depth = m.m_depth;
  grid_size = m.grid_size;

  this->m_host_storage = m.m_host_storage;
  this->m_device_storage = m.m_device_storage;
}

void Mesh_nvgpu::clear() {
  SAFECALL (
      cudaMemset(m_device_storage, 0.0, grid_size),
      "Error in cleaning memory"
  );
}

void Mesh_nvgpu::copy_memory(int direction) {
  cudaMemcpy(m_host_storage, m_device_storage, grid_size * sizeof(float), cudaMemcpyDeviceToHost);
}

void Mesh_nvgpu::clear(int value) {
  SAFECALL (
      cudaMemset(m_device_storage, value, grid_size),
      "Error in cleaning memory"
  );
}

__host__ __device__ float &Mesh_nvgpu::at(unsigned int i, unsigned int j, unsigned int k) {
  return m_device_storage[this->address(i, j, k)];
}

__host__ float &Mesh_nvgpu::ath(unsigned int i, unsigned int j, unsigned int k) {
  return m_host_storage[this->address(i, j, k)];
}

__host__ float Mesh_nvgpu::ath(unsigned int i, unsigned int j, unsigned int k) const {
  return m_host_storage[this->address(i, j, k)];
}


__host__ void Mesh_nvgpu::Free() {
  delete[] m_host_storage;
  if (m_device_storage != nullptr) SAFECALL(cudaFree(m_device_storage), "Error in memory free.");
}

__host__ Mesh_nvgpu::~Mesh_nvgpu() = default;