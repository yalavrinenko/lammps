//
// Created by cheshire on 20.04.17.
//

#ifndef AWPMD_DFT_SPACEMESH_HPP
#define AWPMD_DFT_SPACEMESH_HPP

#include <functional>
#include <memory>

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#ifndef __global__
#define __global__
#endif


class Mesh {
private:
  double *m_host_storage = nullptr;
  double *m_device_storage = nullptr;
  unsigned int grid_size;

  unsigned int inline address(unsigned int i, unsigned int j, unsigned int k) const;

public:
  unsigned int m_width, m_height, m_depth;

  Mesh(unsigned int width, unsigned int height, unsigned int depth);

  Mesh(const Mesh &m);

  Mesh() : m_device_storage(nullptr), m_host_storage(nullptr), grid_size(0) {
  }

  void clear();

  double &at(unsigned int i, unsigned int j, unsigned int k);

  double &ath(unsigned int i, unsigned int j, unsigned int k);

  double ath(unsigned int i, unsigned int j, unsigned int k) const;

  void apply(std::function<void(double &)> function) {
    for (unsigned int i = 0; i < this->m_width; ++i)
      for (unsigned int j = 0; j < this->m_height; ++j)
        for (unsigned int k = 0; k < this->m_depth; k++) {
          function(this->ath(i, j, k));
        }
  }

  void apply(std::function<void(double)> function) const {
    for (unsigned int i = 0; i < this->m_width; ++i)
      for (unsigned int j = 0; j < this->m_height; ++j)
        for (unsigned int k = 0; k < this->m_depth; k++) {
          function(this->ath(i, j, k));
        }
  }

  double reduction();

  void clear(int value);

  double *dev_ptr() { return m_device_storage; }

  double *host_ptr() { return m_host_storage; }

  double const* begin() const { return m_device_storage; }
  double const* end() const { return &m_device_storage[this->count()]; }

  void copy_memory(int direction = 0);

  void Free();

  size_t count() const {
    return grid_size;
  }

  ~Mesh();
};

class Mesh_nvgpu {
private:
  float *m_host_storage = nullptr;
  float *m_device_storage = nullptr;
  unsigned int grid_size;

  __host__ __device__ unsigned int inline address(unsigned int i, unsigned int j, unsigned int k) const;

public:
  unsigned int m_width, m_height, m_depth;

  Mesh_nvgpu(unsigned int width, unsigned int height, unsigned int depth);

  Mesh_nvgpu(const Mesh_nvgpu &m);

  Mesh_nvgpu() : m_device_storage(nullptr), m_host_storage(nullptr), grid_size(0) {
  }

  void clear();

  __host__ __device__ float &at(unsigned int i, unsigned int j, unsigned int k);

  float &ath(unsigned int i, unsigned int j, unsigned int k);

  float ath(unsigned int i, unsigned int j, unsigned int k) const;

  void apply(std::function<void(float &)> function) {
    for (unsigned int i = 0; i < this->m_width; ++i)
      for (unsigned int j = 0; j < this->m_height; ++j)
        for (unsigned int k = 0; k < this->m_depth; k++) {
          function(this->ath(i, j, k));
        }
  }

  void apply(std::function<void(float)> function) const {
    for (unsigned int i = 0; i < this->m_width; ++i)
      for (unsigned int j = 0; j < this->m_height; ++j)
        for (unsigned int k = 0; k < this->m_depth; k++) {
          function(this->ath(i, j, k));
        }
  }

  void clear(int value);

  float *dev_ptr() { return m_device_storage; }

  float *host_ptr() { return m_host_storage; }

  void copy_memory(int direction = 0);

  void Free();

  __host__ __device__ size_t count() const {
    return grid_size;
  }

  ~Mesh_nvgpu();
};

#endif //AWPMD_DFT_SPACEMESH_HPP
