//
// Created by yalavrinenko on 02.08.18.
//

#include "nd4mesh.hpp"
#include "../DataTypes.hpp"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include "gpu_exceptions.h"

nd4mesh::nd4mesh(unsigned int *topology):
    W(topology[0]),
    H(topology[1]),
    D(topology[2]),
    L(topology[3]),
    grid_size(topology[0] * topology[1] * topology[2]){
    alloc_impl();
}

__host__ __device__ unsigned int inline nd4mesh::address(unsigned int i, unsigned int j, unsigned int k){
    return k*(W * H) + i * W + j;
}

float &nd4mesh::at(unsigned int i, unsigned int j, unsigned int k) {
    return this->m_current_ptr[this->address(i,j,k)];
}

void nd4mesh::clear() {
    for (auto dev_ptr : this->m_device_pointers) {
        SAFECALL(
                cudaMemset(dev_ptr, 0.0, grid_size),
                "Error in cleaning memory"
        );
    }
}

__device__ void nd4mesh::change_level(unsigned int level) {
    this->m_current_ptr = this->m_dev_ptrs[level];
}

void nd4mesh::alloc_impl() {
    m_device_pointers.resize(L, nullptr);

    Logger::Info("Create GPU Mesh. Size:", W,'x',H,'x',D,'x',L);
    Logger::Info("Allocate", grid_size * L * sizeof(float) / (1024.0 * 1024.0), "Mb of device memory for mesh.");

    for (auto &ptr : m_device_pointers){
        SAFECALL(
                cudaMalloc(&ptr, sizeof(float)*grid_size),
                "Error in memory allocation for submesh."
                );
    }

    SAFECALL(
            cudaMalloc(&this->m_dev_ptrs, sizeof(float*) * this->L),
            "Error in memory allocation for meshes ptr."
            );

    SAFECALL(
            cudaMemcpy(this->m_dev_ptrs, m_device_pointers.data(), sizeof(float*)*L, cudaMemcpyHostToDevice),
            "Error in ptr copy to device."
            );
}

__global__ void ND4MeshTest(nd4mesh mesh){
    auto idx = threadIdx.x;

    mesh.change_level(idx);
    auto ptr = mesh.current_ptr();

    printf("%d ==> %p\n", idx, ptr);
}

void nd4mesh::free() {
    Logger::Info("nd4mesh deallocation.");
    if (m_dev_ptrs != nullptr) {
        SAFECALL(cudaFree(m_dev_ptrs), "Error in memory free.");
        m_dev_ptrs = nullptr;
        m_current_ptr = nullptr;
    }

    for (auto ptr : m_device_pointers){
        SAFECALL(cudaFree(ptr), "Error in memory free.");
    }

    m_device_pointers.clear();
}
