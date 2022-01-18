//
// Created by yalavrinenko on 02.08.18.
//

#ifndef MDUTILS_LIB_NDMESH_H
#define MDUTILS_LIB_NDMESH_H

#include <vector>

#ifndef __CUDACC__
#define __host__
#define __device__
#endif

/**
 * 4DMesh with dimension W x H x D x L
 * Looks like L meshes with dim W x H x D
 */
class nd4mesh{
public:

    nd4mesh(): m_current_ptr(nullptr), m_dev_ptrs(nullptr){
    }

    explicit nd4mesh(unsigned int topology[4]);

    nd4mesh(nd4mesh const &mesh) = default;

    std::vector<float*> const & dev_ptrs() const{
        return m_device_pointers;
    }

    __host__ __device__ float& at(unsigned int i, unsigned int j, unsigned int k);

    void clear();

    __device__ void change_level(unsigned int level);

    __host__ __device__  unsigned int inline address(unsigned int i, unsigned int j, unsigned int k);

    void free();

    ~nd4mesh() = default;

        __device__ float* current_ptr() const{
        return m_current_ptr;
    }
    ///

    unsigned int W = 0, H = 0, D = 0, L = 0;
    unsigned int grid_size = 0;
private:

    void alloc_impl();

    float* m_current_ptr;
    float** m_dev_ptrs;

    std::vector<float*> m_device_pointers;
};

#endif //MDUTILS_LIB_NDMESH_H
