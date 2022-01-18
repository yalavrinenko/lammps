//
// Created by yalavrinenko on 13.08.18.
//

#ifndef MDUTILS_LIB_CUAPPROXTABLEOBJECT_HPP
#define MDUTILS_LIB_CUAPPROXTABLEOBJECT_HPP

#include "../utils/SpaceMesh.hpp"
#include "../xcfunctionals/IApproximation.hpp"
#include <cuda.h>
#include <cuda_runtime_api.h>

struct ApproximationTableProps{
    float density_range[2] = {0.0f, 100.0f};
    float spin_fraction_range[2] = {-1.0f, 1.0f};
    unsigned mesh_size[2] = {10000, 100};
};

class gpuApproxTablesObject{
public:
    struct TableObject{
        __device__ float value(float rho, float spin_ratio) const;

        __host__ void free(){
            if (m_array)
                cudaFreeArray(m_array);
            if (m_texture)
                cudaDestroyTextureObject(m_texture);
        }

        cudaTextureObject_t m_texture;
        cudaArray_t m_array;
        float2 scale;
    };


    gpuApproxTablesObject(IApproximation *host_approx, IApproximation *dev_approx,
                              ApproximationTableProps const &propertys);

    void free() {
        t_derivatives[0].free();
        t_derivatives[1].free();
        t_kinetic.free();
        t_potential.free();
    }

    __device__ inline bool check(float rho) const {
        return (props.density_range[0] <= rho && rho < props.density_range[1]);
    }

    __device__ inline float density_shift(float rho) const {
      return rho - props.density_range[0];
    }

    IApproximation * origin;
    TableObject t_kinetic, t_potential, t_derivatives[2];
private:
    std::unique_ptr<Mesh_nvgpu>
    build_table(const std::function<float(float, float)>& func) const;

    void bind_texture(std::unique_ptr<Mesh_nvgpu> ptr,
            TableObject &object);

    ApproximationTableProps props;
};

#endif //MDUTILS_LIB_CUAPPROXTABLEOBJECT_HPP
