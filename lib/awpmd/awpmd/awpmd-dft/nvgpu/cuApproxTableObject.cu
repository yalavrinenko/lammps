//
// Created by yalavrinenko on 13.08.18.
//

#include "cuApproxTableObject.hpp"
#include <cuda_runtime.h>
#include <driver_types.h>
#include "../DataTypes.hpp"
#include "gpu_exceptions.h"

gpuApproxTablesObject::gpuApproxTablesObject(IApproximation *host_approx, IApproximation *dev_approx,
                                             ApproximationTableProps const &propertys) :
    origin(dev_approx), props(propertys) {

  bind_texture(build_table([host_approx](float rho, float sup) {
    return host_approx->derives(rho, sup, 1, ElectronSpin::E_UP);
  }), this->t_derivatives[0]);

  bind_texture(build_table([host_approx](float rho, float sup) {
    return host_approx->derives(rho, sup, 1, ElectronSpin::E_DOWN);
  }), this->t_derivatives[1]);

  bind_texture(build_table([host_approx](float rho, float sup) {
    return host_approx->kinetic(rho, sup);
  }), this->t_kinetic);

  bind_texture(build_table([host_approx](float rho, float sup) {
    return host_approx->energy(rho, sup);
  }), this->t_potential);
}

std::unique_ptr<Mesh_nvgpu>
gpuApproxTablesObject::build_table(const std::function<float(float, float)>& func) const {
  auto mesh = std::unique_ptr<Mesh_nvgpu>(new Mesh_nvgpu(props.mesh_size[0] + 1, props.mesh_size[1] + 1, 1));

  float rho_step = (props.density_range[1] - props.density_range[0]) / (props.mesh_size[0]);
  float const &rho_zero = props.density_range[0];

  float sup_step = (props.spin_fraction_range[1] - props.spin_fraction_range[0]) / (props.mesh_size[1]);
  float const &sup_zero = props.spin_fraction_range[0];

  mesh->clear();
  mesh->copy_memory();

  for (int i = 0; i <= props.mesh_size[0]; ++i) {
    auto rho = rho_zero + i * rho_step;
    for (int j = 0; j <= props.mesh_size[1]; ++j) {
      auto sup_f = sup_zero + j * sup_step;
      mesh->ath(i, j, 0) = func(rho, sup_f);
    }
  }

  return mesh;
}

void gpuApproxTablesObject::bind_texture(std::unique_ptr<Mesh_nvgpu> ptr, gpuApproxTablesObject::TableObject &object) {
  cudaArray_t data_storage;

  auto channel_description = cudaCreateChannelDesc<float>();
  SAFECALL(
      cudaMallocArray(&data_storage, &channel_description, props.mesh_size[1] + 1, props.mesh_size[0] + 1,
                      cudaArrayDefault),
      "Error in approximation memory allocation"
  );

  SAFECALL(
      cudaMemcpyToArray(data_storage, 0, 0, ptr->host_ptr(), ptr->count() * sizeof(float), cudaMemcpyHostToDevice),
      "Error in approximation table memory copy."
  );

  ptr->Free();

  auto &tex = object.m_texture;
  object.m_array = data_storage;
  object.scale = {(props.density_range[1] - props.density_range[0]) / static_cast<float>(props.mesh_size[0]),
                  (props.spin_fraction_range[1] - props.spin_fraction_range[0]) / static_cast<float>(props.mesh_size[1])};

  auto resdescr = cudaResourceDesc{};
  resdescr.resType = cudaResourceTypeArray;
  resdescr.res.array.array = data_storage;

  auto texdescr = cudaTextureDesc{};
  texdescr.addressMode[0] = cudaAddressModeClamp;
  texdescr.addressMode[1] = cudaAddressModeClamp;
  texdescr.filterMode = cudaFilterModeLinear;
  texdescr.readMode = cudaReadModeElementType;
  texdescr.normalizedCoords = false;

  SAFECALL(
      cudaCreateTextureObject(&tex, &resdescr, &texdescr, nullptr),
      "Error in texture object creation."
  );
}

__device__ float gpuApproxTablesObject::TableObject::value(float rho, float spin_ratio) const {
  return tex2D<float>(m_texture, (spin_ratio + 1.0f) / scale.y + 0.5f, rho / scale.x + 0.5f);
}
