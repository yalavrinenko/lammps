//
// Created by yalavrinenko on 03.12.2019.
//

#ifndef LAMMPS_MEMORIZED_ADAPTIVE_MESH_CELL_HPP
#define LAMMPS_MEMORIZED_ADAPTIVE_MESH_CELL_HPP
#include "adaptive_mesh_integrator.hpp"

template <typename real_t, typename numeric_t>
struct MemorizedAdaptiveMeshCell: public AdaptiveMeshCell<real_t>{
  using AdaptiveMeshCell<real_t>::AdaptiveMeshCell;
  numeric_t density[3][3][3];
};

template <typename real_t>
struct LinkedMeshCell: public AdaptiveMeshCell<real_t>{
  using AdaptiveMeshCell<real_t>::AdaptiveMeshCell;
  size_t linked_index;
};
#endif //LAMMPS_MEMORIZED_ADAPTIVE_MESH_CELL_HPP
