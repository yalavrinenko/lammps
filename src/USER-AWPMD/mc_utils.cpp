//
// Created by yalavrinenko on 04.03.19.
//

#include "mc_utils.h"

namespace LAMMPS_NS {

  void MC3DVectorSystem::make(size_t size, LAMMPS_NS::mc_stepper &stepper) {
    for (auto i = 0u; i < size; ++i)
      if (_filter(i)) {
        for (auto j = 0; j < 3; ++j)
          stepper.make_shift(src[i][j]);
      }
  }

  std::vector<double> MC3DVectorSystem::pack(size_t size, int *tag) const {
    std::vector<double> buffer;
    buffer.reserve(size + size * 3);
    for (auto i = 0u; i < size; ++i) {
      if (_filter(i)) {
        buffer.push_back(tag[i]);
        for (auto j = 0; j < 3; ++j)
          buffer.push_back(src[i][j]);
      }
    }
    return buffer;
  }

  size_t MC3DVectorSystem::unpack(const double *data, size_t size, std::unordered_map<int, int> const &ghost_map) {
    size_t unpacked = 0;
    auto iter = 0;
    while (iter < size) {
      int tag = (int) data[iter++];
      if (ghost_map.count(tag)) {
        for (auto j = 0u; j < 3; ++j)
          src[ghost_map.at(tag)][j] = data[iter++];
        ++unpacked;
      } else {
        iter += 3;
      }
    }
    return unpacked;
  }

  void MCScalarSystem::make(size_t size, mc_stepper &stepper) {
    for (auto i = 0u; i < size; ++i)
      if (_filter(i)) {
        stepper.make_shift(src[i]);
      }
  }

  std::vector<double> MCScalarSystem::pack(size_t size, int *tag) const {
    std::vector<double> buffer;
    buffer.reserve(size * 2);
    for (auto i = 0u; i < size; ++i) {
      if (_filter(i)) {
        buffer.push_back(tag[i]);
        buffer.push_back(src[i]);
      }
    }
    return buffer;
  }

  size_t MCScalarSystem::unpack(const double *data, size_t size, std::unordered_map<int, int> const &ghost_map) {
    size_t unpacked = 0;
    auto iter = 0;
    while (iter < size) {
      int tag = (int) data[iter++];
      if (ghost_map.count(tag)) {
        src[ghost_map.at(tag)] = data[iter++];
        ++unpacked;
      } else {
        iter += 1;
      }
    }
    return unpacked;
  }
}
