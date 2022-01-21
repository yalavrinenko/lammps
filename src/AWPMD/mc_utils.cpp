//
// Created by yalavrinenko on 04.03.19.
//

#include "mc_utils.h"

namespace LAMMPS_NS {
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
