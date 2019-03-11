//
// Created by yalavrinenko on 04.03.19.
//

#include "mc_utils.h"

namespace LAMMPS_NS{

  void MC3DVectorSystem::make(size_t size, LAMMPS_NS::mc_stepper &stepper) {
    for (auto i = 0u; i < size; ++i)
      if (_filter(i)){
        for (auto j = 0; j < 3; ++j)
          stepper.make_shift(src[i][j]);
      }
  }
}
