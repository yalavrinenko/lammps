//
// Created by yalavrinenko on 21.02.19.
//

#ifndef LAMMPS_MC_UTILS_H
#define LAMMPS_MC_UTILS_H

#include "random_park.h"
namespace LAMMPS_NS{

  struct mc_step{
    double old[3];
    double shift;

    int index;
    int tag;
    bool active;

    enum class mc_type{
      coord = 0,
      vel = 1,
      width = 2,
      pwidth = 3
    } type;
  };

}

#endif //LAMMPS_MC_UTILS_H
