//
// Created by yalavrinenko on 09.04.2020.
//

#ifdef COMPUTE_CLASS
ComputeStyle(ecurrent,ComputeElectricCurrent)
#else
#ifndef LAMMPS_COMPUTE_ELECTRIC_CURRENT_H
#define LAMMPS_COMPUTE_ELECTRIC_CURRENT_H
#include "compute.h"
#include <array>
namespace LAMMPS_NS {
  class ComputeElectricCurrent : public Compute {
  public:
    ComputeElectricCurrent(class LAMMPS *lammps, int i, char **pString);

    void compute_vector() override;
    void init() override;

    ~ComputeElectricCurrent() override ;

  protected:
    std::array<double, 3> compute_cm_velocity() const;
  };
}
#endif //LAMMPS_COMPUTE_ELECTRIC_CURRENT_H
#endif