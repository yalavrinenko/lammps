//
// Created by yalavrinenko on 09.04.2020.
//

#include "compute_electric_current.h"
#include "atom.h"
#include "memory.h"
#include <update.h>

void LAMMPS_NS::ComputeElectricCurrent::compute_vector() {
  //auto mean_v = compute_cm_velocity();
  std::array<double, 3> current{0.0, 0.0, 0.0};
  for (auto i = 0u; i < atom->nlocal; ++i){
    if (atom->mask[i] & groupbit){
      for (auto k: {0, 1, 2}){
        current[k] += atom->v[i][k] * atom->q[i];
      }
    }
  }

  for (auto k: {0, 1, 2}){
    MPI_Allreduce(&current[k], &vector[k], 1, MPI_DOUBLE, MPI_SUM, world);
  }

  invoked_vector = update->ntimestep;
}

LAMMPS_NS::ComputeElectricCurrent::ComputeElectricCurrent(LAMMPS_NS::LAMMPS *lammps, int i, char **pString) : Compute(
    lammps, i, pString) {
  scalar_flag = array_flag = false;

  vector_flag = true;
  size_vector = 3;
  extvector = 1;

  vector = new double[size_vector];
}

LAMMPS_NS::ComputeElectricCurrent::~ComputeElectricCurrent() {
  delete[] vector;
}
std::array<double, 3>
LAMMPS_NS::ComputeElectricCurrent::compute_cm_velocity() const {
  std::array<double, 3> mean{0.0, 0.0, 0.0};
  double M = 0.0;
  for (auto i = 0u; i < atom->nlocal; ++i){
    if (atom->mask[i] & groupbit){
      for (auto k: {0, 1, 2}) mean[k] += atom->v[i][k] * atom->mass[atom->type[i]];
      M += atom->mass[atom->type[i]];
    }
  }
  for (auto k: {0, 1, 2}) mean[k] /= M;
  return mean;
}
void LAMMPS_NS::ComputeElectricCurrent::init() {
}
