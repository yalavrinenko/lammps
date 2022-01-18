//
// Created by yalavrinenko on 13.08.18.
//

#ifndef MDUTILS_LIB_TABLEAPPROXIMATION_HPP
#define MDUTILS_LIB_TABLEAPPROXIMATION_HPP

#include "IApproximation.hpp"

template<class TablesObjectType>
class TableApproximation : public IApproximation {
public:
  __devspec__ explicit TableApproximation(TablesObjectType tables) :
      approx_tables(tables) {
    Type = approx_tables.origin->Type;
  }

  __devspec__ float derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const override {
    return approx_tables.check(rho) ?
           approx_tables.t_derivatives[((ElectronSpin::E_UP == spin) ? 0 : 1)].value(approx_tables.density_shift(rho), spin_ratio) * derivative :
           approx_tables.origin->derives(rho, spin_ratio, derivative, spin);
  }

  __devspec__ float energy(float rho, float spin_ratio) const override {
    return approx_tables.check(rho) ? approx_tables.t_potential.value(approx_tables.density_shift(rho), spin_ratio) :
           approx_tables.origin->energy(rho, spin_ratio);
  }

  __devspec__ float kinetic(float rho, float spin_ratio) const override {
    return approx_tables.check(rho) ? approx_tables.t_kinetic.value(approx_tables.density_shift(rho), spin_ratio) :
           approx_tables.origin->kinetic(rho, spin_ratio);
  }

private:
  TablesObjectType approx_tables;
};

#endif //MDUTILS_LIB_TABLEAPPROXIMATION_HPP
