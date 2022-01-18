//
// Created by yalavrinenko on 09.10.2019.
//

#ifndef DERIVS_CPU_INTEGRATOR_HPP
#define DERIVS_CPU_INTEGRATOR_HPP

#include "adaptive_mesh_integrator.hpp"
#include "numeric.h"
#include <numeric>

template <typename cell_t> class Integrator_cpu {
public:
  template <typename output_value_t, typename FunctionType, typename iterator_t>
  static output_value_t integrate(iterator_t cbegin, iterator_t cend, FunctionType function) {
    return std::accumulate(cbegin, cend, output_value_t{}, [&function](output_value_t sum, cell_t cell){
      double3 a = {cell.begin().x, cell.begin().y, cell.begin().z};
      double3 dh = {cell.end().x - a.x, cell.end().y - a.y, cell.end().z - a.z};
      return sum + IntegrationMethods::simpson_integration<output_value_t>(function, a, dh);
    });
  }
};


#endif //DERIVS_CPU_INTEGRATOR_HPP
