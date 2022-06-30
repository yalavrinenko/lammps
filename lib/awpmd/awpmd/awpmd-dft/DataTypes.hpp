//
// Created by cheshire on 10.03.17.
//

#ifndef AWPMD_DFT_DATATYPES_HPP
#define AWPMD_DFT_DATATYPES_HPP

#include "../../ivutils/include/wavepacket.h"

#include "utils/Logger.hpp"
#include "utils/SpaceMesh.hpp"
#include "xcfunctionals/IApproximation.hpp"
#include <complex>
#include <iostream>

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#ifndef __global__
#define __global__
#endif

#ifndef __VECTOR_TYPES_H__
class double3{
public:
  double x, y, z;
};

class float3 {
public:
  float x, y, z;
};

class uint3 {
public:
  unsigned int x, y, z;
};
#endif

template <class RealType>
struct UnitsSystem{
  RealType Distance2Bohr = 1;
  RealType Hartree2Energy = 1;
};

union ufloat3 {
  class float3 vector{};

  float array[3];

  __host__ __device__ explicit ufloat3(float3 const &vec) {
    vector = vec;
  }

  __host__ __device__ explicit ufloat3(float const *arr) {
    vector.x = arr[0];
    vector.y = arr[1];
    vector.z = arr[2];
  }

  __host__ __device__ ufloat3(float x, float y, float z) {
    vector.x = x;
    vector.y = y;
    vector.z = z;
  }
};

template<class RealType>
struct GaussPacket {

  using fields_type = RealType;

  RealType r[3];
  RealType width;
  ElectronSpin spin{};

  GaussPacket() = default;

  __host__ GaussPacket(class WavePacket const &wp, UnitsSystem<RealType> const &units, ElectronSpin spin) :
      width(wp.get_width() * units.Distance2Bohr),
      r{static_cast<RealType>(wp.get_r()[0]) * units.Distance2Bohr,
        static_cast<RealType>(wp.get_r()[1]) * units.Distance2Bohr,
        static_cast<RealType>(wp.get_r()[2]) * units.Distance2Bohr}, spin{spin} {
  }

  template<typename real_t>
  __host__ GaussPacket(real_t const* coord, real_t width, UnitsSystem<RealType> const &units, ElectronSpin spin) :
      width(width * units.Distance2Bohr),
      r{static_cast<RealType>(coord[0]) * units.Distance2Bohr,
        static_cast<RealType>(coord[1]) * units.Distance2Bohr,
        static_cast<RealType>(coord[2]) * units.Distance2Bohr}, spin{spin} {
  }

  __host__ __device__ void printPacket() const {
    printf("r=(%lf, %lf, %lf) width=%lf\n", r[0], r[1], r[2], width);
  }

  template<typename VectorType3D>
  __host__ __device__ inline RealType overlap_at(VectorType3D const &x) const {
    constexpr auto overlap_prefactor = static_cast<RealType>(0.3299226101861591); //std::pow(3.0 / (2.0 * M_PI), 3.0 / 2.0 );
    return static_cast<RealType>(overlap_prefactor / (width * width * width) * std::exp(
        static_cast<RealType>(-3.0) / (static_cast<RealType>(2.0) * width * width) * rdot(x)));
  }

  template<typename VectorType3D>
  __host__ __device__ inline RealType rdot(VectorType3D const &x) const {
    return (x.x - r[0]) * (x.x - r[0]) + (x.y - r[1]) * (x.y - r[1]) + (x.z - r[2]) * (x.z - r[2]);
  }

  template<typename VectorType3D>
  __host__ __device__ inline RealType dw(VectorType3D const &x) const{
    return overlap_at(x) * (static_cast<RealType>(-3.0) / width + static_cast<RealType>(3.0) / (width * width * width) * rdot(x));
  }

  template <unsigned int index, typename VectorType3D>
  __host__ __device__ inline RealType dr(VectorType3D const &x) const{
    decltype(x.x) v[3] = {x.x, x.y, x.z};
    return static_cast<RealType>(-3.0) / (width * width) * overlap_at(x) * (v[index] - r[index]);
  }

  template <typename VectorType3D>
  __host__ __device__ inline RealType I(VectorType3D const &a, VectorType3D const &b) const{
    constexpr auto overlap_prefactor = static_cast<RealType>(0.125); //std::pow(3.0 / (2.0 * M_PI), 3.0 / 2.0 ) * std::pow(M_PI / 6.0, 3.0 / 2.0);
    constexpr auto erf_factor = static_cast<RealType>(1.224744871391589); //std::sqrt(3.0 / 2.0);

    auto Ix = std::erf(erf_factor/width * (r[0] - a.x)) - std::erf(erf_factor/width * (r[0] - b.x));
    auto Iy = std::erf(erf_factor/width * (r[1] - a.y)) - std::erf(erf_factor/width * (r[1] - b.y));
    auto Iz = std::erf(erf_factor/width * (r[2] - a.z)) - std::erf(erf_factor/width * (r[2] - b.z));

    return overlap_prefactor * Ix * Iy * Iz;
  }

  template<typename VectorType3D>
  __host__ __device__ inline bool in_range(VectorType3D const &begin, VectorType3D const &end) const {
    return (begin.x <= r[0] && r[0] <= end.x) &&
           (begin.y <= r[1] && r[1] <= end.y) &&
           (begin.z <= r[2] && r[2] <= end.z);
  }

  __host__ static inline unsigned derivatives_count() {
    return 4;
  }
};

typedef float (*DerivFunction)(GaussPacket<float> const&, float3 &);

typedef float (*deriv_function)(GaussPacket<float> const&, float3 &);

template<class TValue>
struct MeshSize {
  union {
    struct {
      TValue x, y, z;
    } as_struct;
    TValue as_array[3];
  } size;

  TValue& operator [] (size_t index) { return size.as_array[index]; }
};

template <class MeshType>
class ElectronDensity {
public:
  MeshType total;
  MeshType spinup_fraction;

  ElectronDensity() = default;

  explicit ElectronDensity(MeshSize<unsigned int> const &ms):
          ElectronDensity(ms.size.as_struct.x, ms.size.as_struct.y, ms.size.as_struct.z){
  }

  ElectronDensity(unsigned int width, unsigned int height, unsigned int depth) {
    total = MeshType(width, height, depth);
    spinup_fraction = MeshType(width, height, depth);
  }

  void Free() {
    total.Free();
    spinup_fraction.Free();
  }
};

struct DFTConfig {
  MeshSize<double> mesh_start, mesh_fin;
  MeshSize<unsigned int> mesh_size{};

  MeshSize<float> mesh_step{};

  unsigned int packet_number;

  IApproximation *approximation = nullptr;
  IApproximation *approximation_device = nullptr;

  bool use_adaptive_mesh = true;
  double min_cell = 0.8;
  double max_distance = 2.5;
  size_t force_cell_bins = 10;

  UnitsSystem<float> units{};
  bool calc_force = false;

  size_t nodes = 1;
  size_t node_rank = 0;

  bool use_xc_tables = false;

  DFTConfig() {
    packet_number = 0;
  }
};

#endif //AWPMD_DFT_DATATYPES_HPP
