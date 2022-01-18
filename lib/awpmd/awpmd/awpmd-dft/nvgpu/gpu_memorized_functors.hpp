//
// Created by yalavrinenko on 03.12.2019.
//

#ifndef LAMMPS_GPU_MEMORIZED_FUNCTORS_HPP
#define LAMMPS_GPU_MEMORIZED_FUNCTORS_HPP

#include "../numeric/memorized_adaptive_mesh_cell.hpp"
#include "gpu_functors.hpp"

template<typename _CellType, typename _PacketType>
struct density_functor: public common_functor_env<_CellType, _PacketType>{
  __device__ void operator() (_CellType &c){
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    using RealType = decltype(dh.x);
#pragma unroll
    for (auto i = 0; i < 3; ++i)
#pragma unroll
      for (auto j = 0; j < 3; ++j)
#pragma unroll
        for (auto k = 0; k < 3; ++k){
          float3 r{c.begin_.x + i * dh.x / static_cast<RealType>(2.0),
                   c.begin_.y + j * dh.y / static_cast<RealType>(2.0),
                   c.begin_.z + k * dh.z / static_cast<RealType>(2.0)};
          c.density[i][j][k] = this->rho_function(r);
        }
  }
};

template <typename _CellType, typename _PacketType>
struct memorized_energy_functor: public energy_functor<_CellType, _PacketType>{
private:
  IApproximation *approx = nullptr;

  __device__ NumericType<float, float> xcenergy_function(_CellType const &c, int i, int j, int k) const {
    auto rho_ = c.density[i][j][k];
    auto xc = this->approx->energy(thrust::get<0>(rho_), thrust::get<1>(rho_));
    auto kin = this->approx->kinetic(thrust::get<0>(rho_), thrust::get<1>(rho_));
    return NumericType<float, float>{xc, kin};
  }

  static memorized_energy_functor<_CellType, _PacketType> from_density_functor(density_functor<_CellType, _PacketType> const &dfunc) {
    memorized_energy_functor<_CellType, _PacketType> ef;
    ef.m_spin_up_index = dfunc.m_spin_up_index;
    ef.count = dfunc.count;
    ef.packets = dfunc.packets;
    return ef;
  }

  __device__ NumericType<float, float> operator()(_CellType const c) const{
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    return IntegrationMethods::simpson_integration<NumericType<float, float>>([this, &c](int i, int j, int k){
      return this->xcenergy_function(c, i, j, k);
    }, dh);
  }
};

template <typename _CellType, typename _PacketType>
struct memorized_energy_functor_sic: public energy_functor<_CellType, _PacketType>{
private:
  IApproximation *approx = nullptr;

  __device__ NumericType<float, float> xcenergy_function_sic(_CellType const &c, int i, int j, int k) const {
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    float3 r{c.begin_.x + static_cast<float>(i) * dh.x / 2.0f,
             c.begin_.y + static_cast<float>(j) * dh.y / 2.0f,
             c.begin_.z + static_cast<float>(k) * dh.z / 2.0f};

    float sic_xc = 0.0f, sic_kin = 0.0f;

    for (auto packet = this->packets; packet != &this->packets[this->count]; ++packet) {
      if (packet->rdot(r) < (packet->width * 5.0f) * (packet->width * 5.0f)) {
        auto overlap = packet->overlap_at(r);
        sic_xc += this->approx->energy(overlap, 1);
        sic_kin += this->approx->kinetic(overlap, 1);
      }
    }

    auto rho_ = thrust::get<0>(c.density[i][j][k]);
    auto fraction_v = thrust::get<1>(c.density[i][j][k]);

    auto xc = this->approx->energy(rho_, fraction_v) - sic_xc;
    auto kin = this->approx->kinetic(rho_, fraction_v) - sic_kin;

    return NumericType<float, float>{xc, kin};
  }

  static memorized_energy_functor_sic<_CellType, _PacketType> from_density_functor(density_functor<_CellType, _PacketType> const &dfunc) {
    memorized_energy_functor<_CellType, _PacketType> ef;
    ef.m_spin_up_index = dfunc.m_spin_up_index;
    ef.count = dfunc.count;
    ef.packets = dfunc.packets;
    return ef;
  }

  __device__ NumericType<float, float> operator()(_CellType const c) const{
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    return IntegrationMethods::simpson_integration<NumericType<float, float>>([this, &c](int i, int j, int k){
      return this->xcenergy_function_sic(c, i, j, k);
    }, dh);
  }

};

template <typename _CellType, typename _PacketType>
struct memorized_derivatives_functor: public memorized_energy_functor<_CellType, _PacketType> {
private:
  size_t packet_index{};

  __device__  NumericType<float, float, float, float> derivative(_CellType const &c, int i, int j, int k) const {
    auto rho_ = c.density[i][j][k];
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    float3 r{c.begin_.x + static_cast<float>(i) * dh.x / 2.0f,
             c.begin_.y + static_cast<float>(j) * dh.y / 2.0f,
             c.begin_.z + static_cast<float>(k) * dh.z / 2.0f};

    auto dw = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(), this->packets[packet_index].dw(r),
                                    (packet_index >= this->m_spin_up_index.y) ? ElectronSpin::E_DOWN : ElectronSpin::E_UP);

    auto dx = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(), this->packets[packet_index].dr<0>(r),
                                    (packet_index >= this->m_spin_up_index.y) ? ElectronSpin::E_DOWN : ElectronSpin::E_UP);

    auto dy = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(), this->packets[packet_index].dr<1>(r),
                                    (packet_index >= this->m_spin_up_index.y) ? ElectronSpin::E_DOWN : ElectronSpin::E_UP);

    auto dz = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(), this->packets[packet_index].dr<2>(r),
                                    (packet_index >= this->m_spin_up_index.y) ? ElectronSpin::E_DOWN : ElectronSpin::E_UP);
    return {dx, dy, dz, dw};
  }

  template <typename energy_functor_type>
  static memorized_derivatives_functor<_CellType, _PacketType> from_energy_functor(energy_functor_type const &efunc) {
    memorized_derivatives_functor<_CellType, _PacketType> df;
    df.m_spin_up_index = efunc.m_spin_up_index;
    df.count = efunc.count;
    df.packets = efunc.packets;
    df.approx = efunc.approx;

    return df;
  }

  __device__ NumericType<float, float, float, float> operator() (_CellType const c) const{
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    return IntegrationMethods::simpson_integration<NumericType<float, float, float, float>>([this, &c](int i, int j, int k){
      return this->derivative(c, i, j, k);
    }, dh);

//    return IntegrationMethods::trapz_integration<NumericType<float, float, float, float>>([this, &c](int i, int j, int k){
//      return this->derivative_sic(c, i, j, k);
//    }, dh);
  }
};

template <typename _CellType, typename _PacketType>
struct memorized_derivatives_functor_sic: public memorized_energy_functor_sic<_CellType, _PacketType> {
private:
  size_t packet_index{};

  __device__ NumericType<float, float, float, float> derivative_sic(_CellType const &c, int i, int j, int k) const {
    auto rho_ = c.density[i][j][k];
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    float3 r{c.begin_.x + static_cast<float>(i) * dh.x / 2.0f,
             c.begin_.y + static_cast<float>(j) * dh.y / 2.0f,
             c.begin_.z + static_cast<float>(k) * dh.z / 2.0f};

    auto overlap = this->packets[packet_index].overlap_at(r);

    auto spin = (packet_index >= this->m_spin_up_index.y) ? ElectronSpin::E_DOWN : ElectronSpin::E_UP;

    auto dw = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(),
                                    this->packets[packet_index].dw(r), spin);

    dw -= this->approx->derives(overlap, ((spin == ElectronSpin::E_UP) ? 1 : -1), this->packets[packet_index].dw(r),
                                spin);

    auto dx = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(),
                                    this->packets[packet_index].dr<0>(r), spin);

    auto dy = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(),
                                    this->packets[packet_index].dr<1>(r), spin);

    auto dz = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(),
                                    this->packets[packet_index].dr<2>(r), spin);

    return {dx, dy, dz, dw};
  }

  template <typename energy_functor_type>
  static memorized_derivatives_functor_sic<_CellType, _PacketType> from_energy_functor(energy_functor_type const &efunc) {
    memorized_derivatives_functor<_CellType, _PacketType> df;
    df.m_spin_up_index = efunc.m_spin_up_index;
    df.count = efunc.count;
    df.packets = efunc.packets;
    df.approx = efunc.approx;

    return df;
  }

  __device__ NumericType<float, float, float, float> operator() (_CellType const &c) const{
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    return IntegrationMethods::simpson_integration<NumericType<float, float, float, float>>([this, &c](int i, int j, int k){
      return this->derivative_sic(c, i, j, k);
    }, dh);

//    return IntegrationMethods::trapz_integration<NumericType<float, float, float, float>>([this, &c](int i, int j, int k){
//      return this->derivative_sic(c, i, j, k);
//    }, dh);
  }
};

#endif //LAMMPS_GPU_MEMORIZED_FUNCTORS_HPP
