//
// Created by yalavrinenko on 03.12.2019.
//

#ifndef LAMMPS_GPU_FUNCTORS_HPP
#define LAMMPS_GPU_FUNCTORS_HPP

#include "../numeric/adaptive_mesh_integrator.hpp"
#include "../xcfunctionals/IApproximation.hpp"
#include "numeric_gpu.cuh"

template <typename cell_t, typename packet_t>
struct common_functor_env{
  packet_t const *packets{};
  size_t count{};

  __device__ NumericType<float, float> rho_function(float3 const &r) const {
    float rho_v_up = 0.0f, rho_v_down = 0.0f;
    float fraction_v = 0.0f;
    for (auto packet = this->packets; packet != &this->packets[count]; ++packet) {
      if (packet->rdot(r) < (packet->width * 5.0f) * (packet->width * 5.0f)) {
        if (packet->spin == ElectronSpin::E_UP)
          rho_v_up += packet->overlap_at(r);
        else
          rho_v_down += packet->overlap_at(r);
      }
    }
    fraction_v = (rho_v_up - rho_v_down) / (rho_v_up + rho_v_down);

    return {rho_v_up + rho_v_down, fraction_v};
  }
};

template <typename cell_t, typename packet_t>
struct energy_functor: public common_functor_env<cell_t, packet_t>{
private:
  IApproximation *approx = nullptr;

  __device__ NumericType<float, float> xcenergy_function(float x, float y, float z) const {
    auto rho_ = this->rho_function({x, y, z});
    auto xc = this->approx->energy(thrust::get<0>(rho_), thrust::get<1>(rho_));
    auto kin = this->approx->kinetic(thrust::get<0>(rho_), thrust::get<1>(rho_));
    return NumericType<float, float>{xc, kin};
  }

  __device__ NumericType<float, float> operator()(cell_t const c) const{
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    float3 begin {c.begin_.x, c.begin_.y, c.begin_.z};
    return IntegrationMethods::simpson_integration<NumericType<float, float>>([this](float x, float y, float z){
      return this->xcenergy_function(x, y, z);
    }, begin, dh);
  }
};

template <typename cell_t, typename packet_t>
struct energy_functor_sic: public common_functor_env<cell_t, packet_t>{
  IApproximation *approx = nullptr;

  __device__ NumericType<float, float> xcenergy_function_sic(float x, float y, float z) const {
    float rho_v_up = 0.0f, rho_v_down = 0.0f;
    float fraction_v = 0.0f;
    float3 r{x, y, z};
    float sic_xc = 0, sic_kin = 0;

    for (auto packet = this->packets; packet != &this->packets[this->count]; ++packet) {
      if (packet->rdot(r) < (packet->width * 5.0f) * (packet->width * 5.0f)) {
        auto overlap = packet->overlap_at(r);
        sic_xc += this->approx->energy(overlap, static_cast<int>(packet->spin));
        sic_kin += this->approx->kinetic(overlap, static_cast<int>(packet->spin));

        if (packet->spin == ElectronSpin::E_UP)
          rho_v_up += overlap;
        else
          rho_v_down += overlap;
      }
    }
    auto rho_ = rho_v_up + rho_v_down;
    fraction_v = (rho_v_up - rho_v_down) / (rho_);

    auto xc = this->approx->energy(rho_, fraction_v) - sic_xc;
    auto kin = this->approx->kinetic(rho_, fraction_v) - sic_kin;

    return NumericType<float, float>{xc, kin};
  }

  __device__ NumericType<float, float> operator()(cell_t const c) const{
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    float3 begin {c.begin_.x, c.begin_.y, c.begin_.z};
    return IntegrationMethods::simpson_integration<NumericType<float, float>>([this](float x, float y, float z){
      return this->xcenergy_function_sic(x, y, z);
    }, begin, dh);
  }
};

template <typename cell_t, typename packet_t>
struct derivatives_functor: public energy_functor<cell_t, packet_t> {
private:
  __device__  NumericType<float, float, float, float> derivative(float x, float y, float z, size_t packet_index) const {
    auto rho_ = this->rho_function({x, y, z});
    auto r = float3{x, y, z};
    auto dw = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(), this->packets[packet_index].dw(r),
                                    this->packets[packet_index].spin);

    auto dx = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(), this->packets[packet_index].dr<0>(r),
                                    this->packets[packet_index].spin);

    auto dy = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(), this->packets[packet_index].dr<1>(r),
                                    this->packets[packet_index].spin);

    auto dz = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(), this->packets[packet_index].dr<2>(r),
                                    this->packets[packet_index].spin);
    return {dx, dy, dz, dw};
  }

  __device__ NumericType<float, float, float, float> operator() (cell_t const c) const{
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    float3 begin {c.begin_.x, c.begin_.y, c.begin_.z};
    auto packet_index = c.linked_index;
    return IntegrationMethods::simpson_integration<NumericType<float, float, float, float>>([this, packet_index](float x, float y, float z){
      return this->derivative(x, y, z, packet_index);
    }, begin, dh);
  }

  template <typename energy_functor_type>
  static derivatives_functor<cell_t, packet_t> from_energy_functor(energy_functor_type const &efunc) {
    derivatives_functor<cell_t, packet_t> df;
    df.count = efunc.count;
    df.packets = efunc.packets;
    df.approx = efunc.approx;

    return df;
  }
};

template <typename cell_t, typename packet_t>
struct derivatives_functor_sic: public energy_functor_sic<cell_t, packet_t> {
  __device__ NumericType<float, float, float, float> derivative_sic(float x, float y, float z, size_t packet_index) const {
    auto rho_ = this->rho_function({x, y, z});
    auto r = float3{x, y, z};

    auto overlap = this->packets[packet_index].overlap_at(r);

    auto dw = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(),
                                    this->packets[packet_index].dw(r), this->packets[packet_index].spin);

    dw -= this->approx->derives(overlap, static_cast<float>(this->packets[packet_index].spin), this->packets[packet_index].dw(r),
                                this->packets[packet_index].spin);

    auto dx = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(),
                                    this->packets[packet_index].dr<0>(r), this->packets[packet_index].spin);

    auto dy = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(),
                                    this->packets[packet_index].dr<1>(r), this->packets[packet_index].spin);

    auto dz = this->approx->derives(rho_.template get<0>(), rho_.template get<1>(),
                                    this->packets[packet_index].dr<2>(r), this->packets[packet_index].spin);

    return {dx, dy, dz, dw};
  }

  __device__ NumericType<float, float, float, float> operator() (cell_t const &c) const{
    float3 dh {c.end_.x - c.begin_.x, c.end_.y - c.begin_.y, c.end_.z - c.begin_.z};
    float3 begin {c.begin_.x, c.begin_.y, c.begin_.z};
    auto packet_index = c.linked_index;
    return IntegrationMethods::simpson_integration<NumericType<float, float, float, float>>([this, packet_index](float x, float y, float z){
      return this->derivative_sic(x, y, z, packet_index);
    }, begin, dh);
  }

  template <typename energy_functor_type>
  static derivatives_functor_sic<cell_t, packet_t> from_energy_functor(energy_functor_type const &efunc) {
    derivatives_functor_sic<cell_t, packet_t> df;
    df.count = efunc.count;
    df.packets = efunc.packets;
    df.approx = efunc.approx;

    return df;
  }
};
#endif //LAMMPS_GPU_FUNCTORS_HPP
