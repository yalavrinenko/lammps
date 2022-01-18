//
// Created by yalavrinenko on 01.10.2019.
//

#ifndef DERIVS_NUMERIC_GPU_H
#define DERIVS_NUMERIC_GPU_H
#include <thrust/tuple.h>

template <int index, class TupleType>
struct TupleOperation{
  __device__ static void add(TupleType &t, TupleType const &v){
    thrust::get<index>(t) += thrust::get<index>(v);
    TupleOperation<index-1, TupleType>::add(t, v);
  }

  __device__ static void sub(TupleType &t, TupleType const &v){
    thrust::get<index>(t) -= thrust::get<index>(v);
    TupleOperation<index-1, TupleType>::sub(t, v);
  }

  template <class T>
  __device__ static void mul(TupleType &t, T const &v){
    thrust::get<index>(t) *= v;
    TupleOperation<index-1, TupleType>::mul(t, v);
  }
};

template <class TupleType>
struct TupleOperation<0, TupleType>{
  __device__ static void add(TupleType &t, TupleType const &v){
    thrust::get<0>(t) += thrust::get<0>(v);
  }

  __device__ static void sub(TupleType &t, TupleType const &v){
    thrust::get<0>(t) -= thrust::get<0>(v);
  }

  template <class T>
  __device__ static void mul(TupleType &t, T const &v){
    thrust::get<0>(t) *= v;
  }
};

template <class ... Args>
struct NumericType : public thrust::tuple<Args...>{
  __host__ __device__ NumericType() = default;

  __device__ NumericType(Args... args): thrust::tuple<Args...>(args...){
  }

  __device__ NumericType operator+ (NumericType<Args...> const &v) const{
    NumericType nv{*this};
    TupleOperation<sizeof...(Args)-1, NumericType<Args...>>::add(nv, v);
    return nv;
  }

  __device__ NumericType operator- (NumericType<Args...> const &v) const{
    NumericType nv{*this};
    TupleOperation<sizeof...(Args)-1, NumericType<Args...>>::sub(nv, v);
    return nv;
  }

  template <class T>
  __device__ NumericType operator* (T const &v) const{
    NumericType nv{*this};
    TupleOperation<sizeof...(Args)-1, NumericType<Args...>>::mul(nv, v);
    return nv;
  }

  template <unsigned index>
  __host__ __device__ typename thrust::tuple_element<index, NumericType<Args...>>::type & operator() (){
    return thrust::get<index>(*this);
  }
};

template <class T, class ... Args>
__device__ NumericType<Args...> operator* (T const &v, NumericType<Args...> const &t){
  return t * v;
}

class IntegrationMethods{
public:
  template<class OutputType, class RangeType, class IntegralFunction>
  __device__ static OutputType simpson_integration(IntegralFunction f, RangeType r, RangeType dh) {
    auto g = [&f, &r, &dh](float x, float y) -> OutputType {
      return dh.z / 6.0f * (f(x, y, r.z) + 4.0f * f(x, y, r.z + 0.5f * dh.z) + f(x, y, r.z + dh.z));
    };
    auto h = [&g, &r, &dh](float x) -> OutputType {
      return dh.y / 6.0f * (g(x, r.y) + 4.0f * g(x, r.y + 0.5f * dh.y) + g(x, r.y + dh.y));
    };

    return dh.x / 6.0f * (h(r.x) + 4.0f * h(r.x + 0.5f * dh.x) + h(r.x + dh.x));
  }

  template<class OutputType, class RangeType, class IntegralFunction>
  __device__ static OutputType simpson_integration(IntegralFunction f, RangeType dh) {
    auto g = [&f, &dh](int x, int y) -> OutputType {
      return dh.z / 6.0f * (f(x, y, 0) + 4.0f * f(x, y, 1) + f(x, y, 2));
    };
    auto h = [&g, &dh](int x) -> OutputType {
      return dh.y / 6.0f * (g(x, 0) + 4.0f * g(x, 1) + g(x, 2));
    };

    return dh.x / 6.0f * (h(0) + 4.0f * h(1) + h(2));
  }

  template<class OutputType, class RangeType, class IntegralFunction>
  __device__ static OutputType trapz_integration(IntegralFunction f, RangeType dh) {
    auto g = [&f, &dh](int x, int y) -> OutputType {
      return dh.z / 2.0f * (f(x, y, 0) + f(x, y, 2));
    };
    auto h = [&g, &dh](int x) -> OutputType {
      return dh.y / 2.0f * (g(x, 0) + g(x, 2));
    };

    return dh.x / 2.0f * (h(0) + h(2));
  }

  template<class OutputType, class RangeType, class IntegralFunction>
  __host__ static OutputType simpson_integration_host(IntegralFunction f, RangeType r, RangeType dh) {
    auto g = [&f, &r, &dh](float x, float y) -> OutputType {
      return dh.z / 6.0f * (f(x, y, r.z) + 4.0f * f(x, y, r.z + 0.5f * dh.z) + f(x, y, r.z + dh.z));
    };
    auto h = [&g, &r, &dh](float x) -> OutputType {
      return dh.y / 6.0f * (g(x, r.y) + 4.0f * g(x, r.y + 0.5f * dh.y) + g(x, r.y + dh.y));
    };

    return dh.x / 6.0f * (h(r.x) + 4.0f * h(r.x + 0.5f * dh.x) + h(r.x + dh.x));
  }
};

#endif //DERIVS_NUMERIC_GPU_H
