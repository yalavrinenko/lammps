//
// Created by yalavrinenko on 26.09.2019.
//

#ifndef DERIVS_NUMERIC_H
#define DERIVS_NUMERIC_H

template <int index, class TupleType>
struct TupleOperation{
  static void add(TupleType &t, TupleType const &v){
    std::get<index>(t) += std::get<index>(v);
    TupleOperation<index-1, TupleType>::add(t, v);
  }

  static void sub(TupleType &t, TupleType const &v){
    std::get<index>(t) -= std::get<index>(v);
    TupleOperation<index-1, TupleType>::sub(t, v);
  }

  template <class T>
  static void mul(TupleType &t, T const &v){
    std::get<index>(t) *= v;
    TupleOperation<index-1, TupleType>::mul(t, v);
  }
};

template <class TupleType>
struct TupleOperation<0, TupleType>{
  static void add(TupleType &t, TupleType const &v){
    std::get<0>(t) += std::get<0>(v);
  }

  static void sub(TupleType &t, TupleType const &v){
    std::get<0>(t) -= std::get<0>(v);
  }

  template <class T>
  static void mul(TupleType &t, T const &v){
    std::get<0>(t) *= v;
  }
};

template <class ... Args>
struct NumericType : public std::tuple<Args...>{
  NumericType() = default;

  NumericType(Args... args): std::tuple<Args...>(args...){
  }

  NumericType operator+ (NumericType<Args...> const &v) const{
    NumericType nv{*this};
    TupleOperation<sizeof...(Args)-1, NumericType<Args...>>::add(nv, v);
    return nv;
  }

  NumericType operator- (NumericType<Args...> const &v) const{
    NumericType nv{*this};
    TupleOperation<sizeof...(Args)-1, NumericType<Args...>>::sub(nv, v);
    return nv;
  }

  template <class T>
  NumericType operator* (T const &v) const{
    NumericType nv{*this};
    TupleOperation<sizeof...(Args)-1, NumericType<Args...>>::mul(nv, v);
    return nv;
  }
};

template <class T, class ... Args>
NumericType<Args...> operator* (T const &v, NumericType<Args...> const &t){
  return t * v;
}

class IntegrationMethods{
public:
  template<class OutputType, class IntegralFunction>
   static OutputType simpson_integration(IntegralFunction f, double3 r, double3 dh) {
    auto g = [&f, &r, &dh](double x, double y) -> OutputType {
      return dh.z / 6.0 * (f(x, y, r.z) + 4.0 * f(x, y, r.z + 0.5 * dh.z) + f(x, y, r.z + dh.z));
    };
    auto h = [&g, &r, &dh](double x) -> OutputType {
      return dh.y / 6.0 * (g(x, r.y) + 4.0 * g(x, r.y + 0.5 * dh.y) + g(x, r.y + dh.y));
    };

    return dh.x / 6.0 * (h(r.x) + 4.0 * h(r.x + 0.5 * dh.x) + h(r.x + dh.x));
  }
};

#endif //DERIVS_NUMERIC_H
