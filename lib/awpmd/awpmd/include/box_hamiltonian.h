// Copyright 2014, Nikita Kazeev, JIHT RAS
#ifndef BOX_HAMILTONIAN_H
#define BOX_HAMILTONIAN_H

//#include "./wpmd.h"
#include <assert.h>
#include <complex>
#include <vector>
#include <limits>       // std::numeric_limits

#include "cvector_3.h"
#include "wavepacket.h"
#include "tcpdefs.h"
#include "../../ivutils/include/cerf.h"

#ifdef _WIN32
inline double isnan(double x) { return _isnan(x); }
#else
//inline double isnan(double x) { return std::isnan(x); }
#endif

const double MIN_LOG = log(std::numeric_limits<double>::min());

inline bool isnan(const cdouble &x) {
  return std::isnan(real(x)) || std::isnan(imag(x));
}

inline bool isnan(const cVector_3 &x) {
  for (size_t i = 0; i < 3; ++i) {
    if (isnan(x[i]))
      return true;
  }
  return false;
}

// _check functions check for underflow and if it happens return 0
// To be used with caution, as they don't solve 0*inf uncertainty
inline cdouble cerf_check(const cdouble &z) {
  if (real(z) < sqrt(-MIN_LOG)) {
    return cerf(z);
  } else {
    return real(z) >= 0 ? 1. : -1.;
  }
}

inline cdouble exp_check(const cdouble &z) {
  if (real(z) > MIN_LOG) {
    return exp(z);
  } else {
    return 0.;
  }
}

inline cdouble cerfc(const cdouble &z) {
  return 1. - cerf(z);
}

cdouble packet_norm(
    const cdouble &a,
    const cVector_3 &b);

inline cdouble exped_double_norm(
    const cdouble &a1,
    const cVector_3 &b1,
    const cdouble &a2,
    const cVector_3 &b2) {
  return exp(conj(packet_norm(a1, b1)) + packet_norm(a2, b2));
}

inline cdouble log_double_norm(
    const cdouble &a1,
    const cVector_3 &b1,
    const cdouble &a2,
    const cVector_3 &b2) {
  return conj(packet_norm(a1, b1)) + packet_norm(a2, b2);
}

cdouble norm_deriv_a_real(
    const cdouble a,
    const cVector_3 b);

cdouble norm_deriv_a_imag(
    const cdouble a,
    const cVector_3 b);

cVector_3 norm_deriv_b_real(
    const cdouble a,
    const cVector_3 b);

cVector_3 norm_deriv_b_imag(
    const cdouble a,
    const cVector_3 b);

void get_norm_derivative(
    const cdouble a,
    const cVector_3 b,
    cdouble *a_diff,
    cVector_3 *b_diff);


double factorial(int n);


class HermitePolinomial2 {
protected:
  vector<double> coef;
  int max_n, m_size;

public:
  HermitePolinomial2(int max_order = 10);

  template<class T>
  T hermite(int n, T x, T y);
};


class BoxHamiltonian {
private:
  Vector_3 boundaries; // The box is -boundary; +boundary
  Vector_3 potential_multiplers;

  // f = packet*k*(x-x0)^2
  void set_axis(
      const cVector_3 &b,
      const unsigned char axis,
      cdouble *b_axial,
      cdouble *b_perp1,
      cdouble *b_perp2,
      double *boundary,
      double *force_multiplier) const;

  // Auxiliary variables for projection calculations
  int max_n;
  HermitePolinomial2 hpol;
  vector<double> proj_prefactor[3];
  Vector_3 alpha_k;

  void calc_proj_prefactors();

public:
  BoxHamiltonian(const Vector_3 boundary = Vector_3(0, 0, 0),
                 const Vector_3 potential_multipler = Vector_3(1., 1., 1.),
                 const int max_proj_order = -1);

  // TODO(kazeevn) make private. Use FRIEND_TEST for testing
  cdouble axis_integral(
      const cdouble &a,
      const cVector_3 &b,
      const cdouble &log_norm,
      const unsigned char axis) const;

  // TODO(kazeevn) make private. Use FRIEND_TEST for testing
  void get_axis_integral_derivative(
      const cdouble &a,
      const cVector_3 &b,
      const cdouble &log_norm,
      const unsigned char axis,
      cdouble *a_diff,
      cVector_3 *b_diff) const;

  // TODO(kazeevn) make private. Use FRIEND_TEST for testing
  inline cdouble get_integral(
      const cdouble &a1,
      const cVector_3 &b1,
      const cdouble &a2,
      const cVector_3 &b2) const {
    const cdouble a = conj(a1) + a2;
    const cVector_3 b = conj(b1) + b2;
    const cdouble log_norm = log_double_norm(a1, b1, a2, b2);
    return \
      axis_integral(a, b, log_norm, 0) + \
      axis_integral(a, b, log_norm, 1) + \
      axis_integral(a, b, log_norm, 2);
  }

  // TODO(valuev) Optimize:
  // Cache integral OR have user supply integral OR
  // make get_integral set parameters and get_derivatives to recall them
  // OR use flag
  void get_derivatives(
      const cdouble &a1,
      const cVector_3 &b1,
      const cdouble &a2,
      const cVector_3 &b2,
      cdouble *integral,
      cdouble *a1_real_diff,
      cdouble *a1_imag_diff,
      cVector_3 *b1_real_diff,
      cVector_3 *b1_imag_diff,
      cdouble *a2_real_diff,
      cdouble *a2_imag_diff,
      cVector_3 *b2_real_diff,
      cVector_3 *b2_imag_diff) const;

  cdouble eigen_proj(int n, const WavePacket &wp, int axis);

  ///\en Gets current maximum projection order.
  int get_max_proj_order() const {
    return max_n;
  }

  Vector_3 get_force(const Vector_3 &x, double *dE = NULL) const {
    Vector_3 force;
    if (dE)
      *dE = 0.;
    for (int i = 0; i < 3; i++) {
      if (x[i] < -boundaries[i]) {
        force[i] = -(x[i] + boundaries[i]) * potential_multiplers[i] * 2.;
        if (dE)
          *dE -= force[i] * (x[i] + boundaries[i]) / 2.;
      } else if (x[i] > boundaries[i]) {
        force[i] = -(x[i] - boundaries[i]) * potential_multiplers[i] * 2.;
        if (dE)
          *dE -= force[i] * (x[i] - boundaries[i]) / 2.;
      }
    }
    return force;
  }

  ///\en Gets a vector of coeefficients k[i] in the 3D energy expression sum ki*xi^2
  Vector_3 get_directional_coeffs() const {
    return potential_multiplers;
  }

  ///\en Returns oscillator quants in each direction
  Vector_3 get_eigen_energies() const {
    Vector_3 h_omega;
    for (int i = 0; i < 3; i++)
      h_omega[i] = h_plank * sqrt(2. * potential_multiplers[i] / m_electron);
    return h_omega;
  }

  ///\en Returns the widths of 0th gaussian state in each spatial direction.
  Vector_3 get_eigen_widths() const {
    Vector_3 gamma, eigenE = get_eigen_energies();
    for (int i = 0; i < 3; i++)
      gamma[i] = sqrt(3. / 2 / m_electron / eigenE[i]) * h_plank;
    return gamma;
  }

  double occ_prob_fermi(int axis, int level, int ne, double T) const {
    level++;
    double beta_de = h_plank * sqrt(2. * potential_multiplers[axis] / m_electron) / T;
    double prob = 0.;
    double sign = 1;
    for (int g = 0; g < ne; g++) {
      double prod = sign;
      for (int p = 0; p < g + 1; p++) {
        prod *= exp((ne - p - level) * beta_de) - exp(-level * beta_de);
      }
      sign = -sign;
      prob += prod;
    }
    return prob;
  }

  Vector_3 get_boundaries() const {
    return boundaries;
  }

};

#endif
