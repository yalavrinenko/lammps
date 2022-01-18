// Copyright Nikita Kazeev, Ilya Valuev, Igor Morozov, JIHT RAS 2014
#include "./box_hamiltonian.h"

#ifdef _WIN32
inline int isinf(double x) {
  return !_finite(x);
}
#else
// inline double isinf(double x) { return std::isinf(x); }
#endif

const cdouble I(0, 1);

BoxHamiltonian::BoxHamiltonian(
    const Vector_3 boundary,
    const Vector_3 potential_multipler,
    const int max_proj_order)
    : max_n(max_proj_order), hpol(max_proj_order) {
  boundaries = boundary;
  potential_multiplers = potential_multipler;

  calc_proj_prefactors();
}

void BoxHamiltonian::set_axis(
    const cVector_3& b,
    const unsigned char axis,
    cdouble* b_axial,
    cdouble* b_perp1,
    cdouble* b_perp2,
    double* boundary,
    double* force_multiplier) const {
  *b_axial = b[axis];
  *b_perp1 = b[(axis + 1) % 3];
  *b_perp2 = b[(axis + 2) % 3];
  *boundary = boundaries[axis];
  *force_multiplier = potential_multiplers[axis];
}

cdouble packet_norm(
    const cdouble& a,
    const cVector_3& b) {
  // Assuming Re(a) > 0. Otherwise the packet is not a packet.
  assert(real(a) > 0);
  return 3. / 4. * (log(real(a)) - log(M_PI / 2.)) -
         real(b) * real(b) / 4. / real(a) +
         cdouble(0, 1) *
         (real(b) / 2. / real(a) *
          (real(b) * imag(a) / 2. / real(a) - imag(b)));
}

cdouble norm_deriv_a_real(
    const cdouble a,
    const cVector_3 b) {
  // TODO(kazeevn) use nice vector representation
  return
      (3*real(a)*real(a) - cdouble(0, 2)*imag(a)*
       (real(b[0])*real(b[0]) + real(b[1])*real(b[1]) + real(b[2])*real(b[2])) +
       real(a)*(cdouble(0, 2)*imag(b[0])*real(b[0]) + real(b[0])*real(b[0]) +
                cdouble(0, 2)*imag(b[1])*real(b[1]) +
                real(b[1])*real(b[1]) + cdouble(0, 2)*imag(b[2])*real(b[2]) +
                real(b[2])*real(b[2])))/(4.*real(a)*real(a)*real(a));
}

cdouble norm_deriv_a_imag(
    const cdouble a,
    const cVector_3 b) {
  return cdouble(0, 0.25)*(real(b)*real(b))/real(a)/real(a);
}

cVector_3 norm_deriv_b_real(
    const cdouble a,
    const cVector_3 b) {
  return -(I*imag(b)*real(a)+conj(a)*real(b))/2./real(a)/real(a);
}

cVector_3 norm_deriv_b_imag(
    const cdouble a,
    const cVector_3 b) {
  return cdouble(0, -0.5)*real(b)/real(a);
}

cdouble BoxHamiltonian::axis_integral(
    const cdouble& a,
    const cVector_3& b,
    const cdouble& log_norm,
    const unsigned char axis) const {
  // Returns integral from interaction with a single dimension of the box.
  // Axis - axis in question.
  // log_norm - the result is multiplied by exp(log_norm)
  cdouble b_axial, b_perp1, b_perp2;
  double X0, k;
  set_axis(b, axis, &b_axial, &b_perp1, &b_perp2, &X0, &k);
  const cdouble b_perp_norm = b_perp1*b_perp1 + b_perp2*b_perp2;
  const cdouble general_coeff = k * pow(a, -3.5) * M_PI / 8.;

  const cdouble log_coeff = b_perp_norm/4./a - a*X0*X0 + log_norm;
  const cdouble log_coeff_one = log_coeff - X0*b_axial;
  const cdouble log_coeff_two = log_coeff + X0*b_axial;
  cdouble part_one, part_two;
  part_one = -exp_check(log(b_axial + 2*a*X0) + log_coeff_one);
  part_two = exp_check(log(b_axial - 2*a*X0) + log_coeff_two);

  const cdouble log_coeff_three = (b_perp_norm+b_axial*b_axial)/4./a+log_norm;
  const cdouble part_three =
      sqrt(M_PI)*(b_axial*b_axial + 2*a*(cdouble(1)-2*b_axial*X0+2*a*X0*X0)) *
      (cdouble(1) + cerf_check(b_axial/2./sqrt(a) - sqrt(a) * X0));

  const cdouble cerfc_coeff_four = b_axial/2./sqrt(a) + sqrt(a)*X0;
  cdouble coeff_four;

  // erfc(x) ~ exp(-x^2)
  const cdouble series_log_coeff =
      -cerfc_coeff_four*cerfc_coeff_four + log_coeff_three;
  if (real(cerfc_coeff_four) < sqrt(-MIN_LOG) &&
      real(log_coeff_three) > MIN_LOG) {
      coeff_four = cerfc(cerfc_coeff_four) * exp(log_coeff_three);
  } else {
    coeff_four = exp_check(series_log_coeff)*(
          1/sqrt(M_PI)/cerfc_coeff_four - 0.5/sqrt(M_PI)/
          (cerfc_coeff_four*cerfc_coeff_four*cerfc_coeff_four));
  }
  const cdouble part_four =
      sqrt(M_PI)*(b_axial*b_axial + 2*a*(cdouble(1)+2*X0*(b_axial+a*X0)));
  const cdouble res = general_coeff * (
      2*sqrt(a)*(part_one + part_two) +
      exp(log_coeff_three) * part_three +
          part_four * coeff_four);
  return res;
}

void BoxHamiltonian::get_axis_integral_derivative(
    const cdouble& a,
    const cVector_3& b,
    const cdouble& log_norm,
    const unsigned char axis,
    cdouble* a_diff,
    cVector_3* b_diff) const {
  cdouble b_axial, b_perp1, b_perp2;
  double X0, k;
  set_axis(b, axis, &b_axial, &b_perp1, &b_perp2, &X0, &k);
  const cdouble b_perp_norm = b_perp1*b_perp1 + b_perp2*b_perp2;
  const cdouble b_norm = b.norm2();
  // a_diff
  const cdouble log_coeff_one_two = b_norm / 4. / a + log_norm;
  const cdouble erf_minus = cerf_check((b_axial - 2*a*X0)/(2.*sqrt(a)));
  const cdouble erf_plus = cerf_check((b_axial + 2*a*X0)/(2.*sqrt(a)));
  const cdouble part_one =
      - ((b_axial*b_axial)*b_norm + 24*(a*a*a)*(X0*X0) + 2*a*
      (8*(b_axial*b_axial) + (b_perp1*b_perp1) + (b_perp2*b_perp2) -
       2*b_axial*(b_norm)*X0) + 4*(a*a)*(cdouble(5) - 10*b_axial*X0 + b_norm*(X0*X0)))*
      erf_minus;
  const cdouble part_two =
      (-2*(20*(a*a) + (b_axial*b_axial)*b_norm +
           2*a*(8*(b_axial*b_axial) + b_perp_norm)) -
       8*(a*a)*(6*a + b_norm)*(X0*X0) +
       ((b_axial*b_axial)*(b_norm) + 24*(a*a*a)*(X0*X0) +
        2*a*(8*(b_axial*b_axial) + b_perp_norm + 2*b_axial*b_norm*X0) +
        4*(a*a)*(cdouble(5) + 10*b_axial*X0 + b_norm*(X0*X0))) * erf_plus);
  const cdouble log_coeff_three = b_perp_norm / 4. / a - a*X0*X0 + log_norm;
  const cdouble htrig_coeff = b_axial*X0;
  const cdouble mod_htrig_coeff_plus = htrig_coeff+log_coeff_three;
  const cdouble mod_htrig_coeff_minus = -htrig_coeff+log_coeff_three;
  const cdouble b_2cosh = exp_check(mod_htrig_coeff_plus) +
      exp_check(mod_htrig_coeff_minus);
  const cdouble b_sinh = 0.5*(exp_check(mod_htrig_coeff_plus) -
                              exp_check(mod_htrig_coeff_minus));

  const cdouble part_three = a*(6*a + b_norm)*X0*b_2cosh -
      b_axial*(14*a + b_norm)*b_sinh;

  *a_diff = M_PI / 32. * pow(a, -5.5) * (
      sqrt(M_PI) * exp(log_coeff_one_two) * (part_one + part_two) +
      4 * sqrt(a) * part_three);

  // b_diff[axis]
  const cdouble b_part_one =
      2*b_axial*(6*a + (b_axial*b_axial) + 4*(a*a)*(X0*X0)) +
      ((b_axial*b_axial*b_axial) + 2*a*b_axial*(cdouble(3) - 2*b_axial*X0) +
       4*(a*a)*X0*(cdouble(-2) + b_axial*X0))*erf_minus -
      ((b_axial*b_axial*b_axial) + 4*(a*a)*X0*(cdouble(2) + b_axial*X0) +
       2*a*b_axial*(cdouble(3) + 2*b_axial*X0))*erf_plus;

  const cdouble b_part_two =
      a*b_axial*X0*b_2cosh - (4*a + (b_axial*b_axial))*b_sinh;

  (*b_diff)[axis] = M_PI * pow(a, -4.5) / 16. * (
      exp(log_coeff_one_two) * sqrt(M_PI) * b_part_one -
      4*sqrt(a)*b_part_two);

  // b_diff[perp]
  const cdouble b_perp_part_one = b_axial + 2*a*X0;
  const cdouble b_perp_log_one =
      b_perp_norm / 4. / a - X0*(b_axial + a * X0) + log_norm;
  const cdouble b_perp_part_two = b_axial - 2*a*X0;
  const cdouble b_perp_log_two =
      b_perp_norm / 4. / a + X0*(b_axial - a * X0) + log_norm;

  const cdouble b_perp_part_three =
      2*(2*a + (b_axial*b_axial) + 4*(a*a)*(X0*X0)) +
      ((b_axial*b_axial) + 2*a*(cdouble(1) - 2*b_axial*X0 + 2*a*(X0*X0)))*erf_minus -
      ((b_axial*b_axial) + 2*a*(cdouble(1) + 2*X0*(b_axial + a*X0)))*erf_plus;

  const cdouble deriv_b_perp_common =
      M_PI * pow(a, -4.5) / 16. * (
          -2*sqrt(a)*(exp_check(b_perp_log_one + log(b_perp_part_one)) -
                      exp_check(b_perp_log_two + log(b_perp_part_two))) +
          sqrt(M_PI)*exp_check(log_coeff_one_two + log(b_perp_part_three)));

  (*b_diff)[(axis + 1) % 3] = b_perp1*deriv_b_perp_common;
  (*b_diff)[(axis + 2) % 3] = b_perp2*deriv_b_perp_common;

  *a_diff *= k;
  *b_diff *= k;
}

void BoxHamiltonian::get_derivatives(
    const cdouble& a1,
    const cVector_3& b1,
    const cdouble& a2,
    const cVector_3& b2,
    cdouble* integral,
    cdouble* a1_real_diff,
    cdouble* a1_imag_diff,
    cVector_3* b1_real_diff,
    cVector_3* b1_imag_diff,
    cdouble* a2_real_diff,
    cdouble* a2_imag_diff,
    cVector_3* b2_real_diff,
    cVector_3* b2_imag_diff) const {
  const cdouble a = conj(a1) + a2;
  const cVector_3 b = conj(b1) + b2;

  const cdouble log_norm = log_double_norm(a1, b1, a2, b2);

  *integral = 0;
  *a1_real_diff = 0;
  *a1_imag_diff = 0;
  *b1_real_diff = 0;
  *b1_imag_diff = 0;
  *a2_real_diff = 0;
  *a2_imag_diff = 0;
  *b2_real_diff = 0;
  *b2_imag_diff = 0;

  for (unsigned char axis = 0; axis < 3; ++axis) {
    const cdouble normed_axis_integral = axis_integral(a, b, log_norm, axis);
    *integral += normed_axis_integral;

    cdouble integral_a_derivative;
    cVector_3 integral_b_derivative;
    get_axis_integral_derivative(
        a, b, log_norm, axis, &integral_a_derivative, &integral_b_derivative);
    // We use the fact, that a=a1+a2, thus dJ/da = dJ/da1
    // Also, unnormed_integral is a regular function, so integral_a_derivative
    // is a complex derivative
    *a1_real_diff += norm_deriv_a_real(conj(a1), conj(b1))*normed_axis_integral+
        integral_a_derivative;
    *a1_imag_diff -= norm_deriv_a_imag(conj(a1), conj(b1))*normed_axis_integral+
        I*integral_a_derivative;
    *a2_real_diff +=
        norm_deriv_a_real(a2, b2)*normed_axis_integral + integral_a_derivative;
    *a2_imag_diff += norm_deriv_a_imag(a2, b2)*normed_axis_integral +
        I*integral_a_derivative;

    *b1_real_diff += norm_deriv_b_real(conj(a1), conj(b1))*normed_axis_integral+
        integral_b_derivative;
    *b1_imag_diff -= norm_deriv_b_imag(conj(a1), conj(b1))*normed_axis_integral+
        I*integral_b_derivative;
    *b2_real_diff += norm_deriv_b_real(a2, b2)*normed_axis_integral +
        integral_b_derivative;
    *b2_imag_diff += norm_deriv_b_imag(a2, b2)*normed_axis_integral +
        I*integral_b_derivative;
  }
}


double factorial(int n) {
  assert(n >= 0);
  double res = 1.;
  for (int i = 2; i <= n; ++i)
    res *= i;
  return res;
}


HermitePolinomial2::HermitePolinomial2(int max_order) {
  max_n = max_order;

  if (max_n <= 1)
    return;

  int max_n2 = max_n / 2;

  // Save prefactors for n and m into the table of the size n x (n/2+1)
  m_size = max_n2 + 1;
  coef.resize(max_n*m_size);

  for (int n = 1; n <= max_n; ++n)
    for (int m = 0; m <= n/2; ++m)
      coef[(n-1)*m_size + m] = factorial(n) / (factorial(m)*factorial(n - 2*m));
}


template<class T>
T HermitePolinomial2::hermite(int n, T x, T y) {
  assert(n >= 0 && (n ==0 || (n!=0 && n <= max_n)));

  if (n == 0)
    return 1.;

  double* coef_n = &coef[(n-1)*m_size];
  T y_m = 1.;
  T xx = x*x;
  T res = coef_n[0] * pow(x, n);    // m = 1

  for (int m = 1; m <= n/2; ++m) {
    y_m *= y;                   // y_m = y^m
    res += coef_n[m] * pow(x, n-2*m) * y_m;
  }

  return res;
}


void BoxHamiltonian::calc_proj_prefactors() {
  if(max_n <=0) return;

  for(int j=0; j<3; j++) {
    proj_prefactor[j].resize(max_n + 1);
    alpha_k[j] = pow(2.*potential_multiplers[j]*m_electron/h_sq,1./4);

    for(int n=0; n<=max_n; ++n)
      proj_prefactor[j][n]
        = sqrt( alpha_k[j]*sqrt(M_PI) / (pow(2.,n)*factorial(n)) );
  }
}

cdouble BoxHamiltonian::eigen_proj(int n, const WavePacket& wp, int axis) {
  if(max_n <=0) 
    return 0.;
  cdouble a = conj(wp.a);
  cdouble b = conj(wp.b[axis]);
  double r = wp.get_r()[axis];
  double r2 = r*r;
  cdouble c = cdouble( log(real(a)*(2./M_PI))/4 - r2*real(a),
    - r2*imag(wp.a) + r*imag(wp.b[axis]) );  // c complex conjugate in 1D
  double alpha = alpha_k[axis];
  cdouble d = a + alpha*alpha/2;

  return proj_prefactor[axis][n] / sqrt(d)
    * exp(c + b*b / (4*d))
    * hpol.hermite(n, b*alpha/d, alpha*alpha/d-1.);
}
