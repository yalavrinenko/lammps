/* -*- c++ -*- ----------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 Contributing authors: Andres Jaramillo-Botero, Hai Xiao, Julius Su (Caltech)
------------------------------------------------------------------------- */

#include <cmath>

namespace LAMMPS_NS {

class Constants{
public:
  static inline constexpr double sqrt_pi(){
    return 1.772453850905516027298167483341145182797549456122387128213;
  }
};

inline double ierfoverx1(double x, double *df)
{
  // Computes Erf(x)/x and its first derivative
  auto erf_v = std::erf(x);
  auto f =  erf_v / x;
  *df = (2.0 * std::exp(-x * x) * x / Constants::sqrt_pi()  - erf_v) / (x * x);
  return f;
}

/* ---------------------------------------------------------------------- */

inline void KinElec(double radius, double *eke, double *frc)
{
  *eke += 1.5 / (radius * radius);
  *frc += 3.0 / (radius * radius * radius);
}

/* ---------------------------------------------------------------------- */

inline void ElecNucNuc(double q, double rc, double *ecoul, double *frc)
{
  *ecoul += q / rc;
  *frc += q / (rc * rc);
}

/* ---------------------------------------------------------------------- */

inline void ElecNucElec(double q, double rc, double re1,
                        double *ecoul, double *frc, double *fre1)
{
  double a, arc;
  double coeff_a;

  /* sqrt(2) */
  coeff_a = 1.4142135623730951;

  /* E = -Z/r Erf(a r / re) */
  /* constants: sqrt(2), 2 / sqrt(pi) */
  a = coeff_a / re1;
  arc = a * rc;

  /* Interaction between nuclear point charge and Gaussian electron */
  double E, dEdr, dEdr1, f, df;

  f = ierfoverx1(arc, &df);
  dEdr = q * a * a * df;
  dEdr1 = -q * (a / re1) * (f + arc * df);
  E = -q * a * f;

  *ecoul += E;
  *frc += dEdr;
  *fre1 += dEdr1;
}

/* ---------------------------------------------------------------------- */

inline void ElecElecElec(double rc, double re1, double re2,
                         double *ecoul, double *frc, double *fre1,
                         double *fre2)
{
  double a, arc, re, fre;
  double coeff_a;

  /* sqrt(2) */
  coeff_a = 1.4142135623730951;

  re = sqrt(re1 * re1 + re2 * re2);

  /* constants: sqrt(2), 2 / sqrt(pi) */
  a = coeff_a / re;
  arc = a * rc;

  /* V_elecelec = E * F                              */
  /* where E = -Z/r Erf(a r / re)                    */
  /*       F = (1 - (b s + c s^2) exp(-d s^2))       */
  /* and s = r / re                                  */

  double E, dEdr, dEdr1, dEdr2, f, df;

  f = ierfoverx1(arc, &df);
  dEdr = -a * a * df;
  fre = a * (f + arc * df) / (re * re);
  dEdr1 = fre * re1;
  dEdr2 = fre * re2;

  E = a * f;

  *ecoul += E;
  *frc += dEdr;
  *fre1 += dEdr1;
  *fre2 += dEdr2;
}

/* ---------------------------------------------------------------------- */

inline void ElecCoreNuc(double q, double rc, double re1, double *ecoul, double *frc)
{
  double a, arc;
  double coeff_a;
  double E, dEdr, df, f;

  coeff_a = 1.4142135623730951;    /* sqrt(2) */
  a = coeff_a / re1;
  arc = a * rc;

  f = ierfoverx1(arc, &df);
  dEdr = -q * a * a * df;
  E = q * a * f;

  *ecoul += E;
  *frc += dEdr;
}

/* ---------------------------------------------------------------------- */

inline void ElecCoreCore(double q, double rc, double re1, double re2,
  double *ecoul, double *frc)
{
  double a, arc, re;
  double coeff_a;
  double E, dEdr, f, df;

  coeff_a = 1.4142135623730951;

  re = sqrt(re1 * re1 + re2 * re2);
  a = coeff_a / re;
  arc = a * rc;

  f = ierfoverx1(arc, &df);
  dEdr = -q * a * a * df;
  E = q * a * f;

  *ecoul += E;
  *frc += dEdr;
}

/* ---------------------------------------------------------------------- */

inline void ElecCoreElec(double q, double rc, double re1, double re2,
        double *ecoul, double *frc, double *fre2)
{
  double a,arc, re;
  double coeff_a;
  double E, dEdr, dEdr2, f, df, fre;

  /* sqrt(2) */
  coeff_a = 1.4142135623730951;

  /*
  re1: core size
  re2: electron size
  re3: size of the core, obtained from the electron density function rho(r) of core
  e.g. rho(r) = a1*exp(-((r)/b1)^2), a1 =157.9, b1 = 0.1441 -> re3 = 0.1441 for Si4+
  */

  re = sqrt(re1 * re1 + re2 * re2);

  a = coeff_a / re;
  arc = a * rc;

  f = ierfoverx1(arc, &df);
  E = -q * a * f;
  dEdr = -q * a * df * a;
  fre = q * a * (f + arc * df) / (re * re);
  dEdr2 = fre * re2;

  *ecoul += E;
  *frc -= dEdr;
  *fre2 -= dEdr2;
}

inline void RForce(double dx, double dy, double dz,
                   double rc, double force, double *fx, double *fy, double *fz)
{
  force /= rc;
  *fx = force * dx;
  *fy = force * dy;
  *fz = force * dz;
}

/* ---------------------------------------------------------------------- */

inline void SmallRForce(double dx, double dy, double dz,
                        double rc, double force,
                        double *fx, double *fy, double *fz)
{
  /* Handles case where rc is small to avoid division by zero */

  if (rc > 0.000001){
    force /= rc;
    *fx = force * dx; *fy = force * dy; *fz = force * dz;
  } else {
    if (dx != 0) *fx = force / sqrt(1 + (dy * dy + dz * dz) / (dx * dx));
    else *fx = 0.0;
    if (dy != 0) *fy = force / sqrt(1 + (dx * dx + dz * dz) / (dy * dy));
    else *fy = 0.0;
    if (dz != 0) *fz = force / sqrt(1 + (dx * dx + dy * dy) / (dz * dz));
    else *fz = 0.0;
    //                if (dx < 0) *fx = -*fx;
    //                if (dy < 0) *fy = -*fy;
    //                if (dz < 0) *fz = -*fz;
  }
}

/* ---------------------------------------------------------------------- */

inline double cutoff(double x)
{
  /*  cubic: return x * x * (2.0 * x - 3.0) + 1.0; */
  /*  quintic: return -6 * pow(x, 5) + 15 * pow(x, 4) - 10 * pow(x, 3) + 1; */

  /* Seventh order spline */
  //      return 20 * pow(x, 7) - 70 * pow(x, 6) + 84 * pow(x, 5) - 35 * pow(x, 4) + 1;
  return (((20 * x - 70) * x + 84) * x - 35) * x * x * x * x + 1;
}

/* ---------------------------------------------------------------------- */

inline double dcutoff(double x)
{
  /*  cubic: return (6.0 * x * x - 6.0 * x); */
  /*  quintic: return -30 * pow(x, 4) + 60 * pow(x, 3) - 30 * pow(x, 2); */

  /* Seventh order spline */
  //      return 140 * pow(x, 6) - 420 * pow(x, 5) + 420 * pow(x, 4) - 140 * pow(x, 3);
  return (((140 * x - 420) * x + 420) * x - 140) * x * x * x;
}

}
