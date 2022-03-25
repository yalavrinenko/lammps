/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *****************************************************************************/

/*s****************************************************************************
 * $Log: wpmd.h,v $
 * Revision 1.77  2015/12/04 19:16:06  valuev
 * fixed exchange energy measurement
 *
 * Revision 1.76  2015/12/03 15:02:43  valuev
 * added exchange energies
 *
 * Revision 1.75  2015/11/10 16:29:45  valuev
 * added ion charge
 *
 * Revision 1.74  2015/10/20 15:01:13  valuev
 * switched to width velocity
 *
 * Revision 1.73  2015/02/10 14:08:09  valuev
 * added step function
 *
 * Revision 1.72  2014/07/16 10:32:11  valuev
 * quantum sums for AWPMC
 *
 * Revision 1.71  2014/07/11 18:14:18  valuev
 * added overlap matrix log
 *
 * Revision 1.70  2014/07/08 17:07:21  valuev
 * prepared confined MC simulation with AWPMD nonsplit norm matrix
 *
 * Revision 1.69  2014/07/03 16:59:13  valuev
 * started adapting for confined MC
 *
 * Revision 1.68  2014/06/25 10:53:15  valuev
 * preparing for confined simulations
 *
 * Revision 1.67  2014/04/01 14:07:35  valuev
 * prepared variables for boundary potential
 *
 * Revision 1.66  2013/09/16 22:38:20  kazeev
 * Ported for MKL11, icpc 14
 *
 * Revision 1.65  2013/08/14 02:55:48  valuev
 * more fixes
 *
 * Revision 1.64  2013/07/09 09:22:02  kazeev
 * Fixed compilation for K100 update.
 *
 * Revision 1.63  2013/03/23 09:40:34  valuev
 * new algorithm and fixes for degeneracy constraints
 *
 * Revision 1.62  2013/03/04 16:46:17  valuev
 * added rebuild_wp interface
 *
 * Revision 1.61  2012/09/20 05:38:31  valuev
 * optimizer (test version)
 *
 * Revision 1.60  2012/09/03 15:54:46  valuev
 * fixed sign in WPMD forces
 *
 * Revision 1.59  2012/08/22 11:57:51  valuev
 * changes for KIM compliance
 *
 * Revision 1.58  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.57  2012/04/26 09:21:35  valuev
 * added Yukawa potential (tested)
 *
 * Revision 1.56  2012/04/18 10:36:36  valuev
 * optimizer implementation (not working yet)
 *
 * Revision 1.55  2012/04/17 20:42:31  valuev
 * fixed Psi calculation, started adding e-i Debye screening
 *
 * Revision 1.54  2012/04/03 10:26:41  valuev
 * added variable constraints
 *
 * Revision 1.53  2012/01/23 12:36:17  valuev
 * added wave function overlap calculation
 *
 * Revision 1.52  2011/10/12 16:39:02  valuev
 * restructured  interaction function
 *
 * Revision 1.51  2011/10/08 11:21:14  valuev
 * added external force storage and external force power calculation
 *
 * Revision 1.50  2011/09/18 19:47:18  valuev
 * fixed external force calculation
 *
 * Revision 1.49  2011/07/24 08:56:46  valuev
 * fixed transformation
 *
 * Revision 1.48  2011/07/23 05:06:59  valuev
 * direct inversion of norm matrix
 *
 * Revision 1.47  2011/07/22 19:54:14  valuev
 * fixes for norm matrix
 *
 * Revision 1.46  2011/07/18 16:19:10  valuev
 * divergence test
 *
 * Revision 1.45  2011/06/22 08:31:38  valuev
 * derivative fixes
 *
 * Revision 1.4  2011/06/11 16:53:55  valuev
 * sync with LAMMPS
 *
 * Revision 1.3  2011/06/11 16:45:23  morozov
 * Fixed erf.c for Windows and Unix
 *
 * Revision 1.2  2011/06/10 19:25:17  morozov
 * *** empty log message ***
 *
 * Revision 1.43  2011/06/10 19:20:53  valuev
 * fixes
 *
 * Revision 1.42  2011/06/09 22:55:08  valuev
 * norm matrices
 *
 * Revision 1.41  2011/06/07 19:58:42  valuev
 * corrected partitions
 *
 * Revision 1.40  2011/06/07 17:43:00  valuev
 * added Y derivatives
 *
 * Revision 1.39  2011/06/03 08:13:33  valuev
 * added partitions to account for ghost atoms
 *
 * Revision 1.38  2011/06/02 22:11:17  morozov
 * Compatibility with LAPACK library
 *
 * Revision 1.37  2011/06/01 23:45:35  valuev
 * modified for LAMMPS compatibility
 *
 * Revision 1.36  2011/05/28 17:16:22  valuev
 * fixed template<>, some fixes to UHF
 *
 * Revision 1.35  2011/05/25 21:30:48  morozov
 * Compatibility with ICC 11.1
 *
 * Revision 1.34  2011/05/25 05:23:43  valuev
 * fixed variable transformation for norm matrix
 *
 * Revision 1.33  2011/05/24 19:54:32  valuev
 * fixed sqmatrix::iterator
 *
 * Revision 1.32  2011/05/21 23:06:49  valuev
 * Norm matrix transform to pysical variables
 *
 * Revision 1.31  2011/05/20 21:39:49  valuev
 * separated norm calculation
 *
 * Revision 1.30  2011/05/14 18:56:19  valuev
 * derivative for ee split interactions
 *
 * Revision 1.29  2011/05/05 08:56:02  valuev
 * working split wp version
 *
 * Revision 1.28  2011/05/04 16:48:52  valuev
 * fixed syntax
 *
 * Revision 1.27  2011/05/04 09:04:48  valuev
 * completed wp_split (except for ee forces)
 *
 * Revision 1.26  2011/04/29 03:07:20  valuev
 * new split wp features
 *
 * Revision 1.25  2011/04/22 09:52:49  valuev
 * working on split WP
 *
 * Revision 1.24  2011/04/20 08:43:09  valuev
 * started adding split packet WPMD
 *
 * Revision 1.23  2010/09/03 12:17:48  morozov
 * The order of parameters in Norm matrix is changed to mimic the order of the single WP parameter storage
 *
 * Revision 1.22  2009/08/27 00:01:36  morozov
 * First working MPI equilibration
 *
 * Revision 1.21  2009/04/14 14:44:10  valuev
 * fixed momentum calculation in hartree model, added "fix" constraint and model="hartree"
 *  to parameters
 *
 * Revision 1.20  2009/04/13 17:00:45  morozov
 * Fixed norm-matrix ratio in AWPMC algorithm
 *
 * Revision 1.19  2009/04/06 17:00:28  morozov
 * Fixed Hartree version of WPMC
 *
 * Revision 1.18  2009/04/01 10:06:37  valuev
 * added Hartee factorization to AWPMD
 *
 * Revision 1.17  2009/03/24 16:10:05  morozov
 * Fixed errors in Norm-matrix calculation related to PBC
 *
 * Revision 1.16  2009/03/17 11:40:04  morozov
 * The prefactor of NormMatrix is corrected
 *
 * Revision 1.15  2008/07/23 16:42:12  valuev
 * Added AWPMD Monte-Carlo
 *
 * Revision 1.14  2008/07/23 15:58:32  valuev
 * *** empty log message ***
 *
 * Revision 1.13  2008/07/21 02:23:22  morozov
 * *** empty log message ***
 *
 * Revision 1.12  2008/07/18 18:15:31  morozov
 * *** empty log message ***
 *
 * Revision 1.11  2008/05/29 13:33:05  valuev
 * VASP band structure
 *
 * Revision 1.10  2008/05/14 17:17:26  morozov
 * Passed 2- and 3-electron test. Added Norm matrix.
 *
 * Revision 1.9  2008/05/05 17:29:32  morozov
 * Fixed errors with Hermitian matrix indeces. Redesigned cVector_3 operations.
 *
 * Revision 1.8  2008/04/28 22:16:45  valuev
 * restructured coulomb term
 *
 * Revision 1.7  2008/04/28 09:54:13  valuev
 * corrected summation for Eee part
 *
*******************************************************************************/
# ifndef WPMD_H
# define WPMD_H

/** @file wpmd.h
    @brief Classes for Wave Packet Molecular Dynamics of two component plasma. */


# define AWPMD_VERSION 0.92

# ifndef _USE_MATH_DEFINES
# define _USE_MATH_DEFINES
# endif
# include <complex>
# include <vector>
# include <cmath>
# include "logexc.h"
# include "cvector_3.h"
# include "pairhash.h"
# include "tcpdefs.h"
# include "wavepacket.h"
# include "erf.h"
# include "cerf.h"
# include "math_utils.h"

# include "gradopt.h"


# ifndef DFT_EXTENSION
# define DFT_EXTENSION 0
# endif


# if DFT_EXTENSION
# include "awpmd-dft.hpp"
# include <LDA.hpp>
# include  <ModLDA.hpp>
#include <DataTypes.hpp>
#include "LSDA.hpp"
#include <awpmd-dft-cpu.hpp>
# endif

using namespace std;
#define MKL_Complex16 cdouble
// http://software.intel.com/sites/products/documentation/hpc/mkl/lin/MKL_UG_coding_calls/Using_Complex_Types_in_C_C.htm
# include "lapack_inter.h"

#include "box_hamiltonian.h"

typedef hmatrix<cdouble> chmatrix;

const cdouble i_unit = cdouble(0.,1.), i_unit1 = -i_unit;
const double h_plank2 = h_plank * 2.;



//cdouble ccerf(const cdouble &a);

# if 0
inline cdouble cerf_div(const cdouble &z, const cdouble &c=i_unit){
  cdouble zc = z*c;
  if (zc.imag() > std::numeric_limits<double>::epsilon()){
    if((fabs(real(zc))+fabs(imag(zc)))<1)
      return cerf_div_series(zc)*c;
    else
      return cerf(z*c)/z;
  }

  return std::erf(zc.real() + std::numeric_limits<double>::epsilon()) / (z.real() + std::numeric_limits<double>::epsilon());
}
# endif
inline cdouble cerf_div(const cdouble &z, const cdouble &c = i_unit) {
  cdouble zc = z * c;
  if (zc.imag() > std::numeric_limits<double>::epsilon())
    return cerf(z*c) / z;

  if ((fabs(real(zc)) + fabs(imag(zc))) < 1)
    return cerf_div_series(zc)*c;
  else
    return cerf(z*c) / z;
}



/*
///\en Calculates cerf(c*z)/z
inline cdouble cerf_div(const cdouble &z, const cdouble &c=i_unit){
  if((fabs(real(z))+fabs(imag(z)))<1e-9)
    return c*two_over_sqr_pi;
  else
    return cerf(z*c)/z;
}*/

/*
///\en Calculates 0.5*[cerf(c*z+k)+cerf(cz-k)]/z
inline cdouble cerf_split_div(const cdouble &z, const cdouble &k, const cdouble &c=i_unit){
  if((fabs(real(z))+fabs(imag(z)))<1e-8)
    return c*two_over_sqr_pi*exp(-k*k);
  else
    return 0.5*(cerf(z*c+k)+cerf(z*c-k))/z;
}

///\en Calculates sinh(cz)/z
inline cdouble sinh_div(const cdouble &z, const cdouble &c=i_unit){
  if((fabs(real(z))+fabs(imag(z)))<1e-8)
    return c;
  else
    return sinh(c*z)/z;
}
*/


//e calculates erf(c*z)/z
inline double erf_div(const double &z, double c=1){
  if(fabs(z*c)<1e-9)
    return cerf_div_series(z*c)*c;
  else
    return erf(z*c)/z;
}


/*inline double erf_div(const double &z, double c=1){
  if(fabs(z)<1e-9)
    return c*two_over_sqr_pi;
  else
    return erf(z*c)/z;
}*/

template<class T1, class T2>
struct _myrefpair{
  T1& first;
  T2& second;
  _myrefpair(T1& a, T2 &b):first(a),second(b){}
  _myrefpair operator=(const pair<T1,T2> &other){
    first=other.first;
    second=other.second;
    return *this;
  }
};

template< class T1, class T2>
_myrefpair<T1,T2> _mytie(T1& var1, T2 &var2){
  return _myrefpair<T1,T2>(var1, var2);
}

inline pair<double,double> operator*(const pair<double,double> &right, double left){
  return make_pair(right.first*left,right.second*left);
}

// Auxilary class to handle the normalizing term derivatives
class NormDeriv
{
public:
  cdouble l;    // lambda = (f over a_re)
  double m;     // mu = (f over a_im) / i
  cVector_3 u;  // u = (f over b_re)
  Vector_3 v;   // v = (f over b_im) / i

  NormDeriv() {}
  NormDeriv(const WavePacket& wp) { set(wp); }

  //e Create NormDeriv object and calculate the derivatived for the given WP
  void set(const WavePacket& wp){
    Vector_3 br = real(wp.b), bi = imag(wp.b);
    double ar = real(wp.a), ai = imag(wp.a);
    double i_2ar = 0.5 / ar, ai_ar = ai / ar;

    v = (-i_2ar) * br;
    m = v.norm2();
    u = v * (i_unit1 * ai_ar) - wp.b * i_2ar;
    l = (1.5*i_2ar + m) + cdouble(0.,2.) * ( (br*bi)*i_2ar*i_2ar - m*ai_ar );
  }
};

inline NormDeriv conj(const NormDeriv& src){
  NormDeriv dst;
  dst.l = conj(src.l);
  dst.m = -src.m;
  dst.u = conj(src.u);
  dst.v = - src.v;
  return dst;
}

///\en Auxilary class to handle derivatives of overlaps
class OverlapDeriv{
public:
  WavePacket w1, w2, w12;
  NormDeriv d1, d2;
  cdouble I0, I2;
  cVector_3 I1;
  cdouble bb_4a;
  sqmatrix<cdouble> IDD;


  OverlapDeriv():I0(0),I1(0),IDD(10){}

  void set1(const WavePacket& w1_) {
    w1=w1_;
    d1.set(w1);
    d1=conj(d1);
  }

  //e Create NormDeriv object and calculate the derivatived for the given WP
  void set2(const WavePacket& w2_, const cdouble *I0_=NULL){
    w2=w2_;
    d2.set(w2);
    w12=conj(w1)*w2;
    if(!I0_)
      I0 = w12.integral();
    else
      I0=*I0_;
    I1 = w12.b * (I0 / w12.a / 2);
    bb_4a = w12.b.norm2() / w12.a / 4;
    I2 = I0 * (bb_4a + 1.5) / w12.a;
  }

  cdouble da2_re() const {
    return (d2.l*I0 - I2);
  }

  cdouble da2_im() const {
    return i_unit*(d2.m*I0 - I2);
  }

  cdouble da1_re() const {
    return (d1.l*I0 - I2);
  }

  cdouble da1_im() const {
    return -i_unit*(-d1.m*I0 - I2); //!
  }

  cdouble db2_re(int i) const {
    return d2.u[i]*I0 + I1[i];
  }

  cdouble db2_im(int i) const {
    return i_unit*(d2.v[i]*I0 + I1[i]);
  }

  cdouble db1_re(int i) const {
    return d1.u[i]*I0 + I1[i];
  }

  cdouble db1_im(int i) const {
    return -i_unit*(-d1.v[i]*I0 + I1[i]); //!
  }

  ///\en Calculates  derivative overlap matrix IDD
  void calc_der_overlap(bool self=false, cdouble cc1=0., cdouble c2=0.);


};


class AWPMD {
//protected:
public:
  int ne[2], ni;
  int nwp[2]; ///<\en number of wavepackets (same as ne for unsplit version)
  int nvar[2]; ///<\en full number of dynamic variables for each spin
  int vars_per_wp; ///<\en number of variables per WP, =8 for WPMD, =10 for split WPMD
  chmatrix O[2], Y[2], Te[2], Tei[2];
  smatrix<unsigned char> Oflg[2];  ///<\en equals 0 for non-overlaping packets
  sqmatrix<double> Norm[2];  ///<\en Norm matrix
  vector<WavePacket> wp[2]; ///<\en wave packets for electrons (spin up =0 and down =1)
  vector<double> qe[2];  ///<\en electron charges
  vector<double> qi;  ///<\en ion charges
  vector<Vector_3> xi;  ///<\en ion coordinates
  int pbc; ///<\en pbc flag
  Vector_3 cell; ///<\en cell coordinates (L1,L2,L3)
	Vector_3 center; ///<\en center coordinates, default for PBC-- (L1/2,L2/2,L3/2), otherwise (0,0,0)
  double Lextra; ///<\en width PBC, unset if negative
  double harm_w0_4;
  double w0;
  int calc_ii; ///<\en flag indicating whether to calculate ion-ion interaction
  int calc_ee; ///<\en flag indicating whether to calculate electron-electron interaction
	int calc_ei; ///<\en flag indicating whether to calculate electron-ion interaction
  int norm_needed; ///<\en flag indicating whether to prepare norm matrix data in interaction function

  int screened_ei; ///<\en 1 = use Debye screening in e-i interaction, 0= unscreened
  double kappa_ei; ///<\en Screening length for e-i interaction

  enum {NORM_UNDEFINED, NORM_CALCULATED, NORM_FACTORIZED, NORM_INVERTED};
  int norm_matrix_state[2];
  ///\en bitflags to control the progress of object calculation state 
  enum { CALC_UNDEF=0,CALC_ELECTRONS_SET=0x1, CALC_IONS_SET = 0x2, CALC_NORMS=0x4, CALC_OVERLAPS=0x8, CALC_LMN=0x10, CALC_ENERGIES=0x20 };
  int calc_state;

  // Arrays for temporal data
  chmatrix IDD;  // Second derivatives of the overlap integral (used in Norm matrix)
  vector<cdouble> ID, IDYs;  // First derivatives of the overlap integral (used in Norm matrix)
  vector<int> ipiv;  // The pivot indices array

  recmatrix<cdouble> L[2]; ///<\en overlap derivative matrix for each spin
  recmatrix<cdouble> M[2]; ///<\en 'normalized' overlap derivative matrix for each spin: M=YL

  map<int,int> fixed[2]; ///<\en map (var_id -> flag) listing fixed variables for each spin

	int dft_extension; ///\en< 0 = don't use dft, 1 = add XC energy, 2= add C energy only (not implemented, for AWPMD)
	size_t dft_meshx, dft_meshy, dft_meshz;

  int numprocs;
public:

# if DFT_EXTENSION
    DFTConfig dft_conf;
    XCEnergy *ex;
    double Ee_dft = 0;
# endif

	bool use_box;
  BoxHamiltonian box;



  enum {NONE=0, HARM, FIX, RELAX} constraint;

  ///\em Sets approximation level for quantum problem: \n
  ///    HARTREE Hartree product (no antisymmetrization) \n
  ///    DPRODUCT product of det0*det1 of antisymmetrized functions for spins 0, 1 \n
	///    UHF unrestricted Hartree-Fock (apparently same as DPRODUCT)
  enum APPROX {HARTREE = 0, DPRODUCT = 2, UHF = 1} approx;
  ///\em Sets overlap matrix element to zero if the overlap norm is less than this value
  double ovl_tolerance;
  ///\em Minimal allowed difference between wavepcket overlap norm and 1, overlaps with less difference are normalized 
  double ovl_degeneracy_min;

  double Ee[2]; ///<\en electron kinetic energy for each spin
  double Eei[2]; ///<\en electron-ion energy for each spin
  double Eee, Ew; ///<\en electron-electron energy
  double Eii; ///<\en ion-ion energy
  double Ebord_ion; ///<\en energy due to boundary conditions (confinement border)
  double Eext; ///<\en energy due to external force
  double Ebord; ///<\en energy due to boundary conditions (confinement border)
  double Wext; ///<\en power of external force: Wext = \sum_i (dExt/dq_i)(dq_i/dt)
  double Edk; ///<\en sum of diagonal kinetic energy terms
  double Edc; ///<\en sum of diagonal Coulomb energy terms

  double Ee_exch, Eei_exch, Eee_exch, Eext_exch, Ebord_exch; ///<\en Exchange energies by type (calculated for UHF, Eee_exch stores also XC energy for DFT)
  double Eee_hartree; ///<\en Alternative Hartree ee energy (non-diverging), calculated for UHF
  
  Vector_3 Ext_force; ///<\en External force added to forces in interaction()

  vector<double> Eep[2]; ///<\en per particle electron kinetic energy for each spin
  vector<double> Eeip[2]; ///<\en per particle electron-ion energy for each spin
  vector<double> Eeep[2]; ///<\en per particle electron-electron energy for each spin
  vector<double> Ewp[2]; ///<\en per particle restrain energy for each spin
  vector<double> Eiep; ///<\en per particle ion-electron energy
  vector<double> Eiip; ///<\en per particle ion-ion energy


  vector<double> E_der[2]; ///<\en energy derivative with respect to all WP coordinates (first {a,b}, then changed to physical)
  vector<double> F_extra[2]; ///<\en energy derivative (force) due to external energy contribution: -dVext/dq, in physical representation
  vector<double> F_wall[2]; ///<\en energy derivative (force) due to confinement energy contribution: -dVwall/dq, in physical representation
  vector<double> dq_dt[2]; ///<\en time derivative of all WP variables as calculated using forces and norm matrix

	double norm_det_log[2]; ///<\en logarithm of the split norm-matrix determinant for each particle type
	double ovl_det_log[2]; ///<\en logarithm of the overlap matrix determinant for each particle type


  ///\en \{ Conversion constants that depend on the unit system used (for LAMMPS compatibility).
  ///       Default is GRIDMD units. Change them according to your unit system.
  double me; ///<\en electron mass (LAMMPS: me in the appropriate unit system)
  double one_h; ///<\en inverse of Plancks constant (LAMMPS: conversion [(m*v)/h] to [distance]  )
  double h2_me; ///<\en Plancks constant squared divided by electron mass (LAMMPS: conversion [h^2/(m*r^2)] to [Energy]  )
  double coul_pref; ///<\en Coulomb prefactor (e2 for GRIDMD) (LAMMPS: conversion [q^2/r] to [Energy]  )
  ///    \}
  double mvv2e; ///<\en convert mv^2/2 to energy

  ///\en 0 -- indicates that the inter-partition force should be full, and energy half,\n
  ///    1 -- inter-partition force and energy counts one half (LAMMPS compatibility)
  int newton_pair;

  ///\en How to calculate exchange energy contribution for UHF: 0 = from off diagonal elements, 1= as compared to Hartree in the same configuration 
  int ex_energy_type;

  //int myid; ///<\en id for partitions

  ///\en Partition arrays storing the tags of particles. The initial tags should be >0.
  ///    If the tag stored is <0, then the particle is ghost with -tag.
  ///    partition1[2] is for ions, 0, 1 for each electron spin
  vector<int> partition1[3];
  //vector<int> partition2[3]; ///<\en 2 for ions


  int tag_index(int i, int j) const {
    //return i==j ? -1 : (i>j ? (i-2)*(i-1)/2+j : (j-2)*(j-1)/2+i );
    return i == j ? -1 : (i > j ? (i - 1)*(i) / 2 + j : (j - 1)*(j) / 2 + i);
  }

  //std::map<int, std::map<int, bool>> interaction_tags;
  ///\en Gets ordering for organizing loops without double counting
  /// \return true if the pair should be skipped
  bool skip_by_tag_order(int s1, int ic1, int c1, int s2, int ic2, int c2) const {
    return c1<c2;
# if 0   
    int tag1 = partition1[s1][ic1];
    int tag2 = partition1[s2][ic2];

    //return (interactive_tags.count(std::abs(tag1)) && interactive_tags.at(std::abs(tag1)).count(std::abs(tag2)));
    bool res = true;

    if (tag1 <= 0 && tag2 <= 0) { // all at other partition
      res = tag1 != tag2;
    }
    else if (tag1 > 0 && tag2 > 0) // all at my partition
      res = (c1 < c2); // tag1 < tag2; // take care that tags can be omitted for sequential version
    else
      res = tag2 <= 0;

    //printf("Skip test: (%d,%d,%d)[%d]-(%d,%d,%d)[%d] %s\n", s2, ic2, c2, tag2,s1, ic1, c1, tag1,res ? "skip" : "keep");

    return res;
    //return std::abs(partition1[s][ic]);
# endif
  }

  ///\en 1 -- all my, -1 all other, 2 -- my mixed term, -2 -- other mixed term
  int check_ee(int s1, int icj1, int ick2, int s2, int icj3, int ick4) const {
    int tag[4];
    tag[0] = partition1[s1][icj1];
    tag[1] = partition1[s1][ick2];
    tag[2] = partition1[s2][icj3];
    tag[3] = partition1[s2][ick4];

    int imaxt = 0;
    bool has_neg = false;
    for (int i = 0; i < 4; i++) {
      if (tag[i] <= 0)
        has_neg = true;
      if (std::abs(tag[i]) > std::abs(tag[imaxt]))
        imaxt = i;
    }
    if (!has_neg)
      return 1;
    else if (tag[imaxt] > 0) // dislocated pair is taken by node with the largest own tag
      return 1;
    else
      return -1;
  }


  ///\en 1 -- all my, -1 all other, 2 -- my mixed term, -2 -- other mixed term
  int check_ee(int s1,int icj1,int s2, int ick2) const {
    int tag1 = partition1[s1][icj1];
    int tag2 = partition1[s2][ick2];

    int res = -1;
    if (tag1 <= 0) {
      if (tag2 > 0 && std::abs(tag2) > std::abs(tag1)) // dislocated pair is taken by node with the largest own tag
        return 1;
      else
        return -1;
    }
    else {  // tag1>0
      if (tag2 > 0) // all at this node
        return 1;
      // dislocated pair is taken by node with the largest own tag
      else if (std::abs(tag1) > std::abs(tag2)) // tag2<=0
        return 1;
      else
        return -1;
    }
# if 0


    if (tag1 <= 0 && tag2 <= 0) // all at other partition
      res=  -1;
    else if (tag1 > 0 && tag2 > 0) // all at my partition
      res = 1;
    else {
    
      int atag1 = std::abs(tag1);
      int atag2 = std::abs(tag2);
      int atag_min, atag_max, tag_min;
      if (atag1 < atag2) {
        tag_min = tag1;
        atag_min = atag1;
        atag_max = atag2;
      }
      else {
        tag_min = tag2;
        atag_min = atag2;
        atag_max = atag1;
      }

      int todd = ((atag_max -1)*atag_max/2 + atag_min) % 2 ;
      if (tag_min <= 0) { // 
        if (todd) // odds are skipped
          res =-1;
        else
          res = 1;
      }
      else {
        if (todd) // odds are accepted
          res = 1;
        else
          res =-1;
      }
    }
    //printf("Pair test: (%d)[%d]-(%d)[%d] %s\n", s1, std::abs(tag1), s2, std::abs(tag2), res<0 ? "other" : "my");
    return res;
# endif 
/*
    int c1=(int)(partition1[s1][icj1]>0);
    int c2=(int)(partition1[s2][ick2]>0);
    */
   
/*
    if (c1 != c2) {
      
      int tag1 = abs(partition1[s1][icj1]);
      int tag2 = abs(partition1[s2][ick2]);

      int ind;
      if (s1 == s2) { // pairs form diagonal matrix
        ind = tag_index(tag1 - 1, tag2 - 1);
        if (ind < 0) { // same tags ???
          return 1;
        }
      }
      else // pairs form square matrix
        ind = tag1 + tag2;
      if (ind % numprocs ) { // the first takes it
        if (c1)
          return 1;
        else
          return -1;
      }
      else { // the second takes it
        if (c1)
          return -1;
        else
          return 1;
      }
    }
    if (c1)
      return 1;
    else
      return -1;

*/
    /*
    if (c1 && c2)
      return 1;

    if (c1)
      return 2;

    if (c2)
      return -2;

    return -1;*/

    /*
    int res;
    if(c1!=c2){ // mixed
      int tag1=abs(partition1[s1][icj1]);
      int tag2=abs(partition1[s1][ick2]);
      int num=tag_index(tag1-1,tag2-1);
      if(num<0){ // compare wave packets
        int cmp= s1<2 ?
          wp[s1][icj1].compare(wp[s1][ick2],1e-15) :
          compare_vec(xi[icj1],xi[ick2],1e-15);
        if((cmp>0 && c1) || (cmp<0 && c2))
          res= 2; // my mixed term
        else
          res= -2; // not my term
      }
      else // parity check
        res=num%2 ? 2 : -2;
    }
    else if(c1)
      res=1; // all my
    else
      res=-1; // all other
    return res;*/
  }

  ///\en Returns electron-electron inter-partition multipliers for energy (first) and force (second)
  ///    for 2- electron additive terms (all inter-partition interactions are
  ///    calculated only once based on particle tags)
  ///    If force multiplier is zero, then the term may be omitted (energy will also be zero).
  pair<double, double> check_part1(int s1,int icj1,int s2, int ick2) const {
    int res=check_ee(s1,icj1,s2,ick2);
    if(res==1){ // my term
      //printf(" *\n");
      return make_pair(1.,1.); // all at my partition
    }
    else if(res==-1){
      //printf(" \n");
      return make_pair(0.,0.); // all at other partition
    }
    else if(res == 2){
      //printf(" *\n");
      return make_pair(0.5, 1.0); // my inter-partition
    }
    else if(res==-2){
      //printf(" \n");
      return make_pair(0., newton_pair ? 0.0 : 1. ); // other inter-partition: must add force if newton comm is off
    }
    return make_pair(0.,0.); // nonsense
  }

  ///\en Returns electron-electron inter-partition multipliers for energy (first) and force (second)
  ///    for 4-electron additive terms (all inter-partition interactions are
  ///    calculated only once based on particle tags)
  ///    If force multiplier is zero, then the term may be omitted (energy will also be zero).
  pair<double, double> check_part1(int s1, int icj1, int ick2, int s2, int icj3, int ick4) const {
    int res = check_ee(s1, icj1, ick2, s2, icj3, ick4);
    if (res == 1) { // my term
      //printf(" *\n");
      return make_pair(1., 1.); // all at my partition
    }
    else if (res == -1) {
      //printf(" \n");
      return make_pair(0., 0.); // all at other partition
    }
    else if (res == 2) {
      //printf(" *\n");
      return make_pair(0.5, 1.0); // my inter-partition
    }
    else if (res == -2) {
      //printf(" \n");
      return make_pair(0., newton_pair ? 0.0 : 1.); // other inter-partition: must add force if newton comm is off
    }
    return make_pair(0., 0.); // nonsense
  }

  ///\en Returns elctron-ion inter-partition multipliers for energy (first) and force (second)
  ///    for ion-electron additive terms (all inter-partition interactions are
  ///    calculated only once based on particle tags)
  ///    If force multiplier is zero, then the term may be omitted (energy will also be zero).
  ///    BASED ON ION ATTACHMENT
  pair<double,double> check_part1ei(int s1,int icj1,int ick2, int ion){
    int tagi = partition1[2][ion];
    if (tagi>0) {  // ion's node takes all
      //printf(" *\n");
      return make_pair(1., 1.); // my term
    }
    else {
      //printf(" \n");
      return make_pair(0., 0.); // all at other partition
    }
# if 0
    //printf("%d ",partition1[2][ion]);
    int ci=(int)(partition1[2][ion]>0);

    if(!newton_pair){ // care about mixed terms
      int cee=check_ee(s1,icj1,s1,ick2);
      if((cee==2 || cee==-2) || (ci && cee==-1) || (!ci && cee==1)) // all mixed variants
        make_pair(0., 1. ); // other inter-partition: must add force if newton comm is off
    }
    if(ci){
      //printf(" *\n");
      return make_pair(1.,1.); // my term
    }
    else{
      //printf(" \n");
      return make_pair(0.,0.); // all at other partition
    }
# endif
  }

  ///\en Returns ion-ion inter-partition multipliers for energy (first) and force (second)
  ///    for ion-ion additive terms (all inter-partition interactions are
  ///    calculated only once based on particle tags)
  ///    If force multiplier is zero, then the term may be omitted (energy will also be zero).
  pair<double,double> check_part1ii(int ion1, int ion2){
    return check_part1(2,ion1,2,ion2);
  }



  AWPMD():pbc(0),Lextra(-1),constraint(NONE),newton_pair(1), use_box(false) {
    nwp[0]=nwp[1]=nvar[0]=nvar[1]=ne[0]=ne[1]=ni=0;
    norm_matrix_state[0] = norm_matrix_state[1] = NORM_UNDEFINED;
    ovl_tolerance=0.;
    ovl_degeneracy_min = 0.1;
    approx = HARTREE;

    me=m_electron;
    one_h=1./h_plank;
    h2_me=h_sq/me;
    coul_pref=::coul_pref;

    calc_ee=1;
    calc_ii=0;
		calc_ei=1;
    norm_needed=0;

    screened_ei= 0;
    kappa_ei = 0.;

		norm_det_log[0] = norm_det_log[1] = 0.;
		ovl_det_log[0] = ovl_det_log[1] = 0.;

    ex_energy_type = 1;

		dft_extension = 0;
		dft_meshx = dft_meshy = dft_meshz = 0;
# if DFT_EXTENSION
		ex = NULL;
# endif
    calc_state = CALC_UNDEF;
    numprocs = 1;
  }


	~AWPMD() {
# if DFT_EXTENSION
		if (ex)
			delete ex;
# endif
	}

  ///\en Set the number of workers (processes or threads) for multi-partition processing
  ///    This number shoul be set to obtain correct energies and load balance in multi-partition processing
  void set_numprocs(int numprocs_) {
    numprocs = numprocs_;
  }


  ///\en Set harmonic box. The box starts to be used after this call.
  void set_box(const BoxHamiltonian& box_) {box = box_; use_box=true;}

  ///\en Returns a pointer to harmonic box if it is used, NULL otherwiae
  const BoxHamiltonian* get_box() const {
    if(!use_box)
      return NULL;
    else
      return &box;
  }

  ///\en Copy constructor.
  ///    Copies all structure information and energies, does not copy matrices and temporary data.
  ///    Force/energy recalculation reqires the call to interaction()
  AWPMD(const AWPMD &other){
    *this=other;
  }

  ///\en Copy constructor.
  ///    Copies all structure information and energies, does not copy matrices and temporary data.
  ///    Force/energy recalculation reqires the call to interaction()
  AWPMD &operator=(const AWPMD &other){
    if(this==&other)
      return *this;

    ni       =other.ni       ;
    qi       =other.qi       ;
    xi       =other.xi       ;
    pbc      =other.pbc      ;
    cell     =other.cell     ;
    Lextra   =other.Lextra   ;
    harm_w0_4=other.harm_w0_4;
    w0       =other.w0       ;
    calc_ee  =other.calc_ee  ;
		calc_ee  =other.calc_ei  ;
    calc_ii  =other.calc_ii  ;

    norm_needed          =other.norm_needed          ;
    constraint           =other.constraint           ;
    approx               =other.approx               ;
    ovl_tolerance        =other.ovl_tolerance        ;
    ovl_degeneracy_min   =other.ovl_degeneracy_min   ;
    Eee                  =other.Eee                  ;
    Ew                   =other.Ew                   ;
    Eii                  =other.Eii                  ;
    Eext                 =other.Eext                 ;
    Ebord                =other.Ebord                ;
    Ebord_ion            =other.Ebord_ion            ;
    Wext                 =other.Wext                 ;
    Edk                  =other.Edk                  ;
    Edc                  =other.Edc                  ;

    Ext_force            =other.Ext_force            ;
    Eiep                 =other.Eiep                 ;
    Eiip                 =other.Eiip                 ;

    me                   =other.me                   ;
    one_h                =other.one_h                ;
    h2_me                =other.h2_me                ;
    coul_pref            =other.coul_pref            ;
    newton_pair          =other.newton_pair          ;

    vars_per_wp          =other.vars_per_wp          ;
    

    screened_ei = other.screened_ei;
    kappa_ei = other.kappa_ei;
		
		box = other.box;
		use_box = other.use_box;

    for(int s=0;s<2;s++){
      ne[s]  =other.ne[s]  ;
      nwp[s] =other.nwp[s] ;
      nvar[s]=other.nvar[s];
      wp[s]  =other.wp[s]  ;
      qe[s]  =other.qe[s]  ;

      Ee[s]  =other.Ee[s]  ;
      Eei[s] =other.Eei[s] ;
      Eep[s] =other.Eep[s] ;
      Eeip[s]=other.Eeip[s];
      Eeep[s]=other.Eeep[s];
      Ewp[s] =other.Ewp[s] ;
      partition1[s]=other.partition1[s];

      E_der[s]      =  other.E_der[s];
      F_extra[s]    =  other.F_extra[s];
      F_wall[s]     =  other.F_wall[s];
      dq_dt[s]      =  other.dq_dt[s];
    }
    partition1[2]=partition1[2];
    lin_constraints= other.lin_constraints;

    Ee_exch         =other.Ee_exch         ;
    Eei_exch        =other.Eei_exch        ;
    Eee_exch        =other.Eee_exch        ;
    Eext_exch       =other.Eext_exch       ;
    Ebord_exch      =other.Ebord_exch      ; 

    Eee_hartree     =other.Eee_hartree     ;

    ex_energy_type  = other.ex_energy_type;
		dft_extension   = other.dft_extension;
    calc_state      = other.calc_state;

    numprocs = other.numprocs;
    return *this;
  }

protected:

  vector< vector< vector<double> > > lin_constraints; ///<\en directions in force space (of dimension nvar[s]), along which the motion is suppressed

public:
  int add_force_constraint(int s, int elid, vector<double> &constr_dir, double norm_tolerance = 1.e-10){
    Vector_G vconstr_dir(constr_dir);
    vconstr_dir.normalize();

    int vid=s*ne[0]+elid;
    vector< vector<double> > &constr=lin_constraints[vid]; // constraints for this electron

    double norm=gramm_schmidt_project(constr.begin(), constr.end(), vconstr_dir, vconstr_dir);
    if(norm/vconstr_dir.dim()>norm_tolerance){ // adding if it is not dependent
      constr.push_back(constr_dir);
      return 1;
    }
    else
      return 0;
  }

protected:

  bool apply_force_constraints(int s, int elid, vector<double> &force){
    Vector_G vforce(force);
    int vid=s*ne[0]+elid;
    vector< vector<double> > &constr=lin_constraints[vid]; // constraints for this electron
    gramm_schmidt_project(constr.begin(), constr.end(), vforce, vforce, 0.);
    return true;
  }

  ///\en Translates wp2 to the nearest image postion relative to wp1,
  ///    gets the translation vector
  Vector_3 move_to_image(const WavePacket &wp1, WavePacket &wp2) const {
    Vector_3 r1=wp1.get_r();
    Vector_3 r2=wp2.get_r();
    Vector_3 dr=r2-r1;
    Vector_3 ndr=dr.rcell(cell,pbc); // [0,L)
    ndr=ndr.rcell1(cell,pbc); // [-L/2,L/2)
    ndr-=dr;
    wp2=wp2.translate(ndr);  // wln.b=wln.b+ndr.a
    return ndr;
  }

  //e gets the overlap packet taking PBC into account
  WavePacket pbc_mul(const WavePacket &wp1, const WavePacket &wp2) const {
    if(!pbc)
      return wp1*conj(wp2);
    Vector_3 r1=wp1.get_r();
    Vector_3 r2=wp2.get_r();
    Vector_3 dr=r2-r1; // distance
    Vector_3 drn=dr.rcell1(cell,pbc); // distance within PBC
    Vector_3 rtrans=drn-dr; // new location of wp2 according to PBC (nearest image)
    WavePacket wpn=wp2.translate(rtrans);
    wpn=wp1*(conj(wpn));
    // reducing the result to elementary cell
    //r1=wpn.get_r();
    //r2=r1.rcell(cell,pbc);
    //dr=r2-r1;
    //wpn=wpn.translate(dr);
    return wpn;
  }
  ///\en resizes all internal arrays according to new electrons added
  virtual void resize(int flag);
public:



  ///\en Prepares to setup a new system of particles using \ref add_ion() and add_electron().
  ///    There is no need to call this function when using
  ///    \ref set_electrons() and \ref set_ions() to setup particles.
  virtual void reset(){
    for(int s=0;s<2;s++){
      nwp[s]=ne[s]=nvar[s]=0;
      wp[s].clear();
      qe[s].clear();
      partition1[s].clear();
      //partition2[s].clear();
      fixed[s].clear();
    }
    partition1[2].clear();
    ni=0;
    xi.clear();
    qi.clear();
  }

  //e sets Periodic Boundary Conditions
  //e using bit flags: 0x1 -- PBC along X
  //e                  0x2 -- PBC along Y
  //e                  0x4 -- PBC along Z
  //e cell specifies the lengths of the simulation box in all directions
  //e if PBCs are used, the corresponding coordinates of electrons and ions
  //e in periodic directions must be within the range  [0, cell[per_dir])
  //e @returns 1 if OK
  int set_pbc(const Vector_3P pcell=NULL, int pbc_=0x7);


  ///\en Setup electrons: forms internal wave packet representations.
  ///    If PBCs are used the coords must be within a range [0, cell).
  ///    Default electron mass is AWPMD::me.
  ///    Default (q=NULL )electron charges are -1.
  ///    \a pw_is_vel, if on, indicates that the \a pw array contains velocities (width_momentum[i]/mass) instead of momenta
  virtual int set_electrons(int spin, int n, const Vector_3P x, const Vector_3P v, const double* w, const double* pw, double mass=-1, const double *q=NULL, bool pw_is_vel = false);

  ///\en setup ion charges and coordinates
  /// if PBCs are used the coords must be within a range [0, cell)
  /// \a q0 is a multiplier for charges, or a charge to be set for all ions if q = NULL
  virtual int set_ions(int n, double* q, Vector_3P x, double q0=1);

  ///\en Adds an ion with charge q and position x,
  ///    \return id of the ion starting from 0
  ///    The tags must be nonzero, >0 for the local particle, <0 for ghost particle.
  ///    Unique particle id  is abs(tag).
  ///    Default tag (0) means inserting the current particle id as local particle.
  virtual int add_ion(double q, const Vector_3 &x, int tag=0){
    qi.push_back(q);
    xi.push_back(x);
    ni=(int)xi.size();
    if(tag==0)
      tag=ni;
    partition1[2].push_back(tag);
    return ni-1;
  }

  ///\en Sets flags to calculate specific interaction parts:\n
  ///     \a calc_ee_ -- calculate electron-electron interaction\n
  ///     \a calc_ii_ -- calculate ion-ion interaction\n
	///     \a calc_ei_ -- calculate electron-ion interaction\n
  void set_interaction_parts(int calc_ee_=1, int calc_ii_=0, int calc_ei_=1){
    calc_ee = calc_ee_;
    calc_ii = calc_ii_;
		calc_ei = calc_ei_;
  }


  //e calculates interaction in the system of ni ions + electrons
  //e the electonic subsystem must be previously setup by set_electrons, ionic by set_ions
  //e the iterators are describing ionic system only
  // 0x1 -- give back ion forces
  // 0x2 -- add ion forces to the existing set
  // 0x4 -- calculate derivatives for electronic time step 
  //e if PBCs are used the coords must be within a range [0, cell)
  virtual int interaction(int flag=0, Vector_3P fi=NULL, Vector_3P fe_x=NULL,
                                      Vector_3P fe_p=NULL, double *fe_w=NULL, double *fe_pw=NULL, Vector_2P fe_c=NULL);

  //e same as interaction, but using Hartee factorization (no antisymmetrization)
  virtual int interaction_hartree(int flag=0, Vector_3P fi=NULL, Vector_3P fe_x=NULL,
                                      Vector_3P fe_p=NULL, double *fe_w=NULL, double *fe_pw=NULL, Vector_2P fe_c=NULL);

  ///\en Calculates ion-ion interactions and updates Eii and ion forces if requested. This function
  ///    is called form intaraction() and interaction_hartree if calc_ii is set.

  virtual int interaction_ii(int flag,Vector_3P fi=NULL);

  virtual double interaction_border_ion(int i, double *x, double *f);

  virtual double interaction_border_electron(WavePacket const &packet, double *force, double *erforce, double *ervforce);

  virtual std::pair<double, double>
  interaction_electron_kinetic(WavePacket const &packet, int spin, double *erforce, double *ervfroce);

  virtual std::pair<double, double>
  interaction_electron_kinetic(double width, double pwidth, int spin, double *erforce, double *ervfroce);

  std::pair<double, double> coulomb_cutoff(double r, double cutoff) const{
    if (cutoff < 0)
      return std::make_pair<double, double>(1.0, 0.0);

    auto x = r / cutoff;
    return std::make_pair<double, double>(
        (((20.0 * x - 70.0) * x + 84.0) * x - 35.0) * x * x * x * x + 1, //energy
        (((140.0 * x - 420.0) * x + 420.0) * x - 140.0) * x * x * x / cutoff //denergy
        );
  }

  virtual double interaction_ee_single(WavePacket const &packet_1,
                                       WavePacket const &packet_2,
                                       double **eforce, double **erforce,
                                       double cutoff);

  virtual double interaction_ee_single(double const* coord1, double width1,
                                       double const* coord2, double width2,
                                       double **eforce, double **erforce,
                                       double cutoff);

  virtual double interation_ii_single(int i, int j, double **x, double *q,
                                      double **f, double cutoff);

  virtual double interaction_ei_single(double const *x, double q,
                                       WavePacket const &packet, int spin,
                                       double *f, double *eforce,
                                       double *erforce, double cutoff);

  virtual double interaction_ei_single(double const *x, double q,
                                       double const *packet_r, double width, int spin,
                                       double *f, double *eforce,
                                       double *erforce, double cutoff);
  ///\en adds confinement energies and forces to Ebord_ion and fi if requested. This function
  ///    is called form intaraction() and interaction_hartree if use_box is set.
  int calc_ion_confinement(int flag,Vector_3P fi);

  //e Calculates Norm matrix
  //e The result is saved in AWPMD::Norm[s]
  virtual void norm_matrix(int s);

  ///\en Performs LU-factorization of the Norm matrix
  /// AWPMD::Norm[s] is replaced by the LU matrix
  virtual int norm_factorize(int s);

  ///\en Inverts Norm matrix.
  ///    AWPMD::Norm[s] is replaced by the inverted matrix.
  virtual int norm_invert(int s);

  //e Get the determinant of the norm-matrix for the particles with spin s
  virtual double norm_matrix_det(int s);

  //e Get the determinant logarithm of the norm-matrix for the particles with spin s
  virtual double norm_matrix_detl(int s);

	///\en Get the determinant logarithm of the overlap matrix for the particles with spin s
	///    This value is normalized to the norm matrix dimension.
	virtual double ovl_matrix_detl(int s){
		return 8*ovl_det_log[s];  // eight variables per overlap
	}

  ///\en If \a use_ee_hartree is true, exchange part of Eee is ignored (less divergent at zero overlap)
  double get_energy(bool use_ee_hartree = false);

	///\en Makes timestep \a dt of electronic component: q-> q + (dq_dt)*dt for each variable.
	///    If flag contains 0x10 uses internal variables, \a spin <0 means go through all spins,
	///    the vaector \a dq_dt_ should contain generalized force (for all used spins), if NULL, they are taken from previous \ref interaction(). 
	virtual int step(double dt, int flag =0, int spin =-1,  const vector<double> *dq_dt_=NULL);


	///\en Set parameters of the given wavepacket, parameter \a c is used for Split WPMD only.
  virtual int set_wavepacket(int spin, int wp_id, const Vector_3 *x, const Vector_3 *v, const double *w, const double *pw, const Vector_2 *c=NULL, double mass=-1){
    if(mass<0)
      mass=me;
    wp[spin][wp_id].init(*w,*x,(*v)*mass*one_h,(*pw)*one_h);
    return 1;
  }

  ///\en Gets current electronic coordinates.
  ///    Transforms the momenta to velocity v according to mass setting (-1 means me)
  virtual int get_electrons(int spin, Vector_3P x, Vector_3P v, double* w, double* pw, double mass=-1);

  void set_harm_constr(double w0) {
    constraint = HARM;
    harm_w0_4 = h_sq*9./(8.*m_electron)/(w0*w0*w0*w0);
  }

  void set_fix_constr(double w0_) {
    constraint = FIX;
    w0 = w0_;
  }

  ///\en Prepares force arrays according to \a flag setting for interaction()
  virtual void clear_forces(int flagi,Vector_3P fi, Vector_3P fe_x,
                    Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c=NULL);


  ///\en Creates wave packet acording to the given physical parameters.
  ///    The actual values for the created wp are ajusted according to the existing constraints!
  ///    Default mass (-1) is the electron mass AWPMD::me.
  ///    \a pw_is_vel, if on, indicates that pw contains velocity (width_momentum/mass) instead of momentum
  WavePacket create_wp(const Vector_3 &x, const Vector_3 &v, double w, double pw, double mass=-1, bool pw_is_vel=false);


  ///\en Copies dynamic variables in physical representation into specified pointer, the size of
  ///    the pointer may be checked by calling this function with varptr=NULL\n
  ///    \return number of copied variables if OK, -1 for invalid request.
  ///    \param type: 0=electrons, 1= ions, -1 = all \n
  ///    NOT IMPLEMENTED YET: \n
  ///    \param spin: 0= up, 1=down, -1= all, valid if electrons are selected by type\n
  ///    \param particle_id: id of electron or ion, -1 for all\n
  ///    \param parameter_id ...
  virtual int get_variables(double *varptr=NULL, int type= 0) const {
    int size=-1;
    // calculating size of the requested array
    if(type==0)
      size= nvar[0]+nvar[1];
    else if(type==1)
      size= 3*ni;
    else if(type==-1)
      size= nvar[0]+nvar[1]+3*ni;
    if(!varptr)
      return size;
    // getting values
    if(type==0 || type==-1){ //electrons
      for(int s=0;s<2;s++){
        for(int c=0;c<nwp[s];c++){
          get_wavepacket(s,c,(Vector_3 *)varptr,(Vector_3 *)(varptr+3),varptr+6,varptr+7, vars_per_wp>8 ? (Vector_2 *)(varptr+8): NULL);
          varptr+=vars_per_wp;
        }
      }
    }
    if(type==1 || type==-1){ // ions
      for(int i=0;i<ni;i++){
        *((Vector_3 *)varptr)=xi[i];
        varptr+=3;
      }
    }
    return size;
  }

  ///\en Copies dynamic variables from specified pointer to physical representation, the size of
  ///    the pointer may be checked by calling this function with varptr=NULL\n
  ///    \return number of copied variables if OK, -1 for invalid request.
  ///    parameters are the same as in \ref get_variables
  virtual int set_variables(const double *varptr=NULL, int type= 0){
    int size=get_variables(NULL,type);
    if(!varptr || size<0)
      return size;
    // getting values
    if(type==0 || type==-1){ //electrons
      for(int s=0;s<2;s++){
        for(int c=0;c<nwp[s];c++){
          set_wavepacket(s,c,(Vector_3 *)varptr,(Vector_3 *)(varptr+3),varptr+6,varptr+7, vars_per_wp>8 ? (Vector_2 *)(varptr+8): NULL);
          varptr+=vars_per_wp;
        }
      }
    }
    if(type==1 || type==-1){ // ions
      for(int i=0;i<ni;i++){
        xi[i]=*((Vector_3 *)varptr);
        varptr+=3;
      }
    }
    return size;
  }

  ///\en Get parameters of the given wavepacket, parameter \a c is used for Split WPMD only.
  virtual int get_wavepacket(int spin, int wp_id, Vector_3 *x, Vector_3 *v, double *w, double *pw, Vector_2 *c=NULL, double mass=-1) const {
    if(mass<0)
      mass=me;
    *w=wp[spin][wp_id].get_width();
    *pw=wp[spin][wp_id].get_pwidth()/one_h;
    *x=wp[spin][wp_id].get_r();
    *v=wp[spin][wp_id].get_p()/(mass*one_h);

    /*
    *w=sqrt(3./(4*real(wp[spin][wp_id].a)));
    *pw=-2*(*w)*imag(wp[spin][wp_id].a)/one_h;
    *x=real(wp[spin][wp_id].b)/(2*real(wp[spin][wp_id].a));
    *p=((*pw)*(*x)/(*w) + imag(wp[spin][wp_id].b)/one_h);*/
    return 1;
  }

 
  ///\en If \a use_screening is 1, switches e-i interaction to Yukawa potential with decrement \a kappa_ei,\n
  ///    otherwise uses Coulomb interaction.
  int set_screening(int use_screening, double kappa_ei_=0.){
    screened_ei=use_screening;
    kappa_ei=kappa_ei_;
    return 1;
  }

  ///\en Gets current wave packet number \a Ns for spin variable \a spin (may be 0 or 1).
  ///    The wave packets with ids [0, Ns) may be then accessed by \ref get_wavepacket()
  ///    and \ref set_wavepacket() functions.
  int get_wp_number(int spin) const {
    if(spin<0 || spin >1)
      return 0;
    return nwp[spin];
  }

  ///\en Gets current electron number \a Ns for spin variable \a spin (may be 0 or 1).
  int get_electron_number(int spin) const {
    if(spin<0 || spin >1)
      return 0;
    return ne[spin];
  }

  ///\en Calculates pressure by evaluating forces acting from the box. Electron and ion components a summed up.
  Vector_3  calc_box_pressure() const;
  

  ///\en calculates overlap penalty energy: \sum_ij w0/(1-ovl(i,j))*E0  
  double calc_overlap_penalty(double w0, double E0);
#if DFT_EXTENSION
	void configure_dft(int flag, size_t mesh_nx, size_t mesh_ny, size_t mesh_nz, size_t max_packets,
                       ApproxType atype=ApproxType::T_LDA, IApproximation *custom_approx = nullptr) {
		dft_extension = flag;
//# if DFT_EXTENSION
		if (mesh_nx != dft_meshx || mesh_nz != dft_meshz || mesh_nz != dft_meshz || max_packets !=dft_conf.packet_number) { // same parameters?
			dft_meshx = mesh_nx; dft_meshy = mesh_ny; dft_meshz = mesh_nz;

          if (custom_approx != nullptr)
            dft_conf.approximation = custom_approx;
          else{
            switch (atype){
              case ApproxType::T_LDA:
                dft_conf.approximation = new LDA(0.738558766f, -0.01554534543482745f, 20.4562557f);
                break;
              case ApproxType::T_LDA_2:
                dft_conf.approximation = new ModLDA();
                break;

              case ApproxType::T_LSDA:
                dft_conf.approximation = new LSDA();
                break;

              case ApproxType::T_VOID:
                dft_conf.approximation = new VoidApproximation();
                    break;
            };
          }

			dft_conf.mesh_size.size.as_struct = {(unsigned int)dft_meshx, (unsigned int)dft_meshx, (unsigned int)dft_meshx};
			dft_conf.mesh_start = {-1.8f*(float)cell[0], -1.8f*(float)cell[1], -1.8f*(float)cell[2]};
      dft_conf.mesh_fin = {1.8f*(float)cell[0], 1.8f*(float)cell[1], 1.8f*(float)cell[2]};
			dft_conf.packet_number = max_packets;
			dft_conf.use_adaptive_mesh = true;


			//std::vector<deriv_function> od = DerivsFunction::GetFunctions();
			if (ex)
				delete ex;

			ex =  new XCEnergy_cpu(dft_conf.packet_number, dft_conf);

		}
//# endif
	}
#endif
};



class wpmd_term {
  vector<double> vars;
  vector<int> varmap;
  AWPMD *sys;
public:
  typedef double value_t;
  wpmd_term(AWPMD *sys_):sys(sys_){
    vars.resize((size_t)sys->get_variables(NULL,0));
    sys->get_variables(&vars[0],0);
    //getting the number of relaxed variables
    for(int s=0;s<2;s++){
      for(map<int,int>::iterator it=sys->fixed[s].begin();it!=sys->fixed[s].end();++it){
        if(it->second==2)
          varmap.push_back(it->first+ (s ? sys->nvar[0] :0) );
      }
    }
  }

  template<class it_t>
  int copy_variable(int type, bool write_to_xi, it_t xi){
    if(write_to_xi){
      sys->get_variables(&vars[0],0);
      for(size_t i=0;i<varmap.size();i++)
        *xi++=vars[varmap[i]];
    }
    else{
      for(size_t i=0;i<varmap.size();i++)
        vars[varmap[i]]=*xi++;
      sys->set_variables(&vars[0],0);
    }
    return 1;
  }

  int dimension() const {
    return (int)varmap.size();
  }

  int compute(double *E, double *fi=NULL){
    int flag= fi ? 0x4 : 0;
    int res=sys->interaction(flag);
    if(fi && res>0){
      for(size_t i=0;i<varmap.size();i++){
        int varid=varmap[i];
        int s=0;
        if(varid>sys->nvar[0]){
          s=1;
          varid-=sys->nvar[0];
        }
        fi[i]=-sys->E_der[s][varid];
      }
    }
    if(E)
      *E=sys->get_energy();
    return res;
  }
};

int interaction_relax(AWPMD *sys, double thresh, int flag);



# endif


