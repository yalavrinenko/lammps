/*s***************************************************************************eterm
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *****************************************************************************/

/*s****************************************************************************
 * $Log: wpmd_split.cpp,v $
 * Revision 1.110  2015/12/04 19:16:06  valuev
 * fixed exchange energy measurement
 *
 * Revision 1.109  2015/12/03 15:02:43  valuev
 * added exchange energies
 *
 * Revision 1.108  2015/11/19 08:46:45  valuev
 * common base for AWP
 *
 * Revision 1.107  2015/10/20 15:01:13  valuev
 * switched to width velocity
 *
 * Revision 1.106  2015/10/12 12:59:47  valuev
 * added overlap tolerance for 4-particle integrals
 *
 * Revision 1.105  2015/10/09 19:33:26  valuev
 * modified AWPMC to output the same vars as AWPMD
 *
 * Revision 1.104  2015/04/10 13:18:48  valuev
 * updated state write
 *
 * Revision 1.103  2015/03/13 11:19:43  valuev
 * fixes for ion dynamics
 *
 * Revision 1.102  2015/02/10 14:08:09  valuev
 * added step function
 *
 * Revision 1.101  2014/07/28 10:45:35  valuev
 * added partial WP integration
 *
 * Revision 1.100  2014/07/21 16:55:22  morozov
 * Averaging electron density
 *
 * Revision 1.99  2014/07/21 13:33:17  valuev
 * quantum MC, density for UHF case
 *
 * Revision 1.98  2014/07/11 18:14:18  valuev
 * added overlap matrix log
 *
 * Revision 1.97  2014/07/09 16:22:52  valuev
 * removed reference to h_plank from wavepacket.h
 *
 * Revision 1.96  2014/07/08 17:07:21  valuev
 * prepared confined MC simulation with AWPMD nonsplit norm matrix
 *
 * Revision 1.95  2014/07/03 16:58:42  valuev
 * added overlap matrix det log measurement
 *
 * Revision 1.94  2014/07/03 13:32:02  valuev
 * fixed e-e energy for different spins, generalized 2-center term for vector operations
 *
 * Revision 1.93  2014/07/02 15:26:06  valuev
 * fixed wp order
 *
 * Revision 1.92  2014/07/02 14:59:54  valuev
 * added forgotten yy term
 *
 * Revision 1.91  2014/07/02 09:36:16  valuev
 * added harmonic state projection onto many-electron (antisymmetrized) wave function
 *
 * Revision 1.90  2014/06/27 10:38:33  valuev
 * added second spin to state_io, checked norm_matrix, corrected y derivative for box
 *
 * Revision 1.89  2014/06/25 10:53:15  valuev
 * preparing for confined simulations
 *
 * Revision 1.88  2014/06/23 16:13:58  morozov
 * Added old changes (on August 2013) related to the ion motion. Added box_hamiltonian files. Fixed compilation on VS2005 of cvector_3.h and box_hamiltonian.*.
 *
 * Revision 1.87  2014/04/25 18:31:33  kazeev
 * Added tested integration for box_hamiltonian.
 *
 * Revision 1.86  2014/04/22 02:50:58  kazeev
 * Added tests and implementation for system in a potential box. The implemetations seems to be incorrectly integrated.
 *
 * Revision 1.85  2014/04/21 18:42:41  kazeev
 * Added box_eterm_deriv. Untested.
 *
 * Revision 1.84  2014/04/17 13:51:00  kazeev
 * Fixed gcc compilation. Added BoxHamiltonian.
 *
 * Revision 1.83  2014/04/01 14:07:35  valuev
 * prepared variables for boundary potential
 *
 * Revision 1.82  2013/10/17 12:32:01  morozov
 * Added AWPMD_split::el_density for calculation of |Psi(x)|^2
 *
 * Revision 1.79  2013/08/19 16:08:51  morozov
 * Fixed ion forces
 *
 * Revision 1.78  2013/08/15 10:17:31  valuev
 * overap cutoff fixes
 *
 * Revision 1.77  2013/03/23 09:40:34  valuev
 * new algorithm and fixes for degeneracy constraints
 *
 * Revision 1.76  2013/03/08 10:11:05  kazeev
 * Restored passage through fixed output intervals;
 *
 * Revision 1.75  2012/10/04 16:37:45  valuev
 * sync with KINTECH svn
 *
 * Revision 1.74  2012/10/04 11:42:05  valuev
 * first working version of constrainde dynamics
 *
 * Revision 1.73  2012/09/28 11:25:40  valuev
 * updated degeneracy constraints
 *
 * Revision 1.72  2012/09/20 05:38:31  valuev
 * optimizer (test version)
 *
 * Revision 1.71  2012/09/07 16:41:49  valuev
 * changed N sign
 *
 * Revision 1.70  2012/09/03 15:54:46  valuev
 * fixed sign in WPMD forces
 *
 * Revision 1.69  2012/08/27 16:28:14  valuev
 * fixed a->0 condition
 *
 * Revision 1.68  2012/08/22 12:34:03  valuev
 * Made compatible with icc10
 *
 * Revision 1.67  2012/08/22 11:57:51  valuev
 * changes for KIM compliance
 *
 * Revision 1.66  2012/06/30 12:17:04  valuev
 * fixed constrained norm matrix formation
 *
 * Revision 1.65  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.64  2012/06/15 16:59:00  morozov
 * OpenMP for electron-ion interaction and 2D psi dump. Added switching between dE and dE/dt criteria.
 *
 * Revision 1.63  2012/05/10 15:59:46  morozov
 * Fixed logging on each step. Enabled VTS for minimization (tested for 'Simple' algorithm). Added calculation of closest ions along the trajectory.
 *
 * Revision 1.62  2012/05/10 14:51:29  valuev
 * fixed Yakawa energy and force at zero distance
 *
 * Revision 1.61  2012/05/04 11:08:22  valuev
 * added eigenvalue decomposition norm matrix inversion
 *
 * Revision 1.60  2012/05/03 15:26:56  valuev
 * added fixed variables
 *
 * Revision 1.59  2012/04/29 08:52:29  valuev
 * added power-conditional constraint
 *
 * Revision 1.58  2012/04/26 09:21:35  valuev
 * added Yukawa potential (tested)
 *
 * Revision 1.57  2012/04/17 20:42:31  valuev
 * fixed Psi calculation, started adding e-i Debye screening
 *
 * Revision 1.56  2012/04/03 10:26:41  valuev
 * added variable constraints
 *
 * Revision 1.55  2012/03/30 14:30:55  morozov
 * Added fixed width simulation mode
 *
 * Revision 1.54  2012/03/23 15:50:57  valuev
 * added Psi calculation
 *
 * Revision 1.53  2012/02/01 14:55:53  valuev
 * improved spectral calculation
 *
 * Revision 1.52  2012/01/27 10:54:03  valuev
 * fixed wrong term in the Hartree norm matrix
 *
 * Revision 1.51  2012/01/26 10:03:19  valuev
 * added free norm option
 *
 * Revision 1.50  2012/01/23 12:36:17  valuev
 * added wave function overlap calculation
 *
 * Revision 1.49  2011/10/21 16:29:43  valuev
 * added He init states, fixed calculation for He
 *
 * Revision 1.48  2011/10/12 18:58:48  valuev
 * added polar mode for C
 *
 * Revision 1.47  2011/10/12 16:39:02  valuev
 * restructured  interaction function
 *
 * Revision 1.46  2011/10/11 18:31:50  valuev
 * fixed external force calculation
 *
 * Revision 1.45  2011/10/10 15:52:42  morozov
 * Fixed prefactor for external energy derivatives
 *
 * Revision 1.44  2011/10/10 04:37:22  morozov
 * Fixed indexing at Wext calculation
 *
 * Revision 1.43  2011/10/08 11:23:16  valuev
 * comment on power calculation
 *
 * Revision 1.42  2011/10/08 11:21:14  valuev
 * added external force storage and external force power calculation
 *
 *
*******************************************************************************/
# include <set>
#include "../awpmd-dft/awpmd-dft.hpp"
# include "wpmd_split.h"
# include "math_utils.h"
//# include "erf.h"

#ifdef _OPENMP
#include <omp.h>
#endif


void AWPMD_split::resize(int flag) {
  // find out whether we have a split/no split version: all splits should be 1 for no split
  split_wp = false;
  for (int s1 = 0; s1 < 2; s1++) {
    for (int c1 = 0; c1 < ne[s1]; c1++) {
      if (nspl[s1][c1] > 1) {
        split_wp = true;
        break;
      }
    }
  }
  for (int s = 0; s < 2; s++) {
    wf_norm[s].resize(ne[s]);
    //if(flag&(0x8|0x4)){ //electron forces needed

    //}


    if (flag & (0x8 | 0x4) || norm_needed) { //electron forces or norm matrix needed
      wf_norm_der[s].resize(nvar[s]);
      ovl_der[s].resize(nvar[s]);
      if (approx == HARTREE) { // L and Norm are needed in block form
        Lh[s].resize(nvar[s]);
        if (norm_needed) {
          Normh[s].resize(ne[s]);
          Jh[s].resize(ne[s]);
          rconstr_pivots[s].resize(ne[s]);
          constr_pivots[s].resize(ne[s]);
          for (int i = 0; i < ne[s]; i++)
            Normh[s][i].init(10 * nspl[s][i]);
        }
      } else if (norm_needed) {
        //if(split_wp)
        Norm[s].init(nvar[s]);
        //else
        //Norm[s].init(8*ne[s]);  // single split case
      }
    }

  }
  AWPMD::resize(flag);
}


int
AWPMD_split::add_split(Vector_3 &x, Vector_3 &v, double &w, double &pw, Vector_2 &c, double mass, double q, int tag) {
  if (!spl_add) {
    nspl[s_add].push_back(1);
    ne[s_add]++;
  } else {
    nspl[s_add][ne[s_add] - 1]++; // incrementing the WP number for the last electron
  }
  spl_add++;
  nwp[s_add]++;
  nvar[s_add] += 10;
  wp[s_add].push_back(create_wp(x, v, w, pw, mass));

  if (!c_polar_mode)
    split_c[s_add].push_back(c);
  else { // converting from rho, alpha
    complex<double> cc(polar(c[0], c[1]));
    split_c[s_add].push_back(Vector_2(real(cc), imag(cc)));
  }
  qe[s_add].push_back(q);

  if (tag == 0)
    tag = spl_add;
  partition1[s_add].push_back(tag);

  valid_norms = false;
  return nwp[s_add] - 1; //spl_add-1;
}

int AWPMD_split::set_electrons(int s, int nel, const Vector_3P x, const Vector_3P v, const double *w, const double *pw,
                               const Vector_2 *c, const int *splits, double mass, const double *q, bool pw_is_vel,
                               double q0) {
  if (s < 0 || s > 1)
    return LOGERR(-1, fmt_iv("AWPMD_split.set_electrons: invaid spin setting (%d)!", s), LINFO);
  calc_state &= CALC_IONS_SET;  // refereshes all previously calculated quantities
  calc_state |= CALC_ELECTRONS_SET;
  // calculating the total n
  nvar[s] = 0;
  int n = 0;
  for (int i = 0; i < nel; i++) {
    n += splits[i];
    nvar[s] += 10 * splits[i]; // number of dynamic variables per wp: x[3],p[3],w,pw,c_re,c_im
  }
  nwp[s] = n;

  norm_matrix_state[s] = NORM_UNDEFINED;
  ne[s] = nel;
  wp[s].resize(n);

  split_c[s].resize(n);
  if (!c_polar_mode)
    split_c[s].assign(c, c + n);
  else { // converting from rho, alpha
    for (int i = 0; i < n; i++) {
      complex<double> cc(polar(c[i][0], c[i][1]));
      split_c[s][i][0] = real(cc);
      split_c[s][i][1] = imag(cc);
    }
  }


  nspl[s].resize(nel);
  nspl[s].assign(splits, splits + nel);

  partition1[s].clear();
  for (int i = 0; i < n; i++) {

    /*if(constraint==FIX){
      w[i]=w0;
      pw[i]=0.;
    }

    double rw;
    if(Lextra>0){ // width PBC, keeping the width are within [0,Lextra]
      w[i]=fmod(w[i],Lextra);
      if(w[i]<0) w[i]+=Lextra;
      rw=w[i]; // WP width for energy evaluation is within [0, L/2]
      if(rw > Lextra/2) rw = Lextra - rw;
    }
    else
      rw=w[i];

    wp[s][i].init(rw,x[i],v[i]*m_electron/h_plank,pw[i]/h_plank);*/
    wp[s][i] = create_wp(x[i], v[i], w[i], pw[i], mass, pw_is_vel);
    //printf("%15d %15g\n",i,rw);
    // assign default partition
    partition1[s].push_back(i + 1);
  }

  // assign electronic charge
  if (q)
    qe[s].assign(q, q + nwp[s]);
  else
    qe[s].assign(nwp[s], q0);

  valid_norms = false;
  return 1;
}

int AWPMD_split::get_electrons(int s, Vector_3P x, Vector_3P v, double *w, double *pw, Vector_2 *c, int *splits,
                               double mass) {
  if (s < 0 || s > 1)
    return LOGERR(-1, fmt_iv("AWPMD_split.set_electrons: invaid spin setting (%d)!", s), LINFO);

  if (c) {
    if (!c_polar_mode) {
      for (int i = 0; i < nwp[s]; i++)
        c[i] = split_c[s][i];
    } else { // converting to rho, alpha
      for (int i = 0; i < nwp[s]; i++) {
        complex<double> cc(split_c[s][i][0], split_c[s][i][1]);
        c[i][0] = abs(cc);
        c[i][1] = arg(cc);
      }
    }
  }

  if (mass < 0)
    mass = me;
  for (int i = 0; i < nwp[s]; i++) {
    x[i] = wp[s][i].get_r();
    v[i] = wp[s][i].get_p();
    w[i] = wp[s][i].get_width();
    pw[i] = wp[s][i].get_pwidth();
    //get_wavepacket(s,i,x+i,v+i,w+i,pw+i);
    v[i] /= mass * one_h;
    pw[i] /= one_h;
  }

  return nwp[s];
}

int AWPMD_split::get_splits_spl2(int s, double *c0, double *c1) const {
  if (s < 0 || s > 1)
    return LOGERR(-1, fmt_iv("AWPMD_split.get_splits_spl2: invaid spin setting (%d)!", s), LINFO);

  if (c0 && c1) {
    if (!c_polar_mode) {
      for (int i = 0; i < nwp[s]; i++) {
        c0[i] = split_c[s][i][0];
        c1[i] = split_c[s][i][1];
      }
    } else { // converting to rho, alpha
      for (int i = 0; i < nwp[s]; i++) {
        complex<double> cc(split_c[s][i][0], split_c[s][i][1]);
        c0[i] = abs(cc);
        c1[i] = arg(cc);
      }
    }
  }
  return nwp[s];
}


void AWPMD_split::eterm_deriv(int ic1, int s1, int c1, int j1,
                              int ic2, int s2, int c2, int k2, cdouble pref,
                              const OverlapDeriv &o, cdouble v, cdouble dv_aj_conj,
                              cdouble dv_ak, cVector_3 dv_bj_conj, cVector_3 dv_bk, int external) {
  vector<double> *E_der = external ? AWPMD_split::F_extra : AWPMD_split::E_der;
  cdouble cj(split_c[s1][ic1 + j1][0], split_c[s1][ic1 + j1][1]);
  cdouble ck(split_c[s2][ic2 + k2][0], split_c[s2][ic2 + k2][1]);
  int indw1 = 8 * ic1, indw2 = 8 * ic2;
  int indn1 = (nvar[s1] / 10) * 8 + 2 * ic1, indn2 = (nvar[s2] / 10) * 8 + 2 * ic2;
  cdouble part_jk = conj(cj) * ck;

  int M = 1; //(j1==k2 ? 1 : 2);

  // over a_k_re
  E_der[s2][indw2 + 8 * k2] += M * real(pref * part_jk * (o.da2_re() * v + o.I0 * dv_ak));
  // over a_k_im
  E_der[s2][indw2 + 8 * k2 + 1] += M * real(pref * part_jk * (o.da2_im() * v + i_unit * o.I0 * dv_ak));
  // over a_j_re
  E_der[s1][indw1 + 8 * j1] += M * real(pref * part_jk * (o.da1_re() * v + o.I0 * dv_aj_conj));
  // over a_j_im
  E_der[s1][indw1 + 8 * j1 + 1] += M * real(pref * part_jk * (o.da1_im() * v - i_unit * o.I0 * dv_aj_conj));

  for (int i = 0; i < 3; i++) {
    // over b_k_re
    E_der[s2][indw2 + 8 * k2 + 2 + 2 * i] += M * real(pref * part_jk * (o.db2_re(i) * v + o.I0 * dv_bk[i]));
    // over b_k_im
    E_der[s2][indw2 + 8 * k2 + 2 + 2 * i + 1] +=
        M * real(pref * part_jk * (o.db2_im(i) * v + i_unit * o.I0 * dv_bk[i]));
    // over b_j_re
    E_der[s1][indw1 + 8 * j1 + 2 + 2 * i] += M * real(pref * part_jk * (o.db1_re(i) * v + o.I0 * dv_bj_conj[i]));
    // over b_j_im
    E_der[s1][indw1 + 8 * j1 + 2 + 2 * i + 1] +=
        M * real(pref * part_jk * (o.db1_im(i) * v - i_unit * o.I0 * dv_bj_conj[i]));
  }


  // over ck_re
  E_der[s2][indn2 + 2 * k2] += M * real(pref * conj(cj) * o.I0 * v);
  // over ck_im
  E_der[s2][indn2 + 2 * k2 + 1] += M * real(pref * i_unit * conj(cj) * o.I0 * v);
  // over cj_re
  E_der[s1][indn1 + 2 * j1] += M * real(pref * ck * o.I0 * v);
  // over cj_im
  E_der[s1][indn1 + 2 * j1 + 1] += M * real(-pref * i_unit * ck * o.I0 * v);

  if (norm_mode == NORMALIZE) { // add norm derivative contributions
    double t = -M * real(pref * part_jk * o.I0 * v);
    // nonlocal terms: TODO: make a separate global loop for summation of nonlocal terms
    for (int j = 0; j < nspl[s1][c1]; j++) {
      for (int i = 0; i < 8; i++)
        E_der[s1][indw1 + 8 * j + i] += t * wf_norm_der[s1][indw1 + 8 * j + i];
      E_der[s1][indn1 + 2 * j] += t * wf_norm_der[s1][indn1 + 2 * j];
      E_der[s1][indn1 + 2 * j + 1] += t * wf_norm_der[s1][indn1 + 2 * j + 1];
    }
    for (int k = 0; k < nspl[s2][c2]; k++) {
      for (int i = 0; i < 8; i++)
        E_der[s2][indw2 + 8 * k + i] += t * wf_norm_der[s2][indw2 + 8 * k + i];
      E_der[s2][indn2 + 2 * k] += t * wf_norm_der[s2][indn2 + 2 * k];
      E_der[s2][indn2 + 2 * k + 1] += t * wf_norm_der[s2][indn2 + 2 * k + 1];
    }
  }
}

void AWPMD_split::eterm_deriv(packet_index_info packet_1, packet_index_info packet_2, cdouble pref,
                              const OverlapDeriv &o, cdouble v, cdouble dv_aj_conj,
                              cdouble dv_ak, cVector_3 dv_bj_conj, cVector_3 dv_bk, int external) {
  vector<double> *E_der = external ? AWPMD_split::F_extra : AWPMD_split::E_der;
  if (E_der[packet_1.spin].empty())
    E_der[packet_1.spin].resize( 2 * 8);

  if (E_der[packet_2.spin].empty())
    E_der[packet_2.spin].resize( 2 * 8);
  for (auto &eder : E_der[packet_1.spin])
    eder = 0;
  for (auto &eder : E_der[packet_2.spin])
    eder = 0;
  cdouble part_jk {1.0, 0.0};
  int packet_1_shift = packet_1.index * 8;
  int packet_2_shift = packet_2.index * 8;
  int M = 1;
  // over a_k_re
  E_der[packet_2.spin][packet_2_shift] += M * real(pref * part_jk * (o.da2_re() * v + o.I0 * dv_ak));
  // over a_k_im
  E_der[packet_2.spin][packet_2_shift + 1] += M * real(pref * part_jk * (o.da2_im() * v + i_unit * o.I0 * dv_ak));
  // over a_j_re
  E_der[packet_1.spin][packet_1_shift] += M * real(pref * part_jk * (o.da1_re() * v + o.I0 * dv_aj_conj));
  // over a_j_im
  E_der[packet_1.spin][packet_1_shift + 1] += M * real(pref * part_jk * (o.da1_im() * v - i_unit * o.I0 * dv_aj_conj));

  for (int i = 0; i < 3; i++) {
    // over b_k_re
    E_der[packet_2.spin][packet_2_shift + 2 + 2 * i] += M * real(pref * part_jk * (o.db2_re(i) * v + o.I0 * dv_bk[i]));
    // over b_k_im
    E_der[packet_2.spin][packet_2_shift + 2 + 2 * i + 1] +=
        M * real(pref * part_jk * (o.db2_im(i) * v + i_unit * o.I0 * dv_bk[i]));
    // over b_j_re
    E_der[packet_1.spin][packet_1_shift + 2 + 2 * i] += M * real(pref * part_jk * (o.db1_re(i) * v + o.I0 * dv_bj_conj[i]));
    // over b_j_im
    E_der[packet_1.spin][packet_1_shift + 2 + 2 * i + 1] +=
        M * real(pref * part_jk * (o.db1_im(i) * v - i_unit * o.I0 * dv_bj_conj[i]));
  }
}
void AWPMD_split::box_eterm_deriv(
    const unsigned int ic1,
    const unsigned int s1,
    const unsigned int c1,
    const unsigned int j1,
    const unsigned int ic2,
    const unsigned int s2,
    const unsigned int c2,
    const unsigned int k2,
    const cdouble &pref) {
  cdouble integral, a1_real_diff, a1_imag_diff, a2_real_diff, a2_imag_diff;
  cVector_3 b1_real_diff, b1_imag_diff, b2_real_diff, b2_imag_diff;

  cdouble cj(split_c[s1][ic1 + j1][0], split_c[s1][ic1 + j1][1]);
  cdouble ck(split_c[s2][ic2 + k2][0], split_c[s2][ic2 + k2][1]);
  int indw1 = 8 * ic1, indw2 = 8 * ic2;
  int indn1 = (nvar[s1] / 10) * 8 + 2 * ic1, indn2 = (nvar[s2] / 10) * 8 + 2 * ic2;
  cdouble part_jk = conj(cj) * ck;

  box.get_derivatives(wp[s1][ic1 + j1].a, wp[s1][ic1 + j1].b,
                      wp[s2][ic2 + k2].a, wp[s2][ic2 + k2].b, &integral,
                      &a1_real_diff, &a1_imag_diff,
                      &b1_real_diff, &b1_imag_diff,
                      &a2_real_diff, &a2_imag_diff,
                      &b2_real_diff, &b2_imag_diff);

  vector<double> *E_der_v[2] = { AWPMD_split::F_wall , AWPMD_split::E_der };
  for (int l = 0; l < 2; l++) {

    E_der_v[l][s2][indw2 + 8 * k2] += real(pref * part_jk * a2_real_diff);
    E_der_v[l][s2][indw2 + 8 * k2 + 1] += real(pref * part_jk * a2_imag_diff);
    E_der_v[l][s1][indw1 + 8 * j1] += real(pref * part_jk * a1_real_diff);
    E_der_v[l][s1][indw1 + 8 * j1 + 1] += real(pref * part_jk * a1_imag_diff);

    for (int i = 0; i < 3; i++) {
      E_der_v[l][s2][indw2 + 8 * k2 + 2 + 2 * i] += real(pref * part_jk * b2_real_diff[i]);
      E_der_v[l][s2][indw2 + 8 * k2 + 2 + 2 * i + 1] += real(pref * part_jk * b2_imag_diff[i]);
      E_der_v[l][s1][indw1 + 8 * j1 + 2 + 2 * i] += real(pref * part_jk * b1_real_diff[i]);
      E_der_v[l][s1][indw1 + 8 * j1 + 2 + 2 * i + 1] += real(pref * part_jk * b1_imag_diff[i]);
    }

    E_der_v[l][s2][indn2 + 2 * k2] += real(pref * conj(cj) * integral);
    E_der_v[l][s2][indn2 + 2 * k2 + 1] += real(pref * i_unit * conj(cj) * integral);
    E_der_v[l][s1][indn1 + 2 * j1] += real(pref * ck * integral);
    E_der_v[l][s1][indn1 + 2 * j1 + 1] += real(-pref * i_unit * ck * integral);

    if (norm_mode == NORMALIZE) { // add norm derivative contributions
      double t = -real(pref * part_jk * integral);
      // nonlocal terms: TODO: make a separate global loop for summation of nonlocal terms
      for (int j = 0; j < nspl[s1][c1]; j++) {
        for (int i = 0; i < 8; i++)
          E_der_v[l][s1][indw1 + 8 * j + i] += t * wf_norm_der[s1][indw1 + 8 * j + i];
        E_der_v[l][s1][indn1 + 2 * j] += t * wf_norm_der[s1][indn1 + 2 * j];
        E_der_v[l][s1][indn1 + 2 * j + 1] += t * wf_norm_der[s1][indn1 + 2 * j + 1];
      }
      for (int k = 0; k < nspl[s2][c2]; k++) {
        for (int i = 0; i < 8; i++)
          E_der_v[l][s2][indw2 + 8 * k + i] += t * wf_norm_der[s2][indw2 + 8 * k + i];
        E_der_v[l][s2][indn2 + 2 * k] += t * wf_norm_der[s2][indn2 + 2 * k];
        E_der_v[l][s2][indn2 + 2 * k + 1] += t * wf_norm_der[s2][indn2 + 2 * k + 1];
      }
    }
  }
}

void AWPMD_split::eterm_deriv_omp(vector<double> &E_der1, int ic1, int s1, int c1, int j1,
                                  vector<double> &E_der2, int ic2, int s2, int c2, int k2, cdouble pref,
                                  const OverlapDeriv &o, cdouble v, cdouble dv_aj_conj,
                                  cdouble dv_ak, cVector_3 dv_bj_conj, cVector_3 dv_bk) {
  cdouble cj(split_c[s1][ic1 + j1][0], split_c[s1][ic1 + j1][1]);
  cdouble ck(split_c[s2][ic2 + k2][0], split_c[s2][ic2 + k2][1]);
  int indw1 = 8 * ic1, indw2 = 8 * ic2;
  int indn1 = (nvar[s1] / 10) * 8 + 2 * ic1, indn2 = (nvar[s2] / 10) * 8 + 2 * ic2;
  cdouble part_jk = conj(cj) * ck;

  int M = 1; //(j1==k2 ? 1 : 2);

  // over a_k_re
  E_der2[indw2 + 8 * k2] += M * real(pref * part_jk * (o.da2_re() * v + o.I0 * dv_ak));
  // over a_k_im
  E_der2[indw2 + 8 * k2 + 1] += M * real(pref * part_jk * (o.da2_im() * v + i_unit * o.I0 * dv_ak));
  // over a_j_re
  E_der1[indw1 + 8 * j1] += M * real(pref * part_jk * (o.da1_re() * v + o.I0 * dv_aj_conj));
  // over a_j_im
  E_der1[indw1 + 8 * j1 + 1] += M * real(pref * part_jk * (o.da1_im() * v - i_unit * o.I0 * dv_aj_conj));

  for (int i = 0; i < 3; i++) {
    // over b_k_re
    E_der2[indw2 + 8 * k2 + 2 + 2 * i] += M * real(pref * part_jk * (o.db2_re(i) * v + o.I0 * dv_bk[i]));
    // over b_k_im
    E_der2[indw2 + 8 * k2 + 2 + 2 * i + 1] += M * real(pref * part_jk * (o.db2_im(i) * v + i_unit * o.I0 * dv_bk[i]));
    // over b_j_re
    E_der1[indw1 + 8 * j1 + 2 + 2 * i] += M * real(pref * part_jk * (o.db1_re(i) * v + o.I0 * dv_bj_conj[i]));
    // over b_j_im
    E_der1[indw1 + 8 * j1 + 2 + 2 * i + 1] +=
        M * real(pref * part_jk * (o.db1_im(i) * v - i_unit * o.I0 * dv_bj_conj[i]));
  }


  // over ck_re
  E_der2[indn2 + 2 * k2] += M * real(pref * conj(cj) * o.I0 * v);
  // over ck_im
  E_der2[indn2 + 2 * k2 + 1] += M * real(pref * i_unit * conj(cj) * o.I0 * v);
  // over cj_re
  E_der1[indn1 + 2 * j1] += M * real(pref * ck * o.I0 * v);
  // over cj_im
  E_der1[indn1 + 2 * j1 + 1] += M * real(-pref * i_unit * ck * o.I0 * v);

  if (norm_mode == NORMALIZE) { // add norm derivative contributions
    double t = -M * real(pref * part_jk * o.I0 * v);
    // nonlocal terms: TODO: make a separate global loop for summation of nonlocal terms
    for (int j = 0; j < nspl[s1][c1]; j++) {
      for (int i = 0; i < 8; i++)
        E_der1[indw1 + 8 * j + i] += t * wf_norm_der[s1][indw1 + 8 * j + i];
      E_der1[indn1 + 2 * j] += t * wf_norm_der[s1][indn1 + 2 * j];
      E_der1[indn1 + 2 * j + 1] += t * wf_norm_der[s1][indn1 + 2 * j + 1];
    }
    for (int k = 0; k < nspl[s2][c2]; k++) {
      for (int i = 0; i < 8; i++)
        E_der2[indw2 + 8 * k + i] += t * wf_norm_der[s2][indw2 + 8 * k + i];
      E_der2[indn2 + 2 * k] += t * wf_norm_der[s2][indn2 + 2 * k];
      E_der2[indn2 + 2 * k + 1] += t * wf_norm_der[s2][indn2 + 2 * k + 1];
    }
  }
}


void AWPMD_split::calc_norms(int flag, bool normalize) {
  if (!(calc_state & CALC_NORMS) || flag) { // only default calculation is checked


    if (flag & 0x4 || use_box) { // electron forces requested or box is used 
      for (int s1 = 0; s1 < 2; s1++) { // clearing norm derivatives
        for (int i = 0; i < nvar[s1]; i++) {
          wf_norm_der[s1][i] = 0;
          E_der[s1][i] = 0;
          F_extra[s1][i] = 0;
          F_wall[s1][i] = 0;
          ovl_der[s1][i] = 0;
        }
      }
    }
    if (!split_wp && approx == HARTREE) {
      for (int s1 = 0; s1 < 2; s1++) {
        for (int c1 = 0; c1 < ne[s1]; c1++) {
          wf_norm[s1][c1] = 1.;
        }
      }
      return;
    }

    // calculating block norms and derivatives
    for (int s1 = 0; s1 < 2; s1++) {
      int ic1 = 0; // starting index of the wp for current electron
      int indw1 = 0; // starting index of the electron wp coordinates
      int indn1 = (nvar[s1] / 10) * 8; // starting index of the electron norm coordinates

      for (int c1 = 0; c1 < ne[s1]; c1++) {

        // calculating the block norm
        wf_norm[s1][c1] = 0.;
        for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
          double cj_re = split_c[s1][ic1 + j1][0];
          double cj_im = split_c[s1][ic1 + j1][1];
          cdouble ccj = cdouble(cj_re, -cj_im);

          double part_jj = norm(ccj);
          wf_norm[s1][c1] += part_jj;
          WavePacket &wj = wp[s1][ic1 + j1];
          OverlapDeriv o;
          if (flag & (0x8 | 0x4)) { //electron forces needed
            wf_norm_der[s1][indn1 + 2 * j1] += 2 * cj_re;  // over cj_re
            wf_norm_der[s1][indn1 + 2 * j1 + 1] += 2 * cj_im; // over cj_im
            o.set1(wj);// conjugate: mu -> -mu, v -> -v !!!
          }


          for (int k1 = j1 + 1; k1 < nspl[s1][c1]; k1++) {
            double ck_re = split_c[s1][ic1 + k1][0];
            double ck_im = split_c[s1][ic1 + k1][1];
            cdouble ck = cdouble(ck_re, ck_im);

            WavePacket wk = wp[s1][ic1 + k1];
            if (pbc)
              move_to_image(wj, wk);
            WavePacket wjk = conj(wj) * wk;
            cdouble I0 = wjk.integral();

            cdouble part_jk = ccj * ck;
            wf_norm[s1][c1] += 2 * real(part_jk * I0);


            if (flag & (0x8 | 0x4)) { //electron forces needed
              o.set2(wk, &I0);


              wf_norm_der[s1][indw1 + 8 * k1] += 2 * real(part_jk * o.da2_re());           // over a_k_re
              wf_norm_der[s1][indw1 + 8 * k1 + 1] += 2 * real(part_jk * o.da2_im());    // over a_k_im

              wf_norm_der[s1][indw1 + 8 * j1] += 2 * real(part_jk * o.da1_re());             // over a_j_re
              wf_norm_der[s1][indw1 + 8 * j1 + 1] += 2 * real(part_jk * o.da1_im());    // over a_j_im

              for (int i = 0; i < 3; i++) {
                wf_norm_der[s1][indw1 + 8 * k1 + 2 + 2 * i] += 2 * real(part_jk * o.db2_re(i));           // over b_k_re
                wf_norm_der[s1][indw1 + 8 * k1 + 2 + 2 * i + 1] += 2 * real(part_jk * o.db2_im(i));  // over b_k_im

                wf_norm_der[s1][indw1 + 8 * j1 + 2 + 2 * i] += 2 * real(part_jk * o.db1_re(i));           // over b_j_re
                wf_norm_der[s1][indw1 + 8 * j1 + 2 + 2 * i + 1] += 2 * real(part_jk * o.db1_im(i));  // over b_j_im
              }

              wf_norm_der[s1][indn1 + 2 * j1] += 2 * real(ck * I0);  // over cj_re
              wf_norm_der[s1][indn1 + 2 * j1 + 1] += 2 * real(-i_unit * ck * I0);  // over cj_im
              wf_norm_der[s1][indn1 + 2 * k1] += 2 * real(ccj * I0); // over ck_re
              wf_norm_der[s1][indn1 + 2 * k1 + 1] += 2 * real(i_unit * ccj * I0); // over ck_im


              // overlap derivatives (for norm matrix)
              if (approx != HARTREE) {
                ovl_der[s1][indw1 + 8 * k1] += ccj * o.da2_re();           // over a_k_re
                ovl_der[s1][indw1 + 8 * k1 + 1] += ccj * o.da2_im();    // over a_k_im

                ovl_der[s1][indw1 + 8 * j1] += ck * o.da1_re();             // over a_j_re
                ovl_der[s1][indw1 + 8 * j1 + 1] += ck * o.da1_im();    // over a_j_im

                for (int i = 0; i < 3; i++) {
                  ovl_der[s1][indw1 + 8 * k1 + 2 + 2 * i] += ccj * o.db2_re(i);           // over b_k_re
                  ovl_der[s1][indw1 + 8 * k1 + 2 + 2 * i + 1] += ccj * o.db2_im(i);  // over b_k_im

                  ovl_der[s1][indw1 + 8 * j1 + 2 + 2 * i] += ck * o.db1_re(i);           // over b_j_re
                  ovl_der[s1][indw1 + 8 * j1 + 2 + 2 * i + 1] += ck * o.db1_im(i);  // over b_j_im
                }

                ovl_der[s1][indn1 + 2 * j1] += ck * I0;  // over cj_re
                ovl_der[s1][indn1 + 2 * j1 + 1] += -i_unit * ck * I0;  // over cj_im
                ovl_der[s1][indn1 + 2 * k1] += ccj * I0; // over ck_re
                ovl_der[s1][indn1 + 2 * k1 + 1] += i_unit * ccj * I0; // over ck_im
              }
            }
          } // k1
        }// j1
        if (flag & (0x8 | 0x4)) { //electron forces needed
          // normalizing the norm derivative
          for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
            for (int i = 0; i < 8; i++) // wp parameters
              wf_norm_der[s1][indw1 + 8 * j1 + i] /= 2 * wf_norm[s1][c1];
            // c
            wf_norm_der[s1][indn1 + 2 * j1] /= 2 * wf_norm[s1][c1];
            wf_norm_der[s1][indn1 + 2 * j1 + 1] /= 2 * wf_norm[s1][c1];
          }
        }
        //wf_norm[s1][c1]=1;
        //wf_norm[s1][c1]=sqrt(wf_norm[s1][c1]);
        ic1 += nspl[s1][c1];
        indw1 += 8 * nspl[s1][c1]; // 8 variables in each wavepacket
        indn1 += 2 * nspl[s1][c1]; // 2 variables in each wp norm
      }// c1
    }// s1
    valid_norms = true;
    calc_state |= CALC_NORMS;
  }

  if (normalize) {
    for (int s1 = 0; s1 < 2; s1++) {
      int ic1 = 0; // starting index of the wp for current electron
      for (int c1 = 0; c1 < ne[s1]; c1++) {
        double sqnrm = sqrt(wf_norm[s1][c1]);
        for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
          split_c[s1][ic1 + j1][0] /= sqnrm;
          split_c[s1][ic1 + j1][1] /= sqnrm;
        }
        wf_norm[s1][c1] = 1.;
        ic1 += nspl[s1][c1];
      }
    }
  } //normalize

}

void AWPMD_split::clear_forces(int flag, Vector_3P fi, Vector_3P fe_x,
                               Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c) {
  if (flag & 0x1) {
    for (int i = 0; i < ni; i++)
      fi[i] = Vector_3(0.);
  }
  if (flag & 0x4 && !(flag & 0x10)) { // electron forces requested in physical representation
    int iv1 = 0;
    for (int s1 = 0; s1 < 2; s1++) { // clearing forces
      for (int c1 = 0; c1 < ne[s1]; c1++) {
        for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
          fe_x[iv1 + j1] = Vector_3(0, 0, 0);
          fe_p[iv1 + j1] = Vector_3(0, 0, 0);
          fe_w[iv1 + j1] = 0;
          fe_pw[iv1 + j1] = 0;
          fe_c[iv1 + j1] = Vector_2(0, 0);
        }
        iv1 += nspl[s1][c1];
      }
    }
  }
}

void AWPMD_split::forces2phys() {
  // recalculating derivatives
  int iv1 = 0;
  for (int s1 = 0; s1 < 2; s1++) {
    int ic1 = 0; // starting index of the wp for current electron
    for (int c1 = 0; c1 < ne[s1]; c1++) {
      int indw1 = 8 * ic1;
      int indn1 = (nvar[s1] / 10) * 8 + 2 * ic1;
      for (int k1 = 0; k1 < nspl[s1][c1]; k1++) {
        WavePacket wk = wp[s1][ic1 + k1];
        /*double w=wk.get_width();
        Vector_3 r=wk.get_r();
        double t=3/(2*w*w*w);
        fe_w[ic1+k1]+= t*E_der[s1][indw1+8*k1]+imag(wk.a)*E_der[s1][indw1+8*k1+1]/w;
        fe_pw[ic1+k1]+=E_der[s1][indw1+8*k1+1]/(2*w*h_plank);
        for(int i=0;i<3;i++){
          fe_x[ic1+k1][i]+= -2*real(wk.a)*E_der[s1][indw1+8*k1+2+2*i]-2*imag(wk.a)*E_der[s1][indw1+8*k1+2+2*i+1];
          fe_p[ic1+k1][i]+= (-E_der[s1][indw1+8*k1+2+2*i+1])*(m_electron/h_plank); //*(h_plank/m_electron);
          fe_pw[ic1+k1]+=(r[i]*E_der[s1][indw1+8*k1+2+2*i+1]/w)/h_plank;
          fe_w[ic1+k1]+=2*r[i]*(t*E_der[s1][indw1+8*k1+2+2*i]+imag(wk.a)*E_der[s1][indw1+8*k1+2+2*i+1]/w);
        }*/

        /*
        if(have_force_arrays){ // TODO: update E_Der
          wk.int2phys_der< minus >(E_der[s1].begin()+indw1+8*k1,(double *)&fe_x[iv1+k1],(double *)&fe_p[iv1+k1],&fe_w[iv1+k1],&fe_pw[iv1+k1], 1./one_h);
          if(!c_polar_mode)
            fe_c[iv1+k1]+=-Vector_2(E_der[s1][indn1+2*k1],E_der[s1][indn1+2*k1+1]);
          else{ // converting to polar derivative
            Vector_2 fc(E_der[s1][indn1+2*k1],E_der[s1][indn1+2*k1+1]);
            Vector_2 &c=split_c[s1][ic1+k1];
            double rho=c.norm();
            fe_c[iv1+k1]+=-Vector_2((fc[0]*c[0]+fc[1]*c[1])/rho,-fc[0]*c[1]+fc[1]*c[0]);
          }
        }*/

        vector<double> *F_[3] = {&E_der[s1], &F_extra[s1], &F_wall[s1] };
        for (int i = 0; i < 3; i++) {
          vector<double>::iterator shift_x = F_[i]->begin() + indw1 + 8 * k1;
          vector<double>::iterator shift_p = shift_x + 3;
          vector<double>::iterator shift_w = shift_p + 3;
          vector<double>::iterator shift_pw = shift_w + 1;
          //wk.int2phys_der< eq_minus_second >(shift_x,shift_x,shift_p,shift_w,shift_pw, 1./one_h);
          wk.int2phys_der<eq_second>(shift_x, shift_x, shift_p, shift_w, shift_pw, 1. / one_h);
          // changing sign from dE/dc to force_c
          if (!c_polar_mode) {
            (*(F_[i]))[indn1 + 2 * k1] = (*(F_[i]))[indn1 + 2 * k1]; //-(*(F_[i]))[indn1+2*k1];
            (*(F_[i]))[indn1 + 2 * k1 + 1] = (*(F_[i]))[indn1 + 2 * k1 + 1];//-(*(F_[i]))[indn1+2*k1+1];
          } else { // converting to polar derivative
            Vector_2 fc((*(F_[i]))[indn1 + 2 * k1], (*(F_[i]))[indn1 + 2 * k1 + 1]);
            Vector_2 &c = split_c[s1][ic1 + k1];
            double rho = c.norm();
            (*(F_[i]))[indn1 + 2 * k1] = (fc[0] * c[0] + fc[1] * c[1]) / rho;//-(fc[0]*c[0]+fc[1]*c[1])/rho;
            (*(F_[i]))[indn1 + 2 * k1 + 1] = (-fc[0] * c[1] + fc[1] * c[0]);//-(-fc[0]*c[1]+fc[1]*c[0]);
          }
        }// i
      }// k1
      ic1 += nspl[s1][c1]; // incrementing block1 wp address
      iv1 += nspl[s1][c1]; // incrementing global variable address
    }// c1
  } // s1
}

void AWPMD_split::forces2phsy(packet_index_info packet_index, WavePacket const &packet,
                 double* eforce, double *erforce, double *ervforce){
  auto dx = &E_der[packet_index.spin][packet_index.index * 8];
  auto dp = dx + 3;
  auto dw = dp + 3;
  auto dpw = dw + 1;

  packet.int2phys_der<eq_second>(dx, dx, dp, dw, dpw, 1. / one_h);
  if (eforce)
    for (auto i = 0; i < 3; ++i)
      eforce[i] += dx[i];

  if (erforce)
    *erforce += *dw;

  if (ervforce)
    *ervforce += *dpw;
}

void AWPMD_split::get_el_forces(int flag, Vector_3P fe_x,
                                Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c) {
  bool have_force_arrays = (fe_x && fe_p && fe_w && fe_pw && fe_c);

  if ((flag & 0x4) && have_force_arrays) //need to replace the forces
    clear_forces(0x4, NULL, fe_x, fe_p, fe_w, fe_pw, fe_c);

  // recalculating derivatives
  if (flag & (0x8 | 0x4) && have_force_arrays) { //electron forces needed
    int iv1 = 0;
    for (int s1 = 0; s1 < 2; s1++) {
      int ics = s1 * nvar[0] / 10;
      int ic1 = 0; // starting index of the wp for current electron
      for (int c1 = 0; c1 < ne[s1]; c1++) {
        int indw1 = 8 * ic1;
        int indn1 = (nvar[s1] / 10) * 8 + 2 * ic1;
        for (int k1 = 0; k1 < nspl[s1][c1]; k1++) {
          for (int i = 0; i < 3; i++) // dx/dt
            fe_x[ics][i] = dq_dt[s1][indw1 + 8 * k1 + i];
          for (int i = 3; i < 6; i++) // dp/dt
            fe_p[ics][i - 3] = dq_dt[s1][indw1 + 8 * k1 + i];
          fe_w[ics] = dq_dt[s1][indw1 + 8 * k1 + 6]; // dw/dt
          fe_pw[ics] = dq_dt[s1][indw1 + 8 * k1 + 7]; // dpw/dt
          fe_c[ics][0] = dq_dt[s1][indn1 + 2 * k1]; // c
          fe_c[ics][1] = dq_dt[s1][indn1 + 2 * k1 + 1];
          ics++;  //=nspl[s1][c1];
        }// k1
        ic1 += nspl[s1][c1]; // incrementing block1 wp address
        iv1 += nspl[s1][c1]; // incrementing global variable address
      }// c1
    } // s1
  } // flags
}

//e same as interaction, but using Hartee factorization (no antisymmetrization)
int AWPMD_split::interaction_hartree(int flag, Vector_3P fi, Vector_3P fe_x,
                                     Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c) {

  // resize arrays if needed
  enum APPROX tmp = HARTREE;
  swap(tmp, approx); // do not neeed large matrices
  resize(flag);
  swap(tmp, approx);

  // clearing forces
  clear_forces(flag, fi, fe_x, fe_p, fe_w, fe_pw, fe_c);
  // calculate block norms and (optionally) derivatives
  calc_norms(flag);

  Wext = Eext = Eee = Ew = 0.;
  for (int s1 = 0; s1 < 2; s1++) {
    Ee[s1] = 0.;
    Eei[s1] = 0.;
    int ic1 = 0; // starting index of the wp for current electron

    for (int c1 = 0; c1 < ne[s1]; c1++) {
      // calculating single-electron quantities within block
      double Ee1 = 0., Ew1 = 0., Eei1 = 0.;
      double pref;
      double pref_ei;
      if (norm_mode == NORMALIZE) {  // divide by norms
        pref = -h2_me / (2 * wf_norm[s1][c1]); // ekin
        pref_ei = coul_pref / wf_norm[s1][c1];
      } else {
        pref = -h2_me / 2; // ekin
        pref_ei = coul_pref;
      }

      for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
        cdouble cj(split_c[s1][ic1 + j1][0], split_c[s1][ic1 + j1][1]);
        WavePacket wj = wp[s1][ic1 + j1];

        OverlapDeriv o;
        if (flag & (0x8 | 0x4)) //electron forces needed
          o.set1(wj);

        for (int k1 = j1; k1 < nspl[s1][c1]; k1++) {
          int M1a = (j1 == k1 ? 1 : 2);
          double M1e, M1f;
          _mytie(M1e, M1f) = check_part1(s1, ic1 + j1, s1, ic1 + k1) * M1a;


          cdouble ck(split_c[s1][ic1 + k1][0], split_c[s1][ic1 + k1][1]);

          // electrons kinetic energy
          WavePacket wk = wp[s1][ic1 + k1];
          if (pbc)
            move_to_image(wj, wk);

          WavePacket wjk = conj(wj) * wk;
          cdouble I0 = wjk.integral();
          cdouble part_jk = conj(cj) * ck;

          if (M1f) {
            cVector_3 v1 = conj(wj.b) * wk.a - wk.b * conj(wj.a);
            cdouble v = (v1 * v1) / wjk.a;
            v -= 6. * wk.a * conj(wj.a);
            v /= wjk.a;

            //v=1.;

            // kinetic energy contribution
            Ee[s1] += M1e * real(part_jk * I0 * v) * pref;  // Ejk+Ekj=Ejk+conj(Ejk), j!=k or Ejj

            if (flag & (0x8 | 0x4)) { //electron forces needed
              cVector_3 tv = wk.b * conj(wj.a) - conj(wj.b) * wk.a;
              cdouble ajk2 = wjk.a * wjk.a;
              cdouble ajk3 = ajk2 * wjk.a;
              cdouble dv_aj_conj = -2 * wk.a * (3 * wjk.a * wk.a - tv * wjk.b) / ajk3;
              cdouble dv_ak = -2 * conj(wj.a) * ((3 * wjk.a) * conj(wj.a) + tv * wjk.b) / ajk3;
              cVector_3 dv_bj_conj = (-2 * wk.a / ajk2) * tv;
              cVector_3 dv_bk = (2 * conj(wj.a) / ajk2) * tv;

              /*cdouble dv_aj_conj=0.;
              cdouble dv_ak=0.;
              cVector_3 dv_bj_conj;
              cVector_3 dv_bk;*/

              o.set2(wk, &I0);
              // calculate full derivative of the term pref*conj(cj)*ck*Ijk*vjk/sqrt(nrm(s1)*nrm(s2))
              // denominator must be included in pref
              eterm_deriv(ic1, s1, c1, j1, ic1, s1, c1, k1, M1f * pref, o, v, dv_aj_conj, dv_ak, dv_bj_conj, dv_bk);
            }
          }// M1e

          cVector_3 djk = wjk.b / (2. * wjk.a);
          // e-i energy
          //cdouble sum(0.,0.);
          if (calc_ei) {
            for (int i = 0; i < ni; i++) {  // ions loop
              double M1ie, M1if;
              _mytie(M1ie, M1if) = check_part1ei(s1, ic1 + j1, ic1 + k1, i) * M1a;
              if (!M1if)
                continue;

              cVector_3 gjki = djk - cVector_3(xi[i]);

              if (pbc) // correcting the real part (distance) according to PBC
                gjki = rcell1(gjki, cell, pbc);
              //-Igor- gkli=cVector_3(real(gkli).rcell1(cell,pbc),imag(gkli));

              cdouble ngjki = gjki.norm();
              cdouble sqajk = sqrt(wjk.a);
              //cdouble ttt = cerf_div(ngjki,c);
              cdouble v = cerf_div(ngjki, sqajk);
              double pref_eiq = pref_ei * qi[i] * (qe[s1][ic1 + j1] + qe[s1][ic1 + k1]) * 0.5;
              cdouble dE = pref_eiq * I0 * v * part_jk;

              // ions energy contribution
              Eei[s1] += M1ie * real(dE);
              //sum+=dE;

              cdouble t1, t2;
              bool zero_force = (fabs(real(ngjki)) + fabs(imag(ngjki)) < 1e-10);
              if (flag & (0x8 | 0x4 | 0x3)) { //some derivatives needed
                cdouble arg = ngjki * sqajk;
                t1 = two_over_sqr_pi * exp(-arg * arg);
                t2 = t1 / (2 * sqajk);
                t1 *= sqajk;
              }

              if (flag & 0x3 && !zero_force) {// calculate forces on ions
                cdouble dEw = pref_eiq * I0 * t1 * part_jk;
                dEw = (dE - dEw) / ngjki;
                //Vector_3 dir=-real(gjki);
                //dir.normalize();
                //fi[i]+=M1if*real( dEw)*dir;
                fi[i] += M1if * real(dEw * gjki / ngjki);
              }
              if (flag & (0x8 | 0x4)) { //electron forces needed
                cdouble dv_ak = (zero_force ? 0. : (djk * gjki) * (v - t1) / (wjk.a * ngjki * ngjki)) + t2;
                cVector_3 dv_bk = zero_force ? cVector_3() : -gjki * (v - t1) / (2 * wjk.a * ngjki * ngjki);
                eterm_deriv(ic1, s1, c1, j1, ic1, s1, c1, k1, M1if * pref_eiq, o, v, dv_ak, dv_ak, dv_bk, dv_bk);
              }
            }
          }
          // extra constraint energy
          if (j1 == k1 && constraint == HARM && M1e) {
            cdouble v = conj(wj.a) * wj.a;
            double harm_pref = norm_mode == NORMALIZE ? harm_w0_4 / wf_norm[s1][c1] : harm_w0_4;
            Ew += harm_pref * real(part_jk * v);
            if (flag & (0x8 | 0x4)) { //electron forces needed
              cdouble dv_ak = conj(wj.a);
              cdouble dv_aj_conj = wj.a;
              eterm_deriv(ic1, s1, c1, j1, ic1, s1, c1, k1, harm_pref, o, v, dv_aj_conj,
                          dv_ak, cVector_3(), cVector_3());
            }
          }
# if 1
          // second block
          // e-e interaction
          if (calc_ee) {
            for (int s2 = s1; s2 < 2; s2++) {
              int ic2 = 0; // starting index of the wp for current electron

              for (int c2 = 0; c2 < ne[s2]; ic2 += nspl[s2][c2], c2++) { // incrementing block2 wp address
                if (s1 == s2 && c2 <= c1)
                  continue;

                double pref_ee = norm_mode == NORMALIZE ? coul_pref / (wf_norm[s1][c1] * wf_norm[s2][c2]) : coul_pref;
                double dE = 0.;
                for (int j2 = 0; j2 < nspl[s2][c2]; j2++) {
                  double M2ej, M2fj;
                  _mytie(M2ej, M2fj) = check_part1(s1, ic1 + j1, s2, ic2 + j2);
                  if (!M2fj)
                    continue;


                  cdouble cj2(split_c[s2][ic2 + j2][0], split_c[s2][ic2 + j2][1]);
                  WavePacket &wj2 = wp[s2][ic2 + j2];

                  OverlapDeriv o2;
                  if (flag & (0x8 | 0x4)) //electron forces needed
                    o2.set1(wj2);

                  for (int k2 = j2; k2 < nspl[s2][c2]; k2++) {
                    double M2ek, M2fk;
                    _mytie(M2ek, M2fk) = check_part1(s1, ic1 + k1, s2, ic2 + k2);
                    if (!M2fk)
                      continue;


                    int M2a = (j2 == k2 ? 1 : 2);



                    cdouble ck2(split_c[s2][ic2 + k2][0], split_c[s2][ic2 + k2][1]);

                    WavePacket wk2 = wp[s2][ic2 + k2];
                    if (pbc)
                      move_to_image(wj2, wk2);
                    WavePacket wjk2 = conj(wj2) * wk2;
                    cdouble I02 = wjk2.integral();


                    cdouble part_jk2 = conj(cj2) * ck2;

                    cVector_3 djk2 = wjk2.b / (2 * wjk2.a);
                    cVector_3 ddv = djk - djk2;
                    cdouble dd = ddv.norm();
                    cdouble aa = 1. / sqrt(1. / wjk.a + 1. / wjk2.a);
                    //double ww1=wj.get_width();
                    //double ww2=wj2.get_width();
                    cdouble v = cerf_div(dd, aa);
                    double pref_eeq =
                        pref_ee * (qe[s1][ic1 + j1] + qe[s1][ic1 + k1]) * (qe[s2][ic2 + j2] + qe[s2][ic2 + k2]) * 0.25;
                    cdouble Vj1j2k1k2 = pref_eeq * I0 * I02 * v * part_jk * part_jk2;
                    Eee += M1e * M2ej*M2ek * real(Vj1j2k1k2);
                    //Eee+=(j1==k1 || j2==k2 ? 1: 2)*real(Vj1j2k1k2);


                    cdouble t1, t2;
                    bool zero_force = (fabs(real(dd)) + fabs(imag(dd)) < 1e-10);
                    if (flag & (0x8 | 0x4)) { //electron forces needed
                      cdouble arg = dd * aa;
                      t1 = two_over_sqr_pi * exp(-arg * arg) * aa;
                      t2 = t1 * aa * aa / 2;


                      cdouble dv_ak1 =
                          (zero_force ? 0. : (djk * ddv) * (v - t1) / (wjk.a * dd * dd)) + t2 / (wjk.a * wjk.a);
                      cVector_3 dv_bk1 = zero_force ? cVector_3() : -ddv * (v - t1) / (2 * wjk.a * dd * dd);
                      eterm_deriv(ic1, s1, c1, j1, ic1, s1, c1, k1, M1f * M2fj *M2fk * pref_eeq * I02 * part_jk2,
                                  o, v, dv_ak1, dv_ak1, dv_bk1, dv_bk1);

                      o2.set2(wk2, &I02);
                      cdouble dv_ak2 =
                          (zero_force ? 0. : (-djk2 * ddv) * (v - t1) / (wjk2.a * dd * dd)) + t2 / (wjk2.a * wjk2.a);
                      cVector_3 dv_bk2 = zero_force ? cVector_3() : ddv * (v - t1) / (2 * wjk2.a * dd * dd);
                      eterm_deriv(ic2, s2, c2, j2, ic2, s2, c2, k2, M1f * M2fj *M2fk * pref_eeq * I0 * part_jk,
                                  o2, v, dv_ak2, dv_ak2, dv_bk2, dv_bk2);
                    }

                  }// k2
                } // j2

              } //c2
            }// s2
          }
# endif
        } // k1
      }// j1
      ic1 += nspl[s1][c1]; // incrementing block1 wp address
    }// c1
  } // s1

  // transforming the forces to physical coordinates
  if (flag & (0x8 | 0x4) && !(flag & 0x10))
    get_el_forces((flag & (~0x4)) | 0x8, fe_x, fe_p, fe_w, fe_pw,
                  fe_c); // flag change: electronic forces were cleared already

  Eii = 0.;
  if (calc_ii)
    interaction_ii(flag, fi);
  if (use_box)
    calc_ion_confinement(flag, fi);

  return 1;
}

#define HARTREE_FUNC //enable function from AWPMD for pair interaction

std::pair<double, double>
AWPMD_split::interaction_electron_kinetic(WavePacket const &packet, int spin, double *erforce, double *ervfroce) {
#ifdef HARTREE_FUNC
  return AWPMD::interaction_electron_kinetic(packet, spin, erforce, ervfroce);
#else
  auto pref = -h2_me / 2.0;
  WavePacket const &wk = packet;
  WavePacket const &wj = packet;

  WavePacket wjk = conj(wk) * wj;
  cdouble I0 = wjk.integral();
  cdouble part_jk {1.0, 0.0};

  cVector_3 v1 = conj(packet.b) * packet.a - packet.b * conj(packet.a);
  cdouble v = (v1 * v1) / wjk.a;
  v -= 6. * packet.a * conj(packet.a);
  v /= wjk.a;

  auto ke = packet.get_p().norm2() * (-pref);
  auto energy = real(part_jk * I0 * v) * pref - ke;
  Ee[spin] += energy;

  if (erforce){
    OverlapDeriv o;
    o.set1(wj);
    o.set2(wk, &I0);

    cVector_3 tv = wk.b * conj(wj.a) - conj(wj.b) * wk.a;
    cdouble ajk2 = wjk.a * wjk.a;
    cdouble ajk3 = ajk2 * wjk.a;
    cdouble dv_aj_conj = -2 * wk.a * (3 * wjk.a * wk.a - tv * wjk.b) / ajk3;
    cdouble dv_ak = -2 * conj(wj.a) * ((3 * wjk.a) * conj(wj.a) + tv * wjk.b) / ajk3;
    cVector_3 dv_bj_conj = (-2 * wk.a / ajk2) * tv;
    cVector_3 dv_bk = (2 * conj(wj.a) / ajk2) * tv;

    eterm_deriv({0, spin}, {0, spin}, pref, o, v, dv_aj_conj, dv_ak, dv_bj_conj, dv_bk);
    forces2phsy({0, spin}, packet, nullptr, erforce, ervfroce);
  }
  return energy;
#endif
}

double AWPMD_split::interaction_ee_single(WavePacket const &packet_1,
                                          WavePacket const &packet_2,
                                          double **eforce, double **erforce,
                                          double cutoff) {
#ifdef HARTREE_FUNC
  return AWPMD::interaction_ee_single(packet_1, packet_2, eforce, erforce,
                                      cutoff);
#else
  WavePacket const &wj = packet_1;
  WavePacket const &wk = packet_1;
  WavePacket wjk = conj(wj) * wk;
  cVector_3 djk = wjk.b / (2. * wjk.a);
  double pref_ee = coul_pref;
  cdouble I0 = wjk.integral();
  cdouble part_jk = {1.0, 0.0};

  WavePacket const &wj2 = packet_2;
  WavePacket const &wk2 = packet_2;

  WavePacket wjk2 = conj(wj2) * wk2;
  cdouble I02 = wjk2.integral();
  cdouble part_jk2 = {1.0, 0.0};
  cVector_3 djk2 = wjk2.b / (2 * wjk2.a);
  cVector_3 ddv = djk - djk2;
  cdouble dd = ddv.norm();
  cdouble aa = 1. / sqrt(1. / wjk.a + 1. / wjk2.a);
  cdouble v = cerf_div(dd, aa);

  double pref_eeq =
      pref_ee * (-1.0 + -1.0) * (-1.0 + -1.0) * 0.25;
  cdouble Vj1j2k1k2 = pref_eeq * I0 * I02 * v * part_jk * part_jk2;

  auto energy = real(Vj1j2k1k2);
  if (eforce && erforce){
    OverlapDeriv o;
    o.set1(wj);
    o.set2(wk, &I0);

    OverlapDeriv o2;
    o2.set1(wj2);

    cdouble arg = dd * aa;
    auto t1 = two_over_sqr_pi * exp(-arg * arg) * aa;
    auto t2 = t1 * aa * aa / 2;

    bool zero_force = (fabs(real(dd)) + fabs(imag(dd)) < 1e-10);
    cdouble dv_ak1 =
        (zero_force ? 0. : (djk * ddv) * (v - t1) / (wjk.a * dd * dd)) + t2 / (wjk.a * wjk.a);
    cVector_3 dv_bk1 = zero_force ? cVector_3() : -ddv * (v - t1) / (2 * wjk.a * dd * dd);
    eterm_deriv({0, 0}, {0, 0}, pref_eeq * I02 * part_jk2,
                o, v, dv_ak1, dv_ak1, dv_bk1, dv_bk1);
    forces2phsy({0, 0}, packet_1, eforce[0], erforce[0], nullptr);

    o2.set2(wk2, &I02);
    cdouble dv_ak2 =
        (zero_force ? 0. : (-djk2 * ddv) * (v - t1) / (wjk2.a * dd * dd)) + t2 / (wjk2.a * wjk2.a);
    cVector_3 dv_bk2 = zero_force ? cVector_3() : ddv * (v - t1) / (2 * wjk2.a * dd * dd);
    eterm_deriv({0, 0}, {0, 0}, pref_eeq * I0 * part_jk,
                o2, v, dv_ak2, dv_ak2, dv_bk2, dv_bk2);
    forces2phsy({0, 0}, packet_2, eforce[1], erforce[1], nullptr);
  }
  return energy;
#endif
}

double AWPMD_split::interaction_ei_single(double const *x, double q,
                                          WavePacket const &packet, int spin,
                                          double *f, double *eforce,
                                          double *erforce, double cutoff) {
#ifdef HARTREE_FUNC
  return AWPMD::interaction_ei_single(x, q, packet, spin, f, eforce, erforce,
                                      cutoff);
#else

  WavePacket const &wj = packet;
  WavePacket const &wk = packet;

  WavePacket wjk = conj(wj) * wk;
  cVector_3 djk = wjk.b / (2. * wjk.a);

  cVector_3 gjki = djk - cVector_3(x);
  cdouble ngjki = gjki.norm();
  cdouble sqajk = sqrt(wjk.a);

  cdouble v = cerf_div(ngjki, sqajk);
  auto pref_ei = coul_pref;
  double pref_eiq = pref_ei * q * -1.0;

  cdouble I0 = wjk.integral();
  cdouble part_jk = {1.0, 0.0};
  cdouble dE = pref_eiq * I0 * v * part_jk;

  auto energy = real(dE);
  Eei[spin] += energy;

  if (f != nullptr){
    cdouble arg = ngjki * sqajk;
    auto t1 = two_over_sqr_pi * exp(-arg * arg);
    auto t2 = t1 / (2.0 * sqajk);
    t1 *= sqajk;

    cdouble dEw = pref_eiq * I0 * t1 * part_jk;
    dEw = (dE - dEw) / ngjki;
    auto F = real(-dEw * gjki / ngjki);
    for (auto k = 0; k < 3; ++k)
      f[k] += F[k];

    if (eforce && erforce){
      bool zero_force = (fabs(real(ngjki)) + fabs(imag(ngjki)) < 1e-10);
      cdouble dv_ak = (zero_force ? 0. : (djk * gjki) * (v - t1) / (wjk.a * ngjki * ngjki)) + t2;
      cVector_3 dv_bk = zero_force ? cVector_3() : -gjki * (v - t1) / (2 * wjk.a * ngjki * ngjki);

      OverlapDeriv o;
      o.set1(wj);
      o.set2(wk, &I0);

      eterm_deriv({0, spin}, {0, spin}, pref_eiq, o, v, dv_ak, dv_ak, dv_bk, dv_bk);
      forces2phsy({0, spin}, packet, eforce, erforce, nullptr);
      for (auto k = 0; i < 3; ++k)
        eforce[k] *= -1.0;
    }
  }

  return energy;
#endif
}

///\en Calcualtes the overlap between two electrons taking all split WPs into account.
///    Norms must be pre-calculated.
cdouble AWPMD_split::overlap(int ic1, int s1, int c1, int ic2, int s2, int c2) {
  cdouble sum(0, 0);
  for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
    double cj_re = split_c[s1][ic1 + j1][0];
    double cj_im = split_c[s1][ic1 + j1][1];
    cdouble ccj = cdouble(cj_re, -cj_im);
    WavePacket wj = wp[s1][ic1 + j1];
    for (int k2 = 0; k2 < nspl[s2][c2]; k2++) {
      double ck_re = split_c[s2][ic2 + k2][0];
      double ck_im = split_c[s2][ic2 + k2][1];
      cdouble ck = cdouble(ck_re, ck_im);

      WavePacket wk = wp[s2][ic2 + k2];
      if (pbc)
        move_to_image(wj, wk);
      WavePacket wjk = conj(wj) * wk;
      sum += ccj * ck * wjk.integral();
    }
  }
  if (norm_mode == NORMALIZE)
    return sum / sqrt(wf_norm[s1][c1] * wf_norm[s2][c2]);
  else
    return sum;
}


/// adds the derivatives of Y in the term v*Y[s](c2,c1)
void AWPMD_split::y_deriv(cdouble v, int s, int c2, int c1, int external) {
  //return ;
  vector<double> *E_der = external ? AWPMD_split::F_extra : AWPMD_split::E_der;

  int indn = (nvar[s] / 10) * 8;
  int ic = 0;
  for (int c = 0; c < ne[s]; c++) {
    for (int j = 0; j < nspl[s][c]; j++) {
      for (int k = 0; k < 8; k++)
        E_der[s][8 * ic + k] += real(
            (-M[s](c2, 8 * ic + k) * Y[s](c, c1) - Y[s](c2, c) * conj(M[s](c1, 8 * ic + k))) * v);
      for (int k = 0; k < 2; k++)
        E_der[s][indn + 2 * ic + k] += real(
            (-M[s](c2, indn + 2 * ic + k) * Y[s](c, c1) - Y[s](c2, c) * conj(M[s](c1, indn + 2 * ic + k))) * v);
      ic++;
    }
  }
}

int fnum = 0;

int AWPMD_split::calc_overlaps() {
  if (calc_state & CALC_OVERLAPS) // already done
    return 0;

  if (approx == HARTREE && !split_wp) // not needed for this config
    return 0;
  calc_norms(0);


  for (int s = 0; s < 2; s++) {
    ovl_det_log[s] = 0.;
    int nes = ne[s];
    if (nes == 0) continue;

    int ik = 0;
    for (int k = 0; k < nes; k++) {

      Y[s].set(k, k, 1.);  // Diagonal elements (=1)
      Oflg[s](k, k) = 1;
      int il = 0;
      for (int l = 0/*k+1*/; l < nes; il += nspl[s][l], l++) { // incrementing block2 wp address
        if (l < k + 1)
          continue;
        cdouble Okl = overlap(ik, s, k, il, s, l);
        Y[s].set(k, l, Okl);
        double ovlnorm = norm(Okl);
        Oflg[s](k, l) = ovlnorm > ovl_tolerance;
      }
      ik += nspl[s][k]; // incrementing block1 wp address
    }
    O[s] = Y[s];  // save overlap matrix

# if 0
    // normalizing overlaps
    ik=0;
    for(int k=0;k<nes;k++){
      int il=0;
      for(int l=0;l<nes;il+=nspl[s][l],l++){ // incrementing block2 wp address
        if(l<k+1)
          continue;
        double ovlnorm = norm(Y[s](k,l));
        if(fabs(1.-ovlnorm)<ovl_degeneracy_min){
          Y[s].set(k,l,1. - (1.- Y[s](k,l))/sqrt(norm(1.- Y[s](k,l))/ovl_degeneracy_min) );
        }
      }
      ik+=nspl[s][k]; // incrementing block1 wp address
    }
# endif
    //3. inverting the overlap matrix
    int info = 0;
    if (nes && approx != HARTREE) {
      /*FILE *f1=fopen(fmt("matrO_%d.d",s),"wt");
      fileout(f1,Y[s],"%15g");
      fclose(f1);8*/

      ZPPTRF("L", &nes, Y[s].arr, &info);
      // analyze return code here
      if (info < 0)
        return LOGERR(info, fmt_iv("AWPMD.calc_overlaps: call to ZPTRF failed (exitcode %d)!", info), LINFO);


      cdouble rr = 0.;
      for (int i = 0; i < nes; i++)
        rr += log(Y[s](i, i));
      ovl_det_log[s] += real(rr);
      //if(ovl_det_log<-1.5)
      //ovl_det_log = real(rr);

      ZPPTRI("L", &nes, Y[s].arr, &info);
      if (info < 0)
        return LOGERR(info, fmt_iv("AWPMD.calc_overlaps: call to ZPTRI failed (exitcode %d)!", info), LINFO);


      /*f1=fopen(fmt("matrY_%d.d",s),"wt");
      fileout(f1,Y[s],"%15g");
      fclose(f1);*/
    }
  }
# if 0
  // output the matrices
  FILE *f = fopen(fmt("matr%04d.d",fnum++),"wt");
  if(!f)
    return LOGERR(1,"can't write matr file\n",LINFO);
  fprintf(f,"#1-i 2-j 3-O(0,ij) 4-Y(0,ij) 5-O(1,ij) 6-Y(0,ij)\n");
  for(int i =0; i<ne[0]; i++){
    for(int j=0;j<ne[0];j++){
      fprintf(f,"%d %d %g %g %g %g\n",i,j, abs(O[0](i,j)), abs(Y[0](i,j)), abs(O[1](i,j)), abs(Y[1](i,j)));
    }
    fprintf(f,"\n");
  }
  fclose(f);
# endif
  calc_state |= CALC_OVERLAPS;
  return 1;
}


void AWPMD_split::calc_LMN_s(int s1) {
  if (approx != HARTREE) {
    M[s1].Set(0.);
    L[s1].Set(0.);
    if (norm_needed)
      Norm[s1].Set(0.);
  } else
    Lh[s1].assign(nvar[s1], 0.);

  int indn = (nvar[s1] / 10) * 8;
  int ic1 = 0;
  for (int c1 = 0; c1 < ne[s1]; c1++) {
    // L[s](c1,c1*)=0
    /*for(int j=0;j<nspl[s1][c1];j++){
      for(int i=0;i<10;i++)
        L[s1](c1,10*(ic1+j)+i)=0.;
    }*/

    int ic2 = 0;
    for (int c2 = 0; c2 < ne[s1]; ic2 += nspl[s1][c2], c2++) {
      if (approx == HARTREE && c1 != c2) // only diagonals for Hartree
        continue;
      if (c2 < c1) // taken assymmetry into account
        continue;
      if (!check_overlap(s1,c1,c2))
        continue; // non-overlapping WPs


      if (approx == HARTREE && norm_needed) // initializing the matrices with zero
        Normh[s1][c1].Set(0);


      double sq_norm12 = norm_mode == NORMALIZE ? sqrt(wf_norm[s1][c1] * wf_norm[s1][c2]) : 1.;
      double pref = 1. / (sq_norm12);
      cdouble ovl = 0.;

      // WP blocks:
      for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
        cdouble cj1(split_c[s1][ic1 + j1][0], split_c[s1][ic1 + j1][1]);
        WavePacket wj1 = wp[s1][ic1 + j1];

        OverlapDeriv o12;
        o12.set1(wj1);

        for (int k2 = (c1 == c2 ? j1 : 0); k2 < nspl[s1][c2]; k2++) {
          double M12 = (c1 == c2 && j1 == k2) ? 0.5 : 1.;
          //double M12n= (c1==c2 && j1==k2) ? 1. : 2.;

          cdouble ck2(split_c[s1][ic2 + k2][0], split_c[s1][ic2 + k2][1]);

          // electrons kinetic energy
          WavePacket wk2 = wp[s1][ic2 + k2];
          if (pbc)
            move_to_image(wj1, wk2);

          WavePacket wjk12 = conj(wj1) * wk2;
          cdouble I012 = wjk12.integral();
          cdouble part_jk12 = conj(cj1) * ck2;

          ovl += 2 * M12 * pref * part_jk12 * I012;

          o12.set2(wk2, &I012);


          cdouble der_k[10], der_j[10];
          // over a_k_re
          der_k[0] = part_jk12 * o12.da2_re();
          // over a_k_im
          der_k[1] = part_jk12 * o12.da2_im();
          // over a_j_re
          der_j[0] = part_jk12 * o12.da1_re();
          // over a_j_im
          der_j[1] = part_jk12 * o12.da1_im();

          for (int i = 0; i < 3; i++) {
            // over b_k_re
            der_k[2 + 2 * i] = part_jk12 * o12.db2_re(i);
            // over b_k_im
            der_k[2 + 2 * i + 1] = part_jk12 * o12.db2_im(i);
            // over b_j_re
            der_j[2 + 2 * i] = part_jk12 * o12.db1_re(i);
            // over b_j_im
            der_j[2 + 2 * i + 1] = part_jk12 * o12.db1_im(i);
          }

          // over ck_re
          der_k[8] = conj(cj1) * I012;
          // over ck_im
          der_k[9] = i_unit * conj(cj1) * I012;
          // over cj_re
          der_j[8] = ck2 * I012;
          // over cj_im
          der_j[9] = -i_unit * ck2 * I012;
# if 0
          //if(j1==0){ // add all together instead of adding by parts
            //cdouble t=-O[s1](c1,c2)/pref;   //-part_jk12*I012;
            cdouble t=-part_jk12*I012;
            for(int i=0;i<8;i++){
              der_j[i]+=t*wf_norm_der[s1][8*(ic1+j1)+i];
              der_k[i]+=t*wf_norm_der[s1][8*(ic2+k2)+i];
            }
            der_j[8]+=t*wf_norm_der[s1][indn+2*(c1+j1)];
            der_j[9]+=t*wf_norm_der[s1][indn+2*(c1+j1)+1];
            der_k[8]+=t*wf_norm_der[s1][indn+2*(c2+k2)];
            der_k[9]+=t*wf_norm_der[s1][indn+2*(c2+k2)+1];
          //}
# endif
          if (approx != HARTREE) {
            for (int i = 0; i < 8; i++) {
              L[s1](c1, 8 * (ic2 + k2) + i) += M12 * pref * der_k[i];
              L[s1](c2, 8 * (ic1 + j1) + i) += M12 * pref * conj(der_j[i]);
            }
            for (int i = 0; i < 2; i++) {
              L[s1](c1, indn + 2 * (ic2 + k2) + i) += M12 * pref * der_k[8 + i];
              L[s1](c2, indn + 2 * (ic1 + j1) + i) += M12 * pref * conj(der_j[8 + i]);
            }
          } else { // HARTREE
            for (int i = 0; i < 8; i++) {
              Lh[s1][8 * (ic2 + k2) + i] += M12 * pref * der_k[i];
              Lh[s1][8 * (ic1 + j1) + i] += M12 * pref * conj(der_j[i]);
            }
            for (int i = 0; i < 2; i++) {
              Lh[s1][indn + 2 * (ic2 + k2) + i] += M12 * pref * der_k[8 + i];
              Lh[s1][indn + 2 * (ic1 + j1) + i] += M12 * pref * conj(der_j[8 + i]);
            }
          }
          if (norm_needed) { // filling part of norm matrix
            o12.calc_der_overlap(false, conj(cj1), ck2);
            if (approx != HARTREE) {
              for (int i = 0; i < 10; i++) {  // 10x10
                int indi = i < 8 ? 8 * (ic1 + j1) + i : indn + 2 * (ic1 + j1) + i - 8;
                for (int j = 0; j < 10; j++) {
                  int indj = j < 8 ? 8 * (ic2 + k2) + j : indn + 2 * (ic2 + k2) + j - 8;
                  Norm[s1](indi, indj) = -2 * imag(Y[s1](c2, c1) * o12.IDD((int) i, (int) j)) * pref / one_h;
                }
              }
            } else {
# if 1
              for (int i = 0; i < 10; i++)  // 10x10
                for (int j = 0; j < 10; j++)
                  Normh[s1][c1](10 * j1 + i, 10 * k2 + j) = -2 * imag(o12.IDD((int) i, (int) j)) * pref / one_h;
              //Normh[s1][c1](10*j1+i,10*k2+j)=-2*imag(o12.IDD((int)i,(int)j))*pref/one_h;
# endif
            }
          }
        } // k2
      } //j1

      if (norm_mode == NORMALIZE) { // taking norm derivative into account
        // adding overlap derivative
        for (int k = 0; k < nspl[s1][c1]; k++) { // add all together instead of adding by parts
          cdouble t = -O[s1](c1, c2);   //-part_jk12*I012;
          if (approx != HARTREE) {
            for (int i = 0; i < 8; i++) {
              L[s1](c1, 8 * (ic2 + k) + i) += t * wf_norm_der[s1][8 * (ic2 + k) + i];
              if (c1 != c2)
                L[s1](c2, 8 * (ic1 + k) + i) += conj(t) * wf_norm_der[s1][8 * (ic1 + k) + i];
            }
            for (int i = 0; i < 2; i++) {
              L[s1](c1, indn + 2 * (ic2 + k) + i) += t * wf_norm_der[s1][indn + 2 * (ic2 + k) + i];
              if (c1 != c2)
                L[s1](c2, indn + 2 * (ic1 + k) + i) += conj(t) * wf_norm_der[s1][indn + 2 * (ic1 + k) + i];
            }
          } else { //Hartree
            for (int i = 0; i < 8; i++) {
              Lh[s1][8 * (ic2 + k) + i] += t * wf_norm_der[s1][8 * (ic2 + k) + i];
              if (c1 != c2)
                Lh[s1][8 * (ic1 + k) + i] += conj(t) * wf_norm_der[s1][8 * (ic1 + k) + i];
            }
            for (int i = 0; i < 2; i++) {
              Lh[s1][indn + 2 * (ic2 + k) + i] += t * wf_norm_der[s1][indn + 2 * (ic2 + k) + i];
              if (c1 != c2)
                Lh[s1][indn + 2 * (ic1 + k) + i] += conj(t) * wf_norm_der[s1][indn + 2 * (ic1 + k) + i];
            }
          }
        }
      } // NORMALIZE
    }// c2
    ic1 += nspl[s1][c1]; // incrementing block1 wp address
  }// c1
  if (approx != HARTREE) { // calculating M matrix
    // M=Y*L
    for (int i = 0; i < ne[s1]; i++)  // matrix product
      for (int j = 0; j < nvar[s1]; j++)
        for (int g = 0; g < ne[s1]; g++)
          M[s1](i, j) += Y[s1](i, g) * L[s1](g, j);
  }

  if (norm_needed) {  // filling the rest of norm_matrix
    if (approx != HARTREE) {
      int ic1 = 0;
      for (int c1 = 0; c1 < ne[s1]; c1++) {
        int ic2 = 0;
        for (int c2 = 0; c2 < ne[s1]; ic2 += nspl[s1][c2], c2++) {
          if (c2 < c1)
            continue;
          for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
            for (int j = 0; j < 10; j++) {
              int indj = j < 8 ? 8 * (ic1 + j1) + j : indn + 2 * (ic1 + j1) + j - 8;
              for (int k2 = (c1 == c2 ? j1 : 0); k2 < nspl[s1][c2]; k2++) {
                for (int k = 0; k < 10; k++) {
                  int indk = k < 8 ? 8 * (ic2 + k2) + k : indn + 2 * (ic2 + k2) + k - 8;
                  cdouble y = 0.;
                  for (int g = 0; g < ne[s1]; g++)
                    y += conj(L[s1](g, indj)) * M[s1](g, indk);
                  cdouble elm = -y;
                  if (norm_mode == NORMALIZE)
                    elm += -conj(L[s1](c2, indj)) * wf_norm_der[s1][indk] - L[s1](c1, indk) * wf_norm_der[s1][indj] -
                           O[s1](c1, c2) * wf_norm_der[s1][indk] * wf_norm_der[s1][indj];
                  elm *= Y[s1](c2, c1);
                  Norm[s1](indj, indk) += -2 * imag(elm) / one_h;
                } //k
              }// k2
            } //j
          } // j1
        }// c2
        ic1 += nspl[s1][c1]; // incrementing block1 wp address
      }
# if 0
      // filling the lower triangle according to antisymmetry
      for(int i=0;i<Norm[s1].size;i++){
        for(int j=0;j<i;j++){
          Norm[s1](i,j)=-Norm[s1](j,i);
        }
        Norm[s1](i,i) =  0.; // making exact
      }
# endif
    } else if (norm_mode == NORMALIZE) { // HARTREE AND NORMALIZE
      int ic1 = 0;
#if 1
      for (int c1 = 0; c1 < ne[s1]; c1++) {
        for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
          for (int j = 0; j < 10; j++) {
            int indj = j < 8 ? 8 * (ic1 + j1) + j : indn + 2 * (ic1 + j1) + j - 8;
            for (int k2 = j1; k2 < nspl[s1][c1]; k2++) {
              for (int k = 0; k < 10; k++) {
                int indk = k < 8 ? 8 * (ic1 + k2) + k : indn + 2 * (ic1 + k2) + k - 8;

                //cdouble y=conj(Lh[s1][indj])*Lh[s1][indk];  // WRONG TERM!
                //cdouble elm=-y;

                cdouble elm = -conj(Lh[s1][indj]) * wf_norm_der[s1][indk] - Lh[s1][indk] * wf_norm_der[s1][indj] -
                              wf_norm_der[s1][indk] * wf_norm_der[s1][indj];

                Normh[s1][c1](10 * j1 + j, 10 * k2 + k) += -2 * imag(elm) / one_h;
              } //k
            }// k2
          }//j
        }//j1
        ic1 += nspl[s1][c1]; // incrementing block1 wp address
      }//c1
# endif
      // filling lower triangle: making fully antisymmetric
# if 1
      for (int c = 0; c < ne[s1]; c++) {
        int nvar = 10 * nspl[s1][c];
        for (int j = 0; j < nvar; j++)
          for (int k = j + 1; k < nvar; k++)
            Normh[s1][c](k, j) = -Normh[s1][c](j, k);
        //if(nspl[s1][c]==1){ // special case: adding diagonal elements for arbitraty c, to be able to invert norm matrix
        //  Normh[s1][c](8,8)=Normh[s1][c](9,9)=-1;
        //}
      }
# endif
    }// harttree
    norm_matrix_state[s1] = NORM_CALCULATED;
  } // norm
}

void AWPMD_split::calc_LMN(int flag) {
  if (norm_needed || (flag & (0x8 | 0x4) && approx != HARTREE)) {
    if (!(calc_state & CALC_LMN)) {
      for (int s1 = 0; s1 < 2; s1++) {
        calc_LMN_s(s1);
        norm_matrix_state[s1] = NORM_CALCULATED;
      } // s1
      calc_state |= CALC_LMN;
    }
  } // flag
  if (norm_needed && !(flag & 0x10)) {   // transforming to physical variables
    for (int s1 = 0; s1 < 2; s1++)
      norm_matrix2phys(s1);
# if 0
    // printing norm matrix
    double *arr=Normh[0][0].arr;
    FILE *f=fopen("norm.dat","wt");
    int nn=Normh[0][0].size;
    for(int i=0;i<nn;i++){
      for(int j=0;j<nn;j++){
        double elm=Normh[0][0](i,j);
        if(fabs(elm)<1e-10)
          elm=0.;
        fprintf(f,"%15g ",elm);
      }
      fprintf(f,"\n");
    }
    fclose(f);
# endif
  }

}

void AWPMD_split::norm_matrix(int s) {
  if (approx != HARTREE && split_wp == false) {
    return AWPMD::norm_matrix(s);
  }
  norm_matrix_state[s] = NORM_CALCULATED;
  if (calc_state & CALC_LMN)
    return;
  int tmp_nneeded = norm_needed;
  norm_needed = 1;
  calc_LMN_s(s);
  norm_matrix2phys(s);
  norm_needed = tmp_nneeded;
  calc_state |= CALC_LMN;
}


int AWPMD_split::norm_factorize(int s) {
  if (norm_matrix_state[s] != NORM_CALCULATED)
    norm_matrix(s);
  if (approx == HARTREE) {
    if (split_wp == false) { // assuming simplectic form
      norm_matrix_state[s] = NORM_FACTORIZED;
      return 1;
    }
    ipiv.resize(nvar[s]); // too big ?
    int l = 0;
    norm_det_log[s] = 0.;
    for (int c = 0; c < ne[s]; c++) {
      int dim = (int)Normh[s][c].size, info;
      DGETRF(&dim, &dim, Normh[s][c].arr, &dim, &ipiv[l], &info);
      l += dim;
      if (info < 0)
        return LOGERR(info, fmt_iv("AWPMD_split.norm_factorize: call to DGETRF failed (exitcode %d)!", info), LINFO);
      std::vector<double> sarr(dim);
      for (int i = 0; i < dim; i++)
        sarr[i] = Normh[s][c](i, i);
      std::sort(sarr.begin(), sarr.end());

      for (int i = 0; i < dim; i++)  // assuming degeneracy of 2?
        norm_det_log[s] += log(fabs(sarr[i]));
    }
    norm_matrix_state[s] = NORM_FACTORIZED;
    return 1;
  } else if (split_wp == true)
    return LOGERR(-2, "AWPMD_split.norm_factorize: NOT IMPLEMENTED FOR SPLIT CASE", LINFO);
  else
    return AWPMD::norm_factorize(s);// works for any non-hartree norm mAtrix
}


int AWPMD_split::norm_invert(int s) {
  if (split_wp == true)
    return LOGERR(-2, "AWPMD_split.norm_invert: NOT IMPLEMENTED FOR SPLIT CASE", LINFO);
  return AWPMD::norm_invert(s); // works for any non-hartree norm mAtrix
}


double AWPMD_split::norm_matrix_det(int s) {
  if (split_wp == true) {
    if (approx != HARTREE)
      return LOGERR(-2, "AWPMD_split.norm_matrix_det: NOT IMPLEMENTED FOR SPLIT CASE & UHF", LINFO);
    if (norm_matrix_state[s] != NORM_FACTORIZED) {
      int res = norm_factorize(s);
      if (res < 0)
        return 0.;
    }
    return
        exp(norm_det_log[s]);
  }
  return AWPMD::norm_matrix_det(s); // works for any non-hartree norm mAtrix
}


double AWPMD_split::norm_matrix_detl(int s) {
  if (split_wp == true) {
    if (approx != HARTREE)
      return LOGERR(-2, "AWPMD_split.norm_matrix_detl: NOT IMPLEMENTED FOR SPLIT CASE & UHF", LINFO);
    /*int tnn = norm_needed;
    norm_needed = 1;
    calc_LMN(0);
    norm_needed = tnn;*/
    if (norm_matrix_state[s] != NORM_FACTORIZED) {
      int res = norm_factorize(s);
      if (res < 0)
        return 0.;
    }
    return norm_det_log[s];
  }
  return AWPMD::norm_matrix_detl(s); // works for any non-hartree norm mAtrix
}


int AWPMD_split::calc_norm_forces_simplectic() {
  int iv1 = 0;
  for (int s1 = 0; s1 < 2; s1++) {
    int ic1 = 0; // starting index of the wp for current electron
    for (int c1 = 0; c1 < ne[s1]; c1++) {
      int indw1 = 8 * ic1;
      int indn1 = (nvar[s1] / 10) * 8 + 2 * ic1;
      for (int k1 = 0; k1 < nspl[s1][c1]; k1++) {
        for (int i = 0; i < 3; i++) { // dx/dt= dE/dp
          dq_dt[s1][indw1 + 8 * k1 + i] = E_der[s1][indw1 + 8 * k1 + i];
        }
        for (int i = 3; i < 6; i++) { // dp/dt= -dE/dx
          dq_dt[s1][indw1 + 8 * k1 + i] = -E_der[s1][indw1 + 8 * k1 + i];
        }
        dq_dt[s1][indw1 + 8 * k1 + 6] = E_der[s1][indw1 + 8 * k1 + 6]; // dw/dt= dE/dpw
        dq_dt[s1][indw1 + 8 * k1 + 7] = -E_der[s1][indw1 + 8 * k1 + 7];// dpw/dt= -dE/dw
        dq_dt[s1][indn1 + 2 * k1] = E_der[s1][indn1 + 2 * k1]; // c
        dq_dt[s1][indn1 + 2 * k1 + 1] = -E_der[s1][indn1 + 2 * k1 + 1];
      }// k1
      ic1 += nspl[s1][c1]; // incrementing block1 wp address
      iv1 += nspl[s1][c1]; // incrementing global variable address
    }// c1
  } // s1
  return 1;
}


int AWPMD_split::prepare_constraints(int s, int c) {
  // now getting the constraint Jacobians
  int eid = s * ne[0] + c;
  vector<vector<double> > &constr = lin_constraints[eid]; // constraints for this electron

  int nconstr = (int) constr.size();
  if (!nconstr)
    return 0;
  int dim = 10 * nspl[s][c], ddim = dim - nconstr;

  int info;
  recmatrix<double> cmatr((size_t) nconstr, (size_t) dim);  // must be transpose for Fortran
  for (size_t i = 0; i < constr.size(); i++) {
    for (size_t j = 0; j < constr[i].size(); j++)
      cmatr(i, j) = constr[i][j];
  }

# if 0
  {
  // printing C matrix
  double *arr=cmatr.arr;
  FILE *f=fopen("c1.dat","wt");
  size_t n1=cmatr.sizey;
  size_t n2=cmatr.sizex;
  for(size_t i=0;i<n1;i++){
    for(size_t j=0;j<n2;j++){
      double elm=cmatr(j,i);
      //if(fabs(elm)<1e-10)
        //elm=0.;
      fprintf(f,"%15g ",elm);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  }
# endif

  vector<int> ipiv(nconstr), rpivots(dim);
  // LU decomposition
  DGETRF(&dim, &nconstr, cmatr.arr, &dim, &ipiv[0], &info);
  if (info < 0)
    LOGERR(info, fmt_iv("AWPMD_split.prepare_constraints: call to DGETRF failed (exitcode %d)!", info), LINFO);

  // making pivots
  constr_pivots[s][c].resize(dim);
  rconstr_pivots[s][c].resize(dim);
  for (size_t i = 0; i < (size_t) dim; i++) {
    constr_pivots[s][c][i] = i;
    rconstr_pivots[s][c][i] = i;
  }
  for (int i = 0; i < nconstr; i++) {
    if (ipiv[i] != i + 1)
      swap(constr_pivots[s][c][i], constr_pivots[s][c][ipiv[i] - 1]);
  }
  for (int i = nconstr - 1; i >= 0; i--) {
    if (ipiv[i] != i + 1)
      swap(rconstr_pivots[s][c][i], rconstr_pivots[s][c][ipiv[i] - 1]);
  }


  // making right hand side matrix
  Jh[s][c].init(ddim, nconstr, 1);
  for (int i = 0; i < nconstr; i++) {
    ipiv[i] = i + 1; // new ipiv
    for (int j = 0; j < ddim; j++)
      Jh[s][c](j, i) = -constr[i][constr_pivots[s][c][j +
                                                      nconstr]]; //-cmatr(i,j+nconstr);  // "-" for moving the free vars to rhs
  }


# if 0
  {
  // printing C matrix
  double *arr=cmatr.arr;
  FILE *f=fopen("c2.dat","wt");
  size_t n1=cmatr.sizey;
  size_t n2=cmatr.sizex;
  for(size_t i=0;i<n1;i++){
    for(size_t j=0;j<n2;j++){
      double elm=cmatr(j,i);
      //if(fabs(elm)<1e-10)
        //elm=0.;
      fprintf(f,"%15g ",elm);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  }
# endif

  // solving for Jh
  DGETRS("T", &nconstr, &ddim, cmatr.arr, &dim, &ipiv[0], Jh[s][c].arr, &nconstr, &info);
  // analyze return code here
  if (info < 0)
    return LOGERR(info, fmt_iv("AWPMD_split.prepare_constraints: call to DGETRS failed (exitcode %d)!", info), LINFO);


# if 0
  // printing Jh matrix
  double *arr=Jh[s][c].arr;
  FILE *f=fopen("jh.dat","wt");
  size_t n1=Jh[s][c].sizex;
  size_t n2=Jh[s][c].sizey;
  for(size_t i=0;i<n1;i++){
    for(size_t j=0;j<n2;j++){
      double elm=Jh[s][c](i,j);
      //if(fabs(elm)<1e-10)
        //elm=0.;
      fprintf(f,"%15g ",elm);
    }
    fprintf(f,"\n");
  }
  fclose(f);
# endif
  return 1;
}


# define NORM_MATRIX_INV_ALGORITHM 1  // 0= linear solver, 1= direct matrix inversion, 2= eigenvalue solver
int got_it = 0;
vector<cdouble> ev[2];

// resizes ipiv
void AWPMD_split::prepare_hartree_nm(int s, int c, int indw, int indn, sqmatrix<double> &matr, vector<double> &rhs,
                                     size_t &nconstr, bool apply_lin_constr) {
  //dim=10*nspl[s][c];
  //ipiv.resize(dim);

  nconstr = 0;
  size_t dim = Normh[s][c].size;
  rhs.resize(dim);
  // putting the forces into rhs vector
  for (size_t i = 0; i < dim; i++) {
    int k = (int)i / 10; // wp number within electron c
    int var = i % 10; // variable number
    if (var < 8)
      rhs[i] = E_der[s][indw + 8 * k + var];
    else
      rhs[i] = E_der[s][indn + 2 * k + (var - 8)];
  }


  if (apply_lin_constr) { //checking sizes
    int eid = s * ne[0] + c;
    nconstr = lin_constraints[eid].size();
    if (nconstr && Jh[s][c].sizey != nconstr) {
      LOGMSG(vblWARN, "AWPMD_split::prepare_hartree_nm: Jh matrix does not match the constraints, not prepared?", -1);
      apply_lin_constr = false;
      nconstr = 0;
    }
  }
  if (apply_lin_constr && nconstr) { //inclusion of linear constraints
    matr.init(dim, 1);
    // filling the constrained matrix
    for (size_t i = 0; i < dim; i++) {
      if (rconstr_pivots[s][c][i] < nconstr) { // making pseudo row-column
        matr(i, i) = 2.; // flag
        for (size_t j = i + 1; j < dim; j++)
          matr(i, j) = matr(j, i) = 0.;
        continue;
      }
      matr(i, i) = 0.;

      for (size_t j = i + 1; j < dim; j++) {
        if (rconstr_pivots[s][c][j] < nconstr) {
          matr(i, j) = matr(j, i) = 0.;
          continue;
        }

        double val = Normh[s][c](i, j);
        for (size_t alpha = 0; alpha < nconstr; alpha++) {
          for (size_t beta = 0; beta < nconstr; beta++) {
            val += Normh[s][c](constr_pivots[s][c][alpha], constr_pivots[s][c][beta]) *
                   Jh[s][c](rconstr_pivots[s][c][i] - nconstr, alpha) *
                   Jh[s][c](rconstr_pivots[s][c][j] - nconstr, beta);
          }
          val += Normh[s][c](constr_pivots[s][c][alpha], rconstr_pivots[s][c][j] - nconstr) *
                 Jh[s][c](rconstr_pivots[s][c][i] - nconstr, alpha);
          val += Normh[s][c](rconstr_pivots[s][c][i] - nconstr, constr_pivots[s][c][alpha]) *
                 Jh[s][c](rconstr_pivots[s][c][j] - nconstr, alpha);
        }
        matr(i, j) = val;
        matr(j, i) = -val;
      }
      for (size_t alpha = 0; alpha < nconstr; alpha++)
        rhs[i] += rhs[constr_pivots[s][c][alpha]] * Jh[s][c](rconstr_pivots[s][c][i] - nconstr, alpha);
    }
    for (size_t i = 0; i < dim; i++) { // clearing rhs for fixed vars
      if (rconstr_pivots[s][c][i] < nconstr) // making pseudo row-column
        rhs[i] = 0.;
    }
  } else {
    //matr=Normh[s][c];
    matr.init(dim, 1);
    for (size_t i = 0; i < dim; i++) {
      matr(i, i) = 0.;
      for (size_t j = i + 1; j < dim; j++) {
        double val = Normh[s][c](i, j);
        matr(i, j) = val;
        matr(j, i) = -val;
      }
    }
  }


  if (norm_mode == NORMALIZE && apply_lin_constr && !nconstr) {   // otherwise there is no need to truncate norm matrix
    vector<int> trunc(nspl[s][c], 0);

    // finding wp with maximal c
    double mc = 1e20;
    int imax = 0;
    for (int i = 0; i < nspl[s][c]; i++) {
      //cdouble c(split_c[s][ic+i][0],split_c[s][ic+i][1]);
      //double nc=norm(c);
      double nc = 0., no = 0.;
      // measuring norm matrix there
      for (int j = 0; j < (int) dim; j++) {
        for (int k = 0; k < 10; k++) {
          double dn = matr(10 * i + k, j) * matr(10 * i + k, j);
          no += dn;
          if (k < 8)
            nc += dn;
        }
      }
      if (no < 1e-6)// norm matrix part is too small, truncating
        trunc[i] = 1;
      if (mc > nc && !trunc[i]) {
        mc = nc;
        imax = i;
      }

    }
    if (!trunc[imax])
      trunc[imax] = 2;
    //imax=1;
    // trancating norm matrix there
    for (int j = 0; j < nspl[s][c]; j++) {
      if (!trunc[j])
        continue;
      int k0 = (trunc[j] == 1 ? 0 : 8);
      for (int k = k0; k < 10; k++) {
        for (int i = 0; i < (int) dim; i++)
          matr(10 * j + k, i) = matr(i, 10 * j + k) = 0.;
        matr(10 * j + k, 10 * j + k) = k % 2 ? -3. : 3.; // flag
        rhs[10 * j + k] = 0.;
        /*int var=k%10; // variable number
        if(var<8)
          E_der[s][indw+8*j+var]=0.;
        else
          E_der[s][indn+2*j+(var-8)]=0.;*/
      }
    }
  } // normalize



  // Replace columns and rows for width and p_width with identity matrix
  // if the width is fixed
  if (constraint == FIX) {
    for (int i = 0; i < nspl[s][c]; i++) {
      int idxw = 10 * i + 6, idxpw = idxw + 1;
      for (int j = 0; j < 10 * nspl[s][c]; j++) {
        matr(idxw, j) = 0;
        matr(idxpw, j) = 0;
        matr(j, idxw) = 0;
        matr(j, idxpw) = 0;
      }
      matr(idxw, idxw) = 1.; // flag
      matr(idxpw, idxpw) = 1.;
      // Zero forces for width and p_width if the width is fixed
      rhs[idxw] = 0.;
      rhs[idxpw] = 0.;
    }
  }
  // same for fixed variables
  if (fixed[s].size()) { // some variables are marked as fixed
    for (int j = 0; j < nspl[s][c]; j++) {
      for (int par_id = 0; par_id < 10; par_id++) {
        int var_id;
        if (par_id < 8)
          var_id = indw + 8 * j + par_id;
        else
          var_id = indn + 2 * j + (par_id - 8);
        map<int, int>::iterator it = fixed[s].find(var_id);
        if (it != fixed[s].end() && (it->second & 0x1)) { // making norm matrix diagonal for fixed vars
          int idx = 10 * j + par_id;
          for (int v = 0; v < 10 * nspl[s][c]; v++) {
            matr(idx, v) = 0;
            matr(v, idx) = 0;
          }
          matr(idx, idx) = 1.; // flag
          // Zero forces for fixed variable
          rhs[idx] = 0.;
        }
      }  // parameter
    }// j
  }


# if 0
  // printing norm matrix
  double *arr=matr.arr;
  FILE *f=fopen("normt1.dat","wt");
  int nn=matr.size;
  for(int i=0;i<nn;i++){
    for(int j=0;j<nn;j++){
      double elm=matr(i,j);
      //if(fabs(elm)<1e-10)
        //elm=0.;
      fprintf(f,"%15g ",elm);
    }
    fprintf(f,"* %g\n",rhs[i]);
  }
  fclose(f);
# endif
}

int fix_second = 0;

int AWPMD_split::calc_norm_forces() {
  if (approx == HARTREE) {
    // inverting
    vector<double> rhs;
    for (int s = 0; s < 2; s++) {
      norm_det_log[s] = 0.;
      //int ic=0;
      //int ics=s*nvar[0]/10;
      int indw = 0; // starting index of the electron wp coordinates
      int indn = nwp[s] * 8; // starting index of the electron norm coordinates
      for (int c = 0; c < ne[s]; c++) {
        int info;
        sqmatrix<double> matr;
        size_t nconstr;
        prepare_hartree_nm(s, c, indw, indn, matr, rhs, nconstr, true);
        int dim = (int) matr.size;
        ipiv.resize(dim);


# if NORM_MATRIX_INV_ALGORITHM == 0 || NORM_MATRIX_INV_ALGORITHM == 1  // reduce to upper triangular form

        DGETRF(&dim, &dim, matr.arr, &dim, &ipiv[0], &info); // matrix is transpose !!!
        if (info < 0)
          LOGERR(info, fmt_iv("AWPMD_split.calc_norm_forces: call to DGETRF failed (exitcode %d)!", info), LINFO);


        for (int i = 0; i < dim; i++)
          norm_det_log[s] += log(fabs(matr(i, i)));
# endif

# if NORM_MATRIX_INV_ALGORITHM == 1  // matrix inversion
        int lwork = dim * 64;
        vector<double> work(lwork);
        DGETRI(&dim, matr.arr, &dim, &ipiv[0], &work[0], &lwork, &info);  // matrix is transpose !!
        if (info != 0)
          LOGERR(info, fmt_iv("AWPMD_split.calc_norm_forces: call to DGETRI failed (exitcode %d)!", info), LINFO);

        // multipying the inverse matrix by force
        vector<double> lhs(dim);
        for (int i = 0; i < dim; i++) {
          double sum = 0.;
          for (int j = 0; j < dim; j++)    // double transpose
            sum += matr(i, j) * rhs[j];
          lhs[i] = sum;
        }

# if 0
        // printing norm matrix
        //arr=matr.arr;
        FILE *f=fopen("normi.dat","wt");
        int nn=matr.size;
        for(int i=0;i<nn;i++){
          for(int j=0;j<nn;j++){
            double elm=matr(i,j);  // double transpose
            if(fabs(elm)<1e-10)
              elm=0.;
            fprintf(f,"%15g ",elm);
          }
          fprintf(f,"* %g\n",lhs[i]);
        }
        fclose(f);
# endif


# elif NORM_MATRIX_INV_ALGORITHM == 0 // this will use single equation Nq=f
        int nrhs=1;
        DGETRS("T",&dim,&nrhs,matr.arr,&dim,&ipiv[0],&rhs[0],&dim,&info);
        // analyze return code here
        if(info<0)
          return LOGERR(info,fmt("AWPMD_split.interacton: call to DGETRS failed (exitcode %d)!",info),LINFO);
        vector<double> &lhs=rhs;

# elif NORM_MATRIX_INV_ALGORITHM==2  // using eigenvalue solver

        int &ddim = dim;
        // making hermitian matrix from i*Norm
        int narr=ddim*(ddim+1)/2;
        vector<cdouble> arrp(narr); // matrix in packed storage
        for(int i=0 ;i<ddim;i++){
          arrp[i+(i+1)*i/2]=cdouble(matr(i,i),0);  // truncated matrix parts contain real  1 on diagonal, otherwise 0
          for(int j=i+1;j<ddim;j++)
            arrp[i+(j+1)*j/2]=cdouble(0.,matr(i,j));
        }

        int lwork=2*ddim, lrwork=2*ddim*ddim+5*ddim+1, liwork=5*ddim+3;
        vector<double>  rwork(lrwork);
        vector<int> iwork(liwork);
        vector<double> eigen_val(ddim);
        vector<cdouble> eigen_vect(ddim*ddim), work(lwork);

        ZHPEVD("V","U", &ddim, (MKL_Complex16*)&arrp[0], &eigen_val[0], (MKL_Complex16*)&eigen_vect[0], &ddim,  (MKL_Complex16*)&work[0], &lwork, &rwork[0], &lrwork, &iwork[0], &liwork, &info);
        if(info != 0)
          LOGERR(info,fmt("AWPMD_split.interaction: call to ZHPEVD failed (exitcode %d)!",info),LINFO);
# if 1 // projected eigenvectors
        Vector_G vrhs(rhs);
        vector<double> lhs1(rhs.size(),0.); // solution of equation (dq/dt)
        for(int i=0;i<ddim;i++){ // over eigen vectors
          if(eigen_val[i]>0.) // only negative part of spectrum
            continue;
          vector<double> prj_evectr(ddim),prj_evecti(ddim) ;
          for(int j=0;j<ddim;j++){
            prj_evectr[j]=real(eigen_vect[ddim*i+j]); // fortran ?
            prj_evecti[j]=imag(eigen_vect[ddim*i+j]); // fortran ?
          }

          Vector_G vprj_evectr(prj_evectr), vprj_evecti(prj_evecti);
          double fi=(vprj_evecti*vrhs);
          double fr=(vprj_evectr*vrhs);
          Vector_G vtot=fi*vprj_evectr - fr*vprj_evecti;
          //double force_comp=2*vtot.norm();
          //double cross=vprj_evectr*vprj_evecti;

          for(int j=0;j<ddim;j++)
            lhs1[j]+=2*(fi*prj_evectr[j]-fr*prj_evecti[j])/eigen_val[i];
        } //i
        for(int j=0;j<ddim;j++)
          rhs[j]=lhs1[j];
# else


        // decomposing rhs into eigenvectors
        vector<cdouble> coef(ddim,0.);
        for(int i=0;i<ddim;i++){
          for(int j=0;j<ddim;j++)
            coef[i]+=cdouble(0,rhs[j])*conj(eigen_vect[ddim*i+j]); // fortran ?
        }

        norm_det_log[s] = 0.;
        // filtering eigenvalues
        vector<cdouble> ceval(ddim);
        double lambda_thresh=0.5e-1;
        double lambda_thresh_rel=0.1*fabs(eigen_val[ddim/2+1]); //20*lambda_thresh;
        for(int i=0; i<ddim; i++) {
          norm_det_log[s] += log(fabs( eigen_val[i] ));

          double lambda0=fabs(eigen_val[i])/sqrt(norm(coef[i]));
          double lambda=fabs(eigen_val[i]);
          double angle=0.;
          //if(i==ddim/2-1 || i==ddim/2)
            //angle=M_PI*(1.-erf(4.*(lambda-2*lambda_thresh)/lambda_thresh))/2.;
          ceval[i]=polar(eigen_val[i]/*+lambda_thresh*lambda_thresh/eigen_val[i]*/,angle);

# if 0
          if(fabs(eigen_val[i])<lambda_thresh){
            //eigen_val[i]=eigen_val[i];
            //eigen_val[i]=/*1e20* */  lambda_thresh*(eigen_val[i] <0 ? -1. : 1.);
            /*if(!got_it){
              int j;
              j=ddim/2-1;
              ev[0]=vector<cdouble>(eigen_vect.begin()+ddim*j,eigen_vect.begin()+ddim*j+ddim);
              j=ddim/2;
              ev[1]=vector<cdouble>(eigen_vect.begin()+ddim*j,eigen_vect.begin()+ddim*j+ddim);
            }*/
            got_it=1;
          }
          /*else{
            if(got_it)
              eigen_val[i]=0.;
          }*/
# endif

          //norm_det_log += log(fabs( eigen_val[i] ));
        }


        /*if(got_it){
          if(fabs(eigen_val[ddim/2])>lambda_thresh_rel)
            got_it=0;
        }
        got_it=1;*/
# if 0
        cdouble ceval[2]={ cdouble(eigen_val[ddim/2-1],0),cdouble(eigen_val[ddim/2],0)};

        if(got_it){
          double lambda=fabs(eigen_val[ddim/2-1]);
          double angle=M_PI*(1.-erf(4.*(lambda-2*lambda_thresh)/lambda_thresh))/2.;
          ceval[0]*=polar(1.,angle);
          ceval[1]*=polar(1.,angle);
          //eigen_val[ddim/2-1]=-eigen_val[ddim/2-1];
          //eigen_val[ddim/2]=-eigen_val[ddim/2];
          /*int i=ddim/2-1;
          eigen_val[i]=lambda_thresh*(eigen_val[i] <0 ? -1. : 1.);
          i=ddim/2;
          eigen_val[i]=lambda_thresh*(eigen_val[i] <0 ? -1. : 1.);*/

        }
# endif
        // assembling solutions
        vector<double> rhs(dim);
        for(int i=0;i<ddim;i++){
          rhs[i]=0.;
          for(int j=0;j<ddim;j++){
            //if(j!=ddim/2-1 && j!=ddim/2)
              //continue;
            //if(eigen_val[j]==0.) // filter
            //  continue;
            //if((j==ddim/2-1 || j==ddim/2)){
             // if(ceval[j-(ddim/2-1)]==0.)
               // continue;
              //rhs[i]+=real(coef[j]*eigen_vect[ddim*j+i]/ceval[j-(ddim/2-1)]);
            //}
            //else
              //rhs[i]+=real(coef[j]*eigen_vect[ddim*j+i])/eigen_val[j];
            rhs[i]+=real(coef[j]*eigen_vect[ddim*j+i]/ceval[j]);
          }
        }
# endif // projected eigenvectors
        vector<double> &lhs=rhs;
# endif

        // getting back the forces
        for (int i = 0; i < dim; i++) {
          double val = 0.;
          if (nconstr && rconstr_pivots[s][c][i] < nconstr) { // constrained variable
            size_t ddim = dim - nconstr;
            for (size_t j = 0; j < ddim; j++) // resolving constraints explicitly
              val += lhs[constr_pivots[s][c][nconstr + j]] * Jh[s][c](j, rconstr_pivots[s][c][i]);
          } else
            val = lhs[i];

          int k = i / 10; // wp number within electron c
          int var = i % 10; // variable number
          if (var < 8)
            dq_dt[s][indw + 8 * k + var] = val;
          else
            dq_dt[s][indn + 2 * k + (var - 8)] = val;

        }
        indw += 8 * nspl[s][c]; // 8 variables in each wavepacket
        indn += 2 * nspl[s][c]; // 2 variables in each wp norm
      }//c
    }//s
  }// hartree
  else { // UHF
    if (split_wp)// single split only
      return LOGERR(-2, "AWPMD_split.calc_norm_forces: implemented for HARTREE or non split cases only", LINFO);

    vector<double> rhs;
    for (int s = 0; s < 2; s++) {
      norm_det_log[s] = 0.;
      int info;
      sqmatrix<double> &matr = Norm[s];

      int dim = (int) matr.size;
      if (!dim)   // no electrons in this subsystem
        continue;
      ipiv.resize(dim);
      rhs.resize(dim);
      for (int i = 0; i < dim; i++)
        rhs[i] = E_der[s][i];

# if NORM_MATRIX_INV_ALGORITHM == 0 || NORM_MATRIX_INV_ALGORITHM == 1  // reduce to upper triangular form

      DGETRF(&dim, &dim, matr.arr, &dim, &ipiv[0], &info); // matrix is transpose !!!
      if (info < 0)
        LOGERR(info, fmt_iv("AWPMD_split.calc_norm_forces: call to DGETRF failed (exitcode %d)!", info), LINFO);


      for (int i = 0; i < dim; i++)
        norm_det_log[s] += log(fabs(matr(i, i)));
# endif

# if NORM_MATRIX_INV_ALGORITHM == 1  // matrix inversion
      int lwork = dim * 64;
      vector<double> work(lwork);
      DGETRI(&dim, matr.arr, &dim, &ipiv[0], &work[0], &lwork, &info);  // matrix is transpose !!
      if (info != 0)
        LOGERR(info, fmt_iv("AWPMD_split.calc_norm_forces: call to DGETRI failed (exitcode %d)!", info), LINFO);

      // multipying the inverse matrix by force
      vector<double> lhs(dim);
      for (int i = 0; i < dim; i++) {
        double sum = 0.;
        for (int j = 0; j < dim; j++)    // double transpose
          sum += matr(i, j) * rhs[j];
        dq_dt[s][i] = sum;
      }
      for (int i = dim; i < nvar[s]; i++) // c part is 0
        dq_dt[s][i] = 0;

# if 0
      // printing norm matrix
      //arr=matr.arr;
      FILE *f=fopen("normi.dat","wt");
      int nn=matr.size;
      for(int i=0;i<nn;i++){
        for(int j=0;j<nn;j++){
          double elm=matr(i,j);  // double transpose
          if(fabs(elm)<1e-10)
            elm=0.;
          fprintf(f,"%15g ",elm);
        }
        fprintf(f,"* %g\n",lhs[i]);
      }
      fclose(f);
# endif


# elif NORM_MATRIX_INV_ALGORITHM == 0 // this will use single equation Nq=f
      int nrhs=1;
      DGETRS("T",&dim,&nrhs,matr.arr,&dim,&ipiv[0],&rhs[0],&dim,&info);
      // analyze return code here
      if(info<0)
        return LOGERR(info,fmt("AWPMD_split.interacton: call to DGETRS failed (exitcode %d)!",info),LINFO);
      //vector<double> &lhs=rhs;
      for(int i=0;i<dim;i++) // copying
        dq_dt[s][i] = rhs[i];
      for(int i=dim;i<nvar[s];i++) // c part is 0
        dq_dt[s][i]=0;

# endif
    }// s
  }// UHF
  return 1;
}

/*
void AWPMD_split::calc_ext_power(Vector_3P fe_x,
                           Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c){

  // calculating the power of external force
  // using the fact that dq_i/dt = inv(Norm)*(dHtot/d_qi) = -force_qi
  Wext=0.;
  for(int s=0;s<2;s++){
    int ic=0, ics=s*nvar[0]/10;

    for(int c=0; c<ne[s]; ic+=nspl[s][c++]){
      for(int i=0; i<nspl[s][c]; i++){
        int indw=8*(ic+i);
        int indn=(nvar[s]/10)*8+2*(ic+i);
        for(int j=0;j<3;j++)
          Wext+=-F_extra[s][indw++]*fe_x[ics+ic+i][j];
        for(int j=0;j<3;j++)
          Wext+=-F_extra[s][indw++]*fe_p[ics+ic+i][j];
        Wext+=-F_extra[s][indw++]*fe_w[ics+ic+i];
        Wext+=-F_extra[s][indw++]*fe_pw[ics+ic+i];
        for(int j=0;j<2;j++)
          Wext+=-F_extra[s][indn++]*fe_c[ics+ic+i][j];
      }
    }
  }
}*/


void AWPMD_split::calc_ext_power() {
  // calculating the power of external force
  // using the fact that dq_i/dt = inv(Norm)*(dHtot/d_qi) = -force_qi
  Wext = 0.;
  for (int s = 0; s < 2; s++) {
    int ic = 0, ics = s * nvar[0] / 10;

    for (int c = 0; c < ne[s]; ic += nspl[s][c++]) {
      for (int i = 0; i < nspl[s][c]; i++) {
        int indw = 8 * (ic + i);
        int indn = (nvar[s] / 10) * 8 + 2 * (ic + i);
        for (int j = 0; j < 8; j++) {
          Wext += -F_extra[s][indw] * dq_dt[s][indw];
          indw++;
        }

        for (int j = 0; j < 2; j++) {
          Wext += -F_extra[s][indn] * dq_dt[s][indn];
          indn++;
        }
      }
    }
  }
}

int AWPMD_split::constraint_int_power(Vector_3P fe_x,
                                      Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c, int release_all) {

  // calculating the power of external force
  // using the fact that dq_i/dt = inv(Norm)*(dHtot/d_qi) = -force_qi

  int res = 0; // number of fix changes
  double W_int_thresh = 1e6;
  double decr = 1.; //0.1;

  for (int s = 0; s < 2; s++) {
    int ic = 0, ics = s * nvar[0] / 10;

    for (int c = 0; c < ne[s]; ic += nspl[s][c++]) {
      for (int i = 0; i < nspl[s][c]; i++) {
        int indw = 8 * (ic + i);
        int indn = (nvar[s] / 10) * 8 + 2 * (ic + i);
        double Wi0, Wi1;
        /*       for(int j=0;j<3;j++){
          Wi0=-E_der[s][indw+j]*fe_x[ics+ic+i][j];
          Wi1=-E_der[s][indw+3+j]*fe_p[ics+ic+i][j];
          if(!release_all && fabs(Wi0)+fabs(Wi1)> W_int_thresh*decr ){
            res+=set_var_constraint_by_id(s,indw+j,2);
            res+=set_var_constraint_by_id(s,indw+j+3,2);
            //fe_x[ics+ic+i][j]=0.; // zeroing all conjugate forces
            //fe_p[ics+ic+i][j]=0.;
          }
          else if(release_all || fabs(Wi0)+fabs(Wi1)> W_int_thresh*decr){
            res+=set_var_constraint_by_id(s,indw+j,-2);
            res+=set_var_constraint_by_id(s,indw+j+3,-2);
          }
        } */
        Wi0 = -E_der[s][indw + 7] * fe_w[ics + ic + i];
        Wi1 = -E_der[s][indw + 8] * fe_pw[ics + ic + i];
        if (!release_all && fabs(Wi0) + fabs(Wi1) > W_int_thresh) {
          res += set_var_constraint_by_id(s, indw + 6, 2);
          res += set_var_constraint_by_id(s, indw + 7, 2);
          //fe_w[ics+ic+i]=0.; // zeroing all conjugate forces
          //fe_pw[ics+ic+i]=0.;
        } else if (release_all || fabs(Wi0) + fabs(Wi1) > W_int_thresh * decr) {
          res += set_var_constraint_by_id(s, indw + 6, -2);
          res += set_var_constraint_by_id(s, indw + 7, -2);
        }
/*        Wi0=-E_der[s][indn]*fe_c[ics+ic+i][0];
        Wi1=-E_der[s][indn+1]*fe_c[ics+ic+i][1];
        if(!release_all && fabs(Wi0)+fabs(Wi1)> W_int_thresh ){
          res+=set_var_constraint_by_id(s,indn,2);
          res+=set_var_constraint_by_id(s,indn+1,2);
          //fe_c[ics+ic+i][0]=0.; // zeroing all conjugate forces
          //fe_c[ics+ic+i][1]=0.;
        }
        else if(release_all || fabs(Wi0)+fabs(Wi1)> W_int_thresh*decr){
          res+=set_var_constraint_by_id(s,indn,-2);
          res+=set_var_constraint_by_id(s,indn+1,-2);
        }*/
      }
    }
  }
  return res;
}


///\en Calculates interaction in the system of ni ions + electrons
/// the electonic subsystem must be previously setup by set_electrons, ionic by set_ions
/// 0x1   -- give back ion forces \n
/// 0x2   -- add ion forces to the existing set \n
/// 0x4   -- calculate electronic forces  \n
/// 0x8   -- add electronic forces to the existing arrays \n
/// 0x10  -- calculate internal electronic derivatives only: \n
///          will not update electronic force arrays, which may be NULL, \n
///          the forces may be obtained then using \ref get_el_forces() for all WPs \n
///          or separately for each WP using \ref get_wp_force()
/// if PBCs are used the coords must be within a range [0, cell)
int AWPMD_split::interaction(int flag, Vector_3P fi, Vector_3P fe_x,
                             Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c) {
  //if(approx==HARTREE)
  //return interaction_hartree(flag,fi,fe_x,fe_p,fe_w,fe_pw,fe_c);

  int tmp_nneeded = norm_needed;
  if (flag & 0x20)
    norm_needed = 1;
  //0. resize arrays if needed
  resize(flag);
  bool have_force_arrays = (fe_x && fe_p && fe_w && fe_pw && fe_c);
  //1. clearing forces
  //if(have_force_arrays)
  //clear_forces(flag,fi,fe_x,fe_p,fe_w,fe_pw,fe_c);
  // calculate block norms and (optionally) derivatives
  //if(norm_mode==NORMALIZE)
  if (flag & 0x1)
    clear_forces(0x1, fi, NULL, NULL, NULL, NULL, NULL); // clear ion forces if needed
  calc_norms(flag);

  //2. calculating overlap matrix
  int info = calc_overlaps();
  if (info < 0)
    return LOGERR(info, fmt_iv("AWPMD.interacton: overlap matrix inversion failed!"), LINFO);


# if 1
  // calculating the L, M, N matrices
  if (approx == HARTREE || split_wp)
    calc_LMN(flag); // L, M are correct, N incorrect for UHF
  else {
    calc_LMN(flag);
    if (norm_needed) {
      AWPMD::norm_matrix(0); // replacing N with correct norm matrix from AWPMD
      AWPMD::norm_matrix(1);
    }
  }
# else
  // calculating the L, M, N matrices
  calc_LMN(flag);
# endif
  Vector_3 ndr;

  Ebord_ion = Ebord = Wext = Eext = Edc = Edk = Eee = Ew = 0.;
  Eee_hartree = Ee_exch = Eei_exch = Eee_exch = Eext_exch = Ebord_exch = 0.;
  // BEGIN  main energy loop
  for (int s1 = 0; s1 < 2; s1++) {
    //  BEGIN single particle contribution
    Ee[s1] = Eei[s1] = 0.;
    int ic1 = 0;

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    vector<double> *E_der1 = new vector<double>[nthreads];
    int nvars1 = nvar[s1];
    for(int ith=0; ith<nthreads; ith++)
      E_der1[ith].assign(nvars1,0.);
#endif

    for (int c1 = 0; c1 < ne[s1]; c1++) {

      //double Ee1=0., Ew1=0., Eei1=0.;

      int ic2 = 0;
      for (int c2 = 0; c2 < ne[s1]; ic2 += nspl[s1][c2], c2++) {
        if (!check_overlap(s1,c1,c2))
          continue; // non-overlapping WPs

        if (skip_by_tag_order(s1, ic2, c2, s1,ic1,c1) /*c2<c1*/) // taken as M factor into account
          continue;

        double sq_norm12 = norm_mode == NORMALIZE ? sqrt(wf_norm[s1][c1] * wf_norm[s1][c2]) : 1.;
        double pref_ee = -h2_me / 2;  //-h2_me/(2*sq_norm12); // ekin
        double pref_ei = coul_pref / sq_norm12;
        double pref_norm = 1. / sq_norm12; // for external energy



        cdouble yy;
        if (approx == HARTREE)
          yy = 1.;
        else
          yy = Y[s1](c2, c1);
        cdouble Tc1c2 = 0., Tc1c2e = 0.;

        // only diagonal terms for Hartree, and overlap nonzero
        if ((approx != HARTREE || c2 == c1) && check_overlap(s1,c1,c2)) {
          // WP blocks:
          for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
            cdouble cj1(split_c[s1][ic1 + j1][0], split_c[s1][ic1 + j1][1]);
            WavePacket wj1 = wp[s1][ic1 + j1];

            OverlapDeriv o12;
            if (flag & (0x8 | 0x4)) //electron forces needed
              o12.set1(wj1);

            for (int k2 = (c1 == c2 ? j1 : 0); k2 < nspl[s1][c2]; k2++) {
              int M12 = (c1 == c2 && j1 == k2 ? 1 : 2);
              double M12pe, M12pf;
              _mytie(M12pe, M12pf) = check_part1(s1, ic1 + j1, s1, ic2 + k2) * M12;

              cdouble ck2(split_c[s1][ic2 + k2][0], split_c[s1][ic2 + k2][1]);

              // electrons kinetic energy
              WavePacket wk2 = wp[s1][ic2 + k2];
              if (pbc)
                move_to_image(wj1, wk2);

              WavePacket wjk12 = conj(wj1) * wk2;
              cdouble I012 = wjk12.integral();

              if (norm(I012) < 1e-22) // zero overlap !
                continue;

              cdouble part_jk12 = conj(cj1) * ck2;
              cVector_3 djk12 = wjk12.b / (2. * wjk12.a);

              // kinetic energy contribution
              if (M12pf) {
                cVector_3 v1 = conj(wj1.b) * wk2.a - wk2.b * conj(wj1.a);
                cdouble vext = -cVector_3(Ext_force) * wjk12.b / 2 / wjk12.a; // external energy = b/2a
                // calculate border contribution
                cdouble v = (v1 * v1) / wjk12.a;

                v -= 6. * wk2.a * conj(wj1.a);
                v /= wjk12.a;
                v *= pref_ee;
                v += vext; // += external energy

                cdouble dTc1c2_ext = M12pe * part_jk12 * I012 * vext * pref_norm; // without yy
                cdouble dTc1c2 = M12pe * part_jk12 * I012 * v * pref_norm; // without yy


                double dE = real(dTc1c2 * yy);
                double dEext = real(dTc1c2_ext * yy);
                Eext += dEext;  // counting separately

                cdouble dTc1c2_bord = 0.;
                double dEbord = 0.;
                if (use_box) {
                  dTc1c2_bord = M12 * part_jk12 * box.get_integral(wj1.a, wj1.b, wk2.a, wk2.b) * pref_norm; // without y
                  dEbord = real(dTc1c2_bord * yy);
                  Ebord += dEbord;
                  dE += dEbord;
                }
                //double dE=real(yy);
                //Tc1c2+=1.;
                Ee[s1] += dE;
                Eep[s1][ic1 + j1] += 0.5 * dE; //per particle energy
                Eep[s1][ic2 + k2] += 0.5 * dE;

                // exchange energy
                if (c1 == c2 && (ex_energy_type == 1)) {
                  Ee_exch += dE - real(dTc1c2 + dTc1c2_bord); // E_hartree
                  Ebord_exch = dEbord - real(dTc1c2_bord);
                  Eext_exch = dEext - real(dTc1c2_ext);
                } else {
                  Ee_exch += dE;
                  Ebord_exch = dEbord;
                  Eext_exch = dEext;
                }

                // diagonal term
                if (M12 == 1)
                  Edk += dE;

                if (flag & (0x8 | 0x4)) { //electron forces needed
                  cVector_3 tv = wk2.b * conj(wj1.a) - conj(wj1.b) * wk2.a;
                  cdouble ajk2 = wjk12.a * wjk12.a;
                  cdouble ajk3 = ajk2 * wjk12.a;
                  // external potential contribution
                  cdouble dvext_a = (cVector_3(Ext_force) * wjk12.b) / (2 * ajk2);
                  cVector_3 dvext_b = -cVector_3(Ext_force) / (2 * wjk12.a);
                  // border potential contribution
                  cdouble dvbord_aj_conj = 0.;
                  cdouble dvbord_ak = 0.;
                  cVector_3 dvbord_bj_conj = cVector_3(0.);
                  cVector_3 dvbord_bk = cVector_3(0.);
                  // bord_deriv(dvbord_aj_conj,dvbord_bj_conj,....);

                  cdouble dv_aj_conj = -2 * pref_ee * wk2.a * (3 * wjk12.a * wk2.a - tv * wjk12.b) / ajk3
                                       + dvext_a + dvbord_aj_conj;
                  cdouble dv_ak = -2 * pref_ee * conj(wj1.a) * (3 * wjk12.a * conj(wj1.a) + tv * wjk12.b) / ajk3
                                  + dvext_a + dvbord_ak;
                  cVector_3 dv_bj_conj = -(2 * pref_ee) * (wk2.a / ajk2) * tv + dvext_b + dvbord_bj_conj;
                  cVector_3 dv_bk = (2 * pref_ee) * (conj(wj1.a) / ajk2) * tv + dvext_b + dvbord_bk;

                  o12.set2(wk2, &I012);
                  // calculate full derivative of the term
                  // pref*conj(cj)*ck*Ijk*vjk/sqrt(nrm(s1)*nrm(s2))
                  // denominator must be included in pref
                  eterm_deriv(ic1, s1, c1, j1, ic2, s1, c2, k2, M12pf * pref_norm * yy,
                              o12, v, dv_aj_conj, dv_ak, dv_bj_conj, dv_bk);
                  eterm_deriv(ic1, s1, c1, j1, ic2, s1, c2, k2, M12pf * pref_norm * yy,
                              o12, vext, dvext_a, dvext_a, dvext_b, dvext_b, 1); // external energy term

                  Tc1c2 += M12pf * part_jk12 * I012 * v * pref_norm; // all without Y
                  Tc1c2 += dTc1c2_bord; // contribution from border
                  Tc1c2e += M12pf * part_jk12 * I012 * vext * pref_norm; // external energy
                }
                if (use_box) { // box forces are always calculated (for pressure evaluation)
                  box_eterm_deriv(ic1, s1, c1, j1, ic2, s1, c2, k2, M12pf * pref_norm * yy);
                }

              } // M12pe
              // e-i energy
              double sum_re = 0, sum_im = 0;
              double Eeip_sum = 0;

              if (calc_ei) {
                //printf("ei interaction\n");
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum_re,sum_im,Eeip_sum)
#endif
                for (int i = 0; i < ni; i++) {  // ions loop
#ifdef _OPENMP
                  int ith = omp_get_thread_num();
#endif
                  double M12pie, M12pif;
                  _mytie(M12pie, M12pif) = check_part1ei(s1, ic1 + j1, ic2 + k2, i) * M12;
                  if (!M12pif)
                    continue;

                  cVector_3 gjki = djk12 - cVector_3(xi[i]);

                  if (pbc) // correcting the real part (distance) according to PBC
                    gjki = rcell1(gjki, cell, pbc);
                  //-Igor- gkli=cVector_3(real(gkli).rcell1(cell,pbc),imag(gkli));

                  cdouble ngjki = gjki.norm();
                  cdouble sqajk = sqrt(wjk12.a);
                  cdouble v;
                  cdouble k2_4a, exp_k2_4a, exp_p, exp_m, erf_p, erf_m; // for screened interaction only
                  bool zero_force = (fabs(real(ngjki)) + fabs(imag(ngjki)) < 1e-8);
                  double pref_eiq = pref_ei * qi[i] * (qe[s1][ic1 + j1] + qe[s1][ic2 + k2]) * 0.5;

                  if (!screened_ei) {
                    //cdouble ttt = cerf_div(ngjki,c);
                    v = cerf_div(ngjki, sqajk);
                  } else { // screened e-i interaction
                    // checking for zero contribution at large distances/small a

                    k2_4a = 0.25 * kappa_ei * kappa_ei / wjk12.a;
                    if (!zero_force && abs(k2_4a) > 10.) { // TODO: add full derivative for a->0
                      cdouble x1 = ngjki * sqajk + 0.5 * kappa_ei / sqajk;
                      double dE = abs(pref_eiq * exp(-ngjki * ngjki * wjk12.a) / (6 * x1 * ngjki));
                      if (dE < 1e-4) {
                        printf("-"); // temporary to indicate the critical condition
                        continue;
                      }
                    }


                    exp_k2_4a = exp(k2_4a);
                    exp_p = exp(kappa_ei * ngjki);
                    exp_m = 1. / exp_p;
                    if (zero_force) {
                      erf_p = cerf(0.5 * kappa_ei / sqajk);
                      erf_m = -erf_p;
                      v = sqajk * two_over_sqr_pi + kappa_ei * (1. + erf_p) * exp_k2_4a;
                    } else {
                      erf_p = cerf(ngjki * sqajk + 0.5 * kappa_ei / sqajk) * exp_p;
                      erf_m = cerf(ngjki * sqajk - 0.5 * kappa_ei / sqajk) * exp_m;
                      v = 0.5 * (erf_p + erf_m + exp_p - exp_m) / ngjki;
                      v *= exp_k2_4a;
                    }

                  }


                  cdouble dE = pref_eiq * I012 * v * part_jk12;
                  sum_re += real(M12pie * dE); // sum without Y
                  sum_im += imag(M12pie * dE); // sum without Y
                  dE *= yy;

                  //per particle energy
                  Eeip_sum += 0.25 * real(dE) * M12pie;
                  Eiep[i] += 0.5 * real(dE) * M12pie;

                  if (flag & (0x8 | 0x4 | 0x3)) { //some derivatives or forces  needed
                    if (!screened_ei) {
                      cdouble arg = ngjki * sqajk;
                      cdouble t1 = two_over_sqr_pi * exp(-arg * arg);
                      cdouble t2 = t1 / (2 * sqajk);
                      t1 *= sqajk;

                      if (flag & 0x3 && !zero_force) {// calculate forces on ions
                        // TODO(kazeevn) Implement box
                        cdouble dEw = (pref_eiq * yy) * I012 * t1 * part_jk12;
                        dEw = (dE - dEw) / ngjki;
                        //Vector_3 dir=-real(gjki);
                        //dir.normalize();
                        //fi[i]+=M12pif*real( dEw)*dir;
                        fi[i] += M12pif * real(-dEw * gjki / ngjki);
                      }
                      if (flag & (0x8 | 0x4)) { //electron forces needed
                        cdouble dv_ak = (zero_force ? 0. :
                                         (djk12 * gjki) * (v - t1) / (wjk12.a * ngjki * ngjki)) + t2;
                        cVector_3 dv_bk = zero_force ? cVector_3() :
                                          -gjki * (v - t1) / (2 * wjk12.a * ngjki * ngjki);

#ifdef _OPENMP
                        eterm_deriv_omp(E_der1[ith],ic1,s1,c1,j1,
																				E_der1[ith],ic2,s1,c2,k2,M12pif*pref_eiq*yy,
																				o12,v,dv_ak,dv_ak,dv_bk,dv_bk);
#else
                        eterm_deriv(ic1, s1, c1, j1, ic2, s1, c2, k2, M12pif * pref_eiq * yy,
                                    o12, v, dv_ak, dv_ak, dv_bk, dv_bk);
#endif
                      }
                    } else {  // screened e-i interaction
                      cdouble dv_ak;
                      cVector_3 dv_bk;
                      cdouble sqrpia = sqrt(M_PI) * sqajk;
                      if (zero_force)
                        dv_ak = (1. - 2. * k2_4a) / sqrpia - kappa_ei * (k2_4a / wjk12.a) * (1. + erf_p) * exp_k2_4a;
                      else {
                        cdouble arg_p = ngjki * sqajk + 0.5 * kappa_ei / sqajk;
                        cdouble arg_m = arg_p - kappa_ei / sqajk;
                        cdouble exp2_p = exp(-arg_p * arg_p) * exp_p;
                        cdouble exp2_m = exp(-arg_m * arg_m) * exp_m;
                        cdouble dgg = (djk12 * gjki) / ngjki / ngjki;
                        cdouble cosh_kgg = cosh(kappa_ei * ngjki);
                        dv_ak =
                            v * (dgg - k2_4a) / wjk12.a - 0.5 * exp_k2_4a * (2 * dgg - 1.) * (exp2_p + exp2_m) / sqrpia
                            - kappa_ei * exp_k2_4a *
                              ((exp2_p - exp2_m) / (4 * ngjki * sqrpia) + dgg * (cosh_kgg + 0.5 * (erf_p - erf_m))) /
                              wjk12.a;

                        cdouble dv_g = (exp_k2_4a * (wjk12.a * (exp2_p + exp2_m) / sqrpia +
                                                     kappa_ei * (cosh_kgg + 0.5 * (erf_p - erf_m))) - v) / ngjki;
                        dv_bk = (dv_g / 2 / wjk12.a / ngjki) * gjki;

                        if (flag & 0x3) {// calculate forces on ions
                          cdouble pref = pref_eiq * yy * I012 * part_jk12;
                          Vector_3 dir = -real(gjki);
                          dir.normalize();
                          fi[i] += M12pif * real(pref * dv_g) * dir;
                        }
                      }
                      if (flag & (0x8 | 0x4)) { //electron forces needed
#ifdef _OPENMP
                        eterm_deriv_omp(E_der1[ith],ic1,s1,c1,j1, E_der1[ith],
																				ic2,s1,c2,k2,M12pif*pref_eiq*yy,o12,v,dv_ak,dv_ak,dv_bk,dv_bk);
#else
                        eterm_deriv(ic1, s1, c1, j1, ic2, s1, c2, k2, M12pif * pref_eiq * yy,
                                    o12, v, dv_ak, dv_ak, dv_bk, dv_bk);
#endif
                      }
                    } // screened
                  } // need derivarives
                } // ion loop
              }// calc_ei, ions energy contribution

              cdouble sum(sum_re, sum_im);
              double dEei = real(sum * yy);

              Eei[s1] += dEei;

              // per WP energy
              Eeip[s1][ic1 + j1] += Eeip_sum;
              Eeip[s1][ic2 + k2] += Eeip_sum;

              // exchange energy
              if (c1 == c2 && (ex_energy_type == 1))
                Eei_exch += dEei - real(sum); // subtract E_hartree
              else
                Eei_exch += dEei;


              // diagonal term
              if (M12 == 1)
                Edc += real(sum * yy);;

              if (flag & (0x8 | 0x4)) //electron forces needed
                Tc1c2 += sum;

              // extra constraint energy, will only arrive here when M12p=1
              if ((c1 == c2 && j1 == k2) && constraint == HARM) {
                cdouble v = conj(wj1.a) * wk2.a;
                double pref_harm = norm_mode == NORMALIZE ? harm_w0_4 / wf_norm[s1][c1] : harm_w0_4;
                double dE = M12pe * pref_harm * real(part_jk12 * v);
                Ew += dE;

                Ewp[s1][ic1 + j1] += dE; //per particle energy


                if (flag & (0x8 | 0x4)) { //electron forces needed
                  cdouble dv_ak = conj(wj1.a);
                  cdouble dv_aj_conj = wk2.a;
                  eterm_deriv(ic1, s1, c1, j1, ic1, s1, c1, k2, M12pf * pref_harm,
                              o12, v, dv_aj_conj, dv_ak, cVector_3(), cVector_3());
                }
              }
            }// k2
          }// j1
          if (flag & (0x8 | 0x4) && approx != HARTREE) { //electron forces needed
            // adding Y derivative term for all variables (te and tei terms)
            y_deriv(Tc1c2, s1, c2, c1);
            y_deriv(Tc1c2e, s1, c2, c1, 1); // external energy term

          }
        } // !HARTREE or spins or overlap
        // END single particle contribution

#if 1 // pair by pair sum
        if (calc_ee) { // need to calculate electron-electron interaction
          //printf("ee interaction\n");
          // second block
          // e-e interaction
          if (approx == HARTREE && c1 != c2) // only Vkmkm terms for Hartree
            continue;
          if (/*s1==s2 &&*/ skip_by_tag_order(s1, ic2, c2, s1, ic1, c1)) // pair selection term
            continue;

          for (int s2 = s1; s2 < 2; s2++) {


            int ic3 = 0; // starting index of the wp for current electron
            for (int c3 = 0; c3 < ne[s2]; ic3 += nspl[s2][c3], c3++) { // incrementing block2 wp address
              if (!dft_extension) {
                if (approx == HARTREE && s1 == s2 &&
                    c3 == c1) // [Vkkmn contribution for same spin is 0], also  no Vkkkk terms for Hartree (=)
                  continue;
                if (s1 == s2 && (skip_by_tag_order(s2, ic3, c3, s1, ic1, c1) || c3==c1 ) /*c3 <= c1*/) // and general Coulomb sum: i<j (<) //??
                  continue;
              } else { // for dft extensions we add self interaction as well
                if (s1 == s2 && skip_by_tag_order(s2, ic3, c3, s1, ic1, c1) /*c3 < c1*/) // and general Coulomb sum: i<=j (<) //??
                  continue;
              }

              int ic4 = 0;
              for (int c4 = 0; c4 < ne[s2]; ic4 += nspl[s2][c4], c4++) {
                if (approx == HARTREE && c4 != c3) // only Vkmkm terms for Hartree
                  continue;
                //if(s1==s2 && c4==c2) // Vklmm contribution for same spin is 0 for antisymmetrized approximations, also Vmmmm is 0 for Hartree
                //  continue;
                if (!dft_extension) {
                  if (s1 == s2 && (skip_by_tag_order(s2, ic4, c4, s1, ic2, c2) || c4==c2) /*c4 <= c2*/) // pair selection term: Vklmm
                    continue;
                } else { // for dft extensions we add self interaction as well
                  if (s1 == s2 && skip_by_tag_order(s2, ic4, c4, s1, ic2, c2) /*c4 < c2*/) // pair selection term: Vklmm
                    continue;
                }

                if (/*s1==s2 &&*/ c2 == c1 && skip_by_tag_order(s2, ic4, c4, s2, ic3, c3) /*c4 < c3*/) // pair selection term
                  continue;

                double sq_norm34 = norm_mode == NORMALIZE ? sqrt(wf_norm[s2][c3] * wf_norm[s2][c4]) : 1.;
                double pref_ee = coul_pref / (sq_norm12 * sq_norm34);


                cdouble yy;
                double K = 1.;
                if (approx == HARTREE) {
                  yy = 1.;
                } else {
                  if (s1 == s2) { // same spin antisymmetrized term
                    yy = Y[s1](c2, c1) * Y[s1](c4, c3) - Y[s1](c2, c3) * Y[s1](c4, c1);  // check the order of c
                  } else { // different spin atisymmetrized term
                    //if(approx==UHF)
                    //K=2.;
                    yy = K * Y[s2](c4, c3) * Y[s1](c2, c1);
                  }
                  //yy=1.;
                  //yy=Y[s1](c2,c1)*Y[s1](c4,c3);
                  //yy=Y[s1](c2,c3)*Y[s1](c4,c1);
                  //yy=Y[s1](c2,c1)*Y[s1](c4,c3)-Y[s1](c2,c3)*Y[s1](c4,c1);
                }

                cdouble Vy = 0.;
                // WP blocks: Vc1c3c2c4
                for (int j1 = 0; j1 < nspl[s1][c1]; j1++) {
                  cdouble cj1(split_c[s1][ic1 + j1][0], split_c[s1][ic1 + j1][1]);
                  WavePacket wj1 = wp[s1][ic1 + j1];

                  OverlapDeriv o12, o14;
                  if (flag & (0x8 | 0x4)) { //electron forces needed
                    o12.set1(wj1);
                    o14.set1(wj1);
                  }

                  for (int k2 = (approx == HARTREE ? j1 : 0); k2 < nspl[s1][c2]; k2++) {



                    int M12 = (c1 == c2 && j1 == k2 ? 1 : 2);
                    cdouble ck2(split_c[s1][ic2 + k2][0], split_c[s1][ic2 + k2][1]);
                    WavePacket wk2 = wp[s1][ic2 + k2];
                    if (pbc)
                      move_to_image(wj1, wk2);

                    WavePacket wjk12 = conj(wj1) * wk2;
                    cdouble I012 = wjk12.integral();
                    cdouble part_jk12 = conj(cj1) * ck2;
                    cVector_3 djk12 = wjk12.b / (2. * wjk12.a);

                    if (flag & (0x8 | 0x4)) //electron forces needed
                      o12.set2(wk2, &I012);


                    for (int j3 = 0; j3 < nspl[s2][c3]; j3++) {
# if 0
                      double Mejj, Mfjj;
                      _mytie(Mejj, Mfjj) = check_part1(s1, ic1 + j1, s2, ic3 + j3);
                      double Mekj, Mfkj;
                      _mytie(Mekj, Mfkj) = check_part1(s1, ic2 + k2, s2, ic3 + j3);
                      if (!Mfjj && !Mfkj)
                        continue;
# endif
                      

                      cdouble cj3(split_c[s2][ic3 + j3][0], split_c[s2][ic3 + j3][1]);
                      WavePacket &wj3 = wp[s2][ic3 + j3];

                      OverlapDeriv o34, o32;
                      if (flag & (0x8 | 0x4)) { //electron forces needed
                        o34.set1(wj3);
                        o32.set1(wj3);
                        //o32.set2(wk2);
                      }

                      // 3-2
                      WavePacket wjk32;
                      cdouble I032 = 0., part_jk32;
                      cVector_3 djk32;
                      if (s1 == s2 || approx == UHF) {
                        WavePacket wk2 = wp[s2][ic2 + k2];
                        if (pbc)
                          move_to_image(wj3, wk2);
                        wjk32 = conj(wj3) * wk2;
                        I032 = wjk32.integral();
                        part_jk32 = conj(cj3) * ck2;
                        djk32 = wjk32.b / (2 * wjk32.a);
                        if (flag & (0x8 | 0x4)) //electron forces needed
                          o32.set2(wk2, &I032);
                      }


                      for (int k4 = (approx == HARTREE ? j3 : 0); k4 < nspl[s2][c4]; k4++) {
                        /*double Mekk, Mfkk;
                        _mytie(Mekk, Mfkk) = check_part1(s1, ic2 + k2, s2 , ic4 + k4);
                        double Mejk, Mfjk;
                        _mytie(Mejk, Mfjk) = check_part1(s1, ic1 + j1, s2, ic4 + k4);
                        if (!Mfkk && !Mfjk)
                          continue;*/
                        double Me4, Mf4;
                        _mytie(Me4, Mf4) = check_part1(s1, ic1 + j1, ic2+k2, s2, ic3 + j3, ic4+k4);
                        if (!Me4 && !Mf4)
                          continue; 


                        int M34 = (c3 == c4 && j3 == k4 ? 1 : 2);
                        double M0;
                        if (approx == HARTREE)
                          M0 = M12 * M34;
                        else {
                          M0 = (c1 == c2 && c3 == c4 ? 1
                                                     : 2); // will have exchange term for different pairs instead of M12*M34 factor
                        }
                        if (dft_extension) {
                          if (s1 == s2 && c1 == c3 && c2 == c4) // this is self energy
                            M0 *= 0.5;
                        }


                        cdouble ck4(split_c[s2][ic4 + k4][0], split_c[s2][ic4 + k4][1]);

                        // 3-4
                        WavePacket wk4 = wp[s2][ic4 + k4];
# if 1
                        if (1 /*Mfkk && Mfjj*/) {
                          if (pbc)
                            move_to_image(wj3, wk4);
                          WavePacket wjk34 = conj(wj3) * wk4;
                          cdouble I034 = wjk34.integral();
                          if (norm(I034) > ovl_tolerance && norm(I012) > ovl_tolerance) {
                            double Me = M0 * Me4; // Mejj*Mekk;
                            double Mf = M0 * Mf4; // Mfjj*Mfkk;

                            cdouble part_jk34 = conj(cj3) * ck4;
                            cVector_3 djk34 = wjk34.b / (2 * wjk34.a);

                            cVector_3 ddv = djk12 - djk34;
                            cdouble dd = ddv.norm();
                            cdouble aa = 1. / sqrt(1. / wjk12.a + 1. / wjk34.a);
                            cdouble v = cerf_div(dd, aa);
                            double pref_eeq =
                                pref_ee * (qe[s1][ic1 + j1] + qe[s1][ic2 + k2]) * (qe[s2][ic3 + j3] + qe[s2][ic4 + k4]) *
                                0.25;
                            cdouble Vj1j3k2k4 = pref_eeq * I012 * I034 * v * part_jk12 * part_jk34;  // Vklmn
                            double dE = Me * real(Vj1j3k2k4 * yy);
                            Eee += dE;
                            //if(!_finite(Eee))
                            //  Eee=Eee;


                            if (approx != HARTREE) {
                              // exchange energy
                              if (c1 == c2 && c3 == c4 && ((!dft_extension) || (s1 != s2 || c3 != c1)) &&
                                  (ex_energy_type ==
                                   1)) {  // actually c1, c3 are summation indices, should be different for same spin
                                double dEhartree = real(Vj1j3k2k4) * Me * M12 * M34 / M0;
                                Eee_exch += dE - dEhartree; // subtract E_hartree
                                Eee_hartree += dEhartree; // this one is less divergent at zero overlap
                              }
                              else
                                Eee_exch += dE;
                            }
                            else
                              Eee_hartree += dE;

                            Eeep[s1][ic1 + j1] += 0.25 * dE; // per-particle energy
                            Eeep[s1][ic2 + k2] += 0.25 * dE;
                            Eeep[s2][ic3 + j3] += 0.25 * dE;
                            Eeep[s2][ic4 + k4] += 0.25 * dE;

                            // diagonal term
                            if (M12 == 1 && M34 == 1)
                              Edc += dE;

                            bool zero_force = (fabs(real(dd)) + fabs(imag(dd)) < 1e-10);
                            if (flag & (0x8 | 0x4)) { //electron forces needed
                              cdouble arg = dd * aa;
                              cdouble t1 = two_over_sqr_pi * exp(-arg * arg) * aa;
                              cdouble t2 = t1 * aa * aa / 2;


                              cdouble dv_ak1 = (zero_force ? 0. : (djk12 * ddv) * (v - t1) / (wjk12.a * dd * dd)) +
                                               t2 / (wjk12.a * wjk12.a);
                              cVector_3 dv_bk1 = zero_force ? cVector_3() : -ddv * (v - t1) / (2 * wjk12.a * dd * dd);
                              eterm_deriv(ic1, s1, c1, j1, ic2, s1, c2, k2, Mf * pref_eeq * I034 * part_jk34 * yy,
                                          o12, v, dv_ak1, dv_ak1, dv_bk1, dv_bk1);

                              o34.set2(wk4, &I034);
                              cdouble dv_ak2 = (zero_force ? 0. : (-djk34 * ddv) * (v - t1) / (wjk34.a * dd * dd)) +
                                               t2 / (wjk34.a * wjk34.a);
                              cVector_3 dv_bk2 = zero_force ? cVector_3() : ddv * (v - t1) / (2 * wjk34.a * dd * dd);
                              eterm_deriv(ic3, s2, c3, j3, ic4, s2, c4, k4, Mf * pref_eeq * I012 * part_jk12 * yy,
                                          o34, v, dv_ak2, dv_ak2, dv_bk2, dv_bk2);

                              Vy += Me * Vj1j3k2k4; // all without yy
                            }// force
                          }// ovl_tolerance
                        }
# endif
# if 1
                        // 1-4
                        if (approx != HARTREE && (approx == UHF && s1 == s2) /*&& Mfjk && Mfkj*/) {
                          wk4 = wp[s2][ic4 + k4];
                          if (pbc)
                            move_to_image(wj1, wk4);
                          WavePacket wjk14 = conj(wj1) * wk4;
                          cdouble I014 = wjk14.integral();
                          if (norm(I032) > ovl_tolerance && norm(I014) > ovl_tolerance) {
                            double Me = M0 * Me4; // Mejk*Mekj;
                            double Mf = M0 * Mf4; // Mfjk*Mfkj;


                            cdouble part_jk14 = conj(cj1) * ck4;
                            cVector_3 djk14 = wjk14.b / (2 * wjk14.a);

                            cVector_3 ddv = djk32 - djk14;
                            cdouble dd = ddv.norm();
                            cdouble aa = 1. / sqrt(1. / wjk32.a + 1. / wjk14.a);
                            cdouble v = cerf_div(dd, aa);


                            double pref_eeq = pref_ee * (qe[s1][ic1 + j1] + qe[s2][ic4 + k4]) *
                                              (qe[s2][ic3 + j3] + qe[s1][ic2 + k2]) * 0.25;
                            cdouble Vj1j3k4k2 = pref_eeq * I032 * I014 * v * part_jk32 * part_jk14;  // Vklmn
                            double dE = -Me * real(Vj1j3k4k2 * yy);
                            Eee += dE;
                            //if(!_finite(Eee)){
                            //  v=cerf_div(dd,aa);
                            //}


                            Eee_exch += dE;  // no Hartree part at all here

                            Eeep[s1][ic1 + j1] += 0.25 * dE; // per-particle energy
                            Eeep[s1][ic2 + k2] += 0.25 * dE;
                            Eeep[s2][ic3 + j3] += 0.25 * dE;
                            Eeep[s2][ic4 + k4] += 0.25 * dE;


                            bool zero_force = (fabs(real(dd)) + fabs(imag(dd)) < 1e-10);
                            if (flag & (0x8 | 0x4)) { //electron forces needed
                              cdouble arg = dd * aa;
                              cdouble t1 = two_over_sqr_pi * exp(-arg * arg) * aa;
                              cdouble t2 = t1 * aa * aa / 2;

                              o14.set2(wk4, &I014);
                              cdouble dv_ak1 = (zero_force ? 0. : (-djk14 * ddv) * (v - t1) / (wjk14.a * dd * dd)) +
                                               t2 / (wjk14.a * wjk14.a);
                              cVector_3 dv_bk1 = zero_force ? cVector_3() : ddv * (v - t1) / (2 * wjk14.a * dd * dd);
                              eterm_deriv(ic1, s1, c1, j1, ic4, s2, c4, k4, -Mf * pref_eeq * I032 * part_jk32 * yy,
                                          o14, v, dv_ak1, dv_ak1, dv_bk1, dv_bk1);

                              cdouble dv_ak2 = (zero_force ? 0. : (djk32 * ddv) * (v - t1) / (wjk32.a * dd * dd)) +
                                               t2 / (wjk32.a * wjk32.a);
                              cVector_3 dv_bk2 = zero_force ? cVector_3() : -ddv * (v - t1) / (2 * wjk32.a * dd * dd);
                              eterm_deriv(ic3, s2, c3, j3, ic2, s1, c2, k2, -Mf * pref_eeq * I014 * part_jk14 * yy,
                                          o32, v, dv_ak2, dv_ak2, dv_bk2, dv_bk2);
                              Vy += -Me * Vj1j3k4k2; // all without yy
                            }// force
                          }// ovl_tolerance
                        } // s1==s2
# endif
                      }// k4
                    }// j3
                  } // k2
                } // j1

                // derivative of yy term
                if (approx != HARTREE && flag & (0x8 | 0x4)) { //electron forces needed
                  //yy=Y[s1](c2,c1)*Y[s1](c4,c3)-Y[s1](c2,c3)*Y[s1](c4,c1);
# if 1
                  if (s1 == s2) {
                    // yy
                    y_deriv(Vy * Y[s1](c4, c3), s1, c2, c1);
                    y_deriv(Vy * Y[s1](c2, c1), s1, c4, c3);
                    y_deriv(-Vy * Y[s1](c4, c1), s1, c2, c3);
                    y_deriv(-Vy * Y[s1](c2, c3), s1, c4, c1);
                  } else {
                    y_deriv(K * Vy * Y[s2](c4, c3), s1, c2, c1);
                    y_deriv(K * Vy * Y[s1](c2, c1), s2, c4, c3);
                  }
# endif
                  // y_deriv(Vy*Y[s1](c4,c1),s1,c2,c3);
                  // y_deriv(Vy*Y[s1](c2,c3),s1,c4,c1);
                }
              } // c4
            } // c3
          } // s2
        } // calc_ee: need e-e interaction
# endif
      }// c2
      ic1 += nspl[s1][c1]; // incrementing block1 wp address
    }// c1

#ifdef _OPENMP
    for(int ith=0; ith<nthreads; ith++)
      for(int i=0; i<nvars1; i++)
        E_der[s1][i] += E_der1[ith][i];

    delete[] E_der1;
#endif

  } // s1
  // END main energy loop
# if DFT_EXTENSION
  if (dft_extension) {
    // supply empty vector if forces are not needed
    std::vector<deriv_function> od = (flag & (0x10 | 0x4 | 0x8)) ? DerivsFunction::GetFunctions() : std::vector<deriv_function>();
    auto result = ex->energy(wp[0], wp[1], od);
    double dE = (double)result.energy;
    Eee += dE;
    Eee_exch += dE;
    if (flag & (0x10 | 0x4 | 0x8)) { // forces neeed
      int ind = 0;
      for (int s = 0; s < 2; s++) {
        for (int iwp = 0; iwp < nwp[s]; iwp++, ind++) {
          for(int ider=0;ider<8;ider++) // adding derivatives
            E_der[s][8 * iwp+ider] += (double)result.derivatives.all[ider][ind];
        }
      }
    }
    /// IV: assuming same contributions from all spins!!!
    double dE_kin0 = (double)result.kinetic_energy / 2;
    double dE_kin1 = (double)result.kinetic_energy / 2;
    Ee[0] += dE_kin0;
    Ee[1] += dE_kin1;
    Ee_dft = static_cast<double>(result.kinetic_energy);
  }
# endif

  // transforming the forces to physical coordinates
  if (flag & (0x8 | 0x4) && !(flag & 0x10) || use_box)
    forces2phys();

  //if(flag&(0x8|0x4) && !(flag&0x10))
  //get_el_forces((flag&(~0x4))|0x8,fe_x,fe_p,fe_w,fe_pw,fe_c); // flag change: electronic forces were cleared already

  if (flag & (0x8 | 0x4)) { // need forces
    // multiplying by the inverse of the norm matrix
    if (flag & 0x20 && (approx == HARTREE || !split_wp)) { // now only hartree or 1 split are supported
      calc_norm_forces();
    } else { // assuming simplectic dq/dt: dp_dt= -dE/dq, dq_dt = dE/dp
      calc_norm_forces_simplectic();
    }
  }
  norm_needed = tmp_nneeded;

  Eii = 0.;
  if (calc_ii)
    interaction_ii(flag, fi);


  if (use_box)
    calc_ion_confinement(flag, fi);

  // constraining energy variations
  //if(flag&(0x8|0x4) && have_force_arrays){
  //constraint_int_power(fe_x,fe_p,fe_w,fe_pw,fe_c);
  //}

  // power of external force
  if (flag & (0x8 | 0x4)) {
    calc_ext_power();
  }

  if (flag & (0x8 | 0x4) && have_force_arrays) // fill force arrays from dq_dt
    get_el_forces(flag, fe_x, fe_p, fe_w, fe_pw, fe_c);


  return 1;
}


int AWPMD_split::norm_matrix2phys(int s) {
  // transform norm matrix to the physical variables
  if (approx != HARTREE) {
    //return LOGERR(-1,fmt("AWPMD_split.norm_matrix2phys: not implemented!"),LINFO);
    int ic = 0;
    int dim = (int)Norm[s].size;
    int indn = (nvar[s] / 10) * 8;
    // rotate rows
    for (int i = 0; i < dim; i++) {  // spans all rows the matrix
      int iwp = i < indn ? i / 8 : (i - indn) / 2; // row wave packet index
      // iterator to list all N(i,*) with fixed i
      sqmatrix<double>::iterator mi = Norm[s].fix_first(i, 0);
      for (int jwp = iwp; jwp < nwp[s]; jwp++) { // column wp index
        WavePacket wj = wp[s][jwp];
        wj.int2phys_der<eq_second>(mi + 8 * jwp, mi + 8 * jwp, mi + 8 * jwp + 3, mi + 8 * jwp + 6, mi + 8 * jwp + 7,
                                   h_plank);
        if (c_polar_mode) { // converting to polar derivative
          double &rfc0 = *(mi + indn + 2 * jwp), &rfc1 = *(mi + indn + 2 * jwp + 1);
          Vector_2 fc(rfc0, rfc1);
          Vector_2 &c = split_c[s][jwp];
          double rho = c.norm();
          rfc0 = (fc[0] * c[0] + fc[1] * c[1]) / rho;
          rfc1 = -fc[0] * c[1] + fc[1] * c[0];
        }
      }
    }
    // rotate columns
    for (int j = 0; j < dim; j++) {  // spans all columns of the matrix
      int jwp = j < indn ? j / 8 : (j - indn) / 2; // column packet index
      // iterator to list all N(*,j) with fixed j
      sqmatrix<double>::iterator mj = Norm[s].fix_second(0, j);
      for (int iwp = 0; iwp <= jwp; iwp++) {
        WavePacket wi = wp[s][iwp];
        wi.int2phys_der<eq_second>(mj + 8 * iwp, mj + 8 * iwp, mj + 8 * iwp + 3, mj + 8 * iwp + 6, mj + 8 * iwp + 7,
                                   h_plank);
        if (c_polar_mode) { // converting to polar derivative
          double &rfc0 = *(mj + indn + 2 * iwp), &rfc1 = *(mj + indn + 2 * iwp + 1);
          Vector_2 fc(rfc0, rfc1);
          Vector_2 &c = split_c[s][jwp];
          double rho = c.norm();
          rfc0 = (fc[0] * c[0] + fc[1] * c[1]) / rho;
          rfc1 = -fc[0] * c[1] + fc[1] * c[0];
        }
      }
    }
    // filling the lower triangle according to antisymmetry
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < i; j++) {
        Norm[s](i, j) = -Norm[s](j, i);
      }
      Norm[s](i, i) = 0.; // making exact
    }
  } // hartree
  else {
    int ic = 0;
    for (int c = 0; c < ne[s]; c++) {
      for (int i = 0; i < nspl[s][c]; i++) {
        WavePacket wi = wp[s][ic + i];
        for (int k = 0; k < 10; k++) {
          // iterator to list all N(10*i+k,*) with fixed 10*i+k
          sqmatrix<double>::iterator mi = Normh[s][c].fix_first(10 * i + k, 0);
          for (int j = i;
               j < nspl[s][c]; j++) {  // TO DO: run this loop from i+1 and take simplectic form for (i,i) block
            WavePacket wj = wp[s][ic + j];
            wj.int2phys_der<eq_second>(mi + 10 * j, mi + 10 * j, mi + 10 * j + 3, mi + 10 * j + 6, mi + 10 * j + 7,
                                       h_plank);
            if (c_polar_mode) { // converting to polar derivative
              Vector_2 fc(*(mi + 10 * j + 8), *(mi + 10 * j + 9));
              Vector_2 &c = split_c[s][ic + j];
              double rho = c.norm();
              *(mi + 10 * j + 8) = (fc[0] * c[0] + fc[1] * c[1]) / rho;
              *(mi + 10 * j + 9) = -fc[0] * c[1] + fc[1] * c[0];
            }
          }
        }// finished line of blocks from right
        for (int k = 10 * i;
             k < 10 * nspl[s][c]; k++) { // TO DO: run this loop from 8*i+8 and take simplectic form for (i,i) block
          // iterator to list all N(8i+*,k) by fixed k
          sqmatrix<double>::iterator mi = Normh[s][c].fix_second(10 * i, k);
          wi.int2phys_der<eq_second>(mi, mi, mi + 3, mi + 6, mi + 7, h_plank);
          if (c_polar_mode) { // converting to polar derivative
            Vector_2 fc(*(mi + 8), *(mi + 9));
            Vector_2 &c = split_c[s][ic + i];
            double rho = c.norm();
            *(mi + 8) = (fc[0] * c[0] + fc[1] * c[1]) / rho;
            *(mi + 9) = -fc[0] * c[1] + fc[1] * c[0];
          }

        }// finished line of blocks from left

        for (int k = 0; k < 10; k++) { // filling the lower triangle according to antisymmetry
          for (int j = 10 * i + 10; j < 10 * nspl[s][c]; j++)
            Normh[s][c](j, 10 * i + k) = -Normh[s][c](10 * i + k, j);
        }
      } // i
      ic += nspl[s][c];
    }// c
  } // hartree
  return 1;
}


cdouble AWPMD_split::operator*(AWPMD_split &other) {
  if (approx != HARTREE) {
    LOGMSG(vblWARN, "AWPMD_split::operator*: total WF overlap is not implemented for this approximation!", LINFO);
    return 0.;
  }
  if (!valid_norms)
    calc_norms(0);
  if (!other.valid_norms)
    other.calc_norms(0);

  cdouble ovl = 0.;
  for (int s = 0; s < 2; s++) {
    for (int c1 = 0, ic1 = 0; c1 < ne[s]; ic1 += nspl[s][c1], c1++) {
      for (int c2 = 0, ic2 = 0; c2 < other.ne[s]; ic2 += other.nspl[s][c2], c2++) {
        cdouble sum(0, 0);
        for (int j1 = 0; j1 < nspl[s][c1]; j1++) {
          double cj_re = split_c[s][ic1 + j1][0];
          double cj_im = split_c[s][ic1 + j1][1];
          cdouble ccj = cdouble(cj_re, -cj_im);
          WavePacket wj = wp[s][ic1 + j1];
          for (int k2 = 0; k2 < other.nspl[s][c2]; k2++) {
            double ck_re = other.split_c[s][ic2 + k2][0];
            double ck_im = other.split_c[s][ic2 + k2][1];
            cdouble ck = cdouble(ck_re, ck_im);

            WavePacket wk = other.wp[s][ic2 + k2];
            if (pbc)
              move_to_image(wj, wk);
            WavePacket wjk = conj(wj) * wk;
            sum += ccj * ck * wjk.integral();
          }
        }
        if (norm_mode == NORMALIZE)
          ovl += sum / sqrt(wf_norm[s][c1] * other.wf_norm[s][c2]);
        else
          ovl += sum;
      }
    }
  }
  return ovl;
}


///\en Density term
struct density_term_t {
  typedef cdouble result_t;
  AWPMD_split *awpmds;
  Vector_3 x;
  const iVector_3 &integ_flag;

  density_term_t(AWPMD_split *awpmds_, const Vector_3 &x_, const iVector_3 &integ_flag_) : awpmds(awpmds_), x(x_),
                                                                                           integ_flag(integ_flag_) {}

  result_t operator()(const WavePacket &bra, const WavePacket &ket) const {
    // test: must give exact sum of electrons for all projections
    WavePacket braket = conj(bra) * ket;
    return braket.partial_integrate(x, integ_flag);
    //return braket.integral();
  }

  void add2(cdouble k, const WavePacket &bra, const WavePacket &ket) {}
};


cdouble AWPMD_split::Psi(const Vector_3 &x, const iVector_3 &integ_flag) {
  if (ne[0] > 1 || ne[1] > 1)
    return sqrt(el_density(x, integ_flag));

  if (!valid_norms)
    calc_norms(0);

  cdouble ovl = 0.;
  for (int s = 0; s < 2; s++) {
    for (int c1 = 0, ic1 = 0; c1 < ne[s]; ic1 += nspl[s][c1], c1++) {
      cdouble sum(0, 0);
      for (int j1 = 0; j1 < nspl[s][c1]; j1++) {
        double cj_re = split_c[s][ic1 + j1][0];
        double cj_im = split_c[s][ic1 + j1][1];
        cdouble cj = cdouble(cj_re, cj_im);
        WavePacket wj = wp[s][ic1 + j1];
        sum += cj * wj.partial_integrate(x, integ_flag);
      }
      if (norm_mode == NORMALIZE)
        ovl += sum / sqrt(wf_norm[s][c1]);
      else
        ovl += sum;
    }
  }
  return ovl;
}


double AWPMD_split::el_density(const Vector_3 &x, const iVector_3 &integ_flag) {
  if (approx != HARTREE) {
    if (!valid_norms) {
      calc_norms(0);
      calc_overlaps();
    }
    density_term_t dens(this, x, integ_flag);
    return real(exp_value_2(dens)) / (ne[0] + ne[1]);
  }

  if (!valid_norms)
    calc_norms(0);

  double ovl = 0.;
  for (int s = 0; s < 2; s++) {
    for (int c1 = 0, ic1 = 0; c1 < ne[s]; ic1 += nspl[s][c1], c1++) {
      cdouble sum(0, 0);
      for (int j1 = 0; j1 < nspl[s][c1]; j1++) {
        double cj_re = split_c[s][ic1 + j1][0];
        double cj_im = split_c[s][ic1 + j1][1];
        cdouble cj = cdouble(cj_re, cj_im);
        WavePacket wj = wp[s][ic1 + j1];
        //sum+=cj*wj(x);
        sum += cj * wj.partial_integrate(x, integ_flag);
      }
      if (norm_mode == NORMALIZE)
        ovl += norm(sum) / wf_norm[s][c1];
      else
        ovl += norm(sum);
    }
  }
  return ovl / (ne[0] + ne[1]);
}


int AWPMD_split::set_var_constraint(int s, int electron_id, int wp_id, int param_id, int fix_flag) {
  if (s < 0 || s > 1)
    return -1;
  if (electron_id < 0 || electron_id >= ne[s])
    return -1;
  if (wp_id < 0 || wp_id >= nspl[s][electron_id])
    return -1;
  if (param_id < 0 || param_id >= vars_per_wp)
    return -1;
  //finding full variable id
  int var_id = 0;
  for (int c1 = 0; c1 < electron_id; c1++)
    var_id += nspl[s][c1];
  if (param_id < 8)
    var_id = 8 * (var_id + wp_id) + param_id;
  else
    var_id = 8 * nwp[s] + 2 * (var_id + wp_id) + (param_id - 8);

  set_var_constraint_by_id(s, var_id, fix_flag);


  return 1;
}


int AWPMD_split::set_var_constraint_by_id(int s, int var_id, int fix_flag) {
  int res = 0;
  map<int, int>::iterator it = fixed[s].find(var_id);
  bool found = (it != fixed[s].end());
  if (!fix_flag) { // deleting from map
    if (found) {
      if (it->second)
        res = 1;
      fixed[s].erase(it);
    }
  } else { //inserting or changing
    if (fix_flag < 0) { // removing conditional constraint only if it was previously set
      if (found && it->second == -fix_flag) {
        fixed[s].erase(it);
        res = 1;
      }
    } else {
      if (!found || (found && it->second != fix_flag)) {
        fixed[s][var_id] = fix_flag;
        res = 1;
      }
    }
  }
  return res;
}


//int max_constr=1000;


size_t AWPMD_split::get_constr_count() const {
  size_t count = 0;
  for (int s = 0; s < 2; s++) {
    for (int c = 0; c < ne[s]; c++) {
      int eid = s * ne[0] + c;
      count += lin_constraints[eid].size();
    }
  }
  return count;
}


int AWPMD_split::check_degeneracy_constraint(int action, double prj_tolerance, int max_constr_number) {
  double max_dq_dt_tol = 0. /*1e4*/, max_prj_tol = 1e3; //, min_prj_tol=1e3;
  double max_eval_spread = 1e20;
  int res = 0;

  if (prj_tolerance >= 0) // override default
    max_prj_tol = prj_tolerance;

  for (int s = 0; s < 2; s++) {

    ptr_stor_t<double> st(dq_dt[s]);
    Vector_G vdq_dt(st), vE_der(E_der[s]);
    double max_E_der = fabs(vE_der.maxabscoord());
    double max_dq_dt = fabs(vdq_dt.maxabscoord());
    if (max_dq_dt < max_dq_dt_tol) { // force seems ok
      if (!(action & RELEASE)) // nothing to do with norm matrix
        continue;
    }
    // else getting too large force, investigate exact Norm matrix conditions

    vector<double> rhs;
    int indw = 0; // starting index of the electron wp coordinates
    int indn = nwp[s] * 8; // starting index of the electron norm coordinates
    for (int c = 0; c < ne[s]; c++) {
      int imposed = 0, released = 0;
      int eid = s * ne[0] + c;
      vector<vector<double> > &constr = lin_constraints[eid]; // constraints for this electron
      if (action & RELEASE) {
        released += (int) constr.size();
        constr.clear();
      }
      if (action & IMPOSE || constr.size()) { // no need to invert matrix if there is notning to do with constraints
        constr.clear(); // TODO(valuev): replace this with one by one constraint addition
        int info;
        sqmatrix<double> matr;
        size_t nconstr;
        prepare_hartree_nm(s, c, indw, indn, matr, rhs, nconstr, false); // matrix with no linear constraints (original)
        int dim = (int) matr.size;
        if (max_constr_number < 0)
          max_constr_number = dim - 2;
        int max_constr = min(max_constr_number, dim - 2); // number of constraint should not be less than dim-2

        // making hermitian matrix from i*Norm
        int narr = dim * (dim + 1) / 2;
        vector<cdouble> arrp(narr); // matrix in packed storage
        for (int i = 0; i < dim; i++) {
          arrp[i + (i + 1) * i / 2] = cdouble(matr(i, i),
                                              0);  // truncated matrix parts contain real  1 on diagonal, otherwise 0
          for (int j = i + 1; j < dim; j++)
            arrp[i + (j + 1) * j / 2] = cdouble(0., matr(i, j));
        }

        int lwork = 2 * dim, lrwork = 2 * dim * dim + 5 * dim + 1, liwork = 5 * dim + 3;
        vector<double> rwork(lrwork);
        vector<int> iwork(liwork);
        vector<double> eigen_val(dim);
        vector<cdouble> eigen_vect(dim * dim), work(lwork);

        ZHPEVD("V", "U", &dim, (MKL_Complex16 *) &arrp[0], &eigen_val[0], (MKL_Complex16 *) &eigen_vect[0], &dim,
               (MKL_Complex16 *) &work[0], &lwork, &rwork[0], &lrwork, &iwork[0], &liwork, &info);
        if (info != 0) {
          LOGERR(info, fmt_iv("AWPMD_split.check_degeneracy_constraint: call to ZHPEVD failed (exitcode %d)!", info),
                 LINFO);
          return 0;
        }

        double max_eval = *(max_element(eigen_val.begin(), eigen_val.end()));
        Vector_G vrhs(rhs);
        vector<double> rhs_c = rhs; // constrained force
        //apply_force_constraints(s,c,rhs_c);
        Vector_G vrhs_c(rhs_c);
        vector<double> lhs(rhs.size(), 0.); // solution of equation (dq/dt)
        multimap<double, int> evec_sorter;
        for (int i = 0; i < dim / 2; i++) { // over eigen vectors  (only negative part of spectrum)
          //if(eigen_val[i]>0.) // only negative part of spectrum
          //continue;
          //if(i!=dim/2-1)
          //  continue;
          vector<double> prj_evectr(dim), prj_evecti(dim);
          for (int j = 0; j < dim; j++) {
            prj_evectr[j] = real(eigen_vect[dim * i + j]); // fortran ?
            prj_evecti[j] = imag(eigen_vect[dim * i + j]); // fortran ?
          }

          Vector_G vprj_evectr(prj_evectr), vprj_evecti(prj_evecti);
          double fi = (vprj_evecti * vrhs);
          double fr = (vprj_evectr * vrhs);
          double fi_c = (vprj_evecti * vrhs_c);
          double fr_c = (vprj_evectr * vrhs_c);
          //Vector_G vtot=fi*vprj_evectr - fr*vprj_evecti;
          Vector_G vtot_c = fi_c * vprj_evectr - fr_c * vprj_evecti;
          double force_comp = 2 * vtot_c.norm();
          double prj_scal = fabs(force_comp / eigen_val[i]);
          //double cross=vprj_evectr*vprj_evecti;
          evec_sorter.insert(make_pair(prj_scal, i));

          //for(int j=0;j<dim;j++)
          //lhs[j]+=2*(fi*prj_evectr[j]-fr*prj_evecti[j])/eigen_val[i];
        } //i
        for (multimap<double, int>::reverse_iterator it = evec_sorter.rbegin(); it != evec_sorter.rend(); ++it) {
          double prj_scal = it->first;
          int i = it->second;
          if ((prj_scal > max_prj_tol || fabs(eigen_val[i]) * max_eval_spread < max_eval) &&
              imposed < max_constr) { // need to put constraint
            vector<double> prj_evectr(dim), prj_evecti(dim);
            for (int j = 0; j < dim; j++) {
              prj_evectr[j] = real(eigen_vect[dim * i + j]); // fortran ?
              prj_evecti[j] = imag(eigen_vect[dim * i + j]); // fortran ?
            }
            int res = add_force_constraint(s, c, prj_evectr, 1e-20);
            if (res != 1) {
              LOGERR(-1, fmt_iv("AWPMD_split.check_degeneracy_constraint: adding dependent vector (r)?", 0), LINFO);
              return 0;
            }
            imposed++;
            res = add_force_constraint(s, c, prj_evecti, 1e-20);
            if (res != 1) {
              LOGERR(-1, fmt_iv("AWPMD_split.check_degeneracy_constraint: adding dependent vector (i)?", 0), LINFO);
              return 0;
            }
            imposed++;
          }
        }



        /*
        if(action&RELEASE){ // calculating actual dq/dt projection on a constraint
          Vector_G vlhs(lhs);

          for(size_t i=0;i<constr.size();){
            Vector_G vconstri(constr[i]), vconstri1(constr[i+1]) ;
            double prj=vlhs*vconstri; // ;
            double prj1=vlhs*vconstri1; // projecting
            prj=sqrt(prj*prj+prj1*prj1);
            if(prj<min_prj_tol){ // can remove constraint
              constr.erase(constr.begin()+i); // removing
              constr.erase(constr.begin()+i+1);
              res-=-2;
            }
            else
              i+=2;
          } //i
        }*/
      } // IMPOSE
      if (imposed)
        prepare_constraints(s, c);
      res += (imposed + released);

      indw += 8 * nspl[s][c]; // 8 variables in each wavepacket
      indn += 2 * nspl[s][c]; // 2 variables in each wp norm
    }//c
  }//s
  return res;
}


int AWPMD_split::step(double dt, int flag, int spin, const vector<double> *dq_dt_) {
  int s0 = spin, s1 = spin;
  if (spin < 0 || spin > 2) {
    s0 = 0;
    s1 = 1;
  }
  int k = 0;
  for (int s = s0; s <= s1; s++) {

    const vector<double> *dqp = dq_dt_;
    if (!dqp) {
      dqp = dq_dt + s;
      k = 0;
    }
    const vector<double> &dq = *dqp;
    k += 8 * nwp[s];

    for (int i = 0; i < nwp[s]; i++) {
      split_c[s][i][0] += dq[k++] * dt;
      split_c[s][i][1] += dq[k++] * dt;
    }
  }
  return AWPMD::step(dt, flag, spin, dq_dt_);
}


// axial projection
struct box_projection_term_ax {
  typedef cdouble result_t;
  AWPMD_split *awpmds;
  int num;
  int axis;

  box_projection_term_ax(AWPMD_split *awpmds_, int num_, int axis_) : awpmds(awpmds_), num(num_), axis(axis_) {}

  result_t operator()(const WavePacket &bra, const WavePacket &ket) const {
    // test: must give exact sum of electrons for all projections
    //WavePacket braket = conj(bra)*ket;
    //return braket.integral();

    cdouble left = awpmds->box.eigen_proj(num, bra, axis);
    cdouble right = awpmds->box.eigen_proj(num, ket, axis);

    return conj(right) * left; // the order reflects projection having wp at left: <wp|harm_n>
  }

  void add2(cdouble k, const WavePacket &bra, const WavePacket &ket) {}
};


double AWPMD_split::project_harm_state(int num, int axis) {
  if (axis >= 0 && axis < 3) { // projection onto given axis
    box_projection_term_ax term(this, num, axis);
    return real(exp_value_2(term));
  } else { // sum of all contributing states
    box_projection_term_sum term(this, num, false);
    return real(exp_value_2(term).v[0]);
  }
}





