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
 * $Log: wpmd.cpp,v $
 * Revision 1.65  2015/11/11 18:35:37  valuev
 * removed mass factor in ii
 *
 * Revision 1.64  2015/11/10 16:29:45  valuev
 * added ion charge
 *
 * Revision 1.63  2015/10/20 15:01:13  valuev
 * switched to width velocity
 *
 * Revision 1.62  2015/10/09 19:33:26  valuev
 * modified AWPMC to output the same vars as AWPMD
 *
 * Revision 1.61  2015/03/13 11:19:43  valuev
 * fixes for ion dynamics
 *
 * Revision 1.60  2015/02/10 14:08:09  valuev
 * added step function
 *
 * Revision 1.59  2014/07/11 18:14:18  valuev
 * added overlap matrix log
 *
 * Revision 1.58  2014/07/09 16:22:52  valuev
 * removed reference to h_plank from wavepacket.h
 *
 * Revision 1.57  2014/07/08 17:07:21  valuev
 * prepared confined MC simulation with AWPMD nonsplit norm matrix
 *
 * Revision 1.56  2014/04/01 14:07:35  valuev
 * prepared variables for boundary potential
 *
 * Revision 1.55  2013/08/19 16:08:51  morozov
 * Fixed ion forces
 *
 * Revision 1.54  2013/07/09 09:22:02  kazeev
 * Fixed compilation for K100 update.
 *
 * Revision 1.53  2013/07/09 08:13:38  kazeev
 * Try fixing the callsignature.
 *
 * Revision 1.52  2012/09/20 05:38:31  valuev
 * optimizer (test version)
 *
 * Revision 1.51  2012/09/03 15:54:46  valuev
 * fixed sign in WPMD forces
 *
 * Revision 1.50  2012/06/29 10:50:13  valuev
 * added linear constraints
 *
 * Revision 1.49  2012/04/17 20:42:31  valuev
 * fixed Psi calculation, started adding e-i Debye screening
 *
 * Revision 1.48  2012/03/30 14:30:55  morozov
 * Added fixed width simulation mode
 *
 * Revision 1.47  2011/10/12 16:39:02  valuev
 * restructured  interaction function
 *
 * Revision 1.46  2011/10/08 11:21:14  valuev
 * added external force storage and external force power calculation
 *
 *
*******************************************************************************/
# include "wpmd.h"
# include <array>


// Calculates  derivative overlap matrix IDD
void OverlapDeriv::calc_der_overlap(bool self, cdouble cc1, cdouble c2) {
  cdouble part12 = cc1 * c2;
  cVector_3 I3 = I1 * ((bb_4a + 2.5) / w12.a);
  cdouble I4 = I0 * (bb_4a * (bb_4a + 5.) + 3.75) / w12.a / w12.a;

  // calculate derivatives <(phi_k)'_q_k | (phi_l)'_q_l>:
  IDD.set(0, 0, (I4 - (d1.l + d2.l) * I2 + d1.l * d2.l * I0) * part12);                // over a_k_re and a_l_re
  IDD.set(0, 1, i_unit * (I4 - (d1.l + d2.m) * I2 + d1.l * d2.m * I0) * part12);   // over a_k_re and a_l_im
  if (!self)
    IDD.set(1, 0, i_unit1 * (I4 + (d1.m - d2.l) * I2 - d1.m * d2.l * I0) * part12);  // over a_k_im and a_l_re
  else
    IDD.set(1, 0, conj(IDD(0, 1)));
  IDD.set(1, 1, (I4 + (d1.m - d2.m) * I2 - d1.m * d2.m * I0) * part12);            // over a_k_im and a_l_im

  for (int i = 0; i < 3; i++) {
    IDD.set(0, (i + 1) * 2,
            (-I3[i] + d1.l * I1[i] + d2.u[i] * (d1.l * I0 - I2)) * part12);              // over a_k_re and b_l_re
    IDD.set(0, (i + 1) * 2 + 1,
            i_unit1 * (I3[i] - d1.l * I1[i] + d2.v[i] * (I2 - d1.l * I0)) * part12);   // over a_k_re and b_l_im
    IDD.set(1, (i + 1) * 2,
            i_unit * (I3[i] + d1.m * I1[i] + d2.u[i] * (I2 + d1.m * I0)) * part12); // over a_k_im and b_l_re
    IDD.set(1, (i + 1) * 2 + 1,
            (-I3[i] - d1.m * I1[i] - d2.v[i] * (d1.m * I0 + I2)) * part12);            // over a_k_im and b_l_im
    if (!self) {
      IDD.set((i + 1) * 2, 0,
              (-I3[i] + d2.l * I1[i] + d1.u[i] * (d2.l * I0 - I2)) * part12);              // over b_k_re and a_l_re
      IDD.set((i + 1) * 2 + 1, 0,
              i_unit * (I3[i] - d2.l * I1[i] - d1.v[i] * (I2 - d2.l * I0)) * part12);   // over b_k_im and a_l_re
      IDD.set((i + 1) * 2, 1,
              i_unit1 * (I3[i] - d2.m * I1[i] + d1.u[i] * (I2 - d2.m * I0)) * part12); // over b_k_re and a_l_im
      IDD.set((i + 1) * 2 + 1, 1,
              (-I3[i] + d2.m * I1[i] - d1.v[i] * (d2.m * I0 - I2)) * part12);            // over b_k_im and a_l_im
    } else {
      IDD.set((i + 1) * 2, 0, conj(IDD(0, (i + 1) * 2)));              // over b_k_re and a_l_re
      IDD.set((i + 1) * 2 + 1, 0, conj(IDD(0, (i + 1) * 2 + 1)));   // over b_k_im and a_l_re
      IDD.set((i + 1) * 2, 1, conj(IDD(1, (i + 1) * 2))); // over b_k_re and a_l_im
      IDD.set((i + 1) * 2 + 1, 1, conj(IDD(1, (i + 1) * 2 + 1)));            // over b_k_im and a_l_im
    }

    for (int j = 0; j < 3; j++) {
      if (!self || j >= i) {
        cdouble I2ij = I0 / w12.a *
                       (i == j ? w12.b[i] * w12.b[i] / w12.a / 4 + 0.5
                               : w12.b[i] * w12.b[j] / w12.a / 4);
        // over b_k_re and b_l_re
        IDD.set((j + 1) * 2, (i + 1) * 2, (I2ij + d1.u[i] * I1[j] + d2.u[j] * (I1[i] + d1.u[i] * I0)) * part12);
        // over b_k_re and b_l_im
        IDD.set((j + 1) * 2, (i + 1) * 2 + 1,
                i_unit * (I2ij + d1.u[i] * I1[j] + d2.v[j] * (I1[i] + d1.u[i] * I0)) * part12);
        // over b_k_im and b_l_re
        if (!self || i != j)
          IDD.set((j + 1) * 2 + 1, (i + 1) * 2,
                  i_unit1 * (I2ij - d1.v[i] * I1[j] + d2.u[j] * (I1[i] - d1.v[i] * I0)) * part12);
        else
          IDD.set((j + 1) * 2 + 1, (i + 1) * 2, conj(IDD((i + 1) * 2, (j + 1) * 2 + 1)));
        // over b_k_im and b_l_im
        IDD.set((j + 1) * 2 + 1, (i + 1) * 2 + 1, (I2ij - d1.v[i] * I1[j] + d2.v[j] * (I1[i] - d1.v[i] * I0)) * part12);
      } else { // self && j<i
        // over b_k_re and b_l_re
        IDD.set((j + 1) * 2, (i + 1) * 2, conj(IDD((i + 1) * 2, (j + 1) * 2)));
        // over b_k_re and b_l_im
        IDD.set((j + 1) * 2, (i + 1) * 2 + 1, conj(IDD((i + 1) * 2 + 1, (j + 1) * 2)));
        // over b_k_im and b_l_re
        IDD.set((j + 1) * 2 + 1, (i + 1) * 2, conj(IDD((i + 1) * 2, (j + 1) * 2 + 1)));
        // over b_k_im and b_l_im
        IDD.set((j + 1) * 2 + 1, (i + 1) * 2 + 1, conj(IDD((i + 1) * 2 + 1, (j + 1) * 2 + 1)));
      }
    } // j
  } // i

  if (real(cc1)) { // adding terms for split-packet

    IDD.set(8, 0, c2 * da2_re());   // over c_1_re and a_2_re
    IDD.set(8, 1, c2 * da2_im());   // over c_1_re and a_2_im
    IDD.set(9, 0, -i_unit * c2 * da2_re());   // over c_1_im and a_2_re
    IDD.set(9, 1, -i_unit * c2 * da2_im());   // over c_1_im and a_2_im

    IDD.set(0, 8, cc1 * da1_re());   // over c_2_re and a_1_re
    IDD.set(1, 8, cc1 * da1_im());   // over c_2_re and a_1_im
    IDD.set(0, 9, i_unit * cc1 * da1_re());   // over c_2_im and a_1_re
    IDD.set(1, 9, i_unit * cc1 * da1_im());   // over c_2_im and a_1_im

    for (int i = 0; i < 3; i++) {
      IDD.set(8, 2 + 2 * i, c2 * db2_re(i));   // over c_1_re and b_2_re
      IDD.set(8, 2 + 2 * i + 1, c2 * db2_im(i));   // over c_1_re and b_2_im
      IDD.set(9, 2 + 2 * i, -i_unit * c2 * db2_re(i));   // over c_1_im and b_2_re
      IDD.set(9, 2 + 2 * i + 1, -i_unit * c2 * db2_im(i));   // over c_1_im and b_2_im

      IDD.set(2 + 2 * i, 8, cc1 * db1_re(i));   // over c_2_re and b_1_re
      IDD.set(2 + 2 * i + 1, 8, cc1 * db1_im(i));   // over c_2_re and b_1_im
      IDD.set(2 + 2 * i, 9, i_unit * cc1 * db1_re(i));   // over c_2_im and i_1_re
      IDD.set(2 + 2 * i + 1, 9, i_unit * cc1 * db1_im(i));   // over c_2_im and a_1_im
    }

    IDD.set(8, 8, I0);   // over c_1_re and c_2_re
    IDD.set(8, 9, i_unit * I0);   // over c_1_re and c_2_im
    IDD.set(9, 8, -i_unit * I0);   // over c_1_im and c_2_re
    IDD.set(9, 9, I0);   // over c_1_im and c_2_im
  }

}


WavePacket AWPMD::create_wp(const Vector_3 &x, const Vector_3 &v, double w, double pw, double mass, bool pw_is_vel) {
  if (mass < 0)
    mass = me;
  if (constraint == FIX && w0 > 0) {
    w = w0;
    pw = 0.;
  }

  double rw;
  if (Lextra > 0) { // width PBC, keeping the width are within [0,Lextra]
    w = fmod(w, Lextra);
    if (w < 0) w += Lextra;
    rw = w; // WP width for energy evaluation is within [0, L/2]
    if (rw > Lextra / 2) {
      rw = Lextra - rw;
      //pw = - pw; // reflecting the width momentum
    }
  } else
    rw = w;

  WavePacket wp;
  wp.init(rw, x, v * mass * one_h, pw_is_vel ? pw * mass * one_h : pw * one_h);
  return wp;
}


void AWPMD::resize(int flag) {
  for (int s = 0; s < 2; s++) {
    //0. resizing matrices
    Y[s].init(ne[s], 1);
    O[s].init(ne[s], 1);
    Oflg[s].init(ne[s], 1);
    //Te[s].init(nel,1);
    //Tei[s].init(nel,1);
    Eep[s].assign((size_t) nwp[s], 0);
    Eeip[s].assign((size_t) nwp[s], 0);
    Eeep[s].assign((size_t) nwp[s], 0);
    Ewp[s].assign((size_t) nwp[s], 0);

    if (flag & (0x8 | 0x4 | 0x20) && approx != HARTREE) { //electron forces or norm matrix, L and M are needed
      M[s].init(ne[s], nvar[s]);
      L[s].init(ne[s], nvar[s]);
    }

    if (flag & (0x8 | 0x4) || use_box) { //electron forces needed
      E_der[s].resize(nvar[s]);
      F_extra[s].resize(nvar[s]);
      F_wall[s].resize(nvar[s]);
      dq_dt[s].resize(nvar[s]);
      // resizing constraints
      if (approx == HARTREE)
        lin_constraints.resize(ne[0] + ne[1]); //separate set of constraints for each electron
      else
        lin_constraints.resize(2); //separate set of constraints for each spin

    }

  }
  Eiep.assign((size_t) ni, 0);
  Eiip.assign((size_t) ni, 0);


}


//e sets Periodic Boundary Conditions
//e using bit flags: 0x1 -- PBC along X
//e                  0x2 -- PBC along Y
//e                  0x4 -- PBC along Z
//e cell specifies the lengths of the simulation box in all directions
//e if PBCs are used, the corresponding coordinates of electrons and ions
//e in periodic directions must be within a range  [0, cell[per_dir])
//e @returns 1 if OK
int AWPMD::set_pbc(const Vector_3P pcell, int pbc_) {
  if (!pcell) {
    pbc = 0;
    center = Vector_3();
  } else {
    pbc = pbc_;
    cell = *pcell;
    for (int i = 0; i < 3; i++) {
      if (pbc && (0x1 << i))
        center[i] = cell[i] / 2;
      else
        center[i] = 0.;
    }
  }
# if DFT_EXTENSION
  dft_conf.mesh_start = {-0.9f*(float)cell[0], -0.9f*(float)cell[1], -0.9f*(float)cell[2]};
  dft_conf.mesh_fin = {-0.9f*(float)cell[0], -0.9f*(float)cell[1], -0.9f*(float)cell[2]};
  dft_conf.use_adaptive_mesh = true;

  if(ex) {
    ex->setConfig(dft_conf);
  }
# endif
  return 1;
}

// Setup electrons: forms internal wave packet representations
// If PBCs are used the coords must be within a range [0, cell)
int
AWPMD::set_electrons(int s, int n, const Vector_3P x, const Vector_3P v, const double *w, const double *pw, double mass,
                     const double *q, bool pw_is_vel) {
  if (s < 0 || s > 1)
    return LOGERR(-1, fmt_iv("AWPMD.set_electrons: invaid s setting (%d)!", s), LINFO);

  calc_state &= (CALC_ELECTRONS_SET | CALC_IONS_SET);
  calc_state |= CALC_ELECTRONS_SET;

  norm_matrix_state[s] = NORM_UNDEFINED;
  nwp[s] = ne[s] = n;
  nvar[s] = 8 * n;
  wp[s].resize(n);

  partition1[s].clear();
  for (int i = 0; i < n; i++) {
    wp[s][i] = create_wp(x[i], v[i], w[i], pw[i], mass, pw_is_vel);
    // assign default partition
    partition1[s].push_back(i + 1);
  }

  // assign electronic charge
  if (q)
    qe[s].assign(q, q + nwp[s]);
  else
    qe[s].assign(nwp[s], -1);


  return 1;
}

//e setup ion charges and coordinates
//e if PBCs are used the coords must be within a range [0, cell)
int AWPMD::set_ions(int n, double *q, Vector_3P x, double q0) {
  calc_state &= (CALC_ELECTRONS_SET | CALC_OVERLAPS | CALC_NORMS |
                 CALC_LMN); // makes all previously calculated energies invalid, electronic values remain
  calc_state |= CALC_IONS_SET;
  ni = n;
  qi.resize(n);
  xi.resize(n);
  partition1[2].clear();
  for (int i = 0; i < n; i++) {
    if (q)
      qi[i] = q[i] * q0;
    else
      qi[i] = q0;
    xi[i] = x[i];
    // assign default partition for ions
    partition1[2].push_back(i + 1);
  }

  return 1;
}

double AWPMD::interation_ii_single(int i, int j, double **x, double *q,
                                   double **f, double cutoff) {
  Vector_3 rij{x[i][0] - x[j][0], x[i][1] - x[j][1], x[i][2] - x[j][2]};

  double r = rij.norm();

  if (cutoff > 0 && r > cutoff)
    return 0;

  auto spline_scale = coulomb_cutoff(r, cutoff);

  double energy = coul_pref * q[i] * q[j] / r;
  auto smoothed_energy = energy * spline_scale.first;
  Eii += smoothed_energy;

  if (f) {
    auto denergy = energy / r;
    denergy = denergy * spline_scale.first - energy * spline_scale.second;

    Vector_3 df = denergy * rij / r;

    for (auto k = 0; k < 3; ++k){
      f[i][k] += df[k];
      f[j][k] -= df[k];
    }
  }
  return smoothed_energy;
}


double AWPMD::interaction_border_ion(int i, double *x, double *f) {
  double dE;
  Vector_3 df = box.get_force(*(Vector_3*)x, &dE);
  if (f) // ion forces needed
    for (auto k = 0; k < 3; ++k)
      f[k] += df[k];
  return dE;
}


double AWPMD::interaction_border_electron(WavePacket const &packet, double *force, double *erforce, double *ervforce) {
  double dE{0};
  if (force && erforce && ervforce) {
    cdouble integral;
    cdouble a1_re, a1_im, a2_re, a2_im;
    cVector_3 b1_re, b1_im, b2_re, b2_im;
    box.get_derivatives(packet.a, packet.b, packet.a, packet.b, &integral, &a1_re, &a1_im, &b1_re, &b1_im,
                        &a2_re, &a2_im, &b2_re, &b2_im);

    std::array<double, 8> tmp{2.0 * real(a1_re), 2.0 * real(a1_im),
                              2.0 * real(b1_re[0]), 2.0 * real(b1_im[0]),
                              2.0 * real(b1_re[1]), 2.0 * real(b1_im[1]),
                              2.0 * real(b1_re[2]), 2.0 * real(b1_im[2])};
                             
    auto dx = tmp.begin();
    auto dp = dx + 3;
    auto dw = dp + 3;
    auto pw = dw + 1;

    packet.int2phys_der<eq_second>(dx, dx, dp, dw, pw, 1. / one_h);
    for (auto k = 0u; k < 3; ++k)
      force[k] += -dx[k];
    (*erforce) += *dw;
    (*ervforce) += *pw;
    dE = integral.real();
  } else
    dE = box.get_integral(packet.a, packet.b, packet.a, packet.b).real();
  //Ebord += dE;
  return dE;
}

std::pair<double, double>
AWPMD::interaction_electron_kinetic(WavePacket const &packet, int spin, double *erforce, double *ervfroce) {
  return interaction_electron_kinetic(packet.get_width(), packet.get_pwidth(), spin, erforce, ervfroce);
}

std::pair<double, double>
AWPMD::interaction_electron_kinetic(double w, double pw, int spin, double *erforce, double *ervfroce){
  double pw_eng = pw * pw * h2_me / 2.0;
  Ee[spin] += pw_eng;

  Ew += h2_me * 9. / (8. * w * w);
  auto width_energy = h2_me * 9. / (8. * w * w);

  if (erforce != nullptr)
    *erforce += -2.0 * width_energy / w;

  if (ervfroce != nullptr)
    *ervfroce += h2_me * pw * one_h;
  return {pw_eng, width_energy};
}

double AWPMD::interaction_ee_single(WavePacket const &packet_1,
                                    WavePacket const &packet_2, double **eforce,
                                    double **erforce, double cutoff) {
  return interaction_ee_single(packet_1.get_r().get_ptr(), packet_1.get_width(),
                               packet_2.get_r().get_ptr(), packet_2.get_width(),
                               eforce, erforce, cutoff);
}

double AWPMD::interaction_ee_single(double const* coord1, double w1,
                             double const* coord2, double w2,
                             double **eforce, double **erforce,
                             double cutoff){
  Vector_3 v12{coord2[0] - coord1[0], coord2[1] - coord1[1], coord2[2] - coord1[2]};

  double r12 = v12.normalize();

  if (cutoff > 0 && r12 >= cutoff)
    return 0.0;

  auto scale_splines = coulomb_cutoff(r12, cutoff);

  double wsq = w1 * w1 + w2 * w2;
  double argw = sqrt((2. / 3.) * wsq);

  double Epot = coul_pref * erf_div(r12, 1. / argw);
  auto smoothed_energy = Epot * scale_splines.first;
  Eee += smoothed_energy;

  if (eforce){
    double derf_by_darg = coul_pref / r12 * two_over_sqr_pi * exp(-(r12 / argw) * (r12 / argw));
    double dEpot_by_dr = Epot / r12 - derf_by_darg / argw;

    //scale force
    dEpot_by_dr = dEpot_by_dr * scale_splines.first - Epot * scale_splines.second;

    for (auto k = 0; k < 3; ++k){
      eforce[0][k] -= v12[k] * dEpot_by_dr;
      eforce[1][k] += v12[k] * dEpot_by_dr;
    }

    if (erforce){
      *erforce[0] += -2.0 / 3.0 * derf_by_darg * (r12 / (argw * argw * argw)) * w1 * scale_splines.first;
      *erforce[1] += -2.0 / 3.0 * derf_by_darg * (r12 / (argw * argw * argw)) * w2 * scale_splines.first;
    }
  }

  return smoothed_energy;
}

double AWPMD::interaction_ei_single(double const *x, double q,
                                    WavePacket const &packet, int spin,
                                    double *f, double *eforce, double *erforce,
                                    double cutoff) {
  return interaction_ei_single(x, q, packet.get_r().get_ptr(), packet.get_width(), spin, f, eforce, erforce, cutoff);
}

double AWPMD::interaction_ei_single(double const *x, double q,
                             double const *packet_r, double w1, int spin,
                             double *f, double *eforce,
                             double *erforce, double cutoff){
  double argw = sqr_2_over_3 * w1;

  Vector_3 v12;
  for (auto k = 0; k < 3; ++k)
    v12[k] = x[k] - packet_r[k];

  if (pbc)
    v12 = v12.rcell1(cell, pbc);

  double r12 = v12.normalize();

  if (cutoff > 0 && r12 >= cutoff)
    return 0.0;

  auto scale_splines = coulomb_cutoff(r12, cutoff);

  double cel = -coul_pref * q;
  double Epot = cel * erf_div(r12, 1. / argw);
  auto smoothed_energy = Epot * scale_splines.first;

  Eei[spin] += smoothed_energy;

  if (f) {
    double arg = r12 / argw;
    double dEw = cel * two_over_sqr_pi * exp(-arg * arg) / argw;
    double dEpot = (Epot - dEw) / r12;

    dEpot = dEpot * scale_splines.first - Epot * scale_splines.second;

    for (auto k = 0; k < 3; ++k) {
      f[k] += v12[k] * dEpot;
      if (eforce)
        eforce[k] -= v12[k] * dEpot;
    }

  }

  if (erforce){
    double arg = r12 / argw;
    *erforce += -cel / r12 * two_over_sqr_pi * exp( -arg * arg) * arg / w1 * scale_splines.first;
  }

  return smoothed_energy;
}

//e same as interaction, but using Hartee factorization (no antisymmetrization)
int AWPMD::interaction_hartree(int flag, Vector_3P fi, Vector_3P fe_x,
                               Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c) {

  // 0. resizing the arrays if needed
  enum APPROX tmp = HARTREE;
  swap(tmp, approx); // do not neeed large matrices
  resize(flag);
  swap(tmp, approx);
  //1. clearing forces
  clear_forces(flag, fi, fe_x, fe_p, fe_w, fe_pw, fe_c);

  Ebord = Eext = Eee = Ew = 0.;
  for (int s1 = 0; s1 < 2; s1++) {
    Ee[s1] = 0.;
    Eei[s1] = 0.;
    for (int c1 = 0; c1 < ne[s1]; c1++) {
      // width part
      double w1 = wp[s1][c1].get_width();
      /*double sgn1=1;
      if(Lextra>0){ // width PBC
        if(w1>Lextra-w1){
          w1=-(Lextra-w1); // '-' is to change derivative sign
          sgn1=-1;
        }
      }*/
      Vector_3 r1 = wp[s1][c1].get_r();
      Vector_3 p = wp[s1][c1].get_p() * h_plank;
      Vector_3 pw = wp[s1][c1].get_pwidth() * h_plank;
      // energy contribution
      Ee[s1] += (p.norm2() + pw.norm2()) / (2 * me);
      Ew += h2_me * 9. / (8. * w1 * w1);
      if (constraint == HARM) Ew += harm_w0_4 * w1 * w1;
      // width force contribution
      //double dE=2*Epot/w;
      //if(d->fw1)d->fw1[c1]+=dE;
      //if(fw2 && fw2!=fw1)fw2[c1]+=dE;

      // e-e interaction
      for (int s2 = s1; s2 < 2; s2++) {
        for (int c2 = (s1 == s2 ? c1 + 1 : 0); c2 < ne[s2]; c2++) {
          double w2 = wp[s2][c2].get_width();
          Vector_3 v12 = wp[s2][c2].get_r() - r1;
          // position PBC
          v12 = v12.rcell1(cell, pbc);
          double r12 = v12.normalize();
          /*double sgn2=1; // signs
          if(Lextra>0){ // width PBC
            if(w2>Lextra-w2){
              w2=-(Lextra-w2); // '-' is to change derivative sign
              sgn2=-1;
            }
          }*/
          double wsq = w1 * w1 + w2 * w2;
          double argw = sqrt((2. / 3.) * wsq);

          //double arg=r12/argw;
          //double erfa=erf(arg);
          double Epot = coul_pref * erf_div(r12, 1. / argw);   //erfa/r12;
          Eee += Epot;

          // force contribution
          /*double dEw=coul_pref*two_over_sqr_pi*exp(-arg*arg)/argw;
          double dEpot=(Epot-dEw)/r12;
          if(!d->fixw){
            dEw/=wsq;
            if(d->fw1 && c1>=0){
              d->fw1[c1]+=sgn1*dEw*w1;
            }
            if(d->fw2){
              d->fw2[c2]+=sgn2*dEw*w2;
            }
          }*/
        }
      }
      // e-i interaction
      double wsq1 = w1 * w1;
      double argw = sqr_2_over_3 * w1;
      for (int i = 0; i < ni; i++) {
        Vector_3 v12 = xi[i] - r1;
        // position PBC
        v12 = v12.rcell1(cell, pbc);
        double r12 = v12.normalize();

        //double arg=r12/argw;
        //double erfa=erf(arg);
        double cel = -coul_pref * qi[i]; // electron charge is always -1
        double Epot = cel * erf_div(r12, 1. / argw); //erfa/r12;
        Eei[s1] += Epot;
        //printf("(%g %g %g)- (%g %g %g)\n",r1[0],r1[1],r1[2],xi[i][0],xi[i][1],xi[i][2]);
        //printf("awp(%d,%d:%d)[%g]: %g\n",i,s1,c1,r12,Epot);
        // force contribution
        if (flag & 0x3) {
          double arg = r12 / argw;
          double dEw = cel * two_over_sqr_pi * exp(-arg * arg) / argw;
          double dEpot = (Epot - dEw) / r12;
          fi[i] += v12 * dEpot; // ionic force
        }
        // electron force
        /*if(!d->fixw){
          dEw/=wsq;
          if(d->fw1 && c1>=0){
            d->fw1[c1]+=sgn1*dEw*w1;
          }
        }*/
      }
    }
  }
  if (calc_ii)
    interaction_ii(flag, fi);
  return 1;
}


//e initializes internal buffers for calculations (set_electrons must be called first)
//int init(){}

//e calculates interaction in the system of ni ions + electrons
//e the electonic subsystem must be previously setup by set_electrons, ionic by set_ions
//e the iterators are describing ionic system only
// 0x1 -- give back ion forces
// 0x2 -- add ion forces to the existing set
// 0x4 -- calculate derivatives for electronic time step (NOT IMPLEMENTED)
//e if PBCs are used the coords must be within a range [0, cell)
int AWPMD::interaction(int flag, Vector_3P fi, Vector_3P fe_x,
                       Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c) {
  if (approx == HARTREE)
    return interaction_hartree(flag, fi, fe_x, fe_p, fe_w, fe_pw, fe_c);
  // 0. resizing the arrays if needed
  resize(flag);
  // 1. clearing forces
  clear_forces(flag, fi, fe_x, fe_p, fe_w, fe_pw, fe_c);

  //2. calculating overlap matrix
  for (int s = 0; s < 2; s++) {
    int nes = ne[s];
    if (nes == 0) continue;

    for (int k = 0; k < nes; k++) {
      Y[s].set(k, k, 1.);  // Diagonal elements (=1)
      Oflg[s](k, k) = 1;
      for (int l = k + 1; l < nes; l++) {
        cdouble I0kl = pbc_mul(wp[s][l], wp[s][k]).integral(); // Non-diagonal elements
        Y[s].set(k, l, I0kl);
        Oflg[s](k, l) = (norm(I0kl) > ovl_tolerance);
      }
    }
    O[s] = Y[s];  // save overlap matrix

    //3. inverting the overlap matrix
    int info = 0;
    if (nes) {
      /*FILE *f1=fopen(fmt("matrO_%d.d",s),"wt");
      fileout(f1,Y[s],"%15g");
      fclose(f1);8*/

      ZPPTRF("L", &nes, Y[s].arr, &info);
      // analyze return code here
      if (info < 0)
        return LOGERR(info, fmt_iv("AWPMD.interacton: call to ZPTRF failed (exitcode %d)!", info), LINFO);
      ZPPTRI("L", &nes, Y[s].arr, &info);
      if (info < 0)
        return LOGERR(info, fmt_iv("AWPMD.interacton: call to ZPTRI failed (exitcode %d)!", info), LINFO);


      /*f1=fopen(fmt("matrY_%d.d",s),"wt");
      fileout(f1,Y[s],"%15g");
      fclose(f1);*/
    }

    // Clearing matrices for electronic forces
    if (flag & 0x4) {
      Te[s].Set(0.);
      Tei[s].Set(0.);
    }
  }

  Vector_3 ndr;
  //  calculating single particle contribution
  Wext = 0.;
  Eext = 0.;
  Ebord = 0.;
  for (int s = 0; s < 2; s++) {
    Ee[s] = Eei[s] = 0.;
    for (int k = 0; k < ne[s]; k++) {
      for (int l = k; l < ne[s]; l++) {

        if (!Oflg[s](k, l)) continue; // non-overlapping WPs

        // electrons kinetic energy
        WavePacket wk = wp[s][k];
        WavePacket &wl = wp[s][l];
        if (pbc)
          ndr = move_to_image(wl, wk);

        WavePacket wkl = wl * conj(wk);
        //Vector_3 rrkl=wkl.get_r();
        cVector_3 v1 = wl.b * conj(wk.a) - conj(wk.b) * wl.a;
        cdouble v = (v1 * v1) / wkl.a;
        v -= 6. * conj(wk.a) * wl.a;
        v /= wkl.a;
        cdouble I0kl = O[s](k, l);
        cdouble dE = -h2_me * I0kl * v / 2;
        if (flag & 0x4) // matrix needed only for electronic forces
          Te[s].set(k, l, dE);
        // energy component (trace)
        dE *= Y[s](l, k);
        Ee[s] += (l == k ? 1. : 2.) * real(dE);

        cVector_3 dkl = wkl.b / (2. * wkl.a);

        // e-i energy
        cdouble sum(0., 0.);
        for (int i = 0; i < ni; i++) {  // ions loop
          cVector_3 gkli = dkl - cVector_3(xi[i]);

          if (pbc) // correcting the real part (distance) according to PBC
            gkli = rcell1(gkli, cell, pbc);
          //-Igor- gkli=cVector_3(real(gkli).rcell1(cell,pbc),imag(gkli));

          cdouble ngkli = gkli.norm();
          cdouble c = sqrt(wkl.a);
          //cdouble ttt = cerf_div(ngkli,c);
          dE = -coul_pref * (qi[i]) * I0kl * cerf_div(ngkli, c);

          sum += dE;
          if (flag & 0x3) {// calculate forces on ions
            if (fabs(real(ngkli)) + fabs(imag(ngkli)) > 1e-10) {
              cdouble arg = ngkli * c;
              cdouble dEw = -coul_pref * qi[i] * I0kl * two_over_sqr_pi * exp(-arg * arg) * c;
              dE = (dE - dEw) / ngkli;
              dE *= Y[s](l, k);
              Vector_3 dir = -real(gkli);
              dir.normalize();
              fi[i] += (l == k ? 1. : 2.) * real(dE) * dir;
            }
          }
        }
        dE = sum;
        if (flag & 0x4) // matrix needed only for electronic forces
          Tei[s].set(k, l, dE);
        // energy component (trace)
        dE *= Y[s](l, k);
        Eei[s] += (l == k ? 1. : 2.) * real(dE);
      }
    }
  }

  // calculating e-e interaction
  Eee = Ew = 0.;
  // same spin
  for (int s = 0; s < 2; s++) { // spin
    for (int k = 0; k < ne[s]; k++) {  //c1
      for (int l = k + 1; l < ne[s]; l++) { //c3
        for (int m = k; m < ne[s]; m++) { //c2

          if (Oflg[s](k, m)) {
            WavePacket wkm = pbc_mul(wp[s][m], wp[s][k]);
            cVector_3 dkm = wkm.b / (2 * wkm.a);
            cdouble I0km = O[s](k, m);

            // kl-mn
            for (int n = l; n < ne[s]; n++) { //c4
              if (n <= m || !Oflg[s](l, n)) continue;

              WavePacket wln = pbc_mul(wp[s][n], wp[s][l]);
              if (pbc) // reducing the pair to elementary cell
                ndr = move_to_image(wkm,
                                    wln); // mind the derivative: wln.b+=wln.a*ndr, wln.lz+=-wln.a*ndr^2-i*wln.old_p*ndr;
              //Vector_3 rln=wln.get_r();

              cVector_3 dln = wln.b / (2 * wln.a);
              cdouble dd = (dkm - dln).norm();
              cdouble c = 1. / sqrt(1. / wln.a + 1. / wkm.a);
              cdouble Vklmn = coul_pref * I0km * O[s](l, n) * cerf_div(dd, c);

              //cdouble arge=dkm*dkm*wkm.a+dln*dln*wln.a+wkm.lz+wln.lz;
              //cdouble Vklmn=0.5*coul_pref*M_PI*M_PI*M_PI*exp(arge)*cerf_div(dd,c)/pow(wln.a*wkm.a,3./2.);
              cdouble dE = Vklmn * (Y[s](m, k) * Y[s](n, l) - Y[s](m, l) * Y[s](n, k));
              double rdE = real(dE);
              if (m != k || n != l) // not the same pair
                rdE *= 2;
              Eee += rdE;
            }//n
          }

          if (Oflg[s](l, m)) {
            WavePacket wlm = pbc_mul(wp[s][m], wp[s][l]);
            cVector_3 dlm = wlm.b / (2 * wlm.a);
            cdouble I0lm = O[s](l, m);

            // kl-nm
            for (int n = l; n < ne[s]; n++) {
              if (n <= m || !Oflg[s](k, n)) continue;

              WavePacket wkn = pbc_mul(wp[s][n], wp[s][k]);
              if (pbc) // reducing the pair to elementary cell
                ndr = move_to_image(wlm,
                                    wkn); // mind the derivative: wln.b+=wln.a*ndr, wln.lz+=-wln.a*ndr^2-i*wln.old_p*ndr;

              cVector_3 dkn = wkn.b / (2 * wkn.a);
              cdouble dd = (dkn - dlm).norm();
              cdouble c = 1. / sqrt(1. / wkn.a + 1. / wlm.a);
              cdouble Vklnm = coul_pref * I0lm * O[s](k, n) * cerf_div(dd, c);

              cdouble dE = Vklnm * (Y[s](n, k) * Y[s](m, l) - Y[s](n, l) * Y[s](m, k));
              double rdE = real(dE);
              if (m != k || n != l) // not the same pair
                rdE *= 2;
              Eee += rdE;
            }//n
          }
        }// m
      }// l
    }// k
  }// s

  // different spin
  for (int k = 0; k < ne[0]; k++) { // skm=0  //c1
    for (int l = 0; l < ne[1]; l++) {  // sln=1  //c3
      for (int m = k; m < ne[0]; m++) {         //c2
        if (Oflg[0](k, m)) {
          WavePacket wkm = pbc_mul(wp[0][m], wp[0][k]);
          cVector_3 dkm = wkm.b / (2 * wkm.a);
          cdouble I0km = O[0](k, m);

          for (int n = l; n < ne[1]; n++) { // km-ln   //c4
            if (Oflg[1](n, l)) {
              WavePacket wln = pbc_mul(wp[1][l], wp[1][n]);
              if (pbc) // reducing the pair to elementary cell
                ndr = move_to_image(wkm,
                                    wln); // mind the derivative: wln.b+=wln.a*ndr, wln.lz+=-wln.a*ndr^2-i*wln.old_p*ndr;
              //Vector_3 rln=wln.get_r();

              cVector_3 dln = wln.b / (2 * wln.a);
              cdouble dd = (dkm - dln).norm();
              cdouble c = 1. / sqrt(1. / wln.a + 1. / wkm.a);
              cdouble Vklmn = coul_pref * I0km * wln.integral() * cerf_div(dd, c);

              cdouble dE = Vklmn * Y[0](m, k) * Y[1](n, l);
              // correct version:
              int Mkm = (m == k ? 1 : 2);
              int Mln = (n == l ? 1 : 2);
              double rdE = Mkm * Mln * real(dE); //!!!
              /* old version:
              double rdE=real(dE);
              if(m!=k || n!=l) // not the same pair
                rdE*=2;*/

              Eee += rdE;
            } //if
          } // n
        } //if
      } // m
    }// l
  }// k
  if (calc_ii)
    interaction_ii(flag, fi);
  return 1;
}

///\en calculates \sum_ij w0/(1-ovl(i,j))*E0
double AWPMD::calc_overlap_penalty(double w0, double E0) {
  if (approx == HARTREE) // TODO: add penalty for blocks
    return 0.;
  double sum = 0.;
  for (int s = 0; s < 2; s++) {
    for (int k = 0; k < ne[s]; k++) {
      for (int l = k + 1; l < ne[s]; l++) {
        double ovl = abs(O[s](k, l));
        sum += ovl / (1. - ovl);
      }
    }
  }
  return E0 * w0 * sum;
}


//e Calculates Norm matrix and performs LU-factorization
//e The result is saved in AWPMD::Norm[s]
void AWPMD::norm_matrix(int s) {
  if (approx == HARTREE) {
    norm_matrix_state[s] = NORM_CALCULATED;
    return;
  }
  // Internal variables
  int k, l, i, j, qi, qj;
  int nes = ne[s], nes8 = nes * 8, nnes8 = nes * nes8;

  if (!nes) return;

  // References for frequently used arrays
  sqmatrix<double> &Norms = Norm[s];
  chmatrix &Ys = Y[s];
  smatrix<unsigned char> &Oflgs = Oflg[s];

  // Allocate of vectors and matrices
  Norms.init(nes8, 1);
  IDD.init(nes8, 1);
  if (ID.size() != nnes8)
    ID.resize(nnes8), IDYs.resize(nnes8), ipiv.resize(nes8);

  // Calculate first and second derivatives
  for (k = 0; k < nes; k++) {
    int k8 = k * 8;
    WavePacket &wk = wp[s][k];
    NormDeriv dk(wk);
    dk = conj(dk);  // conjugate: mu -> -mu, v -> -v !!!

    for (l = 0; l < nes; l++) {
      if (!Oflgs(k, l)) continue; // non-overlapping WPs

      int l8 = l * 8;
      WavePacket wl = wp[s][l];
      if (pbc) move_to_image(wk, wl);
      WavePacket wkl = conj(wk) * wl;
      NormDeriv dl(wl);

      // don't normalize O(k,l) here?
      cdouble I0 = O[s](k, l);
      cVector_3 I1 = wkl.b * (I0 / wkl.a / 2);
      cdouble bb_4a = wkl.b.norm2() / wkl.a / 4;
      cdouble I2 = I0 * (bb_4a + 1.5) / wkl.a;

      // calculate derivatives <phi_k | (phi_l)'_q_l>:
      int idx = k + l * nes8;
      if (k != l) {
        ID[idx] = dl.l * I0 - I2;             // over a_l_re
        ID[idx + nes] = i_unit * (dl.m * I0 - I2);    // over a_l_im
        for (i = 0; i < 3; i++) {
          ID[idx + ((i + 1) * 2) * nes] = dl.u[i] * I0 + I1[i];           // over b_l_re
          ID[idx + ((i + 1) * 2 + 1) * nes] = i_unit * (dl.v[i] * I0 + I1[i]);  // over b_l_im
        }
      } else { // k == l
        ID[idx] = i_unit * imag(dl.l);    // over a_l_re
        ID[idx + nes] = i_unit * (dl.m - I2);   // over a_l_im
        for (i = 0; i < 3; i++) {
          ID[idx + ((i + 1) * 2) * nes] = dl.u[i] + I1[i];  // over b_l_re
          ID[idx + ((i + 1) * 2 + 1) * nes] = 0.;               // over b_l_im
        }
      }

      if (k <= l) {
        cVector_3 I3 = I1 * ((bb_4a + 2.5) / wkl.a);
        cdouble I4 = I0 * (bb_4a * (bb_4a + 5.) + 3.75) / wkl.a / wkl.a;

        // calculate derivatives <(phi_k)'_q_k | (phi_l)'_q_l>:
        IDD.set(k8, l8, I4 - (dk.l + dl.l) * I2 + dk.l * dl.l * I0);                // over a_k_re and a_l_re
        IDD.set(k8, l8 + 1, i_unit * (I4 - (dk.l + dl.m) * I2 + dk.l * dl.m * I0));   // over a_k_re and a_l_im
        if (k != l)
          IDD.set(k8 + 1, l8, i_unit1 * (I4 + (dk.m - dl.l) * I2 - dk.m * dl.l * I0));  // over a_k_im and a_l_re
        IDD.set(k8 + 1, l8 + 1, I4 + (dk.m - dl.m) * I2 - dk.m * dl.m * I0);            // over a_k_im and a_l_im

        for (i = 0; i < 3; i++) {
          IDD.set(k8, l8 + (i + 1) * 2,
                  -I3[i] + dk.l * I1[i] + dl.u[i] * (dk.l * I0 - I2));              // over a_k_re and b_l_re
          IDD.set(k8, l8 + (i + 1) * 2 + 1,
                  i_unit1 * (I3[i] - dk.l * I1[i] + dl.v[i] * (I2 - dk.l * I0)));   // over a_k_re and b_l_im
          IDD.set(k8 + 1, l8 + (i + 1) * 2,
                  i_unit * (I3[i] + dk.m * I1[i] + dl.u[i] * (I2 + dk.m * I0))); // over a_k_im and b_l_re
          IDD.set(k8 + 1, l8 + (i + 1) * 2 + 1,
                  -I3[i] - dk.m * I1[i] - dl.v[i] * (dk.m * I0 + I2));            // over a_k_im and b_l_im
          if (k != l) {
            IDD.set(k8 + (i + 1) * 2, l8,
                    -I3[i] + dl.l * I1[i] + dk.u[i] * (dl.l * I0 - I2));              // over b_k_re and a_l_re
            IDD.set(k8 + (i + 1) * 2 + 1, l8,
                    i_unit * (I3[i] - dl.l * I1[i] - dk.v[i] * (I2 - dl.l * I0)));   // over b_k_im and a_l_re
            IDD.set(k8 + (i + 1) * 2, l8 + 1,
                    i_unit1 * (I3[i] - dl.m * I1[i] + dk.u[i] * (I2 - dl.m * I0))); // over b_k_re and a_l_im
            IDD.set(k8 + (i + 1) * 2 + 1, l8 + 1,
                    -I3[i] + dl.m * I1[i] - dk.v[i] * (dl.m * I0 - I2));            // over b_k_im and a_l_im
          }

          for (j = 0; j < 3; j++) {
            cdouble I2ij = I0 / wkl.a *
                           (i == j ? wkl.b[i] * wkl.b[i] / wkl.a / 4 + 0.5
                                   : wkl.b[i] * wkl.b[j] / wkl.a / 4);
            // over b_k_re and b_l_re
            IDD.set(k8 + (j + 1) * 2, l8 + (i + 1) * 2, I2ij + dk.u[i] * I1[j] + dl.u[j] * (I1[i] + dk.u[i] * I0));
            // over b_k_re and b_l_im
            IDD.set(k8 + (j + 1) * 2, l8 + (i + 1) * 2 + 1,
                    i_unit * (I2ij + dk.u[i] * I1[j] + dl.v[j] * (I1[i] + dk.u[i] * I0)));
            // over b_k_im and b_l_re
            if (k != l)
              IDD.set(k8 + (j + 1) * 2 + 1, l8 + (i + 1) * 2,
                      i_unit1 * (I2ij - dk.v[i] * I1[j] + dl.u[j] * (I1[i] - dk.v[i] * I0)));
            // over b_k_im and b_l_im
            IDD.set(k8 + (j + 1) * 2 + 1, l8 + (i + 1) * 2 + 1,
                    I2ij - dk.v[i] * I1[j] + dl.v[j] * (I1[i] - dk.v[i] * I0));
          } // j
        } // i
      } // if(k <= l)
    } // k
  } // l

  // Calculate matrix product IDYs_(k,q_j) = Ys_(k,l) * <phi_l | (phi_j)'_q_j>
  for (qj = 0; qj < nes8; qj++) {
    j = qj / 8;
    int idx = qj * nes;
    for (k = 0; k < nes; k++) {
      cdouble sum = 0.;
      for (l = 0; l < nes; l++)
        if (Oflgs(l, j)) sum += ID[idx + l] * Ys(k, l);
      IDYs[idx + k] = sum;
    }
  }

  // Calculate Norm-matrix
  for (qi = 0; qi < nes8; qi++) {
    i = qi / 8;
    int idxqi = qi * nes;

    Norms(qi, qi) = 0.;  // zero diagonal elements

    for (qj = qi + 1; qj < nes8; qj++) {
      j = qj / 8;
      int idxqj = qj * nes;

      // Calculate matrix product sum = <(phi_i)'_q_i | phi_k> * IDYs_(k,q_j)
      cdouble sum = 0.;
      for (k = 0; k < nes; k++)
        if (Oflgs(i, k))
          sum += IDYs[idxqj + k] * conj(ID[idxqi + k]);

      // Update norm-matrix taking into account its anti-symmetry
      double a = Oflgs(i, j) ?                            // IDD = 0 for non-overlapping WPs
                 h_plank2 * imag((sum - IDD(qi, qj)) * Ys(j, i)) :
                 h_plank2 * imag(sum * Ys(j, i));
      Norms(qi, qj) = a;
      Norms(qj, qi) = -a;
    } // qj
  } // qi

# if 1
  // transform norm matrix to the physical variables
  for (i = 0; i < nes; i++) {
    WavePacket wi = wp[s][i];
    for (k = 0; k < 8; k++) {
      // iterator to list all N(8*i+k,*) with fixed 8*i+k
      sqmatrix<double>::iterator mi = Norms.fix_first(8 * i + k, 0);
      for (j = i; j < nes; j++) {  // TO DO: run this loop from i+1 and take simplectic form for (i,i) block
        WavePacket wj = wp[s][j];
        wj.int2phys_der<eq_second>(mi + 8 * j, mi + 8 * j, mi + 8 * j + 3, mi + 8 * j + 6, mi + 8 * j + 7, h_plank);
      }
    }// finished line of blocks from right
    for (k = 8 * i; k < nes8; k++) { // TO DO: run this loop from 8*i+8 and take simplectic form for (i,i) block
      // iterator to list all N(8i+*,k) by fixed k
      sqmatrix<double>::iterator mi = Norms.fix_second(8 * i, k);
      wi.int2phys_der<eq_second>(mi, mi, mi + 3, mi + 6, mi + 7, h_plank);
    }// finished line of blocks from left

    for (k = 0; k < 8; k++) { // filling the lower triangle according to antisymmetry
      for (j = 8 * i + 8; j < nes8; j++)
        Norms(j, 8 * i + k) = -Norms(8 * i + k, j);
    }
  }
# endif
# if 0

  // transform norm matrix to the physical variables
  for(i=0;i<nes;i++){
    WavePacket wi=wp[s][i];
    for(j=i;j<nes;j++){
      WavePacket wj=wp[s][j];
      for(k=0;k<8;k++){
        // iterator to list all N(8*i+k,*) with fixed 8*i+k
        sqmatrix<double>::iterator mi=Norms.fix_first(8*i+k,8*j);
        wj.int2phys_der< eq_second >(mi,mi,mi+3,mi+6,mi+7);
      }
      for(k=0;k<8;k++){
        // iterator to list all N(8*i+k,*) with fixed 8*i+k
        sqmatrix<double>::iterator mi=Norms.fix_second(8*i,8*j+k);
        wi.int2phys_der< eq_second >(mi,mi,mi+3,mi+6,mi+7);
      }
      if(i!=j){
        for(int k1=0;k1<8;k1++){
          for(int k2=0;k2<8;k2++)
            Norms(8*j+k1,8*i+k2)=-Norms(8*i+k2,8*j+k1);
        }
      }
    }
  }
# endif

  norm_matrix_state[s] = NORM_CALCULATED;
}

///\en Norm matrix LU-factorization
int AWPMD::norm_factorize(int s) {
  if (norm_matrix_state[s] != NORM_CALCULATED)
    norm_matrix(s);
  if (approx == HARTREE) { // assuming simplectic form
    norm_matrix_state[s] = NORM_FACTORIZED;
    return 1;
  }

  int nes8 = ne[s] * 8, info;
  DGETRF(&nes8, &nes8, Norm[s].arr, &nes8, &ipiv[0], &info);
  if (info < 0)
    return LOGERR(info, fmt_iv("AWPMD.norm_factorize: call to DGETRF failed (exitcode %d)!", info), LINFO);

  norm_matrix_state[s] = NORM_FACTORIZED;
  return 1;
}


//e Norm matrix inversion
int AWPMD::norm_invert(int s) {
  if (norm_matrix_state[s] != NORM_FACTORIZED) norm_factorize(s);
  if (approx == HARTREE) { // assuming simplectic form, nothing to do
    norm_matrix_state[s] = NORM_INVERTED;
    return 1;
  }

  int nes8 = ne[s] * 8, info;
  int IDD_size = (int) IDD.get_datasize(nes8);

  DGETRI(&nes8, Norm[s].arr, &nes8, &ipiv[0], (double *) IDD.arr, &IDD_size, &info); // use IDD for work storage
  if (info < 0)
    return LOGERR(info, fmt_iv("AWPMD.norm_invert: call to DGETRI failed (exitcode %d)!", info), LINFO);

  norm_matrix_state[s] = NORM_INVERTED;
  return 1;
}


//e Get the determinant of the norm-matrix for the particles with spin s
double AWPMD::norm_matrix_det(int s) {
  if (approx == HARTREE) { // assuming simplectic form, nothing to do
    return 1.; // simplectic determinant
  }
  double det = 1.;
  int nes8 = ne[s] * 8;

  if (!nes8) return det;
  if (norm_matrix_state[s] != NORM_FACTORIZED)
    norm_factorize(s);

  sqmatrix<double> &Norms = Norm[s];
  for (int i = 0; i < nes8; i++)
    det *= Norms(i, i);  // product of the diagonal elements

  return det;
}


//e Get the determinant logarithm of the norm-matrix for the particles with spin s
double AWPMD::norm_matrix_detl(int s) {
  if (approx == HARTREE) { // assuming simplectic form, nothing to do
    return 0.; // simplectic determinant
  }
  int nes8 = ne[s] * 8;

  norm_det_log[s] = 0.;
  if (!nes8)
    return norm_det_log[s];
  if (norm_matrix_state[s] != NORM_FACTORIZED)
    norm_factorize(s);

  sqmatrix<double> &Norms = Norm[s];

  for (int i = 0; i < nes8; i++)
    norm_det_log[s] += log(fabs(Norms(i, i)));  // product of the diagonal elements

  return norm_det_log[s];
}


double AWPMD::get_energy(bool use_ee_hartree) {
  double res = (use_ee_hartree ? Eee_hartree : Eee) + Ew;
  for (int s = 0; s < 2; s++)
    res += Eei[s] + Ee[s];
  if (calc_ii)
    res += Eii;
  res += Ebord_ion; // electron border energy is included in Ee
  return res;
}


double AWPMD::get_coulomb_energy(bool use_ee_hartree) {
  double res = (use_ee_hartree ? Eee_hartree : Eee);
  for (int s = 0; s < 2; s++)
    res += Eei[s];
  if (calc_ii)
    res += Eii;
  res += Ebord_ion; // electron border energy is included in Ee
  return res;
}


double AWPMD::get_kin_energy() {
  return Ee[0] + Ee[1] + Ew - Ebord - Eext;  // excluding external ebnergies that are included in Ee
}


int AWPMD::step(double dt, int flag /*=0*/, int spin /*=-1*/, const vector<double> *dq_dt_/*=NULL*/) {
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

    if (flag & 0x10) { // working with internal variables
      for (int i = 0; i < nwp[s]; i++) {
        wp[s][i].a += cdouble(dq[k], dq[k + 1]) * dt;
        k += 2;
        wp[s][i].b += cVector_3(cdouble(dq[k], dq[k + 1]),
                                cdouble(dq[k + 2], dq[k + 3]),
                                cdouble(dq[k + 4], dq[k + 5])) * dt;
        k += 6;
      }
    } else { // working with physical variables
      for (int i = 0; i < nwp[s]; i++) {
        Vector_3 x = wp[s][i].get_r();
        Vector_3 p = wp[s][i].get_p() / one_h;
        double w = wp[s][i].get_width();
        double pw = wp[s][i].get_pwidth() / one_h;
        x += Vector_3(dq[k], dq[k + 1], dq[k + 2]) * dt;
        k += 3;
        p += Vector_3(dq[k], dq[k + 1], dq[k + 2]) * dt;
        k += 3;
        w += dq[k++];
        pw += dq[k++];
        wp[s][i].init(w, x, p * one_h, pw * one_h);
      }
    }
  }
  return 1;
}


//e gets current electronic coordinates
int AWPMD::get_electrons(int spin, Vector_3P x, Vector_3P v, double *w, double *pw, double mass) {
  if (spin < 0 || spin > 1)
    return -1; // invalid spin: return LOGERR(-1,fmt("AWPMD.get_electrons: invaid spin setting (%d)!",spin),LINFO);
  if (mass < 0)
    mass = me;
  for (int i = 0; i < ne[spin]; i++) {
    /*w[i]=sqrt(3./(4*real(wp[spin][i].a)));
    pw[i]=-2*w[i]*imag(wp[spin][i].a)/one_h; //pw[i]=-h_plank2*w[i]*imag(wp[spin][i].a);
    x[i]=real(wp[spin][i].b)/(2*real(wp[spin][i].a));
    v[i]=(pw[i]*x[i]/w[i] + imag(wp[spin][i].b)/one_h)/mass; //v[i]=(pw[i]*x[i]/w[i] + h_plank*imag(wp[spin][i].b))/m_electron;*/
    get_wavepacket(spin, i, x + i, v + i, w + i, pw + i);
    v[i] /= mass;
  }
  return 1;
}


void AWPMD::clear_forces(int flag, Vector_3P fi, Vector_3P fe_x,
                         Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c) {
  if (flag & 0x1) {
    for (int i = 0; i < ni; i++)
      fi[i] = Vector_3(0.);
  }
  if (flag & 0x4 && !(flag & 0x10)) { // electron forces requested in physical representation
    for (int s1 = 0; s1 < 2; s1++) { // clearing forces
      for (int c1 = 0; c1 < ne[s1]; c1++) {
        fe_x[c1] = Vector_3(0, 0, 0);
        fe_p[c1] = Vector_3(0, 0, 0);
        fe_w[c1] = 0;
        fe_pw[c1] = 0;
      }
    }
  }
}

int AWPMD::interaction_ii(int flag, Vector_3P fi) {
  Eii = 0.;
  size_t  count = 0;
  for (int i = 0; i < ni; i++) {
    for (int j = i + 1; j < ni; j++) {
      double M12e, M12f;
      _mytie(M12e, M12f) = check_part1ii(i, j);
      if (M12f) {
        Vector_3 rij = xi[i] - xi[j];
        if (pbc)
          rij = rij.rcell1(cell);

        double r = rij.norm();
        double dE = coul_pref * qi[i] * qi[j] / r;
        Eii += M12e * dE;

        Eiip[i] += 0.5 * M12e * dE;
        Eiip[j] += 0.5 * M12e * dE;

        if (flag & 0x3) { // ion forces needed
          Vector_3 df = M12f * dE * rij / (r * r);
          fi[i] += df;
          fi[j] -= df;
        }
      }
    }
  }
  return 1;
}

// adds confinement energies and forces to Ebord_ion and fi if needed
int AWPMD::calc_ion_confinement(int flag, Vector_3P fi) {
  Ebord_ion = 0.;
  for (int i = 0; i < ni; i++) {
    double dE;
    Vector_3 df = box.get_force(xi[i], &dE);
    //df*=100.;   ///m_electron; // normalizing to mass
    //dE*=100.;
    //Eii+=dE;
    Ebord_ion += dE;
    if (flag & 0x3) // ion forces needed
      fi[i] += df;
  }
  return 1;
}

Vector_3  AWPMD::calc_box_pressure() const {
  Vector_3 press;
  if (!use_box)
    return press;

  for (int s = 0; s < 2; s++) {
    for (int ic1 = 0; ic1 < ne[s]; ic1++) {
      int indw1 = 8 * ic1;
      //int indn1 = (nvar[s] / 10) * 8 + 2 * ic1;
      for (int i = 0; i < 3; i++) {
        press[i] += fabs(F_wall[s][indw1 + i]);
      }
    }

  }
  // ion contribution
  for (int i = 0; i < ni; i++) {
    double dE;
    Vector_3 df = box.get_force(xi[i], &dE);
    for (int j = 0; j < 3; j++)
      press[j] += fabs(df[j]);
  }

  Vector_3 floor = 2. * box.get_boundaries();
  for (int i = 0; i < 3; i++) {
    press[i] /= 2 * (floor[(i + 1) % 3] * floor[(i + 2) % 3]);
  }
  return press;
}



int interaction_relax(AWPMD *sys, double thresh, int flag) {
  wpmd_term energy_func(sys);
  md_grad_optimizer<wpmd_term> optimizer(&energy_func);
  optimizer.SetVerbose(0, stdout);
  optimizer.SetStep(0.01, 1.1);
  optimizer.Optimize(1000);
  // preparing normalized forces
  //if(optimizer.forces_valid()){
  //  sys->calc_LMN(flag);
  //  sys->
  // }
  return 1;
}
