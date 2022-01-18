//
// Created by cheshire on 10.03.17.
//
#include "LSDA.hpp"

#define PI_4_3 powf((4.0f * M_PI / 3.0),(-1.0f / 3.0f))

__devspec__ inline float LSDA::F(float kappa) const {
  return (powf(1.0f + kappa, 4.0f / 3.0f) + powf(1.0f - kappa, 4.0f / 3.0f) - 2.0f) / (2.0f * (ex_c_fero - 1.0f));
}

__devspec__ inline float LSDA::dF(float kappa) const {
  return 2.0f / (3.0f * (ex_c_fero - 1.0f)) * ( powf(1.0f + kappa, 1.0f / 3.0f) - powf(1.0f - kappa, 1.0f / 3.0f));
}

/***CHECKED***/
__devspec__ float LSDA::energy(float rho, float spin_ratio) const {
  if (rho < 1e-9f) //NaN in TKin fix. NaN in exc(rho) fix
    return 0.0f;

  float rs = powf(rho, -1.0f / 3.0f) * PI_4_3;

  float ex_para = ex_c_para * (1.0f / rs);
  float ex_fero = ex_para * ex_c_fero;

  return rho * (ex_para + F(spin_ratio) * (ex_fero - ex_para) +
                ((rs >= 1.0f) ? ecorr_para(rs) + F(spin_ratio) * (ecorr_fero(rs) - ecorr_para(rs)) :
                 ecorr_para_hd(rs) + F(spin_ratio) * (ecorr_fero_hd(rs) - ecorr_para_hd(rs))
                ));
}

__devspec__ float LSDA::kinetic(float rho, float spin_ratio) const {
  if (rho < 1e-9f) //NaN in TKin fix. NaN in exc(rho) fix
    return 0.0f;

  float Tkin = 2.862132316206293f *
               (powf(rho * (1.0f - spin_ratio), 5.0f / 3.0f) + powf(rho * (1.0f + spin_ratio), 5.0f / 3.0f));
  return Tkin;
}

__devspec__ float LSDA::derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const {
  if (rho < 1e-9f) //NaN in TKin fix. NaN in exc(rho) fix
    return 0.0f;

  auto Erho = this->energy(rho, spin_ratio);

  float rs = powf(rho, -1.0f / 3.0f) * PI_4_3;
  float drs = -1.0f / 3.0f * PI_4_3 * powf(rho, -4.0f / 3.0f);

  float dspin_ratio_drhp = (spin == ElectronSpin::E_UP) ? (1.0f - spin_ratio) / rho : (-spin_ratio - 1.0f) / rho;

  float Dx = -ex_c_para / (rs * rs) * (1 + F(spin_ratio) * (ex_c_fero - 1)); //dEps_exchange
  float spinDx = dF(spin_ratio) * dspin_ratio_drhp * (ex_c_fero - 1) * ex_c_para * (1.0f / rs);

  float Dc = (
      (rs >= 1.0f) ?
      d_ecorr_para(rs) + F(spin_ratio) * (d_ecorr_fero(rs) - d_ecorr_para(rs)) :
      d_ecorr_para_hd(rs) + F(spin_ratio) * (d_ecorr_fero_hd(rs) - d_ecorr_para_hd(rs))
  ); //dEps_correlation

  float spinDc = dF(spin_ratio) * dspin_ratio_drhp *
                 ((rs >= 1.0f) ? (ecorr_fero(rs) - ecorr_para(rs)) : (ecorr_fero_hd(rs) - ecorr_para_hd(rs)));


  float DTk = 9.540441054020977f * (
      (spin == ElectronSpin::E_UP) ? powf(rho * (1.0f + spin_ratio), 2.0f / 3.0f) :
      powf(rho * (1.0f - spin_ratio), 2.0f / 3.0f));

  return derivative * (
      (Erho / rho + rho * (drs * (Dx + Dc) + spinDx + spinDc)) + DTk);
}
