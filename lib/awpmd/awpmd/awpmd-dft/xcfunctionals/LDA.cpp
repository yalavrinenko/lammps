//
// Created by cheshire on 10.03.17.
//
#include "LDA.hpp"

float LDA::energy(float rho, float spin_ratio) const {
    //rho = rho * powf(this->BohrR, 3);
    float rs = powf((4.0f * (float) M_PI * rho / 3.0f), -1.0f / 3.0f);

    return rho * ( m_a * logf(1.0f + m_b / rs  + m_b / (rs * rs))
                   + (m_cx * 1.0f / rs));
}

float LDA::derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const {
    float rho13 = powf(rho, 1.0f / 3.0f);
    float lnarg = 1.0f + bA_coeff * rho13 + bA_coeff * (rho13 * rho13) / A;
    return ((m_a * logf(lnarg) +
             (rho13 / 3.0f) * (m_a / lnarg) * bA_coeff * (1.0f + 2.0f * rho13 / A)) +
            (-m_cx * 4.0f / 3.0f * rho13)) * derivative;
}

__devspec__ float LDA::kinetic(float rho, float spin_ratio) const {
    if (rho < 1e-9f) //NaN in TKin fix. NaN in exc(rho) fix
        return 0.0f;

    float kF_up = 3.0f * float(M_PI) * float(M_PI) * rho;

    float Tkin = 0.0101f * powf(kF_up, 5.0f / 3.0f);
    return Tkin;
}