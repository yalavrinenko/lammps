//
// Created by cheshire on 10.07.17.
//

#include "ModLDA.hpp"


inline float ModLDA::eps_xc(float rs) const{
    return  ( (rs >= 1.0f) ? -0.9164f / rs - 0.2846f / (1.0f + 1.0529f*sqrtf(rs) + 0.3334f*rs) :
                        -0.9164f / rs - 0.0960f + 0.0622f * logf(rs) - 0.0231f * rs + 0.0040f*rs*logf(rs)
    );
}

float ModLDA::energy(float rho, float spin_ratio) const{
    float rs = powf((4.0f * (float)M_PI * rho / 3.0f), -1.0f / 3.0f);

    return 2.0f * rho  * eps_xc(rs);
}

inline float ModLDA::deps_xc(float rs, float drs) const{
    return (rs >= 1.0f) ?
           0.9164f / (rs * rs) * drs + 0.2846f / powf((1.0f + 1.0529f*sqrtf(rs) + 0.3334f*rs), 2.0f) *
                                              (1.0529f * 0.5f / sqrtf(rs) * drs + 0.3334f * drs):
           0.9164f / (rs * rs) * drs + 0.0622f * 1.0f / rs * drs - 0.0232f * drs + 0.0040f * (drs * logf(rs) + drs);
}

float ModLDA::derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const{
    float rs = powf((4.0f * (float) M_PI * rho / 3.0f), -1.0f / 3.0f);
    float drs = -1.0f / 3.0f * powf((4.0f * (float) M_PI / 3.0f), (-1.0f / 3.0f)) * powf(rho, -4.0f / 3.0f) * derivative;


    return derivative * eps_xc(rs) + rho * deps_xc(rs, drs);
}


__devspec__ float ModLDA::kinetic(float rho, float spin_ratio) const{
    if (rho < 1e-9f) //NaN in TKin fix. NaN in exc(rho) fix
        return 0.0f;

    float kF_up = 3.0f * float(M_PI) * float(M_PI) * rho;

    float Tkin = 0.0101f * powf(kF_up, 5.0f / 3.0f);
    return Tkin;
}