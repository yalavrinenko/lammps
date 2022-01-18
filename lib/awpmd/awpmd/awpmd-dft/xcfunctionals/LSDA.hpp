//
// Created by cheshire on 10.03.17.
//

#ifndef AWPMD_DFT_LSDA_HPP
#define AWPMD_DFT_LSDA_HPP

#define SQR(A) ((A)*(A))

#include "IApproximation.hpp"
#include <cmath>

class LSDA: public IApproximation{
private:
    float ex_c_para;
    float ex_c_fero;
    float K;

public:
    __devspec__ inline float F(float kappa) const;

  __devspec__ inline float dF(float kappa) const;

    /***CHECKED***/
    __devspec__ static inline float ecorr_para(float rs) {
        return -0.1423f / (1.0f + 1.0529f * sqrtf(rs) + 0.3334f * rs);
    }
    __devspec__ static inline float d_ecorr_para(float rs) {
        return -ecorr_para(rs) * ecorr_para(rs) / (-0.1423f) * (0.5f * 1.0529f / sqrtf(rs) + 0.3334f);
    }

    /***CHECKED***/
    __devspec__ static inline float ecorr_fero(float rs) {
        return -0.0843f / (1.0f + 1.3981f * sqrtf(rs) + 0.2611f * rs);
    }
    __devspec__ static inline float d_ecorr_fero(float rs) {
        return -ecorr_fero(rs) * ecorr_fero(rs) / (-0.0843f) * (0.5f * 1.3981f / sqrtf(rs) + 0.2611f);
    }

    /***CHECKED***/
    __devspec__ static inline float ecorr_para_hd(float rs) {
        return 0.0311f * logf(rs) - 0.0480f + 0.0020f * rs * logf(rs) - 0.0116f * rs;
    }
    __devspec__ static inline float d_ecorr_para_hd(float rs) {
        return 0.0311f / rs + 0.0020f * (logf(rs) + 1) - 0.0116f;
    }

    /***CHECKED***/
    __devspec__ static inline float ecorr_fero_hd(float rs) {
        return 0.01555f * logf(rs) - 0.0269f + 0.0007f * rs * logf(rs) - 0.0048f * rs;
    }
    __devspec__ static inline float d_ecorr_fero_hd(float rs) {
        return 0.01555f / rs + 0.0007f * (logf(rs) + 1) - 0.0048f;
    }
public:
    __devspec__ LSDA(){
        ex_c_fero = powf(2.0f, 1.0f / 3.0f);
        ex_c_para = -0.9164f / 2.0f;

        K = powf((3.0f / (4.0f * (float) M_PI)), 1.0f / 3.0f);
        Type = ApproxType::T_LSDA;
    }
    __devspec__ float derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const override ;

    __devspec__ float energy(float rho, float spin_ratio) const override ;

    __devspec__ float kinetic(float rho, float spin_ratio) const override;
};

#endif //AWPMD_DFT_LSDA_HPP
