//
// Created by cheshire on 10.03.17.
//

#ifndef AWPMD_DFT_LDA_HPP
#define AWPMD_DFT_LDA_HPP

#include "IApproximation.hpp"
#include <cmath>

class LDA: public IApproximation{
private:
    float m_cx;
    float m_a;
    float m_b;
    float A;
    float bA_coeff;
public:

    __devspec__ LDA(float cx, float a, float b){
        m_cx = -3.0f / (4.0f * (float)M_PI) * powf(9.0f * (float)M_PI / 4, 1.0f / 3.0f);
        m_a = -0.01554534543482745f;
        m_b = (b );
        A = (powf((4.0f * (float)(M_PI) / 3.0f), -1.0f / 3.0f) );
        bA_coeff = m_b / A;
        Type = ApproxType::T_LDA;
    }

    __devspec__ float derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const override;
    __devspec__ float energy(float rho, float spin_ratio) const override;
    __devspec__ float kinetic(float rho, float spin_ratio) const override;
};

#endif //AWPMD_DFT_LDA_HPP
