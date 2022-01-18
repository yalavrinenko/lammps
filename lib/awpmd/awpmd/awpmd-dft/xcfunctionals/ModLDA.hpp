//
// Created by cheshire on 10.07.17.
//

#ifndef AWPMD_DFT_MODLDA_HPP
#define AWPMD_DFT_MODLDA_HPP

#include "IApproximation.hpp"
#include <cmath>

class ModLDA: public IApproximation{
private:
    __devspec__  inline float eps_xc(float rs) const;
    __devspec__  inline float deps_xc(float rs, float drs) const;
public:

    __devspec__ ModLDA(){
        Type = ApproxType::T_LDA_2;
    }

    __devspec__ float derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const override ;
    __devspec__ float energy(float rho, float spin_ratio) const override ;
    __devspec__ float kinetic(float rho, float spin_ratio) const override;
};

#endif //AWPMD_DFT_MODLDA_HPP
