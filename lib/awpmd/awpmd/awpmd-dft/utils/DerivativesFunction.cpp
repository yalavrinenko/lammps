//
// Created by cheshire on 22.05.17.
//

#include "DerivativesFunction.hpp"

std::vector<DerivFunction> DerivsFunction::GetFunctions(){
    const int NFUNC = 4;
    std::vector<DerivFunction> funcs(NFUNC);

    funcs[0] = &(DerivsFunction::dr<0>); //dx
    funcs[1] = &(DerivsFunction::dr<1>); //dy

    funcs[2] = &(DerivsFunction::dr<2>); //dz
    funcs[3] = &(DerivsFunction::dw); //bw

    return  funcs;
}
