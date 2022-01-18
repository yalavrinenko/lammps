//
// Created by cheshire on 13.03.17.
//

#ifndef AWPMD_DFT_EXTERNALTYPES_HPP
#define AWPMD_DFT_EXTERNALTYPES_HPP

#include <complex>

#ifdef SEPARATE_DEBUG
typedef std::complex<float> cdouble;

class cVector_3{
public:
    cdouble v[3];

    cdouble& operator[] (int index){
        return this->v[index];
    }

    cdouble operator[] (int index) const{
        return this->v[index];
    }
};

class WavePacket{
public:
    cdouble a;
    cVector_3 b;
    cdouble lz;
    WavePacket(){
        throw "Wrong wavepacket";
    }
};
#else
#include <wavepacket.h>
#endif

#endif //AWPMD_DFT_EXTERNALTYPES_HPP
