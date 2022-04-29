//
// Created by cheshire on 21.05.17.
//

#ifndef AWPMD_DFT_IAPPROXIMATION_HPP
#define AWPMD_DFT_IAPPROXIMATION_HPP

#include <string>

#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif

#ifndef __global__
#define __global__
#endif

#define __devspec__ /*__host__*/ __device__
#define __isnan _isnan 

enum class ApproxType{
    T_LDA = 1,
    T_LDA_2 = 3,
    T_LSDA = 2,
    T_VOID = 0,
    T_DUMMY = 4
};

enum class ElectronSpin: int{
  E_UP = 1,
  E_DOWN = -1
};

class IApproximation{
    public:
        ApproxType Type = ApproxType::T_VOID;

        __devspec__ virtual float derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const = 0;
        __devspec__ virtual float energy(float rho, float spin_ratio) const = 0;
        __devspec__ virtual float kinetic(float rho, float spin_ratio) const = 0;

        virtual ~IApproximation() = default;
    };

class VoidApproximation: public  IApproximation{
public:
    __devspec__ VoidApproximation(){
        Type=ApproxType::T_VOID;
    }

    __devspec__ float derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const override {return 0.f;}
    __devspec__ float energy(float rho, float spin_ratio) const override {return 0.f;}
    __devspec__ float kinetic(float rho, float spin_ratio) const override {return 0.f;}
};

class DummyApproximation: public IApproximation{
public:
  __devspec__ DummyApproximation(){
    this->Type = ApproxType::T_DUMMY;
  }

  __devspec__ float derives(float rho, float spin_ratio, float derivative, enum ElectronSpin spin) const override {
    return 2.0f * derivative;
  }

  __devspec__ float energy(float rho, float spin_ratio) const override {
    return rho;
  }

  __devspec__ float kinetic(float rho, float spin_ratio) const override {
    return rho;
  }
};
#endif //AWPMD_DFT_IAPPROXIMATION_HPP
