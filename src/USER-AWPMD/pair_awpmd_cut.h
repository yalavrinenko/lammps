//
// Created by yalavrinenko on 16.03.2020.
//
#ifdef PAIR_CLASS

PairStyle(awpmd/cut,PairAWPMD)

#else
#ifndef LAMMPS_PAIR_AWPMD_CUT_H
#define LAMMPS_PAIR_AWPMD_CUT_H

#include "WavepacketPairCommon.h"
namespace LAMMPS_NS {
  class PairAWPMD: public WavepacketPairCommon {
  public:
    explicit PairAWPMD(class LAMMPS* lmp): WavepacketPairCommon(lmp) {
      throw std::runtime_error("Check pair! AWPMD pair req. tests and full impl.");
    }

    void settings(int i, char **pString) override;

  protected:
    struct awpmd_pair_index {
      int tag{};
      unsigned lmp_index{};
      unsigned wpmd_index{};

      awpmd_pair_index(unsigned l_index, unsigned w_index, int tag) : lmp_index(l_index), wpmd_index(w_index),
                                                                      tag(tag) {}

      bool operator < (awpmd_pair_index const &v) const{
        return this->tag < v.tag;
      }
    };

    using awpmd_ions = std::vector<awpmd_pair_index>;
    using awpmd_electrons = std::map<unsigned, std::vector<awpmd_pair_index>>;
    using awpmd_packets = std::tuple<awpmd_ions, awpmd_electrons>;

    awpmd_packets make_packets() const;
    void init_wpmd(awpmd_ions &ions, awpmd_electrons &electrons);

    awpmd_energies compute_energy_force() override;

    double electron_ke_{};
    double ermscale;
    double width_pbc;
  };
}


#endif //LAMMPS_PAIR_AWPMD_CUT_H
#endif