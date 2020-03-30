//
// Created by yalavrinenko on 16.03.2020.
//

#ifndef LAMMPS_WAVEPACKETPAIRCOMMON_H
#define LAMMPS_WAVEPACKETPAIRCOMMON_H

#include "pair.h"
#include <vector>

class AWPMD_split;
class WavePacket;

namespace LAMMPS_NS{

class WavepacketPairCommon: public Pair {
public:
  explicit WavepacketPairCommon(class LAMMPS *);

  ~WavepacketPairCommon() override;

  void settings(int, char **) override;

  void coeff(int, char **) override;

  void init_style() override;

  double init_one(int, int) override;

  void write_restart(FILE *) override;

  void read_restart(FILE *fp) override;

  void write_restart_settings(FILE *) override;

  void read_restart_settings(FILE *fp) override;

  double memory_usage() override;

  void compute(int i, int i1) override;

  AWPMD_split *awpmd() { return wpmd; }

  std::vector<WavePacket> const& electrons_packets() const {
    return packets;
  }

protected:
  struct awpmd_energies{
    double ke{};
    double ee{};
    double ei{};
    double ii{};
    double ee_w{};

    double sum() const {
      return ke + ee + ei + ii + ee_w;
    }
  };

  //main compute function for pair!
  virtual awpmd_energies compute_energy_force() = 0;

  void allocate();

  void virial_eradius_compute();

  double cut_global{};
  int flexible_pressure_flag{0};
  double **cut = nullptr;

  double interaction_energy_{0};
  awpmd_energies energy_components_{};
  AWPMD_split *wpmd; // solver oybject
  std::vector<WavePacket> packets{};
};

}


#endif //LAMMPS_WAVEPACKETPAIRCOMMON_H
