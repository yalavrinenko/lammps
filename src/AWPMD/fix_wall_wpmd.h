//
// Created by yalavrinenko on 11.02.19.
//
#ifdef FIX_CLASS

FixStyle(wall/wpmd,FixWallWpmd)

#else
#ifndef LAMMPS_FIX_WPMC_WALL_H
#define LAMMPS_FIX_WPMC_WALL_H

#include "fix.h"
#include <box_hamiltonian.h>
#include <memory>
#include <array>
class WavePacket;

namespace LAMMPS_NS {
  class WavepacketPairCommon;
  class FixWallWpmd : public Fix {
  public:
    FixWallWpmd(class LAMMPS *lammps, int i, char **pString);

    ~FixWallWpmd() override;

    int setmask() override;

    double compute_scalar() override;

    double compute_vector(int i) override;

    void post_force(int i) override;

    void setup(int i) override;

  private:
    std::unique_ptr<BoxHamiltonian> construct_box(char **pString, double half_box_length, int pcount);

    double wall_pressure_ = 0;

    double wall_pressure() const;

    double interaction_border_ion(int i, double *x, double *f);
    double interaction_border_electron(WavePacket const &packet, double *rforce, double *erforce, double *ervforce);

    void evaluate_wall_energy(std::vector<WavePacket> const &packets);

    class WavepacketPairCommon* m_pair;

    std::unique_ptr<BoxHamiltonian> box = nullptr;
    double wall_energy = 0;
    std::array<double, 4> wall_pressure_components = {0., 0., 0., 0.} ;
    unsigned int walls_count_ = 3;

    std::array<double, 3> wall_squares{};

    std::array<bool, 3> has_force_{true, true, true};

    bool use_width_force_{false};

    std::vector<class WavePacket> packets;
  };
}

#endif //LAMMPS_FIX_WPMC_WALL_H
#endif
