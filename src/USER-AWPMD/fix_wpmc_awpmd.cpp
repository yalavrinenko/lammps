//
// Created by yalavrinenko on 28.12.18.
//

#include "fix_wpmc_awpmd.h"
#include "atom.h"
#include "force.h"
#include "error.h"
#include "memory.h"
namespace LAMMPS_NS {
  FixWPMCAwpmd::FixWPMCAwpmd(LAMMPS_NS::LAMMPS *lmp, int narg, char **args) : Fix(lmp, narg, args) {
    if (!atom->wavepacket_flag)
      error->all(FLERR, "Fix wpmc/awpmd requires atom style wavepacket");

    vector_flag = 1;
    size_vector = 1;
    global_freq = 1;
    extvector = 0;
    time_depend = 1;

    nevery = 1;

    memory->create(vector_local, 1, "wpmc/awpmd::vector_local");
    memory->create(vector_atom, 1, "wpmc/awpmd::vector_atom");
  }

  void FixWPMCAwpmd::init() {
  }

  double FixWPMCAwpmd::memory_usage() {
    return 0;
  }

  void FixWPMCAwpmd::initial_integrate(int i) {
    static int times = 0;
    std::cout << "WE IN " << __FUNCTION__ << " " << times  << " " << MPI::COMM_WORLD.Get_rank() << std::endl;
    vector_local[0] = 1.0;
    vector_atom[0] = 2.0;
    ++times;
  }

  void FixWPMCAwpmd::post_force(int i) {
    static int times = 0;
    std::cout << "WE IN " << __FUNCTION__ << " " << times  << " " << MPI::COMM_WORLD.Get_rank() << std::endl;
    vector_local[0] = 1.0;
    vector_atom[0] = 2.0;
    ++times;
  }

  FixWPMCAwpmd::~FixWPMCAwpmd() {
    memory->sfree(vector_atom);
    memory->sfree(vector_local);
  }

  double FixWPMCAwpmd::compute_vector(int i) {
    return vector_local[i];
  }
}
