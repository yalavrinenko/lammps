//
// Created by yalavrinenko on 28.12.18.
//

#include "fix_wpmc_awpmd.h"
#include "atom.h"
#include "force.h"
#include "error.h"
#include <systems/interact/TCP/awpmc.h>

namespace LAMMPS_NS {
  FixWPMCAwpmd::FixWPMCAwpmd(LAMMPS_NS::LAMMPS *lmp, int narg, char **args) : Fix(lmp, narg, args) {
    if (!atom->wavepacket_flag)
      error->all(FLERR, "Fix nve/awpmd requires atom style wavepacket");
  }
}
