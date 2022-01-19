//
// Created by yalavrinenko on 19.09.2019.
//

#include "fix_nvt_wpmd.h"
#include "fix_nve_wpmd.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include <modify.h>
#include <group.h>

LAMMPS_NS::FixNVTWpmd::FixNVTWpmd(LAMMPS_NS::LAMMPS *lammps, int argc, char **argv) :
    FixNHWpmd(lammps, argc, argv){
  if (!atom->wavepacket_flag)
    error->all(FLERR,"Fix nve/awpmd requires atom style wavepacket");

  if (!tstat_flag)
    error->all(FLERR,"Temperature control must be used with fix nvt");
  if (pstat_flag)
    error->all(FLERR,"Pressure control can not be used with fix nvt");

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} {} temp",id_temp,group->names[igroup]));
  tcomputeflag = 1;
}