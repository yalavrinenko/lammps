//
// Created by yalavrinenko on 09.02.2022.
//

#include "fix_nph_wpmd.h"
#include <error.h>
#include <group.h>
#include <modify.h>

LAMMPS_NS::FixNPHWpmd::FixNPHWpmd(LAMMPS_NS::LAMMPS *lammps, int i, char **pString) : FixNHWpmd(lammps, i, pString) {
  if (tstat_flag)
    error->all(FLERR,"Temperature control can not be used with fix nph/eff");
  if (!pstat_flag)
    error->all(FLERR,"Pressure control must be used with fix nph/eff");

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  // and thus its KE/temperature contribution should use group all

  id_temp = utils::strdup(std::string(id) + "_temp");
  modify->add_compute(fmt::format("{} all temp",id_temp));
  tcomputeflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  id_press = utils::strdup(std::string(id) + "_press");
  modify->add_compute(fmt::format("{} all pressure {}",id_press, id_temp));
  pcomputeflag = 1;
}
