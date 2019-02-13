//
// Created by yalavrinenko on 11.02.19.
//

#include <force.h>
#include "fix_wpmc_wall.h"
#include "error.h"
#include "pair_awpmd_cut.h"
#include <wpmd_split.h>


LAMMPS_NS::FixWallAwpmd::FixWallAwpmd(LAMMPS_NS::LAMMPS *lammps, int i, char **pString) : Fix(lammps, i, pString) {
  m_pair = dynamic_cast<PairAWPMDCut*>(force->pair);
  if (!m_pair)
    error->all(FLERR, "Fix wall/awpmd require awpmd/cut pair_style.");

  m_pair->awpmd()->w0=force->numeric(FLERR, pString[3]);
  m_pair->awpmd()->set_harm_constr(m_pair->awpmd()->w0);
}

LAMMPS_NS::FixWallAwpmd::~FixWallAwpmd() {
  m_pair->awpmd()->constraint = AWPMD::NONE;
}

int LAMMPS_NS::FixWallAwpmd::setmask() {
  return 0;
}
