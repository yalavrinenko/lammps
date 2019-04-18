//
// Created by yalavrinenko on 11.02.19.
//

#include <force.h>
#include "fix_wpmc_wall.h"
#include "error.h"
#include "pair_awpmd_cut.h"
#include <wpmd_split.h>
#include "domain.h"

LAMMPS_NS::FixWallAwpmd::FixWallAwpmd(LAMMPS_NS::LAMMPS *lammps, int i, char **pString) : Fix(lammps, i, pString) {
  m_pair = dynamic_cast<PairAWPMDCut*>(force->pair);
  if (!m_pair)
    error->all(FLERR, "Fix wall/awpmd require awpmd/cut pair_style.");

  double delx = domain->boxhi[0]-domain->boxlo[0];
  double dely = domain->boxhi[1]-domain->boxlo[1];
  double delz = domain->boxhi[2]-domain->boxlo[2];
  auto half_box_length = 0.5 * MIN(delx, MIN(dely, delz));

  Vector_3 box_size{delx, dely, delz};
  m_pair->awpmd()->set_pbc(&box_size, 0);
  m_pair->awpmd()->set_box(construct_box(pString, half_box_length));
}

LAMMPS_NS::FixWallAwpmd::~FixWallAwpmd() {
  m_pair->awpmd()->use_box = false;
  m_pair->awpmd()->set_pbc(nullptr, 0);
}

int LAMMPS_NS::FixWallAwpmd::setmask() {
  return 0;
}

BoxHamiltonian LAMMPS_NS::FixWallAwpmd::construct_box(char **pString, double half_box_length) {
  auto eigenE = force->numeric(FLERR, pString[4]);
  double prj_ord = force->numeric(FLERR, pString[5]);


  auto floor = half_box_length;
  auto eigenwp = half_box_length / 10.0;

  auto me=force->e_mass;
  auto h2_me=force->hhmrr2e/force->e_mass;
  auto one_h=force->mvh2r;


  if(eigenE>0.){
    eigenwp = sqrt(3./2/me/eigenE) / one_h;
  }
  else
    eigenE = 3./2 * h2_me/(eigenwp*eigenwp);


  double floorYtoX=1., floorZtoX=1., widthYtoX=1., widthZtoX=1.;

  Vector_3 gamma(eigenwp, eigenwp*widthYtoX, eigenwp*widthZtoX), force_k;

  for(int i=0; i<3; ++i){
    force_k[i] = 9./8 * h2_me/(gamma[i]*gamma[i]*gamma[i]*gamma[i]);
  }

  Vector_3 bound(floor,floor*floorYtoX,floor*floorZtoX);
  BoxHamiltonian  box(bound,force_k,(int)prj_ord);

  return box;
}
