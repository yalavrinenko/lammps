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

  m_pair->awpmd()->set_box(construct_box(pString));
}

LAMMPS_NS::FixWallAwpmd::~FixWallAwpmd() {
  m_pair->awpmd()->constraint = AWPMD::NONE;
  m_pair->awpmd()->use_box = false;
}

int LAMMPS_NS::FixWallAwpmd::setmask() {
  return 0;
}

BoxHamiltonian LAMMPS_NS::FixWallAwpmd::construct_box(char **pString) {
  auto eigenwp = force->numeric(FLERR, pString[4]);
  auto eigenE = force->numeric(FLERR, pString[5]);

  auto floor = force->numeric(FLERR, pString[6]);
  double prj_ord = force->numeric(FLERR, pString[7]);

  if(eigenE>0.){
    eigenwp = sqrt(3./2/m_electron/eigenE)*h_plank;
  }
  else
    eigenE = 3./2*h_sq/m_electron/(eigenwp*eigenwp);


  double floorYtoX=1., floorZtoX=1., widthYtoX=1., widthZtoX=1.;

  Vector_3 gamma(eigenwp, eigenwp*widthYtoX, eigenwp*widthZtoX), force_k;

  double gamma_min = VEC_INFTY, gamma_max = 0.;

  for(int i=0; i<3; ++i){
    if(gamma_min>gamma[i])
      gamma_min = gamma[i];
    if(gamma_max<gamma[i])
      gamma_max = gamma[i];
    force_k[i] = 9./8*h_sq/m_electron/(gamma[i]*gamma[i]*gamma[i]*gamma[i]);
  }

  Vector_3 bound(floor,floor*floorYtoX,floor*floorZtoX);
  BoxHamiltonian  box(bound,force_k,(int)prj_ord);
  return box;
}
