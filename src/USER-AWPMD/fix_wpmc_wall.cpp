//
// Created by yalavrinenko on 11.02.19.
//

#include <force.h>
#include "fix_wpmc_wall.h"
#include "error.h"
#include "pair_awpmd_cut.h"
#include <wpmd_split.h>
#include <atom.h>
#include "domain.h"

LAMMPS_NS::FixWallAwpmd::FixWallAwpmd(LAMMPS_NS::LAMMPS *lammps, int i, char **pString) : Fix(lammps, i, pString) {
  double delx = domain->boxhi[0]-domain->boxlo[0];
  double dely = domain->boxhi[1]-domain->boxlo[1];
  double delz = domain->boxhi[2]-domain->boxlo[2];
  auto half_box_length = 0.5 * MIN(delx, MIN(dely, delz));

  Vector_3 box_size{delx, dely, delz};

  m_pair = dynamic_cast<PairAWPMDCut*>(force->pair);
  if (m_pair) {
    m_pair->awpmd()->set_pbc(&box_size, 0);
    m_pair->awpmd()->set_box(*construct_box(pString, half_box_length));
  } else {
    this->box = construct_box(pString, half_box_length);
  }
}

LAMMPS_NS::FixWallAwpmd::~FixWallAwpmd() {
  m_pair->awpmd()->use_box = false;
  m_pair->awpmd()->set_pbc(nullptr, 0);
}

int LAMMPS_NS::FixWallAwpmd::setmask() {
  return LAMMPS_NS::FixConst::PRE_REVERSE;
}

std::unique_ptr<BoxHamiltonian> LAMMPS_NS::FixWallAwpmd::construct_box(char **pString, double half_box_length) {
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
  return std::unique_ptr<BoxHamiltonian>(new BoxHamiltonian(bound,force_k,(int)prj_ord));
}

void LAMMPS_NS::FixWallAwpmd::pre_reverse(int, int) {
  if (box){
    auto const gamma_scale = 0.8660254037844385; //1.0 / (2.0 / std::sqrt(3.0)) EFF packet_width to WPMD packet_width
    auto const one_h=force->mvh2r;

    auto a_coeff = [gamma_scale, one_h](double gamma, double pgamma){
      auto gamma_wpmd = gamma_scale * gamma;
      return std::complex<double>(3.0 / (4.0 * gamma_wpmd * gamma_wpmd), -pgamma * gamma_wpmd / (2.0 * gamma_wpmd) * one_h);
    };

    auto da_coeff = [gamma_scale, one_h](double gamma, double pgamma){

    };

    auto b_coeff = [one_h](std::complex<double> const &a, Vector_3 const &r, Vector_3 const &p) {
      return cVector_3{
          2.0 * a * r[0] + std::complex<double>(0, p[0] * one_h),
          2.0 * a * r[1] + std::complex<double>(0, p[1] * one_h),
          2.0 * a * r[2] + std::complex<double>(0, p[2] * one_h),
      };
    };

    double energy = 0.;
    auto const pref_norm = 1.0;

    for (auto i = 0; i < atom->nlocal; ++i){
      auto m = atom->mass[atom->type[i]];
      auto a = a_coeff(atom->eradius[i], atom->ervel[i] * m);
      auto b = b_coeff(a, Vector_3{atom->x[i][0], atom->x[i][1], atom->x[i][2]},
                       Vector_3{atom->v[i][0] * m, atom->v[i][1] * m, atom->v[i][2] * m});
      energy += box->get_integral(a, b, a, b).real();

      std::complex<double> da_re, da_img, da_dummy_re, da_dummy_img;
      cVector_3 db_re, db_img, db_dummy_re, db_dummy_img;
      std::complex<double> integral = 0;
      box->get_derivatives(a, b, a, b, &integral, &da_re, &da_img, &db_re, &db_img,
          &da_dummy_re, &da_dummy_img, &db_dummy_re, &db_dummy_img);
    }

    wall_energy = energy;
    force->pair->eng_coul += wall_energy;
  }
}

double LAMMPS_NS::FixWallAwpmd::compute_scalar() {
  return wall_energy;
}
