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
#include "neigh_list.h"

LAMMPS_NS::FixWallAwpmd::FixWallAwpmd(LAMMPS_NS::LAMMPS *lammps, int i, char **pString) : Fix(lammps, i, pString) {
  double delx = domain->boxhi[0]-domain->boxlo[0];
  double dely = domain->boxhi[1]-domain->boxlo[1];
  double delz = domain->boxhi[2]-domain->boxlo[2];
  auto half_box_length = 0.5 * MIN(delx, MIN(dely, delz));

  wall_squares = {delz * dely, delx * delz, delx * dely};

  Vector_3 box_size{delx, dely, delz};

  m_pair = dynamic_cast<PairAWPMDCut*>(force->pair);
  this->box = construct_box(pString, half_box_length);

  this->vector_flag = true;
  this->size_vector = 2;
  this->virial_flag = 1;
}

LAMMPS_NS::FixWallAwpmd::~FixWallAwpmd() {
  m_pair->awpmd()->use_box = false;
  m_pair->awpmd()->set_pbc(nullptr, 0);
}

int LAMMPS_NS::FixWallAwpmd::setmask() {
  return LAMMPS_NS::FixConst::POST_FORCE;
}

std::unique_ptr<BoxHamiltonian> LAMMPS_NS::FixWallAwpmd::construct_box(char **pString, double half_box_length) {
  auto box_fraction = force->numeric(FLERR, pString[3]);
  auto eigenE = force->numeric(FLERR, pString[4]);
  double prj_ord = force->numeric(FLERR, pString[5]);

  auto floor = half_box_length;
  auto eigenwp = half_box_length / (box_fraction < 1.0 ? 10.0 : box_fraction);

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

void LAMMPS_NS::FixWallAwpmd::post_force(int i)  {
  wall_energy = 0;
  if (m_pair){
    evaluate_wall_energy(m_pair->electrons_packets());
  } else {
    std::vector<WavePacket> packets(atom->nlocal + atom->nghost);

    auto one_h=force->mvh2r;
    for (auto i = 0; i < atom->nlocal + atom->nghost; ++i){
      if (atom->spin[i] != 0) {
        double width = atom->eradius[i];
        Vector_3 r{atom->x[i][0], atom->x[i][1], atom->x[i][2]}, p{atom->v[i][0], atom->v[i][1], atom->v[i][2]};
        p *= one_h * atom->mass[atom->type[i]];

        double pw = atom->ervel[i];
        pw *= one_h * atom->mass[atom->type[i]];

        packets[i].init(width, r, p, pw);
      }
    }
    evaluate_wall_energy(packets);
  }
  m_pair->eng_coul += wall_energy;
}

double LAMMPS_NS::FixWallAwpmd::compute_scalar() {
  return wall_energy;
}

double LAMMPS_NS::FixWallAwpmd::compute_vector(int i) {
  switch (i){
    case 0: return wall_energy;
    case 1: return wall_pressure();
    default:
      throw std::logic_error("Out of range");
  }
}

double LAMMPS_NS::FixWallAwpmd::interaction_border_ion(int i, double *x, double *f) {
  double dE;
  Vector_3 df = box->get_force(*(Vector_3*)x, &dE);
  if (f) // ion forces needed
    for (auto k = 0; k < 3; ++k)
      f[k] += df[k];
  return dE;
}

double LAMMPS_NS::FixWallAwpmd::interaction_border_electron(WavePacket const &packet, double *rforce, double *erforce,
                                                            double *ervforce) {
  double dE{0};
  if (force && erforce && ervforce) {
    cdouble integral;
    cdouble a1_re, a1_im, a2_re, a2_im;
    cVector_3 b1_re, b1_im, b2_re, b2_im;
    box->get_derivatives(packet.a, packet.b, packet.a, packet.b, &integral, &a1_re, &a1_im, &b1_re, &b1_im,
                        &a2_re, &a2_im, &b2_re, &b2_im);

    std::array<double, 8> tmp{2.0 * real(a1_re), 2.0 * real(a1_im),
                              2.0 * real(b1_re[0]), 2.0 * real(b1_im[0]),
                              2.0 * real(b1_re[1]), 2.0 * real(b1_im[1]),
                              2.0 * real(b1_re[2]), 2.0 * real(b1_im[2])};
    auto dx = tmp.begin();
    auto dp = dx + 3;
    auto dw = dp + 3;
    auto pw = dw + 1;

    packet.int2phys_der<eq_second>(dx, dx, dp, dw, pw, 1. / force->mvh2r);
    for (auto k = 0u; k < 3; ++k)
      rforce[k] += -dx[k];
    (*erforce) += *dw;
    (*ervforce) += *pw;
    dE = integral.real();
  } else
    dE = box->get_integral(packet.a, packet.b, packet.a, packet.b).real();
  //Ebord += dE;
  return dE;
}

void LAMMPS_NS::FixWallAwpmd::evaluate_wall_energy(std::vector<WavePacket> const &packets) {
  auto inum = m_pair->list->inum;
  auto ilist = m_pair->list->ilist;

  wall_pressure_components = {0, 0, 0, 0};

  for (auto ii = 0; ii < inum; ii++) {
    auto i = ilist[ii];
    double f[3] = {0, 0, 0};
    double erf = 0, ervf = 0;
    if (atom->spin[i] == 0) {
      wall_energy += interaction_border_ion(i, atom->x[i], f);
    } else {
      wall_energy += interaction_border_electron(packets[i], f, &erf,
                                                 &ervf);
      atom->erforce[i] += erf;
      atom->ervelforce[i] += ervf;
      wall_pressure_components[3] += std::abs(erf) + std::abs(ervf);
    }

    for (auto k = 0; k < 3; ++k) {
      atom->f[i][k] += f[k];
      wall_pressure_components[k] += std::abs(f[k]);
    }

    double pressure = wall_pressure_components[0] + wall_pressure_components[1] + wall_pressure_components[2];

    double reduce_pressure = 0;
    MPI_Allreduce(&pressure,&reduce_pressure,1,MPI_DOUBLE,MPI_SUM,world);

    double square = 2.0 * (wall_squares[0] + wall_squares[1] + wall_squares[2]);
    wall_pressure_ = reduce_pressure / square * force->nktv2p;
//    virial[0] += f[0]*atom->x[i][0];
//    virial[1] += f[1]*atom->x[i][1];
//    virial[2] += f[2]*atom->x[i][2];
//    virial[3] += f[1]*atom->x[i][0];
//    virial[4] += f[2]*atom->x[i][0];
//    virial[5] += f[2]*atom->x[i][1];
  }
}

double LAMMPS_NS::FixWallAwpmd::wall_pressure() const {
  return wall_pressure_;
}

void LAMMPS_NS::FixWallAwpmd::setup(int i) {
  post_force(i);
}
