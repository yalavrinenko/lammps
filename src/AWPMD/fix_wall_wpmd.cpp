//
// Created by yalavrinenko on 11.02.19.
//

#include "fix_wall_wpmd.h"
#include "WavepacketPairCommon.h"
#include "domain.h"
#include "error.h"
#include "neigh_list.h"
#include <atom.h>
#include <cstring>
#include <force.h>
#include <wpmd_split.h>

LAMMPS_NS::FixWallWpmd::FixWallWpmd(LAMMPS_NS::LAMMPS *lammps, int i,
                                    char **pString)
    : Fix(lammps, i, pString) {
  double delx = domain->boxhi[0] - domain->boxlo[0];
  double dely = domain->boxhi[1] - domain->boxlo[1];
  double delz = domain->boxhi[2] - domain->boxlo[2];
  auto half_box_length = 0.5 * MIN(delx, MIN(dely, delz));

  wall_squares = {delz * dely, delx * delz, delx * dely};

  m_pair = dynamic_cast<WavepacketPairCommon *>(force->pair);
  this->box = construct_box(pString, half_box_length, i);

  this->vector_flag = true;
  this->size_vector = 2;
}

LAMMPS_NS::FixWallWpmd::~FixWallWpmd() {
  m_pair->awpmd()->use_box = false;
  m_pair->awpmd()->set_pbc(nullptr, 0);
}

int LAMMPS_NS::FixWallWpmd::setmask() {
  return LAMMPS_NS::FixConst::POST_FORCE;
}

std::unique_ptr<BoxHamiltonian>
LAMMPS_NS::FixWallWpmd::construct_box(char **pString, double half_box_length,
                                      int pcount) {
  auto numeric = [this](char const* str){
    return utils::numeric(FLERR, str, false, lmp);
  };

  auto eigenE = numeric(pString[3]);

  for (auto i = 3; i < pcount; ++i) {
    if (std::strcmp(pString[i], "box") == 0) {
      auto Lx = numeric(pString[i + 1]);

      half_box_length = 0.5 * Lx;
      wall_squares = {Lx * Lx, Lx * Lx, Lx * Lx};
      i += 1;
    }
    if (std::strcmp(pString[i], "width_force") == 0)
      use_width_force_ = true;
    if (std::strcmp(pString[i], "axes") == 0) {

      has_force_ = {false, false, false};

      auto is_keyword = [](char const* str){
        return std::strcmp(str, "x") == 0
            || std::strcmp(str, "y") == 0
            || std::strcmp(str, "z") == 0;
      };
      auto j = 1;
      while (is_keyword(pString[i + j])){
        has_force_[pString[i + j][0] - 'x'] = true;
        ++j;
      }
      i += j;
    }
  }

  walls_count_ = std::count(has_force_.begin(), has_force_.end(), 1);

  auto floor = half_box_length;

  auto me = force->e_mass;
  auto h2_me = force->hhmrr2e / force->e_mass;
  auto one_h = force->mvh2r;

  double eigenwp = 0.0;

  if (eigenE > 0.) {
    eigenwp = sqrt(3. / 2 / me / eigenE) / one_h;
  }
// else   eigenE = 3. / 2 * h2_me / (eigenwp * eigenwp);

  double floorYtoX = 1., floorZtoX = 1., widthYtoX = 1., widthZtoX = 1.;

  Vector_3 gamma(eigenwp, eigenwp * widthYtoX, eigenwp * widthZtoX), force_k;

  for (int i = 0; i < 3; ++i) {
    force_k[i] = 9. / 8 * h2_me / (gamma[i] * gamma[i] * gamma[i] * gamma[i]) * has_force_[i];
  }

  Vector_3 bound(floor, floor * floorYtoX, floor * floorZtoX);
  auto const PROJ_ORDER_CONST = 10;

  return std::unique_ptr<BoxHamiltonian>(
      new BoxHamiltonian(bound, force_k, PROJ_ORDER_CONST));
}

void LAMMPS_NS::FixWallWpmd::post_force(int flag) {
  wall_energy = 0;
  if (m_pair && !m_pair->electrons_packets().empty()) {
    evaluate_wall_energy(m_pair->electrons_packets());
  } else {
    packets.resize(atom->nlocal + atom->nghost);

    auto one_h = force->mvh2r;
    for (auto i = 0; i < atom->nlocal + atom->nghost; ++i) {
      if (atom->spin[i] != 0) {
        double width = atom->eradius[i];
        Vector_3 r{atom->x[i][0], atom->x[i][1], atom->x[i][2]},
            p{atom->v[i][0], atom->v[i][1], atom->v[i][2]};
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

double LAMMPS_NS::FixWallWpmd::compute_scalar() { return wall_energy; }

double LAMMPS_NS::FixWallWpmd::compute_vector(int i) {
  switch (i) {
    case 0:
      return wall_energy;
    case 1:
      return wall_pressure();
    default:
      throw std::logic_error("Out of range");
  }
}

double LAMMPS_NS::FixWallWpmd::interaction_border_ion(int, double *x,
                                                      double *f) {
  double dE;
  Vector_3 df = box->get_force(*(Vector_3 *) x, &dE);
  if (f) // ion forces needed
    for (auto k = 0; k < 3; ++k)
      f[k] += df[k];
  return dE;
}

double LAMMPS_NS::FixWallWpmd::interaction_border_electron(
    WavePacket const &packet, double *rforce, double *erforce,
    double *ervforce) {
  double dE;
  if (force && erforce && ervforce) {
    cdouble integral;
    cdouble a1_re, a1_im, a2_re, a2_im;
    cVector_3 b1_re, b1_im, b2_re, b2_im;
    box->get_derivatives(packet.a, packet.b, packet.a, packet.b, &integral,
                         &a1_re, &a1_im, &b1_re, &b1_im, &a2_re, &a2_im, &b2_re,
                         &b2_im);

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
      rforce[k] += -dx[k] * has_force_[k];
    (*erforce) += *dw;
    (*ervforce) += *pw;
    dE = integral.real();
  } else
    dE = box->get_integral(packet.a, packet.b, packet.a, packet.b).real();
  // Ebord += dE;
  return dE;
}

void LAMMPS_NS::FixWallWpmd::evaluate_wall_energy(
    std::vector<WavePacket> const &wavepackets) {
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
      wall_energy += interaction_border_electron(wavepackets[i], f, &erf, &ervf);
      atom->erforce[i] += erf;
      atom->ervelforce[i] += ervf;
      wall_pressure_components[3] += std::abs(erf);
    }

    for (auto k = 0; k < 3; ++k) {
      atom->f[i][k] += f[k];
      wall_pressure_components[k] += std::abs(f[k]);
    }
  }

  std::array<double, 4> force_components{0, 0, 0, 0};

  MPI_Allreduce(wall_pressure_components.data(), force_components.data(), 4,
                MPI_DOUBLE, MPI_SUM, world);

  wall_pressure_ = (force_components[0] / (2.0 * wall_squares[0]) +
                    force_components[1] / (2.0 * wall_squares[1]) +
                    force_components[2] / (2.0 * wall_squares[2])) / walls_count_;

  if (use_width_force_)
    wall_pressure_ += force_components[3] /
                      (2.0 * (wall_squares[0] + wall_squares[1] + wall_squares[2]));

  wall_pressure_ = wall_pressure_ * force->nktv2p;
}

double LAMMPS_NS::FixWallWpmd::wall_pressure() const { return wall_pressure_; }

void LAMMPS_NS::FixWallWpmd::setup(int i) { post_force(i); }
