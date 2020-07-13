//
// Created by yalavrinenko on 16.03.2020.
//

#include "pair_wpmd_cut.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "min.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"
#include <wpmd_split.h>

LAMMPS_NS::WavepacketPairCommon::awpmd_energies LAMMPS_NS::PairWPMD::compute_energy_force() {
  auto inum = list->inum;
  auto ilist = list->ilist;
  auto numneigh = list->numneigh;
  auto firstneigh = list->firstneigh;

  wpmd->Eii = 0;
  wpmd->Ebord = wpmd->Eext = wpmd->Eee = wpmd->Ew = 0.;
  wpmd->Ee[0] = wpmd->Ee[1] = 0.;
  wpmd->Eei[0] = wpmd->Eei[1] = 0.;
  wpmd->Ebord_ion = 0;

//    if (packets.size() < atom->nlocal + atom->nghost)
//      packets.resize(atom->nlocal + atom->nghost);

  auto one_h = force->mvh2r;

//    for (auto i = 0; i < atom->nlocal + atom->nghost; ++i){
//      if (atom->spin[i] != 0) {
//        double width = atom->eradius[i];
//        Vector_3 r{atom->x[i][0], atom->x[i][1], atom->x[i][2]},
//        p{atom->v[i][0], atom->v[i][1], atom->v[i][2]}; p *= one_h * atom->mass[atom->type[i]];
//
//        double pw = atom->ervel[i];
//        pw *= one_h * atom->mass[atom->type[i]];
//
//        packets[i].init(width, r, p, pw);
//      }
//    }
  auto pwidth = [one_h, this](size_t index) { return one_h * atom->mass[atom->type[index]] * atom->ervel[index]; };
  awpmd_energies interaction_energy{};

  for (auto ii = 0; ii < inum; ii++) {
    auto i = ilist[ii];

    if (atom->spin[i] != 0) {
      double ke, ee_w;
      std::tie(ke, ee_w) = wpmd->AWPMD::interaction_electron_kinetic(atom->eradius[i], pwidth(i), atom->spin[i] + 1,
                                                                     &atom->erforce[i], &atom->ervelforce[i]);
      interaction_energy.ke += ke;
      interaction_energy.ee_w += ee_w;
    }

    for (auto jj = 0; jj < numneigh[i]; jj++) {
      auto j = firstneigh[i][jj];
      j &= NEIGHMASK;

      double cutoff = (is_pbc_) ? std::sqrt(cutsq[atom->type[i]][atom->type[j]]) : -1;

      if (atom->spin[i] == 0 && atom->spin[j] == 0 && wpmd->calc_ii) {
        interaction_energy.ii += wpmd->interation_ii_single(i, j, atom->x, atom->q, atom->f, cutoff);
      } else if (atom->spin[i] != 0 && atom->spin[j] != 0 && wpmd->calc_ee) {
        double *eforce_ptrs[]{atom->f[i], atom->f[j]};
        double *erforce_ptrs[]{&atom->erforce[i], &atom->erforce[j]};
        interaction_energy.ee += wpmd->AWPMD::interaction_ee_single(
            atom->x[i], atom->eradius[i], atom->x[j], atom->eradius[j], eforce_ptrs, erforce_ptrs, cutoff);
      } else if (atom->spin[i] == 0 && atom->spin[j] != 0 && wpmd->calc_ei) {
        interaction_energy.ei +=
            wpmd->AWPMD::interaction_ei_single(atom->x[i], atom->q[i], atom->x[j], atom->eradius[j], atom->spin[j],
                                               atom->f[i], atom->f[j], &atom->erforce[j], cutoff);
      } else if (atom->spin[i] != 0 && atom->spin[j] == 0 && wpmd->calc_ei) {
        interaction_energy.ei +=
            wpmd->AWPMD::interaction_ei_single(atom->x[j], atom->q[j], atom->x[i], atom->eradius[i], atom->spin[i],
                                               atom->f[j], atom->f[i], &atom->erforce[i], cutoff);
      }
    }
  }

  return interaction_energy;
}
void LAMMPS_NS::PairWPMD::settings(int argc, char **pString) {
  WavepacketPairCommon::settings(argc, pString);

  is_pbc_ = true;
  for (auto i = 0; i < 3; ++i)
    if (domain->boundary[i][0] != 0 || domain->boundary[i][1] != 0) {
      is_pbc_ = false;
      break;
    }
}
