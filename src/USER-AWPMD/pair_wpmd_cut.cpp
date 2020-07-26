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

class neighlist_wrapper {
public:
  class interaction_pairs_iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = std::pair<std::size_t, std::size_t>;
    using pointer = value_type const *;
    using reference = value_type const &;

    explicit interaction_pairs_iterator(LAMMPS_NS::NeighList const *list): list_{list} {
    }

    bool operator == (neighlist_wrapper::interaction_pairs_iterator const &it) const{
      return list_ == it.list_ && i_index_ == it.i_index_ && j_index_ == it.j_index_;
    }
    bool operator != (neighlist_wrapper::interaction_pairs_iterator const &it) const{
      return !(*this == it);
    }

    neighlist_wrapper::interaction_pairs_iterator & operator++ () {
      next(); return *this;
    }

    neighlist_wrapper::interaction_pairs_iterator operator++ (int) {
      auto tmp = *this; ++(*this); return tmp;
    }

    value_type operator* () const{
      return atom_pair();
    }

    value_type operator-> () const {
      return atom_pair();
    }
  private:

    interaction_pairs_iterator(size_t i_index, size_t j_index, LAMMPS_NS::NeighList const *list):
      i_index_{i_index}, j_index_{j_index}, list_{list}{
    }

    std::size_t i_particle_index() const {
      return list_->ilist[i_index_];
    }

    value_type atom_pair() const {
      return {i_particle_index(), list_->firstneigh[i_particle_index()][j_index_] & NEIGHMASK};
    }

    void next(){
      if (list_ != nullptr) {
        if (i_index_ >= list_->inum) {
          set_end();
          return;
        }

        if (j_index_ >= list_->numneigh[i_particle_index()]) {
          ++i_index_; j_index_ = 0;
        }

        if (i_index_ < list_->inum) {
          ++j_index_;
        } else {
          set_end();
        }
      }
    }

    void set_end() {
      list_ = nullptr;
    }

    size_t i_index_ {0};
    size_t j_index_ {0};
    LAMMPS_NS::NeighList const *list_ = nullptr;
  };

  explicit neighlist_wrapper(LAMMPS_NS::NeighList const* list): list_{list}{
  }

  interaction_pairs_iterator begin() const {
    return interaction_pairs_iterator(list_);
  }

  interaction_pairs_iterator end() const {
    return interaction_pairs_iterator(nullptr);
  }

private:
  LAMMPS_NS::NeighList const* list_;
};

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

  auto one_h = force->mvh2r;

  auto pwidth = [one_h, this](size_t index) { return one_h * atom->mass[atom->type[index]] * atom->ervel[index]; };
  awpmd_energies interaction_energy{};

  neighlist_wrapper nlist_wrapper{list};

  for (auto const &atom_pair : nlist_wrapper){
    std::cout << atom_pair.first << "\t" << atom_pair.second << std::endl;
  }

  std::cout << "END PAIR" << std::endl;

//  for (auto ii = 0; ii < inum; ii++) {
//    auto i = ilist[ii];
//
//    if (atom->spin[i] != 0) {
//      double ke, ee_w;
//      std::tie(ke, ee_w) = wpmd->AWPMD::interaction_electron_kinetic(atom->eradius[i], pwidth(i), atom->spin[i] + 1,
//                                                                     &atom->erforce[i], &atom->ervelforce[i]);
//      interaction_energy.ke += ke;
//      interaction_energy.ee_w += ee_w;
//    }
//
//    for (auto jj = 0; jj < numneigh[i]; jj++) {
//      auto j = firstneigh[i][jj];
//      j &= NEIGHMASK;
//
//      double cutoff = (is_pbc_) ? std::sqrt(cutsq[atom->type[i]][atom->type[j]]) : -1;
//
//      if (atom->spin[i] == 0 && atom->spin[j] == 0 && wpmd->calc_ii) {
//        interaction_energy.ii += wpmd->interation_ii_single(i, j, atom->x, atom->q, atom->f, cutoff);
//      } else if (atom->spin[i] != 0 && atom->spin[j] != 0 && wpmd->calc_ee) {
//        double *eforce_ptrs[]{atom->f[i], atom->f[j]};
//        double *erforce_ptrs[]{&atom->erforce[i], &atom->erforce[j]};
//        interaction_energy.ee += wpmd->AWPMD::interaction_ee_single(
//            atom->x[i], atom->eradius[i], atom->x[j], atom->eradius[j], eforce_ptrs, erforce_ptrs, cutoff);
//      } else if (atom->spin[i] == 0 && atom->spin[j] != 0 && wpmd->calc_ei) {
//        interaction_energy.ei +=
//            wpmd->AWPMD::interaction_ei_single(atom->x[i], atom->q[i], atom->x[j], atom->eradius[j], atom->spin[j],
//                                               atom->f[i], atom->f[j], &atom->erforce[j], cutoff);
//      } else if (atom->spin[i] != 0 && atom->spin[j] == 0 && wpmd->calc_ei) {
//        interaction_energy.ei +=
//            wpmd->AWPMD::interaction_ei_single(atom->x[j], atom->q[j], atom->x[i], atom->eradius[i], atom->spin[i],
//                                               atom->f[j], atom->f[i], &atom->erforce[i], cutoff);
//      }
//    }
//  }

  return interaction_energy;
}
void LAMMPS_NS::PairWPMD::settings(int argc, char **pString) {
  WavepacketPairCommon::settings(argc, pString);

  is_pbc_ = true;
  for (auto & i : domain->boundary)
    if (i[0] != 0 || i[1] != 0){
      is_pbc_ = false;
      break;
    }
}
