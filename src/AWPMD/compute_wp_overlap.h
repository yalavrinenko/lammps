//
// Created by Igor Morozov on 12/04/2025
// Based on compute_rdf.h (Copyright (2003) Sandia Corporation)
//

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(wp_overlap,ComputeWPOverlap);
// clang-format on
#else

#ifndef LMP_COMPUTE_WP_OVERLAP_H
#define LMP_COMPUTE_WP_OVERLAP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeWPOverlap : public Compute {
 public:
  ComputeWPOverlap(class LAMMPS *, int, char **);
  ~ComputeWPOverlap() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;

 private:
  int nbin;                // # of bins
  int cutflag;             // user cutoff flag
  int npairs;              // # of pairs
  double delr, delrinv;    // bin width and its inverse
  double cutoff_user;      // user-specified cutoff
  double mycutneigh;       // user-specified cutoff + neighbor skin
  int ***wppair;          // map 2 type pair to WP pair for each histo
  int **nwppair;          // # of histograms for each type pair
  int *ilo, *ihi, *jlo, *jhi;
  double **hist;       // histogram bins
  double **histall;    // summed histogram bins across all procs

  double wp_bound_width;      // IM: width of a bound electron

  int *typecount;
  int *icount, *jcount;
  int *duplicates;

  class NeighList *list;    // half neighbor list
  void init_norm();
  bigint natoms_old;
};

}    // namespace LAMMPS_NS

#endif
#endif
