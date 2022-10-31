//
// Created by yalavrinenko on 29.11.2019.
//

#include "compute_norm_matrix.h"
#include "wpmd_split.h"
#include <cstring>
#include <domain.h>
#include <error.h>
#include <atom.h>
#include <force.h>
#include <update.h>
#include <DataTypes.hpp>


double ComputeNormAwpmd::compute_scalar() {
  
  if (invoked_scalar != update->ntimestep) {

    scalar = 0;
    double result = wppair->awpmd()->norm_matrix_detl(0) + wppair->awpmd()->norm_matrix_detl(1) ;
    MPI_Allreduce(&result, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

    invoked_scalar = update->ntimestep;
  }

  return scalar;
}

ComputeNormAwpmd::ComputeNormAwpmd(LAMMPS_NS::LAMMPS *lmp, int argc, char **argv) : Compute(lmp, argc, argv) {
  //array_flag = true;
  scalar_flag = true;
  wppair = dynamic_cast<LAMMPS_NS::PairAWPMD *>(force->pair);
  if(!wppair)
    error->all(FLERR, "Norm matrix compute requires pair_awpmd");


  //
# if 0
  auto argv_index = 3;
  auto nbins = std::stol(argv[argv_index]);

  auto is_par_equal = [&argv](int index, char const* value){
    return std::strcmp(argv[index], value) == 0;
  };

  std::array<double, 3> begin{};

  while (argv_index < argc){
    if (is_par_equal(argv_index, "axis")){
      ++argv_index;
      for (auto j = 0; j < 3; ++j, ++argv_index){
        if (argv[argv_index][0] == '-'){
          nbins_[j] = nbins;
          vary_axis_[j] = true;
        } else {
          axis_[j] = std::stod(argv[argv_index]);
        }
      }
    } else if (is_par_equal(argv_index, "region")){
      ++argv_index;
      if (is_par_equal(argv_index, "block")){
        ++argv_index;
        region_ = this->domain->get_region_by_id(argv[argv_index]);
        if (region_) {
          begin = {region_->extent_xlo, region_->extent_ylo, region_->extent_zlo};
          L_ = {region_->extent_xhi - region_->extent_xlo,
               region_->extent_yhi - region_->extent_ylo,
               region_->extent_zhi - region_->extent_zlo};
        } else {
          throw std::runtime_error("Box not found");
        }
      } else if (is_par_equal(argv_index, "cbox")){
        ++argv_index;
        for (auto j = 0; j < 3; ++j, ++argv_index){
          L_[j] = std::stod(argv[argv_index]);
          begin[j] = -L_[j] / 2;
        }
        --argv_index;
      }
    } else if (is_par_equal(argv_index, "scalef")){
      ++argv_index;
      scalef_ = std::stod(argv[argv_index]);
    } else if (is_par_equal(argv_index, "center")) {
      use_center_ = true;
      ++argv_index;
    } else {
      ++argv_index;
    }
  }

  for (auto &l : L_) {
    l *= scalef_;
  }
  for (auto &b : begin) { b *= scalef_; }

  std::array<double, 3> delta{L_[0] / nbins_[0], L_[1] / nbins_[1], L_[2] / nbins_[2]};

  for (auto i = 0; i < 3; ++i)
    if (!vary_axis_[i]) {
      begin[i] = axis_[i] - delta[i] / 2;
    }

  std::cout << "Electron density profile info ...\n";
  for (auto i = 0; i < 3; ++i){
    std::cout << "\tCompute density along axe " << i << ". Start point: " << begin[i] << " step delta: "
      << delta[i] << " bins for axe " << nbins_[i] << std::endl;
  }

  plane_info info{begin[0], begin[1], begin[2],
                  delta[0],delta[1], delta[2],
                  nbins_[0],
                  nbins_[1],
                  nbins_[2]};

  create_cell_list(info);

  size_array_cols = 4;
  size_array_rows = cells_.size();

  memory->create(array,size_array_rows,size_array_cols,"denprof:array");

  config_.use_adaptive_mesh = true;
  config_.max_distance = 4.0;
  config_.min_cell = 0.7;

  for (auto i : {0, 1, 2}){
    config_.mesh_start[i] = -L_[i] / 2.0;
    config_.mesh_fin[i] = L_[i] / 2.0;
  }

  config_.approximation = new VoidApproximation();
# endif 
}





ComputeNormAwpmd::~ComputeNormAwpmd() {
  //memory->destroy(array);
}
