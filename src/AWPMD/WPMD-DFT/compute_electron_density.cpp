//
// Created by yalavrinenko on 29.11.2019.
//

#include "compute_electron_density.h"
#include <cstring>
#include <region.h>
#include <domain.h>
#include <error.h>
#include <atom.h>
#include <force.h>
#include <update.h>
#include <DataTypes.hpp>

struct mesh_stepper {
  double init;
  double value;
  double dvalue;

  mesh_stepper(double init, double dvalue) : init{init - dvalue}, value{init}, dvalue{dvalue} {
  }

  double next() {
    value += dvalue;
    return value;
  }

  void reset() { value = init; }
};

double ComputeDensityAwpmd::compute_scalar() {
  if (invoked_scalar != update->ntimestep) {
    make_packets();
    std::vector<XCEnergy_cpu::cell_density> inner_cell(1);

    double3 begin{-L_[0] / 2.0, -L_[1] / 2.0, -L_[2] / 2.0};
    double3 end{L_[0] / 2.0, L_[1] / 2.0, L_[2] / 2.0};
    AdaptiveMeshCell<double> range{AdaptiveMeshCell<double>::MeshPoint{begin},
                                   AdaptiveMeshCell<double>::MeshPoint{end}};
    inner_cell[0].cell = range;

    double result = xc_energy_->build_density_map(packets_, inner_cell, use_center_);
    result = result / (force->angstrom * force->angstrom * force->angstrom) * 1e+24;

    scalar = 0;
    MPI_Allreduce(&result, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

    invoked_scalar = update->ntimestep;
  }

  return scalar;
}

ComputeDensityAwpmd::ComputeDensityAwpmd(LAMMPS_NS::LAMMPS *lmp, int argc, char **argv) : Compute(lmp, argc, argv) {
  array_flag = true;
  scalar_flag = true;

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
        auto domain_index = this->domain->find_region(argv[argv_index]);
        if (domain_index != -1) {
          region_ = this->domain->regions[domain_index];
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
}

void ComputeDensityAwpmd::create_cell_list(plane_info const &info) {
  mesh_stepper s[3] = {mesh_stepper(info.x, info.dx),
                       mesh_stepper(info.y, info.dy),
                       mesh_stepper(info.z, info.dz)};

  s[0].reset();
  for (auto i = 0; i < info.Nx; ++i) {
    double x = s[0].next();
    s[1].reset();
    for (auto j = 0; j < info.Ny; ++j) {
      double y = s[1].next();
      s[2].reset();
      for (auto k = 0; k < info.Nz; ++k) {
        double z = s[2].next();
        double3 begin{x, y, z};
        double3 end{begin.x + info.dx, begin.y + info.dy, begin.z + info.dz};

        XCEnergy_cpu::cell_density cell{};
        AdaptiveMeshCell<double> range{AdaptiveMeshCell<double>::MeshPoint{begin}, AdaptiveMeshCell<double>::MeshPoint{end}};
        cell.cell = range;
        cells_.emplace_back(cell);
      }
    }
  }
}

void ComputeDensityAwpmd::compute_array() {
  if (invoked_array != update->ntimestep) {
    make_packets();
    xc_energy_->build_density_map(packets_, cells_, use_center_);

    auto row = 0;
    std::vector<double> tmp_density(cells_.size());
    for (auto const &cell : cells_) {
      double cell_[3] = {(cell.cell.begin().x + cell.cell.end().x) * 0.5,
                         (cell.cell.begin().y + cell.cell.end().y) * 0.5,
                         (cell.cell.begin().z + cell.cell.end().z) * 0.5};
      for (auto k = 0; k < 3; ++k)
        array[row][k] = cell_[k];
      tmp_density[row] = cell.rho / (force->angstrom * force->angstrom * force->angstrom) * 1e+24;
      ++row;
    }

    std::vector<double> recv_buffer(cells_.size(), 0);
    MPI_Allreduce(tmp_density.data(), recv_buffer.data(), tmp_density.size(), MPI_DOUBLE, MPI_SUM, world);
    for (row = 0; row < recv_buffer.size(); ++row)
      array[row][3] = recv_buffer[row];

    invoked_array = update->ntimestep;
  }

  //array[row][3] = cell.rho / (force->angstrom * force->angstrom * force->angstrom) * 1e+24;
}

void ComputeDensityAwpmd::init() {
  xc_energy_  = std::unique_ptr<XCEnergy_cpu>(new XCEnergy_cpu(1, config_)); //old compiler support
}

void ComputeDensityAwpmd::make_packets() {
  packets_.clear();

  auto one_h=force->mvh2r;
  for (auto i = 0; i < atom->nlocal; ++i){
    if (atom->mask[i] & groupbit){
      packets_.emplace_back(atom->x[i], atom->eradius[i], ElectronSpin(atom->spin[i]), false);
    }
  }
}

ComputeDensityAwpmd::~ComputeDensityAwpmd() {
  memory->destroy(array);
}
