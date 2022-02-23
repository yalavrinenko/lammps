//
// Created by yalavrinenko on 09.10.2019.
//

#ifndef DERIVS_ADAPTIVE_MESH_INTEGRATOR_HPP
#define DERIVS_ADAPTIVE_MESH_INTEGRATOR_HPP

#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include "../utils/Logger.hpp"

template<typename real_t>
class value_range {
public:
  class value_range_iterator{
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = real_t const;
    using difference_type = std::ptrdiff_t;
    using pointer = real_t;
    using reference = real_t;

    value_range_iterator(real_t begin, real_t step, size_t init_step) : begin_{std::move(begin)}, dh_{std::move(step)},
                                                                        stepper_{init_step} {
    }

    bool operator == (value_range_iterator const &rhs) const {
      return this->stepper_ == rhs.stepper_;
    }
    bool operator != (value_range_iterator const &rhs) const {
      return !(*this == rhs);
    }

    value_range_iterator& operator++ (){
      ++stepper_; return *this;
    }
    value_range_iterator operator++ (int){
      auto copy = *this; ++stepper_; return copy;
    }

    reference operator* () const { return begin_ + dh_ * stepper_; }
  private:
    real_t const begin_;
    real_t const dh_;
    size_t stepper_;

    real_t const base_precision_ = std::numeric_limits<real_t>::epsilon();
  };

  value_range(real_t begin, real_t end, real_t dh) :
      begin_{std::move(begin)}, end_{std::move(end)}, dh_{std::move(dh)} {
  }

  value_range_iterator begin() const {
    return value_range_iterator(begin_, dh_, 0);
  }

  value_range_iterator end() const {
    return value_range_iterator(begin_, dh_, static_cast<size_t>((end_ - begin_) / dh_));
  }

private:
  real_t begin_, end_, dh_;
};

template<typename real_t>
struct AdaptiveMeshCell {
  struct MeshPoint {
    real_t x, y, z;

    MeshPoint() = default;

    template<typename src_range_t>
    explicit MeshPoint(src_range_t const &r): x(r.x), y(r.y), z(r.z) {}

    MeshPoint(real_t x, real_t y, real_t z) : x(x), y(y), z(z) {}

    MeshPoint operator-(MeshPoint const &p) const { return MeshPoint{x - p.x, y - p.y, z - p.z}; }

    MeshPoint operator+(MeshPoint const &p) const { return MeshPoint{x + p.x, y + p.y, z + p.z}; }

    MeshPoint operator/(real_t v) const { return MeshPoint{x / v, y / v, z / v}; }
  };

  typedef MeshPoint RangeType;
  MeshPoint begin_, end_;

  AdaptiveMeshCell() = default;

  AdaptiveMeshCell(RangeType begin, RangeType end) : begin_{std::move(begin)}, end_{std::move(end)} {}

  inline MeshPoint const &begin() const { return begin_; }

  inline MeshPoint const &end() const { return end_; }

  inline MeshPoint dh() const { return end_ - begin_; }

  inline MeshPoint middle() const { return begin_ + this->dh() / static_cast<real_t>(2.0); }
};

template<typename cell_t, typename integration_unit_t>
class AdaptiveMeshIntegrator {
public:
  struct WorldTopology {
    size_t nodes{1};
    size_t me{0};

    WorldTopology() = default;

    WorldTopology(size_t nodes, size_t me) : nodes(nodes), me(me) {}
  };

  template<typename test_functor_t>
  void
  refine_mesh(typename cell_t::RangeType begin, typename cell_t::RangeType end, test_functor_t const &checker,
              WorldTopology) {
    auto iter = 0;
    clear();

    square_block_decomposition(begin, end);

    while (iter < cells_.size()) {
      auto a = cells_[iter].begin();
      auto b = cells_[iter].end();

      typename cell_t::RangeType dh{b.x - a.x, b.y - a.y, b.z - a.z};

      using RealType = decltype(dh.x);
      if (checker(a, b)) {
        for (auto i = 0; i < 2; ++i)
          for (auto j = 0; j < 2; ++j)
            for (auto k = 0; k < 2; ++k) {
              typename cell_t::RangeType h{static_cast<RealType>(dh.x / 2.0),
                                           static_cast<RealType>(dh.y / 2.0),
                                           static_cast<RealType>(dh.z / 2.0)};
              typename cell_t::RangeType start_subrange{a.x + i * h.x, a.y + j * h.y, a.z + k * h.z};
              typename cell_t::RangeType end_subrange{a.x + (i + 1) * h.x, a.y + (j + 1) * h.y,
                                                      a.z + (k + 1) * h.z};
              if (i == 0 && j == 0 && k == 0)
                cells_[iter] = std::move(cell_t{start_subrange, end_subrange});
              else {
                cells_.emplace_back(std::move(cell_t{start_subrange, end_subrange}));
                if (cells_.size() > 1e+7){
                  Logger::Error("Adaptive mesh overflow.");
                  Logger::Error("Last cell: begin (", start_subrange.x, start_subrange.y, start_subrange.z, ") end (",
                                end_subrange.x, end_subrange.y, end_subrange.z, ") size (",
                                h.x, h.y, h.z, ")");
                  throw std::runtime_error("Adaptive mesh overflow.");
                }
              }
            }
      } else
        ++iter;

    }

    cells_begin_ = cells_.begin();
    cells_end_ = cells_.end();
  }

  template<class real_t, typename mesh_size_t>
  void refine_mesh(typename cell_t::RangeType begin, typename cell_t::RangeType end, mesh_size_t const &step) {
    if (!cells_.empty())
      return;

    for (auto const &x : value_range<real_t>(begin.x, end.x, step.x))
      for (auto const &y : value_range<real_t>(begin.y, end.y, step.y))
        for (auto const &z : value_range<real_t>(begin.z, end.z, step.z)){
          typename cell_t::RangeType start_subrange{x - step.x / 2, y - step.y / 2, z - step.z / 2};
          typename cell_t::RangeType end_subrange{x + step.x / 2, y + step.y / 2, z + step.z / 2};
          cells_.emplace_back(start_subrange, end_subrange);
        }

    cells_begin_ = cells_.begin();
    cells_end_ = cells_.end();
  }

  void refine_linked_mesh(typename cell_t::RangeType begin, typename cell_t::RangeType end, size_t dim_bins,
                          size_t linked_key) {
    using Range = typename cell_t::RangeType;
    auto min_dh = (end.x - begin.x) / dim_bins;

    Range dh{min_dh, min_dh, min_dh};
    for (auto xf = 0; xf < dim_bins; ++xf)
      for (auto yf = 0; yf < dim_bins; ++yf)
        for (auto zf = 0; zf < dim_bins; ++zf) {

          Range fbegin = begin + Range{xf * min_dh, yf * min_dh, zf * min_dh};
          Range fend = fbegin + Range{min_dh, min_dh, min_dh};

          cells_.emplace_back(fbegin, fend);
          cells_.back().linked_index = linked_key;
        }

    cells_begin_ = cells_.begin();
    cells_end_ = cells_.end();
  }

  template<typename output_t, typename ... integration_args_t>
  output_t integrate(integration_unit_t &engine, integration_args_t &&... args) {
    return engine.template integrate<output_t>(cells_begin_, cells_end_, std::forward<integration_args_t>(args)...);
  }
# if 0
  template<typename ... foreach_args_t>
  void for_each(integration_unit_t &engine, foreach_args_t &&... args) {
    engine.template for_each(cells_begin_, cells_end_, std::forward<foreach_args_t>(args)...);
  }
# endif

  std::vector<cell_t> const &mesh_cells() const { return cells_; }

  void clear() {
    cells_.clear();
  }

  using cell_iterator = typename std::vector<cell_t>::const_iterator;

protected:
  void square_block_decomposition(typename cell_t::RangeType begin, typename cell_t::RangeType end) {
    auto &a = begin;
    auto &b = end;

    typename cell_t::RangeType dh{b.x - a.x, b.y - a.y, b.z - a.z};
    auto min_dh = std::min(dh.x, std::min(dh.y, dh.z));

    for (auto x = a.x; x < b.x - min_dh * 0.1; x += min_dh)
      for (auto y = a.y; y < b.y - min_dh * 0.1; y += min_dh)
        for (auto z = a.z; z < b.z - min_dh * 0.1; z += min_dh) {
          typename cell_t::RangeType begin_{x, y, z};
          typename cell_t::RangeType end_{x + min_dh, y + min_dh, z + min_dh};
          cells_.emplace_back(std::move(cell_t{begin_, end_}));
        }
  }

  std::vector<cell_t> cells_{};
  typename std::vector<cell_t>::const_iterator cells_begin_{};
  typename std::vector<cell_t>::const_iterator cells_end_{};
};

#endif //DERIVS_ADAPTIVE_MESH_INTEGRATOR_HPP
