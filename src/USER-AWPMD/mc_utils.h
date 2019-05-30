//
// Created by yalavrinenko on 21.02.19.
//

#ifndef LAMMPS_MC_UTILS_H
#define LAMMPS_MC_UTILS_H

#include "random_park.h"
#include <mcarlo.h>
#include <memory>
#include <functional>
#include <vector>
#include <array>
#include "atom.h"
#include <iterator>
#include <unordered_map>

namespace LAMMPS_NS{

  enum class stepper_type{
    ion_r = 0,
    electron_r = 1,
    electron_p = 2,
    electron_w = 3,
    electron_pw = 4,
    electron_c0 = 5, //Not impl
    electron_c1 = 6, //Not impl
    ion_p = 7
  };

  struct mc_stepper;

  class MCSystem{
  public:
    using filter_func = std::function<bool(int)>;

    explicit MCSystem(filter_func &&filter): _filter(std::move(filter)){
    }

    virtual void save(size_t size) = 0;
    virtual void restore(size_t size) = 0;
    virtual void make(size_t size, mc_stepper& stepper) = 0;
    virtual std::vector<double> pack(size_t size, int *tag) const = 0;
    virtual size_t object_size() const = 0;
    virtual size_t unpack(const double *data, size_t size, std::unordered_map<int, int> const &ghost_map) = 0;
  protected:
    filter_func _filter;
  };

  class MC3DVectorSystem: public MCSystem{
  public:
    void save(size_t size) override {
      storage.clear();
      storage.reserve(size);
      for (auto i = 0; i < size; ++i)
        if (_filter(i)){
          storage.push_back({src[i][0], src[i][1], src[i][2]});
        }
      storage.shrink_to_fit();
    }

    void restore(size_t size) override {
      auto insert_iterator = storage.begin();
      for (auto i = 0; i < size; ++i)
        if (_filter(i)){
          src[i][0] = (*insert_iterator)[0];
          src[i][1] = (*insert_iterator)[1];
          src[i][2] = (*insert_iterator)[2];
          ++insert_iterator;
        }
    }

    void make(size_t size, mc_stepper &stepper) override;

    std::vector<double> pack(size_t size, int *tag) const override;

    size_t object_size() const override {
      return 4;
    }

    size_t unpack(const double *data, size_t size, std::unordered_map<int, int> const &ghost_map) override;

  public:
    MC3DVectorSystem(double** &source, filter_func &&filter):
        src(source), MCSystem(std::forward<filter_func>(filter)){
    }

  private:
    double** &src;
    std::vector<std::array<double, 3>> storage;
  };

  class MCScalarSystem: public MCSystem{
  public:
    void save(size_t size) override {
      storage.clear();
      storage.reserve(size);

      for (auto i = 0; i < size; ++i)
        if (_filter(i)){
          storage.push_back(src[i]);
        }
      storage.shrink_to_fit();
    }

    void restore(size_t size) override {
      auto insert_iterator = storage.begin();
      for (auto i = 0; i < size; ++i)
        if (_filter(i)){
          src[i] = *insert_iterator;
          ++insert_iterator;
        }
    }

    size_t unpack(const double *data, size_t size, std::unordered_map<int, int> const &ghost_map) override;

    std::vector<double> pack(size_t size, int *tag) const override;

    void make(size_t size, mc_stepper &stepper) override;

    size_t object_size() const override {
      return 2;
    }

  public:
    MCScalarSystem(double* &source, filter_func &&filter):
        src(source), MCSystem(std::forward<filter_func>(filter)){
    }

  private:
    double* &src;
    std::vector<double> storage;
  };

  struct mc_stepper{
    MonteCarlo engine;
    double max_shift{0.0};
    stepper_type type;

    std::unique_ptr<MCSystem> system;
    std::unique_ptr<RanPark> random;

    explicit mc_stepper(LAMMPS *lmp, stepper_type type = stepper_type::ion_r, size_t seed = 42, size_t engine_seed = 42)
        : type(type), engine(engine_seed) {
      random = std::unique_ptr<RanPark>(new RanPark(lmp, seed));
    }

    void assign_subsystem(std::unique_ptr<MCSystem> &&ssystem){
      system = std::move(ssystem);
    }

    void save(size_t system_size){
      system->save(system_size);
    }

    void restore(size_t system_size){
      system->restore(system_size);
    }

    void make(size_t system_size){
      system->make(system_size, *this);
    }

    std::vector<double> pack(size_t system_size, int *tag){
      return system->pack(system_size, tag);
    }

    size_t unpack(double const *data, size_t size, std::unordered_map<int, int> const &ghost_map){
      return system->unpack(data, size, ghost_map);
    }

    inline double shift(){
      return 2.0 * max_shift * (random->uniform() - 0.5);
    }

    template <class TValue>
    void make_shift(TValue& v){
      auto dx = shift();
      v += dx;
      if (type == stepper_type::electron_w && v < 0)
        v = std::abs(v);
    }

    void adjust(){
      max_shift *= engine.adjust();
    }
  };

  class MCStepperSet{
  public:
    using stepper_container = std::vector<mc_stepper>;

    stepper_container& all() {
      return steppers;
    }

    mc_stepper& current() {
      return steppers[current_stepper];
    }

    void next(){
      current_stepper = (current_stepper + 1) % steppers.size();
    }

    template <class ... Targs>
    mc_stepper& add(Targs ... args){
      steppers.emplace_back(args...);
      return steppers.back();
    }

    size_t count() const {
      return steppers.size();
    }

    mc_stepper& get(size_t index){
      return steppers[index];
    }

    size_t current_stepped_id() const{
      return current_stepper;
    }
  private:
    stepper_container steppers;
    size_t current_stepper{0};
  };

}

#endif //LAMMPS_MC_UTILS_H
