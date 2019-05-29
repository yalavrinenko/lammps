//
// Created by yalavrinenko on 29.05.19.
//

#ifndef LAMMPS_TIMEMETRICS_H
#define LAMMPS_TIMEMETRICS_H
#include <chrono>
#include <iostream>

class TimeMetrics{
public:
  TimeMetrics(){
    make_tick("Init");
  }

  template <class tick_type = decltype(std::chrono::high_resolution_clock::now())>
  struct timepoint{
    std::string label;
    tick_type point;

    timepoint(std::string const &l, tick_type const &tick): label(l), point(tick){}
  };

  void make_tick(std::string const &label){
    points.emplace_back(label, std::chrono::high_resolution_clock::now());
  }

  void summary() const {
    for (auto i = 1; i < points.size(); ++i){
      std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(points[i].point - points[i-1].point);
      std::cout << points[i-1].label << "-" << points[i].label << ":\t" << duration.count() << " ms." << std::endl;
    }
  }

  ~TimeMetrics(){
    summary();
  }
private:
  std::vector<timepoint<>> points;
};

#endif //LAMMPS_TIMEMETRICS_H
