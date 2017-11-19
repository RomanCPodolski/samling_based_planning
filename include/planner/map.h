/*
 * map.h
 * Copyright (C) 2017 Roman C. Podolski <maito:roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef MAP_H
#define MAP_H

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "planner/typedefs.h"

DECLARE_string(map);

namespace planner {

// TODO(roman): add structs for waypoints and fenceposts and define how to
// implicit cast them in required types like double* or matlab*

// static const double earth_radius = 6371000.0;
// using boost::math::constants::two_pi;

// TODO(roman): wrap std::pair<double, double> as a type Point

/*! \class Map
 *  \brief Brief class description
 *  \author Roman C. Podolski - <roman.podolski@tum.de>
 *
 *  Detailed description
 */
class Map {
 private:
  std::vector<double> _waypoints_x;
  std::vector<double> _waypoints_y;
  std::vector<double> _fenceposts_x;
  std::vector<double> _fenceposts_y;
  std::pair<double, double> _datum;

 public:
  explicit Map(std::string path = FLAGS_map);
  virtual ~Map() {}
  void load(std::string path = FLAGS_map);
  void parse_from_stream(std::istream *is);
  void save(std::string path = FLAGS_map);

  // TODO(roman): manage access to private members for the planner
  std::vector<std::pair<double, double>> get_way_points() const {
    auto x = get_waypoints_x();
    auto y = get_waypoints_y();
    CHECK_EQ(x.size(), y.size())
        << "way point coordinates are not of equal size";
    std::vector<std::pair<double, double>> waypoints;
    std::pair<double, double> fp;
    for (int i = 0; i < x.size(); i++) {
      fp = std::pair<double, double>(x[i], y[i]);
      waypoints.push_back(fp);
    }

    return waypoints;
  }

  std::vector<std::pair<double, double>> get_fenceposts() const {
    auto x = get_fenceposts_x();
    auto y = get_fenceposts_y();
    CHECK_EQ(x.size(), y.size())
        << "way point coordinates are not of equal size";
    std::vector<std::pair<double, double>> fenceposts;
    std::pair<double, double> fp;
    for (int i = 0; i < x.size(); i++) {
      fp = std::pair<double, double>(x[i], y[i]);
      fenceposts.push_back(fp);
    }

    return fenceposts;
  }

  const std::vector<double> get_waypoints_x() const { return _waypoints_x; }
  const std::vector<double> get_waypoints_y() const { return _waypoints_y; }
  const std::vector<double> get_fenceposts_x() const { return _fenceposts_x; }
  const std::vector<double> get_fenceposts_y() const { return _fenceposts_y; }

  void add_waypoint(const std::pair<double, double> &wp) {
    _waypoints_x.push_back(wp.first);
    _waypoints_y.push_back(wp.second);
  }

  void add_fencepost(const std::pair<double, double> &fp) {
    _fenceposts_x.push_back(fp.first);
    _fenceposts_y.push_back(fp.second);
  }
};

}  // namespace planner

#endif /* !MAP_H */
