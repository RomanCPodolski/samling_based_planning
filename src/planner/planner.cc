/*
 * parallel_planner.cc
 * Copyright (C) 2017 Roman C. Podolski <maito:roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#include "planner/planner.h"
#include <glog/logging.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include "planner/baseframe.h"
#include "planner/maneuver.h"
#include "planner/map.h"

DEFINE_double(max_manouver_offset, 4, "maximum offset from baseframe");
DEFINE_uint64(number_manouvers, 30, "how many manouvers should be considered");

/// TODO(roman): add debug output using glog
namespace planner {

typedef bg::model::box<point> box;
typedef std::pair<box, std::shared_ptr<polygon>> value;
typedef bgi::rtree<value, bgi::quadratic<16>> Rtree;

Planner::Planner(std::shared_ptr<Baseframe> const &base)
    : _m_max{static_cast<double>(FLAGS_number_manouvers)},
      _q_max{static_cast<double>(FLAGS_max_manouver_offset)},
      _baseframe{base},
      _rtree{new Rtree()},
      _last_known_s{0.0},
      _maneuvers(),
      _maneuver(_baseframe->P(0.0), _baseframe->theta(0.0), 0.0, 0.0, 0.0, 0.0,
                _baseframe, _rtree) {}

Planner::Planner(const Map &map)
    : Planner(std::shared_ptr<Baseframe>(new Baseframe(map))) {}

Planner::Planner(const std::vector<std::pair<double, double>> &wp)
    : Planner(std::shared_ptr<Baseframe>(new Baseframe(wp))) {}

Planner::Planner(const std::vector<point> &wp)
    : Planner(std::shared_ptr<Baseframe>(new Baseframe(wp))) {}

const Maneuver &Planner::update(const point position, const double heading,
                                const double velocity,
                                const double dt,  // TODO(roman): braucht es dt?
                                std::vector<polygon> obstacles) {
  _maneuvers.clear();

  // add obstacles to spatial index
  // TODO(roman): move to own function
  for (auto obstacle : obstacles) {
    box b = bg::return_envelope<box>(obstacle);
    auto p = std::make_shared<polygon>(obstacle);
    _rtree->insert(std::make_pair(b, p));
  }

  std::pair<double, double> p =
      _baseframe->localize(position, velocity, heading, dt, _last_known_s);
  _last_known_s = p.first;      // arc-length at current position
  const double q_i = p.second;  // lateral offset at current position

  for (double q_f = -_q_max; q_f <= _q_max; q_f += (2 * _q_max / _m_max)) {
    Maneuver m = Maneuver(position, heading, velocity, q_i, q_f, _last_known_s,
                          _baseframe, _rtree);

    if (m.driveable()) {  // only keep reasonable maneuvers
      _maneuvers.push_back(m);
    }
  }

  // throw exception on this
  CHECK(!_maneuvers.empty());

  std::map<double, Maneuver> candidates;
  std::map<double, Maneuver> collisions_paths;

  // collision check
  for (auto &m : _maneuvers) {
    if (m.collision()) {
      collisions_paths.insert(std::make_pair(m.collision(), m));
    } else {  // calculate path costs
      double cost =
          m.J(_maneuver, std::make_shared<std::vector<Maneuver>>(_maneuvers));
      candidates.insert(std::make_pair(cost, m));
    }
  }

  _maneuver = candidates.empty() ? collisions_paths.rbegin()->second
                                 : candidates.begin()->second;

  return _maneuver;
}
}  // namespace planner
