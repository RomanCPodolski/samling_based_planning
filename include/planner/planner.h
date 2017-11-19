/*
 * parallel_planner.h
 * Copyright (C) 2017 Roman C. Podolski <maito:roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef PLANNER_PLANNER_H
#define PLANNER_PLANNER_H

#include <glog/logging.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <string>
#include <utility>
#include <vector>
#include "planner/baseframe.h"
#include "planner/maneuver.h"
#include "planner/map.h"
#include "planner/typedefs.h"

DECLARE_double(max_steering_angle);
DECLARE_double(max_manouver_offset);
DECLARE_double(granularity);
DECLARE_uint64(number_manouvers);

namespace planner {

/*!
 * \class Planner
 * \brief A local path planner for autonomous driving.
 * \see Maneuver, Baseframe
 *
 * \author Roman C. Podolski - <roman.podolski@tum.de>
 *
 * This is a implementation of the path planner introduced in
 * [[Chu2012]]. This path planner is suited for Off-road autonomous
 * driving with static obstacle avoidance. It is assumed that a set of precise
 * waypoints is provided. The planner itself will compute a collision-free path
 * that minimizes the costfunction Maneuver::J. If no collision free path can be
 * found, the longest path without collision is choosen. The waypoints provide
 * the base frame of a curvilinear coordinate system to generate path candidates
 * for autonomous vehicle path planning. Each path is converted to a Carthesian
 * coordinate system and evaluated using obstacle data. To select the optimal
 * path, the priority of each path is determined by considering the weigted sum
 * Maneuver::J of the path safety cost Maneuver::C_s, path smoothness
 * Maneuver::C_k, and path consistency Maneuver::C_c.
 *
 * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
 */
class Planner {
 public:
  /*!
   * \brief Constructs a planner from waypoints loaded in a Map.
   * \see Map, Maneuver, Baseframe
   *
   * Sets Planner::_q_max, Planner::_m_max according to the values specified
   * from the cmd-flags. Calls the Baseframe::Baseframe(Map) constructor and
   * stores a pointer to the object in a shared pointer. The last known vehicle
   * position (Planner::_last_known_s) is assumed to be at the beginnig of the
   * Baseframe. The last/first maneuver (Planner::_maneuver) is assumed to be a
   * maneuver with offset 0.0 from the baseframe and the heading equal to the
   * tangent-angle of the baseframe.
   */
  explicit Planner(const Map& map);
  explicit Planner(const std::vector<std::pair<double, double>>& wp);
  explicit Planner(const std::vector<point>& wp);
  explicit Planner(std::shared_ptr<Baseframe> const& baseframe);

  /*! \brief Compute a path that minimizes the cost function Maneuver::J.
   *
   * \param position position of vehicle in Carthesian coordinates given in the
   * local coordinatesystem. X and Y coordinates given in [m].
   * \param heading orientation of the vehicle in [rad].
   * \param velocity forward velocity of the vehicle in [m/s].
   * \param dt time difference since last update in [s].
   * \param obstacles Polygons represent obstacles to avoid. All dimensions in
   * Carthesian coordinates in [m].
   *
   * TODO
   *
   */
  const Maneuver& update(const point position, const double heading,
                         const double velocity, const double dt,
                         std::vector<polygon> obstacles);

  /*!
   * \brief set the last known location on the Baseframe.
   * \param s coordinate to set the last known position to
   *
   * mainly exists for testability reasons.
   */
  void set_last_known_s(const double s) { _last_known_s = s; }

  const Baseframe& baseframe() const noexcept { return *_baseframe; }

  const std::vector<Maneuver>& maneuvers() const noexcept { return _maneuvers; }

 private:
  std::shared_ptr<Baseframe> _baseframe;  ///< Baseframe for path computation.
  const double _m_max;                    ///< quanity of paths to generate.
  const double _q_max;   ///< maximal offset of a path from the Baseframe [m].
  double _last_known_s;  ///< last known location on the Baseframe [m].

  /*!
   * Spatial index for obstacle collision detection. See
   * [boost::geometry::index] for implementation details and [2] for details on
   * R-Trees.
   *
   * [boost::geometry::index]:http://www.boost.org/doc/libs/1_65_1/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html
   * [2]: https://dl.acm.org/citation.cfm?id=602266
   */
  std::shared_ptr<Rtree> _rtree;
  std::map<double, Maneuver> _candidates;
  std::map<double, Maneuver> _collisions_paths;
  std::vector<Maneuver> _maneuvers;
  Maneuver _maneuver;  ///< Maneuver computed by latest call to Planner::update.
};
}  // namespace planner

#endif  // PLANNER_PLANNER_H_
