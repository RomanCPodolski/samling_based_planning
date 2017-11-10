/*
 * playgound_main.cc
 * Copyright (C) 2017 roman <roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "planner/io/beamer.h"
#include "planner/map.h"
#include "planner/planner.h"

using boost::math::constants::pi;
using boost::math::constants::two_pi;

int main(int argc, char *argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  // planner::Map map("data/simple.csv");

  // std::vector<std::pair<double, double>> waypoints_a{
  //{0, 10},  {10, 10}, {20, 10}, {30, 10}, {40, 10}, {50, 10},
  //{60, 10}, {70, 20}, {80, 30}, {90, 40}, {100, 50}};

  // std::vector<std::pair<double, double>> waypoints{
  //{0, 0},  {10, 0}, {20, 0}, {30, 0}, {40, 0}, {50, 0},
  //{60, 0}, {70, 0}, {80, 0}, {90, 0}, {100, 0}};

  std::ofstream tex("obstacles.tex");
  planner::io::beamer_mapper mapper(tex);

  std::vector<std::pair<double, double>> waypoints{{0, 0},  {10, 0}, {20, 0},
                                                   {30, 0}, {31, 0}, {32, 0}};
  planner::Planner p(waypoints);

  planner::pose pos{{5, 2.0}, pi<double>() / 16};

  // create some obstacles
  std::vector<planner::polygon> obstacles;
  for (int i = 0; i < 40; i += 10) {
    planner::polygon p1;
    for (double j = 0; j < two_pi<double>(); j += pi<double>() / 4.0) {
      double x = i + cos(j);
      double y = 5 + sin(j);
      p1.outer().push_back(planner::point(x, y));
    }

    planner::polygon p2;
    for (double j = 0; j < two_pi<double>(); j += pi<double>() / 4.0) {
      double x = i + cos(j);
      double y = -5 + sin(j);
      p2.outer().push_back(planner::point(x, y));
    }
    obstacles.push_back(p1);
    obstacles.push_back(p2);
  }

  p.set_last_known_s(0.0);

  auto my = p.update(pos.position, pos.theta, 10, 1, obstacles);

  // mapper.map(waypoints);

  // for (auto m : p.maneuvers()) {
  // if (m.collision()) {
  // mapper.map(m, "color=red");
  //} else {
  // mapper.map(m, "color=black");
  //}
  //}
  // mapper.map(my, "color=green");
  // for (const auto &o : obstacles) {
  // mapper.map(planner::bg::return_envelope<planner::box>(o), "color=red");
  //}

  for (const auto &pi : my.path()) {
    mapper.map(p.baseframe());
    mapper.map(my, "color=green");

    for (const auto &o : obstacles) {
      mapper.map(o);
    }

    for (const auto &pipi : my.path()) {
      mapper.map(pipi);
    }
    mapper.map(pi, "color=red");
    mapper.next_frame();
  }
}
