/*
 * Copyright (C) 2017 roman <roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <boost/geometry.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/random/normal_distribution.hpp>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include "planner/car.h"
#include "planner/io/beamer.h"
#include "planner/io/tex.h"
#include "planner/planner.h"

DEFINE_string(filename, "simulation.tex", "name of file");
DEFINE_double(road_width, 3.0, "Road with");
DEFINE_double(cone_dist, 5.0, "Cone dist");
DEFINE_double(cone_radius, 0.15, "Cone dist");
DEFINE_double(cone_sigma, 0.5, "cone placement standart deviation");
DEFINE_double(velocity, 2, "velocity of car");
DEFINE_double(frequency, 1, "update frequency");
DEFINE_double(position_sigma, 0.01, "standart deviation of the position");
DEFINE_double(heading_sigma, 0.01, "standart deviation of the heading");

using boost::math::constants::pi;
using boost::math::constants::half_pi;
using boost::math::constants::third_pi;
using boost::math::constants::sixth_pi;
using boost::math::constants::two_pi;
using boost::random::normal_distribution;

int main(int argc, char *argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  std::vector<std::pair<double, double>> waypoints{
      {0, 0},    {10, 0},   {20, 0},   {30, 10},  {40, 20},  {40, 30},
      {40, 40},  {40, 50},  {40, 60},  {20, 70},  {10, 80},  {0, 80},
      {-10, 70}, {-20, 70}, {-30, 75}, {-40, 80}, {-50, 80}, {-60, 70},
      {-60, 60}, {-60, 50}, {-50, 40}, {-45, 30}, {-40, 20}, {-30, 10},
      {-20, 0},  {-10, 0},  {0, 0}};
  double road_width = FLAGS_road_width;
  double cone_radius = FLAGS_cone_radius;
  double cone_dist = FLAGS_cone_dist;
  double cone_sigma = FLAGS_cone_sigma;

  std::ofstream tex(FLAGS_filename);
  planner::io::beamer_mapper mapper(tex, "xmin=-70,xmax=50,ymin=-10,ymax=90");

  planner::Planner jack(waypoints);

  std::vector<planner::polygon> obstacles;
  std::mt19937 rnd;

  double x_left, y_left, x_rigth, y_rigth, x1, x2, y1, y2;
  for (int j = 1; j < jack.baseframe().length() - 1; j += cone_dist) {
    auto point_1 = jack.baseframe().P(j - 0.1);
    auto point_2 = jack.baseframe().P(j + 0.1);

    Eigen::Vector2d v_1(point_1.x(), point_1.y());
    Eigen::Vector2d v_2(point_2.x(), point_2.y());

    Eigen::Vector2d n = v_2 - v_1;
    n.normalize();

    double x_left = jack.baseframe().P(j).x() + road_width * n[1];
    double y_left = jack.baseframe().P(j).y() - road_width * n[0];
    double x_rigth = jack.baseframe().P(j).x() - road_width * n[1];
    double y_rigth = jack.baseframe().P(j).y() + road_width * n[0];

    normal_distribution<> dist_x_left(x_left, cone_sigma);
    normal_distribution<> dist_y_left(y_left, cone_sigma);
    normal_distribution<> dist_x_right(x_rigth, cone_sigma);
    normal_distribution<> dist_y_right(y_rigth, cone_sigma);

    x_left = dist_x_left(rnd);
    y_left = dist_y_left(rnd);
    x_rigth = dist_x_right(rnd);
    y_rigth = dist_y_right(rnd);

    planner::polygon p1, p2;
    for (double i = 0; i < two_pi<double>(); i += sixth_pi<double>()) {
      x1 = x_left + cone_radius * cos(i);
      y1 = y_left + cone_radius * sin(i);
      x2 = x_rigth + cone_radius * cos(i);
      y2 = y_rigth + cone_radius * sin(i);
      p1.outer().push_back(planner::point(x1, y1));
      p2.outer().push_back(planner::point(x2, y2));
    }
    obstacles.push_back(p1);
    obstacles.push_back(p2);
  }

  planner::Car car;

  double x = waypoints[0].first, y = waypoints[0].second,
         heading = sixth_pi<double>();

  double s = 0.0;

  std::vector<planner::pose> optimal_path;

  while (s <= jack.baseframe().length() - 1) {
    std::cout << s << '\n';

    for (auto o : obstacles) {
      mapper.map(o);
    }

    mapper.map(jack.baseframe(), "color=blue, dotted");

    planner::point position(x, y);
    planner::pose pos{position, heading, 0, 0};

    auto m_opt = jack.update(pos.position, pos.theta, FLAGS_velocity,
                             1.0 / FLAGS_frequency, obstacles);
    for (auto m : jack.maneuvers()) {
      if (m.collision()) {
        mapper.map(m, "color=red");
      } else {
        mapper.map(m, "color=black");
      }
    }
    mapper.map(m_opt, "color=green");
    mapper.map(optimal_path, "color=green");

    mapper.map(car.obb(pos));
    normal_distribution<> n_postion_x(m_opt.path()[1].x(),
                                      FLAGS_position_sigma);
    normal_distribution<> n_postion_y(m_opt.path()[1].y(),
                                      FLAGS_position_sigma);
    normal_distribution<> n_heading(m_opt.path()[1].theta, FLAGS_heading_sigma);
    optimal_path.push_back(m_opt.path()[1]);
    x = n_postion_x(rnd);
    y = n_postion_y(rnd);
    heading = n_heading(rnd);
    s = m_opt.path()[1].s;
    mapper.next_frame();
  }

  // std::ofstream bla("last_frame.tex");
  // planner::io::tex_mapper tex_m(bla);

  for (auto o : obstacles) {
    // tex_m.map(o);
    mapper.map(o);
  }

  mapper.map(jack.baseframe(), "color=blue, dotted");
  mapper.map(optimal_path, "color=green");
}
