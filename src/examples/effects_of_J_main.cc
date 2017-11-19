#include <gflags/gflags.h>
#include <boost/math/constants/constants.hpp>
#include <boost/random/normal_distribution.hpp>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <random>
#include <utility>
#include <vector>
#include "planner/io/tex.h"
#include "planner/planner.h"
#include "planner/pose.h"

using boost::math::constants::pi;
using boost::math::constants::two_pi;
using boost::math::constants::sixth_pi;
using boost::random::normal_distribution;

DEFINE_string(filename, "effects_of_J.tex", "name of file");
DEFINE_double(road_width, 4.0, "Road with");
DEFINE_double(cone_dist, 5.0, "Cone dist");
DEFINE_double(cone_radius, 0.15, "Cone dist");
DEFINE_double(cone_sigma, 0.5, "cone placement standart deviation");
DEFINE_double(velocity, 2, "velocity of car");
DEFINE_double(frequency, 1, "update frequency");

int main(int argc, char *argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  std::ofstream ofs(FLAGS_filename);
  planner::io::tex_mapper mapper(ofs);

  std::mt19937 rnd;

  std::vector<std::pair<double, double>> waypoints;
  for (double x = 0; x < 30; ++x) {
    double y = 0.8 * exp(x / 8.0);
    waypoints.push_back(std::make_pair(x, y));
    // std::cout << x << ' ' << y << std::endl;
  }
  planner::Planner plan(waypoints);

  planner::pose pos{planner::point(7, 3), 0, 10, 0};

  double road_width = FLAGS_road_width;
  double cone_radius = FLAGS_cone_radius;
  double cone_dist = FLAGS_cone_dist;
  double cone_sigma = FLAGS_cone_sigma;

  std::vector<planner::polygon> obstacles = {};

  double x_left, y_left, x_rigth, y_rigth, x1, x2, y1, y2;
  for (int j = 1; j < plan.baseframe().length() - 1; j += cone_dist) {
    auto point_1 = plan.baseframe().P(j - 0.1);
    auto point_2 = plan.baseframe().P(j + 0.1);

    Eigen::Vector2d v_1(point_1.x(), point_1.y());
    Eigen::Vector2d v_2(point_2.x(), point_2.y());

    Eigen::Vector2d n = v_2 - v_1;
    n.normalize();

    double x_left = plan.baseframe().P(j).x() + road_width * n[1];
    double y_left = plan.baseframe().P(j).y() - road_width * n[0];
    double x_rigth = plan.baseframe().P(j).x() - road_width * n[1];
    double y_rigth = plan.baseframe().P(j).y() + road_width * n[0];

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

  // mapper.map(plan.baseframe());
  mapper.map(pos);

  plan.set_last_known_s(pos.s);
  plan.update(planner::point(0.0, -2.0), 0, 5, 1, obstacles);
  auto m_opt = plan.update(pos.position, pos.theta, 5, 1, obstacles);

  for (auto o : obstacles) {
    mapper.map(o);
  }

  for (auto m : plan.maneuvers()) {
    if (m.collision()) {
      mapper.map(m, "color=red");
    } else {
      mapper.map(m, "color=black");
    }
  }

  mapper.map(m_opt, "color=green");

  return 0;
}
