/*
 * baseframe.cc
 * Copyright (C) 2017 roman <roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#include "planner/baseframe.h"
#include <gflags/gflags.h>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include "planner/map.h"

namespace planner {

using boost::math::quadrature::trapezoidal;
using boost::math::cubic_b_spline;
using boost::math::tools::bisect;
using boost::math::tools::eps_tolerance;
using boost::math::tools::newton_raphson_iterate;

Baseframe::Baseframe(const Map &map) {
  const auto vx = map.get_waypoints_x();
  const auto vy = map.get_waypoints_y();
  init(vx, vy);
}

Baseframe::Baseframe(const std::vector<std::pair<double, double>> &wp) {
  std::vector<double> vx, vy;
  for (const auto &p : wp) {
    vx.push_back(p.first);
    vy.push_back(p.second);
  }
  init(vx, vy);
}

Baseframe::Baseframe(const std::vector<point> &wp) {
  std::vector<double> vx, vy;
  for (const auto &p : wp) {
    vx.push_back(p.x());
    vy.push_back(p.y());
  }
  init(vx, vy);
}

void Baseframe::init(const std::vector<double> vx,
                     const std::vector<double> vy) {
  CHECK_EQ(vx.size(), vy.size());
  _N = vx.size() - 1;
  _x = cubic_b_spline<double>(vx.begin(), vx.end(), 0, 1.0);
  _y = cubic_b_spline<double>(vy.begin(), vy.end(), 0, 1.0);
}

point Baseframe::Q(const double t) const {
  if (t < 0 || t > _N)
    throw std::domain_error("t = " + std::to_string(t) +
                            " not within the domain [0,N]");
  return point(_x(t), _y(t));
}

point Baseframe::Q_prime(const double t) const {
  if (t < 0 || t > _N)
    throw std::domain_error("t = " + std::to_string(t) +
                            " not within the domain [0,N]");
  return point(_x.prime(t), _y.prime(t));
}

point Baseframe::Q_pprime(const double t) const {
  const double ddt_x =
      derivative(t, 0.001, [this](const double t) { return _x.prime(t); });
  const double ddt_y =
      derivative(t, 0.001, [this](const double t) { return _y.prime(t); });
  return point(ddt_x, ddt_y);
}

double Baseframe::theta_Q(const double t) const {
  auto p = Q_prime(t);
  // return atan2(p.y(), p.x());
  return atan(p.y() / p.x());
}

static std::map<double, double> _A_lut{{0, 0}};
static std::vector<double> _segment_lut;
double Baseframe::A(const double t) const {
  if (t < 0 || t > _N)
    throw std::domain_error("t = " + std::to_string(t) +
                            " not within the domain [0,N]");
  if (_A_lut.find(t) != _A_lut.end()) {
    return _A_lut[t];
  }

  double result = 0.0;
  if (t == 0) return result;

  auto f = [this](double t) {
    return sqrt(pow(_x.prime(t), 2) + pow(_y.prime(t), 2));
  };

  result = trapezoidal(f, 0.0, t);

  _A_lut[t] = result;

  return result;
}

static std::map<double, double> _inv_A_lut{{0, 0}};
double Baseframe::inv_A(const double s) const {
  if (s < 0 || s > length())
    throw std::domain_error("s = " + std::to_string(s) +
                            " not within the domain [0,L]");

  if (_segment_lut.empty()) {
    for (int t = 0; t < _N; ++t) {
      _segment_lut.push_back(A(t));
    }
  }

  if (_inv_A_lut.find(s) != _inv_A_lut.end()) {
    return _inv_A_lut[s];
  }

  std::uintmax_t iterations = 1000;
  auto low = std::lower_bound(_segment_lut.begin(), _segment_lut.end(), s);
  const double segment = low - _segment_lut.begin();
  double result = 0.0;
  result = bisect([this, s](double t) { return A(t) - s; }, segment - 1,
                  segment, eps_tolerance<float>(), iterations)
               .first;

  _inv_A_lut[s] = result;

  return result;
}

static std::map<double, double> _curvature_lut;
double Baseframe::curvature(const double s) {
  if (_curvature_lut.find(s) != _curvature_lut.end()) {
    return _curvature_lut[s];
  }

  const auto d_P = P_prime(s);
  double d_x = d_P.x();
  double d_y = d_P.y();

  if (d_x == 0.0 && d_y == 0.0) {  // avoid division by 0.0
    d_x = std::numeric_limits<double>::denorm_min();
    d_y = std::numeric_limits<double>::denorm_min();
  }

  const auto dd_P = P_pprime(s);
  const double dd_x = dd_P.x();
  const double dd_y = dd_P.y();

  const double kappa =
      (d_x * dd_y - dd_x * d_y) / pow((pow(d_x, 2) + pow(d_y, 2)), 3.0 / 2.0);

  _curvature_lut[s] = kappa;

  return kappa;
}

// TODO(roman): check for numerical problems
std::pair<double, double> Baseframe::localize(const point position,
                                              const double v,
                                              const double heading,
                                              const double dt,
                                              const double last_known_s) const {
  auto D_prime = [ this, x_0 = position.x(), y_0 = position.y() ](double s) {
    auto p = P(s);
    double x = p.x();
    double y = p.y();

    auto p_prime = P_prime(s);
    double x_prime = p_prime.x();
    double y_prime = p_prime.y();

    auto p_pprime = P_pprime(s);
    double x_pprime = p_pprime.x();
    double y_pprime = p_pprime.y();

    double dist_prime = 2 * (x - x_0) * x_prime + 2 * (y - y_0) * y_prime;
    double dist_pprime = 2 * ((x - x_0) * x_pprime + pow(x_prime, 2) +
                              (y - y_0) * y_pprime + pow(y_prime, 2));

    return std::pair<double, double>(dist_prime, dist_pprime);
  };

  std::uintmax_t iterations = 1000;
  double minimum_extra_offset = 10;
  double lower_barrier = last_known_s - minimum_extra_offset;
  if (lower_barrier <= 0.0) lower_barrier = 0.0;

  double initial_guess = last_known_s + v * dt;
  if (initial_guess > length()) initial_guess = length();
  double upper_barrier = initial_guess + minimum_extra_offset;
  if (upper_barrier > length()) upper_barrier = length();

  auto result = newton_raphson_iterate(
      D_prime, lower_barrier, initial_guess, upper_barrier,
      std::numeric_limits<double>::digits / 2, iterations);

  auto p = P(result);

  return std::pair<double, double>(result, bg::distance(p, position));
}

}  // namespace planner
