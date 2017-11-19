/*
 * manouver.cc
 * Copyright (C) 2017 romancpodolski <mailto:roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#include "planner/maneuver.h"
#include <gflags/gflags.h>
#include <boost/format.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include "planner/baseframe.h"
#include "planner/car.h"

DEFINE_double(min_manouver_length, 20.0, "Minimal length of manouver in [m]");
DEFINE_double(manouver_speed_gain, 1.0,
              "Gain of manouver length from vehicle speed");
DEFINE_double(weight_safety, 1.0, "weigthing factor for the safety cost");
DEFINE_double(weight_smoothness, 1.0,
              "weigthing factor for the smoothness cost");
DEFINE_double(weight_consistency, 1.0,
              "weigthing factor for the consistency cost");
DEFINE_double(granularity, 1.0, "Stepsize of the path generation");
DEFINE_double(collision_standart_deviation, 1.0,
              "standart deviation of risk of collision");
DEFINE_double(max_curvature, 0.5, "standart deviation of risk of collision");

namespace planner {

Maneuver::Maneuver(const point position, const double heading,
                   const double velocity, const double q_i, const double q_f,
                   const double s_i,
                   const std::shared_ptr<Baseframe> &baseframe,
                   const std::shared_ptr<Rtree> &rtree)
    : _q_f{q_f},
      _s_i{s_i},
      _baseframe{baseframe},
      _rtree{rtree},
      _theta{heading - baseframe->theta(s_i)},
      _v{velocity},
      _position{position},
      _heading{heading},
      _collision_checked{false},
      _drivable{true},
      _collision_length{0.0} {
  const double delta_s_f =
      FLAGS_manouver_speed_gain * velocity + FLAGS_min_manouver_length;
  s_f(_s_i + delta_s_f);

  Eigen::Matrix4d A;
  Eigen::Vector4d b, c;
  A << 0, 0, 0, 1,                                         // end first row
      pow(delta_s_f, 3), pow(delta_s_f, 2), delta_s_f, 1,  // second row
      0, 0, 1, 0,                                          // third
      3 * pow(delta_s_f, 2), 2 * delta_s_f, 1, 0;          // fourth
  c << q_i, q_f, tan(_theta), 0.0;
  b = A.inverse() * c;
  _a = b(0), _b = b(1), _c = b(2), _d = b(3);
}

double Maneuver::q(double s) const {
  double result = 0.0;
  double delta_s = s - _s_i;
  if (_s_i <= s && s < _s_f) {
    result = _a * pow(delta_s, 3) + _b * pow(delta_s, 2) + _c * delta_s + _d;
  } else if (_s_f <= s) {
    result = _q_f;
  } else {
    // result = _d;
    throw std::domain_error("q is not defined for s < " + std::to_string(_s_i));
  }
  return result;
}

double Maneuver::dq(double s) const {
  double result = 0.0;
  double delta_s = s - _s_i;

  if (_s_i <= s && s < _s_f) {
    result = 3 * _a * pow(delta_s, 2) + 2 * _b * delta_s + _c;
  } else if (_s_f <= s) {
    result = 0.0;  // this is redundant - but makes the code more readable
  } else {
    // result = 0.0;  // this is redundant - but makes the code more readable
    throw std::domain_error("dq/ds is not defined for s < " +
                            std::to_string(_s_i));
  }
  return result;
}

double Maneuver::ddq(double s) const {
  double result = 0.0;
  double delta_s = s - _s_i;

  if (_s_i <= s && s < _s_f) {
    result = 6 * _a * delta_s + 2 * _b;
  } else if (_s_f <= s) {
    result = 0.0;  // this is redundant - but makes the code more readable
  } else {
    // result = 0.0;  // this is redundant - but makes the code more readable
    throw std::domain_error("d^2q/ds^2 is not defined for s < " +
                            std::to_string(_s_i));
  }
  return result;
}

double Maneuver::Q(const double s) const {
  return sqrt(pow(dq(s), 2) + pow((1 - q(s) * _baseframe->curvature(s)), 2));
}

double Maneuver::S(const double s) const {
  return sign(1 - q(s) * _baseframe->curvature(s));
}

double Maneuver::curvature(const double s) const {
  double K_b = _baseframe->curvature(s);

  double k =
      (S(s) / Q(s)) *
      (K_b + ((1 - q(s) * K_b) * ddq(s) + K_b * pow(dq(s), 2)) / pow(Q(s), 2));

  return k;
}

/*! \struct push_back_state_and_arc_length
 *  \brief Brief struct description
 *
 *  Detailed description
 */
struct push_back_state_and_arc_length {
  std::vector<pose> &_path; /*!< Description */

  push_back_state_and_arc_length(std::vector<pose> &path) : _path(path) {}

  void operator()(const state_type &x, double s) {
    _path.push_back(pose{point(x[0], x[1]), x[2], s});
  }
};

std::vector<pose> Maneuver::path() {
  // TODO(roman): make sure that the path does not extend L without setting s_f
  // to L
  if (_path.empty()) {
    double ds = FLAGS_granularity;

    state_type x_0 = {_position.x(), _position.y(), _heading};

    auto system = [this](const state_type &x, state_type &dxds,
                         const double s) {
      double K_b{_baseframe->curvature(s)};
      double K{curvature(s)};

      if (fabs(q(s)) >= fabs(1 / K_b) || fabs(K) > FLAGS_max_curvature) {
        _drivable = false;
      }

      double Q_s{Q(s)};
      double theta{x[2]};

      dxds[0] = Q_s * cos(theta);
      dxds[1] = Q_s * sin(theta);
      dxds[2] = Q_s * K;
    };

    boost::numeric::odeint::euler<state_type> stepper;

    boost::numeric::odeint::integrate_const(
        stepper, system, x_0, s_i(), s_f(), ds,
        push_back_state_and_arc_length(_path));
  }
  return _path;
}

double Maneuver::J(const Maneuver &previous,
                   const std::shared_ptr<std::vector<Maneuver>> maneuvers,
                   double w_s, double w_k, double w_c) const {
  return w_s * C_s(maneuvers) + w_k * C_k() + w_c * C_c(previous);
}

double Maneuver::C_s(const std::shared_ptr<std::vector<Maneuver>> maneuvers,
                     const double sigma) const {
  double result = 0.0;
  normal g(q_f(), sigma);

  for (auto &m : *maneuvers) {
    result += static_cast<bool>(m.collision()) * pdf(g, m.q_f());
  }

  return result;
}

double Maneuver::C_k() const {
  auto l = [this](const double s) { return pow(curvature(s), 2) * Q(s); };
  return trapezoidal(l, s_i(), s_f());
}

double Maneuver::C_c(const Maneuver &previous) const {
  const double s_1 = s_i();
  const double s_2 = s_f();
  // TODO(roman): handle this better!
  if (s_1 > previous.s_f()) {
    throw std::domain_error("Maneuvers do not overlap!");
  }
  auto l = [this, &previous](double s) {
    point a(s, q(s));
    point b(s, previous.q(s));
    return bg::distance(a, b);
  };
  double result = 1 / (s_2 - s_1) * trapezoidal(l, s_1, s_2);
  return result;
}

double Maneuver::collision() {
  // if a collision check was already performed, return the result
  if (_collision_checked) return _collision_length;
  _collision_length = 0.0;
  _collision_checked = true;

  Car car;

  for (const auto &p : path()) {
    // search spatial index for possible collisions
    std::vector<value> collision_risk;
    _rtree->query(bgi::intersects(car.aabb(p)),
                  std::back_inserter(collision_risk));
    // perfom the actual collision check.
    for (const auto &risk : collision_risk) {
      polygon obstacle{*risk.second};
      if (bg::intersects(car.obb(p), obstacle)) {  // collision detected
        _collision_length = p.s;  // return length at collison and stop loop.
        return _collision_length;
      }
    }
  }
  return _collision_length;  // return 0.0 - no collision detected
}

std::ostream &operator<<(std::ostream &os, const Maneuver &m) {
  os << "\\addplot coordinates {\n";
  for (auto point : m.path()) {
    os << boost::format("(%d,%d)\n") % point.x() % point.y();
  }
  os << "};\n";
  return os;
}

}  // namespace planner
