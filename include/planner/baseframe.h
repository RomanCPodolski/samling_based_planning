/*
 * baseframe.h
 * Copyright (C) 2017 roman <roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef PLANNER_BASEFRAME_H_
#define PLANNER_BASEFRAME_H_

#include <glog/logging.h>
#include <boost/math/constants/constants.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <cmath>
#include <utility>
#include <vector>
#include "planner/map.h"
#include "planner/typedefs.h"

namespace planner {

/*!
 * TODO(roman): add doc
 * code snipped taken from
 * http://www.boost.org/doc/libs/1_65_0/libs/multiprecision/doc/html/boost_multiprecision/tut/floats/fp_eg/nd.html
 */
template <typename value_type, typename function_type>
inline const value_type derivative(const value_type x, const value_type dx,
                                   function_type func) {
  // Compute d/dx[func(*first)] using a three-point
  // central difference rule of O(dx^6).

  const value_type dx1 = dx;
  const value_type dx2 = dx1 * 2;
  const value_type dx3 = dx1 * 3;

  const value_type m1 = (func(x + dx1) - func(x - dx1)) / 2;
  const value_type m2 = (func(x + dx2) - func(x - dx2)) / 4;
  const value_type m3 = (func(x + dx3) - func(x - dx3)) / 6;

  const value_type fifteen_m1 = 15 * m1;
  const value_type six_m2 = 6 * m2;
  const value_type ten_dx1 = 10 * dx1;

  return ((fifteen_m1 - six_m2) + m3) / ten_dx1;
}

using boost::math::cubic_b_spline;
using boost::math::constants::pi;

/*! \class Baseframe
 *  \brief Baseframe for the curvilinear coordinate system
 *  \author Roman C. Podolski - <roman.podolski@tum.de>
 *  \see Planner, Maneuver, [B-spline], [Wang2002], [WangKearney2002]
 *
 * This provides the base for the curvilinear coordinatesystem in which the
 * Planner generates Maneuvers. A parametric 2D-cubic-B-spline is employed to
 * model the geometry of the road. It is assumed that a
 * precise set of waypoints is given, trough which the baseframe passes.
 * The baseframe is parametrizied by \f$t\f$ where \f$t \in [0,N]\f$ and \f$N\f$
 * is the number of waypoints, or the arc-length \f$s\f$ where \f$s \in [0,L]\f$
 * and \f$L\f$ is the lenght of the entire baseframe [[Wang2002]].
 *
 * [Wang2002]:http://homepage.divms.uiowa.edu/~kearney/pubs/CurvesAndSurfacesArcLength.pdf
 * [WangKearney2002]:http://homepage.divms.uiowa.edu/~kearney/pubs/CurvesAndSufacesClosestPoint.pdf
 * [B-spline]: https://en.wikipedia.org/wiki/B-spline
 */
class Baseframe {
 public:
  /*!
   * \brief Construct baseframe from waypoints stored in a Map.
   * \see Map, [B-spline], [boost::math::cubic_b_spline]
   *
   * Generates a 2D-B-spline from a set of waypoints stored in
   * a map. The generated spline is implemented using
   * [boost::math::cubic_b_spline].
   *
   * [boost::math::cubic_b_spline]:http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/interpolate/cubic_b.html.
   * [B-spline]: https://en.wikipedia.org/wiki/B-spline
   */
  explicit Baseframe(const Map &map);
  explicit Baseframe(const std::vector<std::pair<double, double>> &wp);
  explicit Baseframe(const std::vector<point> &wp);

  /*!
   * \brief curvature of the baseframe at arc-length \f$s\f$.
   * \param s where \f$s \in [0,L]\f$ and \f$L\f$ is the length of the baseframe
   * in [m].
   * \see P_prime, P_pprime, [Curvature]
   * \returns curvature of the baseframe at arc-length \f$s\f$.
   * \throws std::domain_error if \f$s \notin [0,L] \f$.
   *
   * \f[ \kappa_b = \frac{x'y'' - x''y'}{(x'^2+y'^2)^\frac{3}{2}}  \f]
   *
   * The curvature of the baseframe at arc-length \f$s\f$ where \f$s \in
   * [0,L]\f$ and \f$L\f$ is the length of the baseframe
   * in [m], can be calulated from the first and second derivatices of the
   * cubic-B-spline. For a definition of curvature in this sense please refere
   * to [Curvature].
   *
   * [Curvature]: https://en.wikipedia.org/wiki/Curvature
   */
  double curvature(const double s);

  /*! \brief calculates the tangent angle in [rad] of baseframe at arc-length
   * \f$s\f$.
   * \param s where \f$s \in [0,L]\f$ and \f$L\f$ is the length of the baseframe
   * in [m].
   * \see P_prime, theta_P, theta_Q
   * \returns the tangent angle in [rad] of baseframe at arc-length \f$s\f$.
   * \throws std::domain_error if \f$s \notin [0,L] \f$.
   *
   * \f[\tan(\theta) = \frac{y'}{x'} \rightarrow \theta =
   * \mathrm{atan2}(y',x')\f]
   *
   * this can be derived from \f$x' = \cos(\theta)\f$ and \f$y' =
   * \sin(\theta)\f$. The use of atan2 gurantees the right sign of the angle.
   *
   * this function is a name alias for theta_P.
   */
  double theta(const double s) const { return theta_P(s); }

  /*! \brief getter for length of the Baseframe in [m].
   * \returns length of the Baseframe in [m].
   *
   * \see A, _N, [Wang2002]
   *
   * This is equal to calculating \f[A(N)\f]
   *
   * where \f$A(\cdot)\f$ is the non-linear transformation from the parameter
   * space of
   * \f$t\f$ to \f$s\f$, and N is the number of waypoints and therefore the last
   * waypoint. This calculation yields the arc-length at the last waypoint and
   * therefore the length of the entire curve.
   *
   * [Wang2002]:http://homepage.divms.uiowa.edu/~kearney/pubs/CurvesAndSurfacesArcLength.pdf
   */
  double length() const { return A(_N); }

  /*! \brief Localization of a vehicle in the curvilinear-coordinatesystem.
   * \see [WangKearney2002],
   * [boost::math::tools::roots],
   * [boost::math::tools::roots::newton_raphson_iterate], [boost::geometry]
   * \returns std::pair containing the position of the vehicle in curvilinear
   * coordinates - first is arc-length \f$s\f$ along the baseframe and second is
   * the offset from the baseframe \f$q\f$.
   *
   * \param position position of the vehicle in Carthesian coordinates [m].
   * \param velocity of the vehicle [m/s].
   * \param heading of the vehicle [rad].
   * \param delta_t time difference from the last call to Planner::update [s].
   * \param last_known_s the last known vehicle location on the baseframe [m].
   *
   * Localizing a vehicle in the curvilinear coordinatesystem is equal to
   * minimizing the squared distance between the vehilce posistion and the
   * baseframe w.r.t the arc-length \f$s\f$.
   *
   * \f[\min_{s \in [0,L]} d(s) = (x_v - x(t))^2 + (y_v - y(t))^2\f]
   *
   * this is equal to finding the root of the function
   *
   *  \f[\frac{\mathrm{d}}{\mathrm{d}s}d(s) = 2(x_v - x(t))x'(t) + 2(y_v -
   * y(t))y'(t)\f]
   *
   *
   * This can be archieved with any root finding algorithm.
   * [boost::math::tools::roots] implements a number of possible choises.
   * However, because fast convergence is a criterium for the localisation,
   * [Newton-Raphson] is the algorithm of choise in this scenario. The
   * performance of the root-finding algortihm is highly-depentend on the
   * initial guess for the root and the boundaries in which the algorihmn
   * searches for the root. Therefore this implementation tries to guess where
   * the vehicle might be in curvilinear coordinates, based on its current
   * velocity, last known position in cuvilinear coordinates and time past since
   * last update.
   *
   * \todo This implementation can be improoved by providing a choise for the
   * rootfinding algorithm and/or implementing a robust and realtime capable
   * solution as described in
   * [Xu2008](http://ieeexplore.ieee.org/document/4722219/).
   *
   * \todo improve the initial guess for the rootfinding algortihm - delta_t can
   * propably omitted.
   *
   * [WangKearney2002]:http://homepage.divms.uiowa.edu/~kearney/pubs/CurvesAndSufacesClosestPoint.pdf
   * [boost::math::tools::roots]:http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/roots.html
   * [Newton-Raphson]:https://en.wikipedia.org/wiki/Newton%27s_method
   * [Xu2008]:http://ieeexplore.ieee.org/document/4722219/
   * [boost::math::tools::roots::newton_raphson_iterate]:http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/roots/roots_deriv.html
   * [boost::geometry]:http://www.boost.org/doc/libs/1_65_1/libs/geometry/doc/html/geometry/reference/algorithms/distance/distance_2.html
   */
  std::pair<double, double> localize(const point position,
                                     const double velocity,
                                     const double heading, const double delta_t,
                                     const double last_known_s) const;
  /*! \brief 2D-cubic-B-spline parametrized by \f$s\f$.
   * \see Q, inv_A, [Wang2002]
   * \param s where \f$s \in [0,L]\f$ and \f$L\f$ is the total length of the
   * baseframe in [m].
   * \throws std::domain_error if \f$s \notin [0,L] \f$.
   * \returns A point on the Baseframe in Carthesian coordinates [m].
   *
   * A parametric representation of a cubic spline curve.
   *
   * \f[ P(s) = Q(A^{-1}) = (x(A^{-1}),y(A^{-1}))\f]
   *
   * where \f$s \in [0,L]\f$ and \f$L\f$ is the total length of the baseframe in
   * [m].
   *
   * [Wang2002]:http://homepage.divms.uiowa.edu/~kearney/pubs/CurvesAndSurfacesArcLength.pdf
   */
  point P(const double s) const { return Q(inv_A(s)); }

  /*!
   * \brief 2D-cubic-B-spline parametrized by t.
   * \see P, [B-Spline]
   *
   * \param t where \f$t \in [0,N]\f$ and \f$N\f$ is the number of
   * waypoints.
   * \returns A point on the Baseframe in Carthesian coordinates [m].
   * \throws std::domain_error if \f$t \notin [0,N] \f$.
   *
   * A parametric representation of a cubic spline curve.
   *
   * \f[ Q(t) = (x(t), y(t)) \f]
   *
   * where \f$t \in [0,N]\f$ and \f$N\f$ is the number of waypoints.
   *
   * [B-spline]: https://en.wikipedia.org/wiki/B-spline
   */
  point Q(const double t) const;

  /*!
   * \brief The first derivative of \f$Q(t)\f$ w.r.t. \f$t\f$.
   * \see P_prime, [boost::math::cubic_b_spline]
   * \param t where \f$t \in [0,N]\f$ and \f$N\f$ is the number of
   * waypoints.
   * \returns the first derivatives of x(t) and y(t) w.r.t. t.
   * \throws std::domain_error if \f$t \notin [0,N] \f$.
   *
   * \f[\frac{\mathrm{d}}{\mathrm{d}t}Q(t)=\left(\frac{\mathrm{d}}{\mathrm{d}t}x(t),\frac{\mathrm{d}}{\mathrm{d}t}y(t)\right)
   * = (x'(t), y'(t)) \f]
   *
   * where \f$t \in [0,N]\f$ and \f$N\f$ is the number of waypoints. Thanks to
   * [boost::math::cubic_b_spline] the derivates can be directly returned from
   * the underlying datastructure. This is the tangent of \f$Q(t)\f$ at \f$t\f$.
   *
   * [boost::math::cubic_b_spline]:http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/interpolate/cubic_b.html
   */
  point Q_prime(const double t) const;

  /*!
   * \brief The second derivative of \f$Q(t)\f$ w.r.t. \f$t\f$.
   * \see P_pprime, curvature
   *
   * \param t where \f$t \in [0,N]\f$ and \f$N\f$ is the number of
   * waypoints.
   * \returns the second derivatives of x(t) and y(t) w.r.t. t.
   * \throws std::domain_error if \f$t \notin [0,N] \f$.
   *
   * \f[\frac{\mathrm{d}^2}{\mathrm{d}t^2}Q(t)=\left(\frac{\mathrm{d}^2}{\mathrm{d}t^2}x(t),\frac{\mathrm{d}^2}{\mathrm{d}^2t}y(t)\right)
   * = (x''(t), y''(t))\f]
   * where \f$t \in [0,N]\f$ and \f$N\f$ is the number of waypoints. The
   * second derivate is computed by using the three point interpolation
   * implemented in planner::derivative.
   */
  point Q_pprime(const double t) const;

  /*!  \brief The first derivative of \f$P(s)\f$ w.r.t. \f$s\f$.
   * \see Q_prime, inv_A
   * \param s where \f$s \in [0,L]\f$ and \f$L\f$ is the total length of the
   * baseframe in [m].
   * \throws std::domain_error if \f$s \notin [0,L] \f$.
   * \returns the first derivate of \f$P(s)\f$ w.r.t s.
   *
   * The derivate can be easily calculated by transforming the parameter domain
   * to \f$t\f$.
   *
   * \f[\frac{\mathrm{d}}{\mathrm{d}s}P(s)=
   * \frac{\mathrm{d}}{\mathrm{d}t}Q(A^{-1}(s))\f]
   *
   * where \f$s \in [0,L]\f$ and \f$L\f$ is the total length of the baseframe in
   * [m]. This is the tangent of the baseframe at arc-length \f$s\f$.
   */
  point P_prime(const double s) const { return Q_prime(inv_A(s)); }

  /*!  \brief The second derivative of \f$P(s)\f$ w.r.t. \f$s\f$.
   * \see Q_pprime, inv_A
   * \param s where \f$s \in [0,L]\f$ and \f$L\f$ is the total length of the
   * baseframe in [m].
   * \throws std::domain_error if \f$s \notin [0,L] \f$.
   * \returns the second derivate of \f$P(s)\f$ w.r.t s.
   *
   * The derivate can be easily calculated by transforming the parameter domain
   * to \f$t\f$.
   *
   * \f[\frac{\mathrm{d}^2}{\mathrm{d}s^2}P(s)=
   * \frac{\mathrm{d}^2}{\mathrm{d}t^2}Q(A^{-1}(s))\f]
   *
   * where \f$s \in [0,L]\f$ and \f$L\f$ is the total length of the baseframe in
   * [m].
   */
  point P_pprime(const double s) const { return Q_pprime(inv_A(s)); }

  /*! \brief calculates the tangent angle in [rad] of baseframe at arc-length
   * \f$s\f$.
   * \param s where \f$s \in [0,L]\f$ and \f$L\f$ is the length of the baseframe
   * in [m].
   * \see P_prime, theta, theta_Q, inv_A
   * \returns the tangent angle in [rad] of baseframe at arc-length \f$s\f$.
   * \throws std::domain_error if \f$s \notin [0,L] \f$.
   *
   * \f[\tan(\theta) = \frac{y'}{x'} \rightarrow \theta =
   * \mathrm{atan2}(y',x')\f]
   *
   * this can be derived from \f$x' = \cos(\theta)\f$ and \f$y' =
   * \sin(\theta)\f$. The use of atan2 gurantees the right sign of the angle.
   *
   * In this implementation the parameter domain is transformed to \f$t\f$ and
   * the tangent angle is calculated from \f$Q(t)\f$.
   */
  double theta_P(const double s) const { return theta_Q(inv_A(s)); }

  /*! \brief calculates the tangent angle in [rad] of baseframe at \f$t\f$.
   * \param t where \f$t \in [0,N]\f$ and \f$N\f$ is the number of waypoints.
   * \see Q_prime, Q_pprime
   * \returns the tangent angle in [rad] of baseframe at arc-length \f$s\f$.
   * \throws std::domain_error if \f$t \notin [0,N] \f$.
   *
   * \f[\tan(\theta) = \frac{y'}{x'} \rightarrow \theta =
   * \mathrm{atan2}(y',x')\f]
   *
   * this can be derived from \f$x' = \cos(\theta)\f$ and \f$y' =
   * \sin(\theta)\f$. The use of atan2 gurantees the right sign of the angle.
   */
  double theta_Q(const double t) const;

 private:
  void init(std::vector<double> vx, std::vector<double> vy);
  /*! \brief \f$ A: t \in [0,N] \rightarrow  s \in [0,L] \f$.
   * \see [Wang2002], [boost::math::quadrature::trapezoidal], inv_A.
   * \param t where \f$t \in [0,N]\f$ and \f$N\f$ is the number of waypoints.
   * \returns s where \f$s \in [0,L]\f$ and \f$L\f$ is the length of the
   * baseframe.
   * \throws std::domain_error if \f$t \notin [0,N] \f$
   *
   * \f[
   * s = A(t) = \int^t_{t_0} \sqrt{((x'(t))^2 + (y'(t))^2)}\mathrm{d}t
   * \f]
   *
   * The transformation function \f$ A: t \in [0,N] \rightarrow  s \in [0,L] \f$
   * maps from the parameter domain of waypoints to the parameter domain of the
   * arc-length of the Baseframe. The complexity of this function is
   * \f$\mathcal{O}(n\log(n))\f$ in computation time. This is archieved by
   * saving previously calculated results in a static lookup table
   * planner::_A_lut to avoid unnecessary computation. The numerical
   * integration is implemented using [boost::math::quadrature::trapezoidal].
   *
   * [Wang2002]:http://homepage.divms.uiowa.edu/~kearney/pubs/CurvesAndSurfacesArcLength.pdf
   * [boost::math::quadrature::trapezoidal]:http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/quadrature/trapezoidal.html
   */
  double A(const double t) const;

  /*! \brief \f$ A^{-1}: s \in [0,L] \rightarrow  t \in [0,N] \f$.
   * \see A, [Wang2002], [boost::math::tools::roots], [bisection method].
   * \param s where \f$s \in [0,L]\f$ and \f$L\f$ is the length of the
   * baseframe.
   * \returns t where \f$t \in [0,N]\f$ and \f$N\f$ is the number of waypoints.
   * \throws std::domain_error if \f$t \notin [0,N] \f$.
   *
   * \f[ t = A^{-1}(s) \f]
   *
   * The transformation function \f$ A^{-1}: s \in [0,L] \rightarrow t \in [0,N]
   * \f$
   * maps from the parameter domain of the arc-length of the Baseframe to the
   * parameter domain of waypoints. This is the inverse function of A.
   * The complexity of this function is \f$\mathcal{O}(n\log(n))\f$ in
   * computation time. This is archieved by
   * saving previously calculated results in a static lookup table
   * planner::_A_inv_lut to avoid unnecessary computation.
   * Since an analytical solution of the inverse function is nor needed or
   * feasible, a numerical root finding algorithm from
   * [boost::math::tools::roots] is used. Because of its robustness, the
   * [bisection method] was choosen.
   *
   * [Wang2002]:http://homepage.divms.uiowa.edu/~kearney/pubs/CurvesAndSurfacesArcLength.pdf
   * [boost::math::tools::roots]:http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/roots/roots_noderiv.html
   * [bisection method]:https://en.wikipedia.org/wiki/Bisection_method
   */
  double inv_A(const double s) const;

  /*!
   * \brief \f$x(t)\f$ where \f$t \in [0,N]\f$
   * \see Q, [boost::math::cubic_b_spline]
   *
   * \f$x(t)\f$ in function \f$Q(t) = (x(t), y(t))\f$
   * [boost::math::cubic_b_spline]:http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/interpolate/cubic_b.html.
   */
  cubic_b_spline<double> _x;
  /*!
   * \brief \f$y(t)\f$ where \f$t \in [0,N]\f$
   * \see Q, [boost::math::cubic_b_spline]
   *
   * \f$y(t)\f$ in function \f$Q(t) = (x(t), y(t))\f$
   * [boost::math::cubic_b_spline]:http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/interpolate/cubic_b.html.
   */
  cubic_b_spline<double> _y;

  /*!
   * \brief number of waypoints \f$N\f$
   */
  int _N;
};
}  // namespace planner

#endif  // PLANNER_BASEFRAME_H_
