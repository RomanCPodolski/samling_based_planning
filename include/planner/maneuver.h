/*
 * maneuver.h
 * Copyright (C) 2017 roman <maito:roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef MANEUVER_H
#define MANEUVER_H

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <algorithm>
#include <boost/geometry.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <exception>
#include <memory>
#include "planner/baseframe.h"
#include "planner/pose.h"
#include "planner/typedefs.h"

DECLARE_double(min_manouver_length);
DECLARE_double(manouver_speed_gain);
DECLARE_double(granularity);
DECLARE_double(collision_standart_deviation);
DECLARE_double(max_curvature);
DECLARE_double(weight_safety);
DECLARE_double(weight_smoothness);
DECLARE_double(weight_consistency);

namespace planner {

using boost::math::sign;
using boost::math::quadrature::trapezoidal;
using boost::geometry::distance;
using boost::math::normal;

/*! \class Maneuver
 *  \brief A Maneuver describes a path in curvilinear- and carthesian
 * coordinates.
 *  \see Planner, Baseframe, [Chu2012]
 *  \author Roman C. Podolski - <roman.podolski@tum.de>
 *
 *  A Maneuver is a cubic-polynomnial in curvilinear coordinates of form
 *
 *  \f[
 *  q(s) = a(s - s_i)^3 + b(s - s_i)^2 + c(s - s_i) + d
 *  \f]
 *
 *  that satisfies the following boundary conditions
 *
 *  \f{eqnarray*}
 *  q(s_i) &=& q_i  \\
 *  q(s_f) &=& q_f  \\
 *  \frac{\mathrm{d}}{\mathrm{d}s}q(s_i) &=& \tan(\Delta \theta)  \\
 *  \frac{\mathrm{d}}{\mathrm{d}s}q(s_i) &=& 0  \\
 *  \f}
 *
 *  with the following constants:
 *  * \f$s_i \in [0,L] \f$ is the initial arc-length along the baseframe [m].
 *  * \f$s_f \in [0,L] \wedge s_f > s_i \f$ is the final arc-length [m].
 *  * \f$q_i \in [-q_\mathrm{max}, q_\mathrm{max}]\f$ initial offset from
 * Baseframe
 * [m].
 *  * \f$q_f \in [-q_\mathrm{max}, q_\mathrm{max}]\f$ final offset from
 * Baseframe
 * [m].
 *  * \f$ \Delta \theta \f$ difference between Baseframe::theta and vehicle
 * heading [rad].
 *
 *  By forward integration of a suitable vehicle model, a discrete path in
 *  carthesian coordinates can be calculated from a maneuver.
 *
 *  [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
 */
class Maneuver {
 public:
  /*!
   * Constructs a Manouver by calculating the coefficients of the polynomnial.
   *
   * \see Baseframe, [Eigen::Dense], [boost::geometry::index::rtree], [Chu2012]
   *
   * \param q_i initial offset \f$ q_i [m] \f$.
   * \param q_f final offset \f$ q_f [m] \f$.
   * \param s_i inital arc-length \f$ s_i [m] \f$.
   * \param v vehicle speed \f$ v [m/s] \f$.
   * \param theta difference between vehilce heading and tangent \f$\theta\f$
   * [rad].
   * \param baseframe observer pointer to Baseframe of this maneuver.
   * \param rtree observer pointer to Spatial index used for collision
   * detection, see [boost::geometry::index::rtree].
   * \param obstacles observer pointer to obstacles used for collision
   * detection.
   *
   * The Constructor calculates the coefficients for the cubic-polynomnial
   *
   *  \f[ q(s) = a(s - s_i)^3 + b(s - s_i)^2 + c(s - s_i) + d \f]
   *
   *  that satisfies the following boundary conditions
   *
   *  \f{eqnarray*}
   *  q(s_i) &=& q_i  \\
   *  q(s_f) &=& q_f  \\
   *  \frac{\mathrm{d}}{\mathrm{d}s}q(s_i) &=& \tan(\Delta \theta)  \\
   *  \frac{\mathrm{d}}{\mathrm{d}s}q(s_i) &=& 0  \\
   *  \f}
   *
   *  with the following constants:
   *  * \f$s_i \in [0,L] \f$ is the initial arc-length along the baseframe [m].
   *  * \f$s_f \in [0,L] \wedge s_f > s_i \f$ is the final arc-length [m].
   *  * \f$q_i \in [-q_\mathrm{max}, q_\mathrm{max}]\f$ initial offset from
   * Baseframe [m].
   *  * \f$q_f \in [-q_\mathrm{max}, q_\mathrm{max}]\f$ final offset from
   * Baseframe [m].
   *  * \f$ \Delta \theta \f$ difference between Baseframe::theta and vehicle
   * heading [rad].
   *
   *
   * Finding the cofficients that satisfy the boundary conditions is simply a
   * matther of solving the linear system where \f$ \Delta s = s_f - s_i \f$.
   *
   * \f[
   * \begin{bmatrix}
   * 0 & 0 & 0 & 1 \\
   * \Delta s^3 & \Delta s^2 & \Delta s & 1 \\
   * 0 & 0 & 1 & 0 \\
   * 3\Delta s^2 & 2\Delta s & 1 & 0 \\
   * \end{bmatrix} \begin{bmatrix} a \\ b \\ c \\ d \end{bmatrix} =
   * \begin{bmatrix} q_i \\ q_f \\ \tan(\Delta \theta) \\ 0.0 \end{bmatrix}
   * \f]
   *
   * Because \f$s_i < s_f\f$ and therefore \f$ \Delta s > 0 \f$ this matrix
   * positive definte and therefore can be inverted. The implementation uses
   * [Eigen::Dense].
   *
   * \f[ \bm{A}\bm{b}= \bm{c} \rightarrow \bm{c} = \bm{A}^{-1}\bm{c} \f]
   *
   * The Baseframe arc-length at the end of the maneuver \f$ s_f \f$ is w.r.t.
   * the vehicle speed \f$v\f$, maneuver length speed gain \f$k_v\f$, and
   * minimum maneuver length \f$\Delta s_\mathrm{min}\f$.
   *
   * \f[ s_f = s_i + \underbrace{k_v v + \Delta s_\mathrm{min}}_{\Delta s} \f]
   *
   * The design parameters \f$\Delta s_\mathrm{min}\f$ and \f$k_v\f$ are
   * implemented as command line flags using [gflags].
   *
   * To set the maneuver speed gain \f$k_v\f$ from the command line, use the
   * command line flag `--maneuver_speed_gain=value`, where value is a floating
   * point number of double precision. Default is 1.0.
   *
   * To set the minimum maneuver length \f$\Delta s_\mathrm{min} [m]\f$ from the
   * command line, use the command line flag `--min_maneuver_length=value`,
   * where value is a floating point number of double precision. Default is
   * 10.0.
   *
   * [boost::geometry::index::rtree]:http://www.boost.org/doc/libs/1_65_1/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html
   * [Eigen::Dense]:https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   * [gflags]:gflags.github.io/gflags/
   */
  Maneuver(const point position, const double heading, const double velocity,
           const double q_i, const double q_f, const double s_i,
           std::shared_ptr<Baseframe> const &baseframe,
           std::shared_ptr<Rtree> const &rtree);

  /*! \brief Curvature of the Maneuver at arc-length of Baseframe s.
   * \see Baseframe::curvature, S, Q, q, dq, ddq, [Curvature], [Werling2010],
   * [Barfoot2004]
   * \param s \f$ \in [0,L] \f$ where \f$L\f$ is the total length of the
   * Basefame.
   * \throws std::domain_error if \f$s \notin [0,L]\f$
   *
   * The curvature of the continous path at Baseframe arc-length \f$s\f$ can be
   * computed as
   *
   * \f[
   * \kappa =
   * \frac{S}{Q}\left(\frac{(1-q\kappa_b)\frac{\mathrm{d}^2q}{\mathrm{d}s^2}+\kappa_s\left(\frac{\mathrm{d}q}{\mathrm{d}s}\right)^2}{Q^2}\right)
   * \f]
   *
   * for a full derivation of this, please refere to [Barfoot2004] and
   * [Werling2010].
   *
   * [Curvature]:https://en.wikipedia.org/wiki/Curvature
   * [Werling2010]:https://www.researchgate.net/publication/224156269_Optimal_Trajectory_Generation_for_Dynamic_Street_Scenarios_in_a_Frenet_Frame
   * [Barfoot2004]:http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.362.8916&rep=rep1&type=pdf
   */
  double curvature(const double s) const;

  /*! \brief getter for the final arc-length along the Baseframe \f$ s_f \f$.
   */
  double s_f() const noexcept { return _s_f; }

  /*! \brief getter for the final arc-length along the Baseframe \f$ s_f \f$.
   * \todo chech if this makes sense, rather check that the path is only
   * generated until the end of the baseframe
   */
  void s_f(const double s) {
    _s_f = s > _baseframe->length() ? _baseframe->length() : s;
  }

  /*! \brief getter for the initial arc-length along the Baseframe \f$ s_i \f$
   * [m].
   */
  double s_i() const noexcept { return _s_i; }

  /*! \brief setter for the initial arc-length along the Baseframe \f$ s_i \f$
   * [m].
   * \todo rather throw an exception than adjust s
   */
  void s_i(const double s) {
    _s_i = s > _baseframe->length() ? _baseframe->length() : s;
  }

  /*! \brief getter for the initial offset from the Baseframe \f$ q_i \f$ [m].
   * this is equal to _d
   */
  double q_i() const noexcept { return _d; }

  /*! \brief getter for the final offset from the Baseframe \f$ q_f \f$ [m].
   */
  double q_f() const noexcept { return _q_f; }

  /*! \brief getter difference between vehicle heading and Baseframe::theta \f$
   * \Delta \theta \f$ [rad].
   */
  double theta() const noexcept { return _theta; }

  /*! \brief indicates if this maneuver is driveable for a non-holonomic vehicle
   *
   * This is true if, and only if
   * \f[\left\{q(s) < \kappa_b(s)^{-1} \wedge |\kappa| \leq
   * \kappa_\mathrm{max} \; | \; \forall s \in [s_i, s_f]\right\}\f]
   */
  bool driveable() {
    if (_path.empty()) path();
    return _drivable;
  }

  /*! \brief Collision checking function \f$ \Phi: m \in \mathcal{M} \rightarrow
   * s \in
   * [0,L] \f$
   * \see C_s, J, _rtree, Planner::_rtree, Planner::update
   * [LaValleChap5.],[boost::geometry],[boost::geometry::index::rtree],[boost::geomerty::strategy::transform]
   * \returns length at which collsion was detected, 0.0 if no collison appears
   * [m].
   *
   * This performs the collision check for all discrete poses in this maneuver.
   * Since naive collision checking involves checking all distances to all
   * stored obstacles it is not feasible from a performance point of view.
   * This implements a more sophisticated approach using both, axis aligned
   * bounding boxes and object orientated bounding boxes [LaValleChap5.].
   *
   * The load is reduced by using a Rtree [boost::geometry::index::rtree] as
   * spatial index for efficient searching. The index itself is populated in
   * Planner::update and availabe in Manouver via Maneuver::_rtree. For every
   * pose in a discrete path, first the object aligned bounding box is computed
   * using [boost::geometry] and [boost::geomerty::strategy::transform], from
   * which the axis aligned bounding box can be obtaind via it's envelope.
   *
   * This axis aligned bounding box is used to search the rtree for overlapping
   * bounding boxes (that indicate the risk of collision). For all risks an
   * actual collison check by checking if the obstacle aligned bounding boxes
   * overlap, can be performed. If a collision is detected, the lenght \f$s\f$
   * along the Baseframe at which the collsion would occur is returned.
   *
   * [LaValleChap5.]:http://planning.cs.uiuc.edu/ch5.pdf
   * [boost::geometry]:http://www.boost.org/doc/libs/1_65_1/libs/geometry/doc/html/geometry/quickstart.html
   * [boost::geometry::index::rtree]:http://www.boost.org/doc/libs/1_65_1/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html
   * [boost::geomerty::strategy::transform]:http://www.boost.org/doc/libs/1_65_1/libs/geometry/doc/html/geometry/reference/algorithms/transform/transform_3_with_strategy.html
   */
  double collision();

  /*! \brief linear combination of cost functions with weighting factors.
   * \see C_s, C_k, C_c,[Chu2012]
   * \returns path cost for this maneuver.
   * \throws std::domain_error see Maneuver::C_s
   *
   * The overall cost of a maneuver is a linear combination with weigting
   * factors.
   * The weigths function as design parameters and influence the driving style
   * of the vehicle.
   * They can be set from the commandline using the following [gflags].
   * * `--weight_saefty=value` sets the weight \f$w_s\f$ of the safety cost
   * \f$C_s\f$.
   * * `--weight_smoothness=value` sets the weight \f$w_s\f$ of the smoothness
   * cost \f$C_\kappa\f$.
   * * `--weight_consistency=value` sets the weight \f$w_s\f$ of the consistency
   * cost \f$C_c\f$.
   *
   *  \f[ J = w_sC_s + w_\kappa C_\kappa + w_cC_c \f]
   *
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   * [gflags]:gflags.github.io/gflags/
   */
  double J(const Maneuver &previous,
           const std::shared_ptr<std::vector<Maneuver>> maneuvers,
           double w_s = FLAGS_weight_safety,
           double w_k = FLAGS_weight_smoothness,
           double w_c = FLAGS_weight_consistency) const;

  /*!
   * \brief Safety cost function \f$ C_s \f$.
   * \see J, [boost::math::normal_distribution],[Chu2012],[LaValleChap5.]
   * \param maneuver all maneuver that are considered in this planning cycle.
   * \param sigma standard deviation \f$ \sigma \f$ of the collision risk normal
   * distribution. This can be set from the command line using the [gflag]
   * `--collision_standart_deviation=value`. Default is 1.0.
   * \return the saefty cost for this maneuver.
   *
   * The safety-cost of a maneuver is equal to it's risk of collison. This can
   * be calculated by
   *
   * \f[
   * \sum_{i = 1}^N \Phi(m_i)g(q_{i,f};\mu = q_f, \sigma^2)
   * \f]
   *
   * where:
   * * \f$\Phi: m \rightarrow \{0,1\} \f$ where \f$ m \in \mathcal{M} \f$ and
   * \f$
   * \mathcal{M} \f$ is the set of all maneuvers in this planning cycle. This
   * function performs the collision check. See [LaValleChap5.],
   * Maneuver::collision.
   * * \f$g(x; \mu, \sigma) \f$ is the propability density function of the
   * normal distribution \f$\mathcal{N}(\mu,\sigma)\f$.
   * * \f$m_i \in \mathcal{M}\f$ is a maneuver of this planning cycle and \f$
   * \mathcal{M} \f$ is the set of all maneuvers in this planning cycle.
   * * \f$ q_{i,f} \f$ is the final offset from the Baseframe of the maneuver
   * \f$m_i\f$.
   * * \f$ q_f\f$ is the final offset from the baseframe of this maneuver.
   *
   * This implementation uses [boost::math::normal_distribution].
   *
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   * [boost::math::normal_distribution]:http://www.boost.org/doc/libs/1_65_1/libs/math/doc/html/math_toolkit/dist_ref/dists/normal_dist.html
   * [gflag]:gflags.github.io/gflags/
   * [LaValleChap5.]:http://planning.cs.uiuc.edu/ch5.pdf
   */
  double C_s(std::shared_ptr<std::vector<Maneuver>> maneuvers,
             const double sigma = FLAGS_collision_standart_deviation) const;

  /*! \brief Smoothness cost function \f$ C_\kappa \f$.
   * \see J, [boost::math::trapezoidal],[Chu2012]
   * \return the smoothness cost for this maneuver.
   *
   * The smoothess of the path can be described as the integral of it's
   * curvature over the length of the path. Since the true length of the path is
   * unknown, and not needed, the integration can be done by integrating ofer
   * the lenght of the baseframe that the maneuver covers and projecting the
   * curvature onto the baseframe using Q.
   *
   * \f[
   * C_\kappa = \int \kappa^2(s)\; \mathrm{d}s_p = \int \kappa^2(s)Q(s)\;
   * \mathrm{d}s
   * \f]
   *
   * The numerical integration is implemented using [boost::math::trapezoidal].
   *
   * [boost::math::trapezoidal]:http://www.boost.org/doc/libs/1_65_0/libs/math/doc/html/math_toolkit/quadrature/trapezoidal.html
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   */
  double C_k() const;

  /*! \brief Consistency cost function \f$C_c\f$.
   * \see J, [boost::math::trapezoidal],[Chu2012]
   *
   * \param previous reference to the previous maneuver, choosen by the Planner.
   * \return the consistency cost for this maneuver.
   * \throws std::domain_error if \f$ s_1 \geq s_2 \vee s_1 \notin [0,L] \vee
   * s_2
   * \notin [0,L] \f$
   *
   * The consistency cost of a path influences the driving dynamic by
   * calculating the cost for "change'' that a maneuver introduces w.r.t. to
   * the last maneuver that was choosen by the planner. The cost can be
   * calculated by integrating the distance \f$|q - q_p|\f$ for the overlapping
   * sections of the maneuvers.
   *
   * \f[
   * C_c = \frac{1}{s_2 - s_1} \int^{s_2}_{s_1} |q(s) - q_p(s)| \;\mathrm{d}s
   * \f]
   *
   * where
   * * \f$s_1\f$ is the Baseframe arc-length at the start of the overlapping
   * section
   * * \f$s_2\f$ is the Baseframe arc-length at the end of the overlapping
   * section
   * * \f$q: s \in [0, L] \rightarrow q \in [-q_\mathrm{max},q_\mathrm{max}]\f$
   * is the offset function q for this maneuver.
   * * \f$q_p: s \in [0, L] \rightarrow q_p \in
   * [-q_\mathrm{max},q_\mathrm{max}]\f$
   * is the offset function q for the previous maneuver.
   *
   * The numerical integration is implemented using [boost::math::trapezoidal].
   *
   * Note: the definition in [Chu2012] of this function uses the euclidiean
   * distance in the Carthesian coordinatesystem to calculate the safetycost,
   * but since the transformation between the curvilinear and the Carthesian
   * coordinatesystem is linear (see path and gen_path) the distance in
   * curvilinear and Carthesian coordinates express the same relation. The
   * implementation introduced here is far simpler, more performant and exact,
   * since this implementation allows the use of trapezoidial quadrature while
   * with the discrete path explicit euler would have been used.
   *
   * [boost::math::trapezoidal]:http://www.boost.org/doc/libs/1_65_0/libs/math/doc/html/math_toolkit/quadrature/trapezoidal.html
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   */
  double C_c(const Maneuver &previous) const;

  /*! \brief calculates a discrete path in carthesian coordinates as a set of
   * vehicle poses.
   * \see [Chu2012], [boost::numeric::odeint]
   * \returns a discrete path as a set of vehilce poses.
   *
   * Vehicle motion in Carthesian-coordinates can be described with the
   * following set of non-linear differential equations.
   * \f[
   * \dot{x}=v\cos\theta, \qquad
   * \dot{y}=v\sin\theta, \qquad
   * \dot{\theta}=v\kappa
   * \f]
   *
   * where:
   * * \f$v\f$ is the forward velocity of the vehicle [m/s].
   * * \f$\theta\f$ is the heading of the vehicle [rad].
   * * \f$\kappa\f$ is the curvature of the path.
   *
   * To describe this set of ODEs not w.r.t. time \f$t\f$ but with respect to
   * arc-length along Baseframe \f$s\f$, the following substitution is made.
   *
   * \f[
   * v = SQ\frac{\mathrm{d}s}{\mathrm{d}t}
   * \f]
   *
   * This yieds
   *
   * \f[
   * \frac{\mathrm{d}x}{\mathrm{d}s}=Q\cos\theta, \qquad
   * \frac{\mathrm{d}y}{\mathrm{d}s}=Q\sin\theta, \qquad
   * \frac{\mathrm{d}\theta}{\mathrm{d}s}=Q\kappa
   * \f]
   *
   * Now any numerical integration stepper, as defined in
   * [boost::numeric::odeint], can be used to solve this set of ODEs by forward
   * integration. The discrete step size of the numerical integrator is a design
   * parameter and can be set from the command line using the [gflag]
   * `--path-granularity=value`.
   * Default is 1.0[m]
   *
   * During path generation the relation \f$\left\{q(s) < \kappa_b(s)^{-1}
   * \wedge |\kappa| \leq
   * \kappa_\mathrm{max} \; | \; \forall s \in [s_i, s_f]\right\}\f$ is checked
   * for violations at the discrete poses. A path that violates this relation,
   * is marked as not-driveable.
   *
   * The path generation is implemented as a lazy evaluation. It is only done
   * once, afterwards this function becomes just a getter for _path.
   *
   * \todo The vehilce model used to describe vehilce motion is a very simple
   * one. It is suitable for proofing the concept of this planner, but may not
   * be suitable for high performance racing. To improve this planner a set of
   * ODEs that model the movement of the SAE-Vehicle in a more accurate way
   * would be needed. This could be the topic of a final-years project,
   * Bachelors- or Masters-Thesis in mechanical engineering. The model can be
   * implemented in this function (or better - improove the overall design by
   * moving this to it's own class `PathGenerator`).
   *
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   * [boost::numeric::odeint]:http://www.boost.org/doc/libs/1_65_1/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/integrate_functions.html
   * [gflag]:gflags.github.io/gflags/
   */
  std::vector<pose> path();
  std::vector<pose> path() const { return _path; };

 private:
  /*!
   * \brief \f[
   * Q = \sqrt{\left(\frac{\mathrm{d}q}{\mathrm{d}s}\right)^2+(1-q\kappa_b)^2}
   * \f]
   *
   * \see curvature, path, C_k,[Chu2012]
   * \throws std::domain_error if \f$ s \notin [0,L] \f$
   *
   * helper function.
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   */
  double Q(const double s) const;
  /*!  \f[
   * S = \mathrm{sgn}(1 - q\kappa_b)
   * \f]
   * \see curvature,[Chu2012],[boost::math::sign]
   * \throws std::domain_error if \f$ s \notin [0,L] \f$
   *
   * helper function.
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   * [boost::math::sign]:http://www.boost.org/doc/libs/1_65_1/boost/math/special_functions/sign.hpp
   */
  double S(const double s) const;

  /*! \brief Calculate the offset from the Baseframe at s \f$q: s \in
   [s_i,\infty] \rightarrow q \in [q_i, q_f]\f$
   * \see _a, _b,_c, _d, dq, ddq, curvature, Maneuver, [Chu2012]
   * \throws std::domain_error if \f$ s < s_i \f$.
   * \returns the offset \f$q\f$ from the Baseframe at arc-length \f$s\f$.
   *
   * A maneuver in curvilinear coordinates can be described as a cubic
   polynomnial
   *
   * \f[
        q(s) = \begin{cases} a(s - s_i)^3 + b(s - s_i)^2 + c(s - s_i) + d,
   &\text{ if } s_i \leq s < s_f\\
                q_f, &\text{ if } s_f \leq s
        \end{cases}
   * \f]
   *
   * that satisfies the following boundary conditions
   *
   *  \f{eqnarray*}
   *  q(s_i) &=& q_i  \\
   *  q(s_f) &=& q_f  \\
   *  \frac{\mathrm{d}}{\mathrm{d}s}q(s_i) &=& \tan(\Delta \theta)  \\
   *  \frac{\mathrm{d}}{\mathrm{d}s}q(s_i) &=& 0  \\
   *  \f}
   *
   *  with the following constants:
   *  * \f$s_i \in [0,L] \f$ is the initial arc-length along the baseframe [m].
   *  * \f$s_f \in [0,L] \wedge s_f > s_i \f$ is the final arc-length [m].
   *  * \f$q_i \in [-q_\mathrm{max}, q_\mathrm{max}]\f$ initial offset from
   * Baseframe [m].
   *  * \f$q_f \in [-q_\mathrm{max}, q_\mathrm{max}]\f$ final offset from
   * Baseframe [m].
   *  * \f$ \Delta \theta \f$ difference between Baseframe::theta and vehicle
   * heading [rad].
   *
   * See the constructor Maneuver::Maneuver for details on calculating the
   coefficients.
   *
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   */
  double q(double s) const;
  /*!
   * \see _a, _b,_c, _d, dq, ddq, curvature, Maneuver, [Chu2012]
   * \throws std::domain_error if \f$ s < s_i \f$.
   * \returns the first derivate of Maneuver::q
   *
   * Can be analytically computed using
   *
   *
   * \f[
        \frac{\mathrm{d}}{\mathrm{d}s}q(s) = \begin{cases}
                3a(s - s_i)^2 + 2b(s - s_i) + c, &\text{ if } s_i \leq s < s_f\\
                0, &\text{ if } s_f \leq s
        \end{cases}
   * \f]
   *
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   */
  double dq(double s) const;
  /*!
   * \see _a, _b,_c, _d, dq, ddq, curvature, Maneuver, [Chu2012]
   * \throws std::domain_error if \f$ s < s_i \f$.
   * \returns the second derivate of Maneuver::q
   *
   * Can be analytically computed using
   *
   * \f[
   \frac{\mathrm{d^2}}{\mathrm{d}s^2}q(s) = \begin{cases}
           6a(s - s_i) + 2b, &\text{ if } s_i \leq s < s_f\\
           0, &\text{ if } s_f \leq s
   \end{cases}
   * \f]
   *
   * [Chu2012]: http://ieeexplore.ieee.org/document/6203588/
   */
  double ddq(double s) const;

  double _q_f;  ///< final offset from Baseframe \f$q_f\f$ [m].
  double _s_i;  ///< initial arc-length along Baseframe \f$s_i\f$ [m].
  double _s_f;  ///< final arc-length along Baseframe \f$s_f\f$ [m].

  /*! \brief difference between vehicle heading and tangent angle of baseframe
   * \f$ \Delta \theta \f$
   * */
  double _theta;
  double _v;  ///< forward velocity of the vehicle \f$v\f$ [m/s]

  /*! \brief coefficient a
   * \see q, dq, ddq, Maneuver::Maneuver
   * */
  double _a;
  /*! \brief coefficient b
   * \see q, dq, ddq, Maneuver::Maneuver
   * */
  double _b;
  /*! \brief coefficient c
   * \see q, dq, ddq, Maneuver::Maneuver
   * */
  double _c;
  /*! \brief coefficient d
   * \see q, dq, ddq, Maneuver::Maneuver
   * */
  double _d;

  point _position;  ///< position of the vehicle, used for Maneuver::path
  double _heading;  ///< heading of the vehicle, used in Maneuver::path
  std::vector<pose> _path;  ///< discrete path as vector of poses

  /*!
   * \brief Baseframe for this maneuver
   * \see Baseframe, curvature, Q, S, path
   *
   * This is an observer-pointer of a Baseframe that was constructed at another
   * place - normally in the Planner.
   */
  std::shared_ptr<Baseframe> _baseframe;
  /*!
   * \brief Spatial index observer for collision checking
   * \see [boost::geometry::index::rtree], collision
   *
   * This is an observer-pointer of a Rtree that was constructed at another
   * place - normally in the Planner.
   *
   * [boost::geometry::index::rtree]:http://www.boost.org/doc/libs/1_65_1/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html
   */
  std::shared_ptr<Rtree> _rtree;

  /*! \brief flag indicates if collision was already calculated.
   * \see collision
   */
  bool _collision_checked;

  /*! \brief stores length where collision occured - 0.0 for non collision.
   * \see collision
   */
  double _collision_length;
  /*! \brief flag if this path is reasonable for a non-holonomic vehicle with
   * the maximal curvature of \f$k_\mathrm{max}\f$
   * \see driveable, path
   *
   * This is true if, and only if
   * \f[\left\{q(s) < \kappa_b(s)^{-1} \wedge |\kappa| \leq
   * \kappa_\mathrm{max} \; | \; \forall s \in [s_i, s_f]\right\}\f]
   */
  bool _drivable;

  friend std::ostream &operator<<(std::ostream &os, const Maneuver &m);
};

} /* namespace planner */

#endif /* !MANEUVER_H */
