#ifndef POSE_H
#define POSE_H
#include "planner/typedefs.h"
namespace planner {

/*! \struct pose
 *  \brief A vehicle pose
 *
 *  The pose of a vehicle consits of its position and orientation in Carthesian
 * space.
 *  The pose can also store the position of the vehicle in curvinilear
 * coordinates, if known.
 */
struct pose {
  point position; /*!< position of the vehicle in cartesian coordinates */
  double theta;   /*!< heading of the vehicle*/
  double s;       /*!< arc-length s along the baseframe in [m] */
  double q;       /*!< offset from the baseframe in [m] */

  /*! \brief getter for the x-coordinate of the position in
   * Carthesian-coordinates [m].
   */
  double x() const { return position.x(); }
  /*! \brief getter for the y-coordinate of the position in
   * Carthesian-coordinates [m].
   */
  double y() const { return position.y(); }
};
}  // namespace planner

#endif /* POSE_H */
