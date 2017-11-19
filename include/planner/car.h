#ifndef CAR_OBB_H
#define CAR_OBB_H

#include <gflags/gflags.h>
#include "planner/pose.h"
#include "planner/typedefs.h"

DECLARE_double(car_length);
DECLARE_double(car_width);

namespace planner {

class Car {
 private:
  polygon m_car;  // polygon that represents the dimensons of the car.

 public:
  Car() {
    double width = FLAGS_car_width;
    double length = FLAGS_car_length;
    m_car.outer().push_back(point(-length / 2.0, -width / 2));
    m_car.outer().push_back(point(+length / 2.0, -width / 2));
    m_car.outer().push_back(point(+length / 2.0, +width / 2));
    m_car.outer().push_back(point(-length / 2.0, +width / 2));
    m_car.outer().push_back(
        point(-length / 2.0, -width / 2.0));  // close the polygon
  }

  polygon obb(const pose &p) {
    polygon oobb, tmp;  // object oriented bounding box
    // transform vehicle bounding box to match pose
    // roatate will perform a clockwise rotation around the origin, but the
    // heading is defined anti-clockwise
    bg::strategy::transform::rotate_transformer<bg::radian, double, 2, 2>
        rotate(-p.theta);
    bg::strategy::transform::translate_transformer<double, 2, 2> translate(
        p.x(), p.y());
    bg::transform(m_car, tmp, rotate);    // rotate first...
    bg::transform(tmp, oobb, translate);  // ...then translate
    return oobb;
  }

  box aabb(const pose &p) {
    return bg::return_envelope<box>(obb(p));  // axis aligned bounding box
  }

  virtual ~Car(){};
};
}  // namespace planner

#endif /* CAR_OBB_H */
