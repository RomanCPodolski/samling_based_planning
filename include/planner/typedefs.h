#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <vector>

namespace planner {
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::d2::point_xy<double> point;
typedef bg::model::polygon<point> polygon;

typedef std::vector<double> state_type;  ///< state used in Maneuver::path
typedef bg::model::box<point> box;       ///< aabb in Maneuver::collision
typedef std::pair<box, std::shared_ptr<polygon>> value;  ///< see Rtree
typedef bgi::rtree<value, bgi::quadratic<16>> Rtree;     ///< Spatial index

}  // namespace planner
#endif /* TYPEDEFS_H */
