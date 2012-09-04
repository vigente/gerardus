// Boost.Polygon library voronoi_utils.hpp header file

//          Copyright Andrii Sydorchuk 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_POLYGON_VORONOI_UTILS
#define BOOST_POLYGON_VORONOI_UTILS

#include <cmath>
#include <stack>
#include <vector>

#include "voronoi_diagram.hpp"

namespace boost {
namespace polygon {

// Bounding rectangle data structure. Contains coordinates
// of the bottom left and upper right corners of rectangle.
template <typename T>
class bounding_rectangle {
public:
  typedef T coordinate_type;

  bounding_rectangle() : x_min_(1), x_max_(0) {}

  bounding_rectangle(coordinate_type x1, coordinate_type y1,
                     coordinate_type x2, coordinate_type y2) {
    x_min_ = (std::min)(x1, x2);
    y_min_ = (std::min)(y1, y2);
    x_max_ = (std::max)(x1, x2);
    y_max_ = (std::max)(y1, y2);
  }

  // Update bounding rectangle.
  void update(coordinate_type x, coordinate_type y) {
    if (x_min_ > x_max_) {
      x_min_ = x_max_ = x;
      y_min_ = y_max_ = y;
    } else {
      x_min_ = (std::min)(x_min_, x);
      y_min_ = (std::min)(y_min_, y);
      x_max_ = (std::max)(x_max_, x);
      y_max_ = (std::max)(y_max_, y);
    }
  }

  bool is_empty() const {
    return x_min_ > x_max_;
  }

  void clear() {
    x_min_ = 1;
    x_max_ = 0;
  }

  bool contains(coordinate_type x, coordinate_type y) const {
    return x > x_min_ && x < x_max_ && y > y_min_ && y < y_max_;
  }

  // Return the x-coordinate of the bottom left point of the rectangle.
  coordinate_type x_min() const {
    return x_min_;
  }

  // Return the y-coordinate of the bottom left point of the rectangle.
  coordinate_type y_min() const {
    return y_min_;
  }

  // Return the x-coordinate of the upper right point of the rectangle.
  coordinate_type x_max() const {
    return x_max_;
  }

  // Return the y-coordinate of the upper right point of the rectangle.
  coordinate_type y_max() const {
    return y_max_;
  }

private:
  coordinate_type x_min_;
  coordinate_type y_min_;
  coordinate_type x_max_;
  coordinate_type y_max_;
};

template <typename fpt>
struct voronoi_utils_traits {
  typedef fpt coordinate_type;
  typedef detail::point_2d<coordinate_type> point_type;
  typedef std::vector<point_type> point_set_type;
  typedef bounding_rectangle<coordinate_type> brect_type;
  typedef struct {
    template <typename CT>
      coordinate_type operator()(const CT& value) const {
        return static_cast<coordinate_type>(value);
    }
  } ctype_converter_type;
};

// Voronoi output post-processing tools.
template <typename T, typename TRAITS = voronoi_utils_traits<T> >
class voronoi_utils {
public:
  typedef typename TRAITS::coordinate_type coordinate_type;
  typedef typename TRAITS::point_type point_type;
  typedef typename TRAITS::point_set_type point_set_type;
  typedef typename TRAITS::brect_type brect_type;

  // Get scaled by a factor bounding rectangle.
  template <typename CT>
  static brect_type scale(const bounding_rectangle<CT> &brect,
      coordinate_type factor) {
    coordinate_type x_len = to_fpt(brect.x_max()) - to_fpt(brect.x_min());
    coordinate_type y_len = to_fpt(brect.y_max()) - to_fpt(brect.y_min());
    coordinate_type x_mid = to_fpt(brect.x_max()) + to_fpt(brect.x_min());
    coordinate_type y_mid = to_fpt(brect.y_max()) + to_fpt(brect.y_min());
    x_mid *= to_fpt(0.5);
    y_mid *= to_fpt(0.5);
    coordinate_type offset = (std::max)(x_len, y_len) * factor * to_fpt(0.5);
    brect_type new_brect(x_mid - offset, y_mid - offset,
                         x_mid + offset, y_mid + offset);
    return new_brect;
  }

  // Discretizes finite Voronoi edge.
  // Discretization absolute error is defined by max_error value.
  template <typename CT>
  static void discretize(const voronoi_edge<CT> &edge,
      coordinate_type max_error, point_set_type &discretization) {
    // Don't process infinite edges.
    if (!edge.is_finite())
      return;

    // Retrieve the cell to the left of the current edge.
    const typename voronoi_edge<CT>::voronoi_cell_type *cell1 = edge.cell();
    // Retrieve the cell to the right of the current edge.
    const typename voronoi_edge<CT>::voronoi_cell_type *cell2 =
        edge.twin()->cell();

    discretization.push_back(get_point(edge.vertex0()->vertex()));
    discretization.push_back(get_point(edge.vertex1()->vertex()));
    if (edge.is_curved()) {
      bool flag = cell1->contains_segment();
      // point1 - site point;
      point_type point1 = flag ?
          get_point(cell2->point0()) : get_point(cell1->point0());
      // point2 - start-point of the segment site;
      point_type point2 = flag ?
          get_point(cell1->point0()) : get_point(cell2->point0());
      // point3 - endpoint of the segment site;
      point_type point3 = flag ?
          get_point(cell1->point1()) : get_point(cell2->point1());
      fill_intermediate_points(
          point1, point2, point3, max_error, discretization);
    }
  }

  // Clip a linear Voronoi edge with a given rectangle.
  template <typename CT>
  static void clip(const voronoi_edge<CT> &edge,
      const brect_type &brect, point_set_type &clipped_edge) {
    // Don't process curved edges.
    if (edge.is_curved())
      return;

    if (edge.is_finite()) {
      clip(edge.vertex0()->vertex(), edge.vertex1()->vertex(),
          brect, clipped_edge);
    } else {
      const typename voronoi_edge<CT>::voronoi_cell_type *cell1 = edge.cell();
      const typename voronoi_edge<CT>::voronoi_cell_type *cell2 =
          edge.twin()->cell();
      point_type point1 = get_point(cell1->point0());
      point_type point2 = get_point(cell2->point0());
      if (point1 == point2) {
        if (cell1->contains_segment())
          point1 = get_point(cell1->point1());
        else
          point2 = get_point(cell2->point1());
      }

      if (edge.vertex0()) {
        // Ray edge.
        point_type origin = get_point(edge.vertex0()->vertex());
        point_type direction(point1.y() - point2.y(), point2.x() - point1.x());
        if (brect.contains(origin.x(), origin.y()))
          clipped_edge.push_back(origin);
        intersect(origin, direction, RAY, brect, clipped_edge);
      } else if (edge.vertex1()) {
        // Ray edge.
        point_type origin = get_point(edge.vertex1()->vertex());
        point_type direction(point2.y() - point1.y(), point1.x() - point2.x());
        if (brect.contains(origin.x(), origin.y()))
          clipped_edge.push_back(origin);
        intersect(origin, direction, RAY, brect, clipped_edge);
        // Keep correct ordering.
        (std::reverse)(clipped_edge.begin(), clipped_edge.end());
      } else if (cell1->contains_point() && cell2->contains_point()) {
        // Primary line edge.
        point_type origin((point1.x() + point2.x()) * 0.5, (point1.y() + point2.y()) * 0.5);
        point_type direction(point1.y() - point2.y(), point2.x() - point1.x());
        intersect(origin, direction, LINE, brect, clipped_edge);
      } else {
        point_type origin = cell1->contains_point() ? point1 : point2;
        point_type direction(point1.y() - point2.y(), point2.x() - point1.x());
        intersect(origin, direction, LINE, brect, clipped_edge);
      }
    }
  }

  // Clip an edge with the given rectangle.
  template <typename PointType>
  static void clip(const PointType &p1, const PointType &p2,
      const brect_type &brect, point_set_type &clipped_edge) {
    point_type start = get_point(p1);
    point_type end = get_point(p2);
    point_type direction(end.x() - start.x(), end.y() - start.y());
    if (brect.contains(start.x(), start.y()))
      clipped_edge.push_back(start);
    if (p1 != p2)
      intersect(start, direction, SEGMENT, brect, clipped_edge);
    if (brect.contains(end.x(), end.y()))
      clipped_edge.push_back(end);
  }

private:
  typedef typename TRAITS::ctype_converter_type ctype_converter_type;

  voronoi_utils();

  // There are three different types of linear primitive:
  //   1) SEGMENT - has two finite endpoints;
  //   2) RAY - has one finite and one infinite endpoint;
  //   3) LINE - has two infinite endpoints.
  enum EdgeType {
    SEGMENT = 0,
    RAY = 1,
    LINE = 2
  };

  template <typename P>
  static point_type get_point(const P &point) {
    coordinate_type x = to_fpt(point.x());
    coordinate_type y = to_fpt(point.y());
    return point_type(x, y);
  }

  static void intersect(
      const point_type &origin, const point_type &direction,
      EdgeType edge_type, const brect_type &brect,
      point_set_type &clipped_edge) {
    // Find the center of the rectangle.
    coordinate_type center_x = (brect.x_min() + brect.x_max()) * to_fpt(0.5);
    coordinate_type center_y = (brect.y_min() + brect.y_max()) * to_fpt(0.5);

    // Find the half-diagonal vector of the rectangle.
    coordinate_type len_x = brect.x_max() - center_x;
    coordinate_type len_y = brect.y_max() - center_y;

    // Find the vector between the origin and the center of the
    // rectangle.
    coordinate_type diff_x = origin.x() - center_x;
    coordinate_type diff_y = origin.y() - center_y;

    // Find the vector perpendicular to the direction vector.
    coordinate_type perp_x = direction.y();
    coordinate_type perp_y = -direction.x();

    // Projection of the vector between the origin and the center of
    // the rectangle on the line perpendicular to the direction vector.
    coordinate_type lexpr = magnitude(perp_x * diff_x + perp_y * diff_y);

    // Maximum projection among any of the half-diagonals of the
    // rectangle on the line perpendicular to the direction vector.
    coordinate_type rexpr =
        magnitude(perp_x * len_x) + magnitude(perp_y * len_y);

    // Intersection check. Compare projections.
    if (lexpr > rexpr)
      return;

    // Intersection parameters. If fT[i]_used is true than:
    // origin + fT[i] * direction = intersection point, i = 0 .. 1.
    // For different edge types next fT values are acceptable:
    // Segment: [0; 1];
    // Ray: [0; infinity];
    // Line: [-infinity; infinity].
    bool fT0_used = false;
    bool fT1_used = false;
    coordinate_type fT0 = 0;
    coordinate_type fT1 = 0;

    // Check for the intersections with the lines
    // going through the sides of the bounding rectangle.
    clip_line(+direction.x(), -diff_x - len_x, fT0_used, fT1_used, fT0, fT1);
    clip_line(-direction.x(), +diff_x - len_x, fT0_used, fT1_used, fT0, fT1);
    clip_line(+direction.y(), -diff_y - len_y, fT0_used, fT1_used, fT0, fT1);
    clip_line(-direction.y(), +diff_y - len_y, fT0_used, fT1_used, fT0, fT1);

    // Update intersections vector.
    if (fT0_used && check_extent(fT0, edge_type))
      clipped_edge.push_back(point_type(
          origin.x() + fT0 * direction.x(), origin.y() + fT0 * direction.y()));
    if (fT1_used && fT0 != fT1 && check_extent(fT1, edge_type))
      clipped_edge.push_back(point_type(
          origin.x() + fT1 * direction.x(), origin.y() + fT1 * direction.y()));
  }

  // Find intermediate points of the parabola.
  // Parabola is a locus of points equidistant from the point and segment
  // sites. intermediate_points should contain two initial endpoints
  // of the edge (Voronoi vertices). Intermediate points are inserted
  // between the given two endpoints.
  // Max_dist is the maximum distance allowed between parabola and line
  // segments that discretize it.
  static void fill_intermediate_points(
      point_type point_site, point_type segment_site_start,
      point_type segment_site_end, coordinate_type max_dist,
      point_set_type &intermediate_points) {
    // Apply the linear transformation to move start point of the
    // segment to the point with coordinates (0, 0) and the direction
    // of the segment to coincide the positive direction of the x-axis.
    coordinate_type segm_vec_x =
        segment_site_end.x() - segment_site_start.x();
    coordinate_type segm_vec_y =
        segment_site_end.y() - segment_site_start.y();
    coordinate_type sqr_segment_length =
        segm_vec_x * segm_vec_x + segm_vec_y * segm_vec_y;

    // Compute x-coordinates of the endpoints of the edge
    // in the transformed space.
    coordinate_type projection_start = sqr_segment_length *
        get_point_projection(
            intermediate_points[0], segment_site_start, segment_site_end);
    coordinate_type projection_end = sqr_segment_length *
        get_point_projection(
            intermediate_points[1], segment_site_start, segment_site_end);

    // Compute parabola parameters in the transformed space.
    // Parabola has next representation:
    // f(x) = ((x-rot_x)^2 + rot_y^2) / (2.0*rot_y).
    coordinate_type point_vec_x = point_site.x() - segment_site_start.x();
    coordinate_type point_vec_y = point_site.y() - segment_site_start.y();
    coordinate_type rot_x =
        segm_vec_x * point_vec_x + segm_vec_y * point_vec_y;
    coordinate_type rot_y =
        segm_vec_x * point_vec_y - segm_vec_y * point_vec_x;

    // Save the last point.
    point_type last_point = intermediate_points[1];
    intermediate_points.pop_back();

    // Use stack to avoid recursion.
    std::stack<coordinate_type> point_stack;
    point_stack.push(projection_end);
    coordinate_type cur_x = projection_start;
    coordinate_type cur_y = parabola_y(cur_x, rot_x, rot_y);

    // Adjust max_dist parameter in the transformed space.
    max_dist *= max_dist * sqr_segment_length;

    while (!point_stack.empty()) {
      coordinate_type new_x = point_stack.top();
      coordinate_type new_y = parabola_y(new_x, rot_x, rot_y);

      // Compute coordinates of the point of the parabola that is
      // furthest from the current line segment.
      coordinate_type mid_x =
          (new_y - cur_y) / (new_x - cur_x) * rot_y + rot_x;
      coordinate_type mid_y = parabola_y(mid_x, rot_x, rot_y);

      // Compute maximum distance between the given parabolic arc
      // and line segment that discretize it.
      coordinate_type dist = (new_y - cur_y) * (mid_x - cur_x) -
          (new_x - cur_x) * (mid_y - cur_y);
      dist = dist * dist / ((new_y - cur_y) * (new_y - cur_y) +
          (new_x - cur_x) * (new_x - cur_x));
      if (dist <= max_dist) {
        // Distance between parabola and line segment is
        // not greater than max_dist.
        point_stack.pop();
        coordinate_type inter_x = (segm_vec_x * new_x - segm_vec_y * new_y) /
            sqr_segment_length + segment_site_start.x();
        coordinate_type inter_y = (segm_vec_x * new_y + segm_vec_y * new_x) /
            sqr_segment_length + segment_site_start.y();
        intermediate_points.push_back(point_type(inter_x, inter_y));
        cur_x = new_x;
        cur_y = new_y;
      } else {
        point_stack.push(mid_x);
      }
    }

    // Update last point.
    intermediate_points.back() = last_point;
  }

  // Compute y(x) = ((x - a) * (x - a) + b * b) / (2 * b).
  static coordinate_type parabola_y(
      coordinate_type x, coordinate_type a, coordinate_type b) {
    return ((x - a) * (x - a) + b * b) / (to_fpt(2.0) * b);
  }

  // Check whether extent is compatible with the edge type.
  static bool check_extent(coordinate_type extent, EdgeType etype) {
    switch (etype) {
      case SEGMENT: return extent >= to_fpt(0.0) && extent <= to_fpt(1.0);
      case RAY: return extent >= to_fpt(0.0);
      case LINE: return true;
    }
    return true;
  }

  // Compute the absolute value.
  static inline coordinate_type magnitude(coordinate_type value) {
    return (value >= to_fpt(0.0)) ? value : -value;
  }

  // Find fT min and max values: fT = numer / denom.
  static void clip_line(
      coordinate_type denom, coordinate_type numer,
      bool &fT0_used, bool &fT1_used,
      coordinate_type &fT0, coordinate_type &fT1) {
    if (denom > to_fpt(0.0)) {
      if (fT1_used && numer > denom * fT1)
        return;
      if (!fT0_used || numer > denom * fT0) {
        fT0_used = true;
        fT0 = numer / denom;
      }
    } else if (denom < to_fpt(0.0)) {
      if (fT0_used && numer > denom * fT0)
        return;
      if (!fT1_used || numer > denom * fT1) {
        fT1_used = true;
        fT1 = numer / denom;
      }
    }
  }

  // Get normalized length of the distance between:
  //   1) point projection onto the segment;
  //   2) start point of the segment.
  // Return this length divided by the segment length.
  // This is made to avoid sqrt computation during transformation from
  // the initial space to the transformed one and vice versa.
  // Assumption is made that projection of the point lies
  // between the start-point and endpoint of the segment.
  static coordinate_type get_point_projection(
      const point_type &point,
      const point_type &segment_start,
      const point_type &segment_end) {
    coordinate_type segment_vec_x = segment_end.x() - segment_start.x();
    coordinate_type segment_vec_y = segment_end.y() - segment_start.y();
    coordinate_type point_vec_x = point.x() - segment_start.x();
    coordinate_type point_vec_y = point.y() - segment_start.y();
    coordinate_type sqr_segment_length =
        segment_vec_x * segment_vec_x + segment_vec_y * segment_vec_y;
    coordinate_type vec_dot =
        segment_vec_x * point_vec_x + segment_vec_y * point_vec_y;
    return vec_dot / sqr_segment_length;
  }

  template <typename CT>
  static coordinate_type to_fpt(const CT& value) {
    static ctype_converter_type converter;
    return converter(value);
  }
};
}  // polygon
}  // boost

#endif  // BOOST_POLYGON_VORONOI_UTILS
