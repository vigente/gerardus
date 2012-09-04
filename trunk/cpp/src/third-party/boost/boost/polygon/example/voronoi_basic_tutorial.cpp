// Boost.Polygon library voronoi_basic_tutorial.cpp file

//          Copyright Andrii Sydorchuk 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

#include <cstdio>
#include <vector>

#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/voronoi_utils.hpp>
using namespace boost::polygon;

struct Point {
  int a;
  int b;
  Point (int x, int y) : a(x), b(y) {}
};

struct Segment {
  Point p0;
  Point p1;
  Segment (int x1, int y1, int x2, int y2) : p0(x1, y1), p1(x2, y2) {}
};

namespace boost {
namespace polygon {

template <>
struct geometry_concept<Point> { typedef point_concept type; };
  
template <>
struct point_traits<Point> {
  typedef int coordinate_type;

  static inline coordinate_type get(const Point& point, orientation_2d orient) {
    return (orient == HORIZONTAL) ? point.a : point.b;
  }
};

template <>
struct geometry_concept<Segment> { typedef segment_concept type; };

template <>
struct segment_traits<Segment> {
  typedef int coordinate_type;
  typedef Point point_type;
    
  static inline point_type get(const Segment& segment, direction_1d dir) {
    return dir.to_int() ? segment.p1 : segment.p0;
  }
};
}  // polygon
}  // boost

// Traversing Voronoi edges using edge iterator.
int iterate_primary_edges1(const voronoi_diagram<double> &vd) {
  int result = 0;
  for (voronoi_diagram<double>::const_edge_iterator it = vd.edges().begin();
       it != vd.edges().end(); ++it) {
    if (it->is_primary())
      ++result;
  }
  return result;
}

// Traversing Voronoi edges using cell iterator.
int iterate_primary_edges2(const voronoi_diagram<double> &vd) {
  int result = 0;
  for (voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin();
       it != vd.cells().end(); ++it) {
    const voronoi_diagram<double>::cell_type &cell = *it;
    const voronoi_diagram<double>::edge_type *edge = cell.incident_edge();
    // This is convenient way to iterate edges around Voronoi cell.
    do {
      if (edge->is_primary())
        ++result;
      edge = edge->next();
    } while (edge != cell.incident_edge());
  }
  return result;
}

// Traversing Voronoi edges using vertex iterator.
// As opposite to the above two functions this one will not iterate through edges
// without finite endpoints and will iterate only once through edges with single
// finite endpoint.
int iterate_primary_edges3(const voronoi_diagram<double> &vd) {
  int result = 0;
  for (voronoi_diagram<double>::const_vertex_iterator it = vd.vertices().begin();
       it != vd.vertices().end(); ++it) {
    const voronoi_diagram<double>::vertex_type &vertex = *it;
    const voronoi_diagram<double>::edge_type *edge = vertex.incident_edge();
    // This is convenient way to iterate edges around Voronoi vertex.
    do {
      if (edge->is_primary())
        ++result;
      edge = edge->rot_next();
    } while (edge != vertex.incident_edge());
  }
  return result;
}

// Prototype of a function that renders segments.
void draw_segment(double x1, double y1, double x2, double y2) {
  printf("Rendering segment: ");
  printf("%7.3f %7.3f %7.3f %7.3f\n", x1, y1, x2, y2);
}

void render_diagram(const voronoi_diagram<double> &vd,
                    const voronoi_utils<double>::brect_type &bbox) {
  const std::size_t VISITED = 1;
  for (voronoi_diagram<double>::const_edge_iterator it = vd.edges().begin();
       it != vd.edges().end(); ++it) {
    // We use color to mark visited edges.
    it->color(VISITED);
    // Don't render the same edge twice.
    if (it->twin()->color()) continue;
    voronoi_utils<double>::point_set_type polyline;
    if (it->is_linear())
      voronoi_utils<double>::clip(*it, bbox, polyline);
    else
      // Parabolic edges are always finite.
      voronoi_utils<double>::discretize(*it, 1E-1, polyline);
    // Note: discretized edges may also lie outside of the bbox.
    // So user might do additional clipping before rendering each such edge.  
    for (std::size_t i = 1; i < polyline.size(); ++i)
      draw_segment(polyline[i-1].x(), polyline[i-1].y(),
                   polyline[i].x(), polyline[i].y());
  }
}

int main() {
  // Preparing Input Geometries.
  std::vector<Point> points;
  points.push_back(Point(0, 0));
  points.push_back(Point(1, 6));
  std::vector<Segment> segments;
  segments.push_back(Segment(-4, 5, 5, -1));
  segments.push_back(Segment(3, -11, 13, -1));

  // Construction of the Voronoi Diagram.
  voronoi_diagram<double> vd;
  construct_voronoi(points.begin(), points.end(), segments.begin(), segments.end(), &vd);

  // Traversing Voronoi Graph.
  {
    printf("Traversing Voronoi graph.\n");
    printf("Number of visited primary edges using edge iterator: %d\n", iterate_primary_edges1(vd));
    printf("Number of visited primary edges using cell iterator: %d\n", iterate_primary_edges2(vd));
    printf("Number of visited primary edges using vertex iterator: %d\n", iterate_primary_edges3(vd));
    printf("\n");
  }

  // Using color member of the Voronoi primitives..
  {
    for (voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin();
         it != vd.cells().end(); ++it) {
      const voronoi_diagram<double>::cell_type &cell = *it;
      const voronoi_diagram<double>::edge_type *edge = cell.incident_edge();
      std::size_t count = 0;
      do {
        ++count;
        edge = edge->next();
      } while (edge != cell.incident_edge());
      cell.color(count);
    }

    // Count the average number of edges.
    double total = 0;
    for (voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin();
         it != vd.cells().end(); ++it) {
      total += it->color();
    }
    total /= vd.cells().size();
    printf("The average number of edges per Voronoi cell is equal to: %3.1f\n", total);
    printf("\n");
  } 

  // Rendering Voronoi diagram.
  {
    // Construct clipping bounding rectangle.
    bounding_rectangle<double> bbox;
    for (std::vector<Point>::iterator it = points.begin(); it != points.end(); ++it)
      bbox.update(it->a, it->b);
    for (std::vector<Segment>::iterator it = segments.begin(); it != segments.end(); ++it) {
      bbox.update(it->p0.a, it->p0.b);
      bbox.update(it->p1.a, it->p1.b);
    }
    // Add 10% offset to the bounding rectangle.
    bbox = voronoi_utils<double>::scale(bbox, 1.1);
    // Render Voronoi diagram.
    render_diagram(vd, bbox);
  }

  return 0;
}
