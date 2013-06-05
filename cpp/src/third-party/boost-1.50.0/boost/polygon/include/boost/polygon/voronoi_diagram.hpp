// Boost.Polygon library voronoi_diagram.hpp header file

//          Copyright Andrii Sydorchuk 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_POLYGON_VORONOI_DIAGRAM
#define BOOST_POLYGON_VORONOI_DIAGRAM

#include <vector>

#include "detail/voronoi_ctypes.hpp"
#include "detail/voronoi_structures.hpp"

namespace boost {
namespace polygon {

// Forward declarations.
template <typename T>
class voronoi_edge;

// Represents Voronoi cell.
// Data members:
//   1) pointer to the incident edge
//   2) site inside cell
//   3) mutable color member
// Cell may contain point or segment site inside.
template <typename T>
class voronoi_cell {
public:
  typedef T coordinate_type;
  typedef detail::point_2d<coordinate_type> point_type;
  typedef std::size_t color_type;
  typedef voronoi_edge<coordinate_type> voronoi_edge_type;

  voronoi_cell(const point_type &p1, voronoi_edge_type *edge) :
      point0_(p1),
      point1_(p1),
      color_(0),
      incident_edge_(edge) {}

  voronoi_cell(const point_type &p1,
               const point_type &p2,
               voronoi_edge_type *edge) :
      point0_(p1),
      point1_(p2),
      color_(0),
      incident_edge_(edge) {}

  // Returns true if the cell contains point site, false else.
  bool contains_point() const { return point0_ == point1_; }

  // Returns true if the cell contains segment site, false else.
  bool contains_segment() const { return point0_ != point1_; }

  // Degenerate cells don't have any incident edges.
  bool is_degenerate() const { return incident_edge_ == NULL; }

  // Returns site point in case cell contains point site,
  // the first endpoint of the segment site else.
  const point_type &point0() const { return point0_; }

  // Returns site point in case cell contains point site,
  // the second endpoint of the segment site else.
  const point_type &point1() const { return point1_; }

  color_type color() const { return color_; }
  void color(const color_type& color) const { color_ = color; }

  voronoi_edge_type *incident_edge() { return incident_edge_; }
  const voronoi_edge_type *incident_edge() const { return incident_edge_; }
  void incident_edge(voronoi_edge_type *e) { incident_edge_ = e; }

private:
  point_type point0_;
  point_type point1_;
  mutable color_type color_;
  voronoi_edge_type *incident_edge_;
};

// Represents Voronoi vertex.
// Data members:
//   1) vertex point itself
//   2) pointer to the incident edge
//   3) mutable color member
template <typename T>
class voronoi_vertex {
public:
  typedef T coordinate_type;
  typedef detail::point_2d<T> point_type;
  typedef std::size_t color_type;
  typedef voronoi_edge<coordinate_type> voronoi_edge_type;

  voronoi_vertex(const point_type &vertex, voronoi_edge_type *edge) :
      vertex_(vertex),
      color_(0),
      incident_edge_(edge) {}

  const point_type &vertex() const { return vertex_; }

  bool is_degenerate() const { return incident_edge_ == NULL; }

  color_type color() const { return color_; }
  void color(const color_type& color) const { color_ = color; }

  voronoi_edge_type *incident_edge() { return incident_edge_; }
  const voronoi_edge_type *incident_edge() const { return incident_edge_; }
  void incident_edge(voronoi_edge_type *e) { incident_edge_ = e; }

private:
  point_type vertex_;
  mutable color_type color_;
  voronoi_edge_type *incident_edge_;
};

// Half-edge data structure. Represents Voronoi edge.
// Data members:
//   1) pointer to the corresponding cell
//   2) pointer to the vertex that is the starting
//      point of the half-edge
//   3) pointer to the twin edge
//   4) pointer to the CCW next edge
//   5) pointer to the CCW prev edge
//   6) mutable color member
template <typename T>
class voronoi_edge {
public:
  typedef T coordinate_type;
  typedef voronoi_cell<coordinate_type> voronoi_cell_type;
  typedef voronoi_vertex<coordinate_type> voronoi_vertex_type;
  typedef voronoi_edge<coordinate_type> voronoi_edge_type;
  typedef std::size_t color_type;

  voronoi_edge() :
      cell_(NULL),
      vertex_(NULL),
      twin_(NULL),
      next_(NULL),
      prev_(NULL),
      color_(0) {}

  voronoi_cell_type *cell() { return cell_; }
  const voronoi_cell_type *cell() const { return cell_; }
  void cell(voronoi_cell_type *c) { cell_ = c; }

  voronoi_vertex_type *vertex0() { return vertex_; }
  const voronoi_vertex_type *vertex0() const { return vertex_; }
  void vertex0(voronoi_vertex_type *v) { vertex_ = v; }

  voronoi_vertex_type *vertex1() { return twin_->vertex0(); }
  const voronoi_vertex_type *vertex1() const { return twin_->vertex0(); }
  void vertex1(voronoi_vertex_type *v) { twin_->vertex0(v); }

  voronoi_edge_type *twin() { return twin_; }
  const voronoi_edge_type *twin() const { return twin_; }
  void twin(voronoi_edge_type *e) { twin_ = e; }

  voronoi_edge_type *next() { return next_; }
  const voronoi_edge_type *next() const { return next_; }
  void next(voronoi_edge_type *e) { next_ = e; }

  voronoi_edge_type *prev() { return prev_; }
  const voronoi_edge_type *prev() const { return prev_; }
  void prev(voronoi_edge_type *e) { prev_ = e; }

  // Returns a pointer to the rotation next edge
  // over the starting point of the half-edge.
  voronoi_edge_type *rot_next() {
    return (vertex_) ? prev_->twin() : NULL;
  }
  const voronoi_edge_type *rot_next() const {
    return (vertex_) ? prev_->twin() : NULL;
  }

  // Returns a pointer to the rotation prev edge
  // over the starting point of the half-edge.
  voronoi_edge_type *rot_prev() {
    return (vertex_) ? twin_->next() : NULL;
  }
  const voronoi_edge_type *rot_prev() const {
    return (vertex_) ? twin_->next() : NULL;
  }

  // Returns true if the edge is finite (segment, parabolic arc).
  // Returns false if the edge is infinite (ray, line).
  bool is_finite() const { return vertex0() && vertex1(); }

  // Returns true if the edge is linear (segment, ray, line).
  // Returns false if the edge is curved (parabolic arc).
  bool is_linear() const {
    if (!is_primary())
      return true;
    return !(cell()->contains_segment() ^ twin()->cell()->contains_segment());
  }

  // Returns true if the edge is curved (parabolic arc).
  // Returns false if the edge is linear (segment, ray, line).
  bool is_curved() const {
    if (!is_primary())
      return false;
    return (cell()->contains_segment() ^ twin()->cell()->contains_segment());
  }

  // Returns false if edge goes through the endpoint of the segment.
  // Returns true else.
  bool is_primary() const {
    bool flag1 = cell_->contains_segment();
    bool flag2 = twin_->cell()->contains_segment();
    if (flag1 && !flag2) {
      return cell_->point0() != twin_->cell()->point0() &&
             cell_->point1() != twin_->cell()->point0();
    }
    if (!flag1 && flag2) {
      return twin_->cell()->point0() != cell_->point0() &&
             twin_->cell()->point1() != cell_->point0();
    }
    return true;
  }

  color_type color() const { return color_; }
  void color(const color_type& color) const { color_ = color; }

private:
  voronoi_cell_type *cell_;
  voronoi_vertex_type *vertex_;
  voronoi_edge_type *twin_;
  voronoi_edge_type *next_;
  voronoi_edge_type *prev_;
  mutable color_type color_;
};

template <typename T>
struct voronoi_diagram_traits {
  typedef T coordinate_type;
  typedef struct {
    template <typename CT>
    coordinate_type operator()(const CT& value) {
      return static_cast<coordinate_type>(value);
    }
  } ctype_converter_type;
  typedef detail::point_2d<coordinate_type> point_type;
  typedef voronoi_cell<coordinate_type> cell_type;
  typedef voronoi_vertex<coordinate_type> vertex_type;
  typedef voronoi_edge<coordinate_type> edge_type;
  typedef class {
  public:
    enum { ULPS = 128 };
    bool operator()(const point_type &v1, const point_type &v2) const {
      return (ulp_cmp(v1.x(), v2.x(), ULPS) ==
              detail::ulp_comparison<T>::EQUAL) &&
             (ulp_cmp(v1.y(), v2.y(), ULPS) ==
              detail::ulp_comparison<T>::EQUAL);
    }
  private:
    typename detail::ulp_comparison<T> ulp_cmp;
  } vertex_equality_predicate_type;
};

// Voronoi output data structure.
// CCW ordering is used on the faces perimeter and around the vertices.
template <typename T, typename TRAITS = voronoi_diagram_traits<T> >
class voronoi_diagram {
public:
  typedef typename TRAITS::coordinate_type coordinate_type;
  typedef typename TRAITS::point_type point_type;
  typedef typename TRAITS::cell_type cell_type;
  typedef typename TRAITS::vertex_type vertex_type;
  typedef typename TRAITS::edge_type edge_type;

  typedef std::vector<cell_type> cell_container_type;
  typedef typename cell_container_type::iterator cell_iterator;
  typedef typename cell_container_type::const_iterator const_cell_iterator;

  typedef std::vector<vertex_type> vertex_container_type;
  typedef typename vertex_container_type::iterator vertex_iterator;
  typedef typename vertex_container_type::const_iterator const_vertex_iterator;

  typedef std::vector<edge_type> edge_container_type;
  typedef typename edge_container_type::iterator edge_iterator;
  typedef typename edge_container_type::const_iterator const_edge_iterator;

  // This builder class is mainly used to hide from the user methods that
  // construct the Voronoi diagram.
  class voronoi_diagram_builder {
  public:
    void vd(voronoi_diagram *vd) {
      vd_ = vd;
    }

    bool done() {
      return vd_ == NULL;
    }

    void reserve(int num_sites) {
      vd_->reserve(num_sites);
    }

    template <typename SEvent>
    void process_single_site(const SEvent &site) {
      vd_->process_single_site(site);
    }

    template <typename SEvent>
    std::pair<void*, void*> insert_new_edge(
        const SEvent &site1, const SEvent &site2) {
      return vd_->insert_new_edge(site1, site2);
    }

    template <typename SEvent, typename CEvent>
    std::pair<void *, void *> insert_new_edge(
        const SEvent &site1, const SEvent &site3, const CEvent &circle,
        void *data12, void *data23) {
      return vd_->insert_new_edge(site1, site3, circle, data12, data23);
    }

    void build() {
      vd_->build();
      vd_ = NULL;
    }

  private:
    voronoi_diagram *vd_;
  };

  voronoi_diagram() {
    builder_.vd(&(*this));
  }

  void clear() {
    cells_.clear();
    vertices_.clear();
    edges_.clear();
    builder_.vd(&(*this));
  }

  const cell_container_type &cells() const {
    return cells_;
  }

  const vertex_container_type &vertices() const {
    return vertices_;
  }

  const edge_container_type &edges() const {
    return edges_;
  }

  unsigned int num_cells() const {
    return cells_.size();
  }

  unsigned int num_edges() const {
    return edges_.size() >> 1;
  }

  unsigned int num_vertices() const {
    return vertices_.size();
  }

  voronoi_diagram_builder *builder() {
    if (builder_.done()) {
      return NULL;
    } else {
      return &builder_;
    }
  }

private:
  typedef typename TRAITS::ctype_converter_type ctype_converter_type;
  typedef typename TRAITS::vertex_equality_predicate_type
    vertex_equality_predicate_type;

  friend class voronoi_diagram_builder;

  void reserve(int num_sites) {
    cells_.reserve(num_sites);
    vertices_.reserve(num_sites << 1);
    edges_.reserve((num_sites << 2) + (num_sites << 1));
  }

  // Update the Voronoi output in case of a single point input.
  template <typename SEvent>
  void process_single_site(const SEvent &site) {
    // Update cell records.
    point_type p = prepare_point(site.point0());
    cells_.push_back(cell_type(p, NULL));
  }

  // Insert a new half-edge into the output data structure.
  // Takes as input left and right sites that form a new bisector.
  // Returns a pair of pointers to a new half-edges.
  template <typename SEvent>
  std::pair<void*, void*> insert_new_edge(
      const SEvent &site1, const SEvent &site2) {
    // Get sites' indexes.
    int site_index1 = site1.sorted_index();
    int site_index2 = site2.sorted_index();

    // Create a new half-edge that belongs to the first site.
    edges_.push_back(edge_type());
    edge_type &edge1 = edges_.back();

    // Create a new half-edge that belongs to the second site.
    edges_.push_back(edge_type());
    edge_type &edge2 = edges_.back();

    // Add the initial cell during the first edge insertion.
    if (cells_.empty()) {
      process_single_site(site1);
    }

    // The second site represents a new site during site event
    // processing. Add a new cell to the cell records.
    cells_.push_back(cell_type(prepare_point(site2.point0()),
                               prepare_point(site2.point1()),
                               NULL));

    // Set up pointers to cells.
    edge1.cell(&cells_[site_index1]);
    edge2.cell(&cells_[site_index2]);

    // Set up twin pointers.
    edge1.twin(&edge2);
    edge2.twin(&edge1);

    // Return a pointer to the new half-edge.
    return std::make_pair(&edge1, &edge2);
  }

  // Insert a new half-edge into the output data structure with the
  // start at the point where two previously added half-edges intersect.
  // Takes as input two sites that create a new bisector, circle event
  // that corresponds to the intersection point of the two old half-edges,
  // pointers to those half-edges. Half-edges' direction goes out of the
  // new Voronoi vertex point. Returns a pair of pointers to a new half-edges.
  template <typename SEvent, typename CEvent>
  std::pair<void *, void *> insert_new_edge(
      const SEvent &site1, const SEvent &site3, const CEvent &circle,
      void *data12, void *data23) {
    edge_type *edge12 = static_cast<edge_type*>(data12);
    edge_type *edge23 = static_cast<edge_type*>(data23);

    // Add a new Voronoi vertex.
    vertices_.push_back(vertex_type(prepare_point(circle), NULL));
    vertex_type &new_vertex = vertices_.back();

    // Update vertex pointers of the old edges.
    edge12->vertex0(&new_vertex);
    edge23->vertex0(&new_vertex);

    // Add a new half-edge.
    edges_.push_back(edge_type());
    edge_type &new_edge1 = edges_.back();
    new_edge1.cell(&cells_[site1.sorted_index()]);

    // Add a new half-edge.
    edges_.push_back(edge_type());
    edge_type &new_edge2 = edges_.back();
    new_edge2.cell(&cells_[site3.sorted_index()]);

    // Update twin pointers.
    new_edge1.twin(&new_edge2);
    new_edge2.twin(&new_edge1);

    // Update vertex pointer.
    new_edge2.vertex0(&new_vertex);

    // Update Voronoi prev/next pointers.
    edge12->prev(&new_edge1);
    new_edge1.next(edge12);
    edge12->twin()->next(edge23);
    edge23->prev(edge12->twin());
    edge23->twin()->next(&new_edge2);
    new_edge2.prev(edge23->twin());

    // Return a pointer to the new half-edge.
    return std::make_pair(&new_edge1, &new_edge2);
  }

  void build() {
    // Remove degenerate edges.
    edge_iterator last_edge = edges_.begin();
    for (edge_iterator it = edges_.begin(); it != edges_.end(); it += 2) {
      const vertex_type *v1 = it->vertex0();
      const vertex_type *v2 = it->vertex1();
      if (v1 && v2 && vertex_equality_predicate_(v1->vertex(), v2->vertex())) {
        remove_edge(&(*it));
      } else {
        if (it != last_edge) {
          edge_type *e1 = &(*last_edge = *it);
          edge_type *e2 = &(*(last_edge + 1) = *(it + 1));
          e1->twin(e2);
          e2->twin(e1);
          if (e1->prev()) {
            e1->prev()->next(e1);
            e2->next()->prev(e2);
          }
          if (e2->prev()) {
            e1->next()->prev(e1);
            e2->prev()->next(e2);
          }
        }
        last_edge += 2;
      }
    }
    edges_.erase(last_edge, edges_.end());

    // Set up incident edge pointers for cells and vertices.
    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      it->cell()->incident_edge(&(*it));
      if (it->vertex0()) {
        it->vertex0()->incident_edge(&(*it));
      }
    }

    // Remove degenerate vertices.
    vertex_iterator last_vertex = vertices_.begin();
    for (vertex_iterator it = vertices_.begin(); it != vertices_.end(); ++it) {
      if (it->incident_edge()) {
        if (it != last_vertex) {
          *last_vertex = *it;
          vertex_type *v = &(*last_vertex);
          edge_type *e = v->incident_edge();
          do {
            e->vertex0(v);
            e = e->rot_next();
          } while (e != v->incident_edge());
        }
        ++last_vertex;
      }
    }
    vertices_.erase(last_vertex, vertices_.end());

    // Set up next/prev pointers for infinite edges.
    if (vertices_.empty()) {
      if (!edges_.empty()) {
        // Update prev/next pointers for the line edges.
        edge_iterator edge_it = edges_.begin();
        edge_type *edge1 = &(*edge_it);
        edge1->next(edge1);
        edge1->prev(edge1);
        ++edge_it;
        edge1 = &(*edge_it);
        ++edge_it;

        while (edge_it != edges_.end()) {
          edge_type *edge2 = &(*edge_it);
          ++edge_it;

          edge1->next(edge2);
          edge1->prev(edge2);
          edge2->next(edge1);
          edge2->prev(edge1);

          edge1 = &(*edge_it);
          ++edge_it;
        }

        edge1->next(edge1);
        edge1->prev(edge1);
      }
    } else {
      // Update prev/next pointers for the ray edges.
      for (cell_iterator cell_it = cells_.begin();
         cell_it != cells_.end(); ++cell_it) {
        if (cell_it->is_degenerate())
          continue;
        // Move to the previous edge while
        // it is possible in the CW direction.
        edge_type *left_edge = cell_it->incident_edge();
        while (left_edge->prev() != NULL) {
          left_edge = left_edge->prev();
          // Terminate if this is not a boundary cell.
          if (left_edge == cell_it->incident_edge())
            break;
        }

        if (left_edge->prev() != NULL)
          continue;

        edge_type *right_edge = cell_it->incident_edge();
        while (right_edge->next() != NULL)
          right_edge = right_edge->next();
        left_edge->prev(right_edge);
        right_edge->next(left_edge);
      }
    }
  }

  template <typename P>
  point_type prepare_point(const P& p) {
    coordinate_type nx = convert_(p.x());
    coordinate_type ny = convert_(p.y());
    return point_type(nx, ny);
  }

  // Remove degenerate edge.
  void remove_edge(edge_type *edge) {
    // Update the endpoints of the incident edges to the second vertex.
    vertex_type *vertex = edge->vertex0();
    edge_type *updated_edge = edge->twin()->rot_next();
    while (updated_edge != edge->twin()) {
      updated_edge->vertex0(vertex);
      updated_edge = updated_edge->rot_next();
    }

    edge_type *edge1 = edge;
    edge_type *edge2 = edge->twin();

    edge_type *edge1_rot_prev = edge1->rot_prev();
    edge_type *edge1_rot_next = edge1->rot_next();

    edge_type *edge2_rot_prev = edge2->rot_prev();
    edge_type *edge2_rot_next = edge2->rot_next();

    // Update prev/next pointers for the incident edges.
    edge1_rot_next->twin()->next(edge2_rot_prev);
    edge2_rot_prev->prev(edge1_rot_next->twin());
    edge1_rot_prev->prev(edge2_rot_next->twin());
    edge2_rot_next->twin()->next(edge1_rot_prev);
  }

  cell_container_type cells_;
  vertex_container_type vertices_;
  edge_container_type edges_;

  voronoi_diagram_builder builder_;
  ctype_converter_type convert_;
  vertex_equality_predicate_type vertex_equality_predicate_;

  // Disallow copy constructor and operator=
  voronoi_diagram(const voronoi_diagram&);
  void operator=(const voronoi_diagram&);
};
}  // polygon
}  // boost

#endif  // BOOST_POLYGON_VORONOI_DIAGRAM
