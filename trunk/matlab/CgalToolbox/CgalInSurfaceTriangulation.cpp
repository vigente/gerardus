/* CgalInSurfaceTriangulation.cpp
 *
 * CGAL_INSURFACETRIANGULATION  Find whether a point is inside or outside a
 * closed surface
 *
 *   This function evaluates whether one or more points belong inside a
 *   closed surface. First, we check whether the point is on the surface
 *   itself (in that case, it's considered inside). If not, a ray is
 *   projected from the point and the intersections with the surface are
 *   counted. An odd number means that the point is inside. This approach
 *   fails if the point is not on the surface, but the ray lies on the
 *   surface or crosses a vertex, because this spans many arbitrary
 *   intersections. To solve this problem, a few rays are used for each
 *   point, and the majority vote decides whether it's inside or outside.
 *
 * ISIN = cgal_insurftri(TRI, X, XI)
 *
 *   TRI is a 3-column matrix. Each row represents the indices of the tree
 *   vertices that form a triangle. TRI as a whole represents the closed
 *   surface.
 *
 *   X is a 3-column matrix. Each row represents the Cartesian coordinates
 *   of a vertex on the surface, indexed by TRI values.
 *
 *   XI is a 3-column matrix. Each row represents the Carterian coordinates
 *   of a point for which we want to find whether it's inside or outside the
 *   closed surface.
 *
 *   ISIN is a boolean vector with one element per point in XI. True means
 *   that the corresponding point is inside the closed surface (or on the
 *   surface), and false means that it's outside.
 *
 * ISIN = cgal_insurftri(TRI, X, XI, DIRECTIONS, TOL)
 *
 *   DIRECTIONS is a 3-column matrix. Each row represents a vector with a
 *   ray direction. By default, 
 *
 *          DIRECTIONS=[ 1.0,  0.0,  0.0; ...
 *                      -1.0,  1.0,  1.0; ...
 *                      -1.0, -1.0, -1.0]
 *
 *   This default can fail with regular voxels, as rays may cross vertices.
 *   A good practical alternative is to use a few random directions, e.g.
 *
 *          DIRECTIONS=rand(4, 3);
 *
 *   TOL is a scalar with the distance tolerance. Points at distance <= TOL
 *   are considered to be on the surface, and thus "inside". By default,
 *   TOL=1e-15.
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012 University of Oxford
  * Version: 0.1.0
  * $Rev$
  * $Date$
  *
  * University of Oxford means the Chancellor, Masters and Scholars of
  * the University of Oxford, having an administrative office at
  * Wellington Square, Oxford OX1 2JD, UK. 
  *
  * This file is part of Gerardus.
  *
  * This program is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
  *
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details. The offer of this
  * program under the terms of the License is subject to the License
  * being interpreted in accordance with English Law and subject to any
  * action against the University of Oxford being under the jurisdiction
  * of the English Courts.
  *
  * You should have received a copy of the GNU General Public License
  * along with this program.  If not, see
  * <http://www.gnu.org/licenses/>.
  */

#ifndef CGALINSURFACETRIANGULATION
#define CGALINSURFACETRIANGULATION

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>

/* Gerardus headers */
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double>            K;

typedef CGAL::Ray_3<K>                            Ray;
typedef CGAL::Line_3<K>                           Line;
typedef CGAL::Point_3<K>                          Point;
typedef CGAL::Direction_3<K>                      Direction;
typedef CGAL::Triangle_3<K>                       Triangle;

typedef std::list<Triangle>::iterator             Iterator;
typedef CGAL::AABB_triangle_primitive<K,Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive>           AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits>     Tree;

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->SetMatlabArgumentsPointer(nrhs, prhs);

  // check that we have at least a filter name and input image
  matlabImport->CheckNumberOfArguments(3, 5);

  // interface to deal with outputs to Matlab
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->SetMatlabArgumentsPointer(nlhs, plhs);

  // check number of outputs the user is asking for
  matlabExport->CheckNumberOfArguments(0, 1);

  // if any of the inputs is empty, the output is empty too
  if (mxIsEmpty(prhs[0]) || mxIsEmpty(prhs[1]) || mxIsEmpty(prhs[2])) {
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    return;
  }

  // if user provides the ray directions, read them
  std::vector<Direction> direction; // ray query directions
  direction.push_back(Direction(1.0, 0.0, 0.0));    // default directions if 
  direction.push_back(Direction(-1.0, 1.0, 1.0));   // not provided by the
  direction.push_back(Direction(-1.0, -1.0, -1.0)); // user
  direction = matlabImport->GetVectorOfStaticVector3Argument<Direction>(3, "DIR", direction);

  // distance tolerance value
  double tol = matlabImport->GetScalarArgument<double>(4, "TRI", 1e-15);

  // point coordinates with NaN values in case there's a problem reading them
  Point def(mxGetNaN(), mxGetNaN(), mxGetNaN());

  // get size of input matrix
  mwSize nrowsTri = mxGetM(prhs[0]);
  mwSize nrowsXi = mxGetM(prhs[2]);
  mwSize ncolsTri = mxGetN(prhs[0]);
  mwSize ncolsX = mxGetN(prhs[1]);
  mwSize ncolsXi = mxGetN(prhs[2]);
  if ((ncolsTri != 3) || (ncolsX != 3) || (ncolsXi != 3)) {
    mexErrMsgTxt("All input arguments must have 3 columns");
  }

  // initialise output
  plhs[0] = mxCreateNumericMatrix(nrowsXi, 1, mxLOGICAL_CLASS, mxREAL);
  if (plhs[0] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output");
  }
  
  // pointer to the output
  bool *isin = (bool *)mxGetData(plhs[0]);
  if (isin == NULL) {
    mexErrMsgTxt("Cannot get pointer to allocated output");
  }

  // read triangular mesh from function
  std::list<Triangle> triangles;
  mwIndex v0, v1, v2; // indices of the 3 vertices of each triangle
  Point x0, x1, x2; // coordinates of the 3 vertices of each triangle
  for (mwIndex i = 0; i < nrowsTri; ++i) {

    // get indices of the 3 vertices of each triangle. These indices
    // follow Matlab's convention v0 = 1, 2, ..., n
    v0 = matlabImport->GetScalarArgument<mwIndex>(0, i, 0, "TRI", mxGetNaN());
    v1 = matlabImport->GetScalarArgument<mwIndex>(0, i, 1, "TRI", mxGetNaN());
    v2 = matlabImport->GetScalarArgument<mwIndex>(0, i, 2, "TRI", mxGetNaN());
    if (mxIsNaN(v0) || mxIsNaN(v1) || mxIsNaN(v2)) {
      mexErrMsgTxt("Parameter TRI: Vertex index is NaN");
    }
    
    // get coordinates of the 3 vertices (substracting 1 so that
    // indices follow the C++ convention 0, 1, ..., n-1)
    x0 = matlabImport->GetStaticVector3Argument<Point>(1, v0 - 1, "X", def);
    x1 = matlabImport->GetStaticVector3Argument<Point>(1, v1 - 1, "X", def);
    x2 = matlabImport->GetStaticVector3Argument<Point>(1, v2 - 1, "X", def);

    // add triangle to the list of triangles in the surface
    triangles.push_back(Triangle(x0, x1, x2));
  }

  // construct AABB tree
  Tree tree(triangles.begin(),triangles.end());

  // construct internal data structure to accelerate distance queries
  if (!tree.accelerate_distance_queries()) {
    mexErrMsgTxt("Not enough memory to accelerate distance queries");
  }

  // loop every point that is tested to see whether it's inside or
  // outside the surface
  Point xi; // test point coordinates
  double d; // distance from test point to surface
  for (mwIndex i = 0; i < nrowsXi; ++i) {

    // get point coordinates to be tested
    xi = matlabImport->GetStaticVector3Argument<Point>(2, i, "XI", def);

    // minimum distance from the test point to the surface
    d = tree.squared_distance(xi);

    // if the test point is close enough to the surface, we consider
    // it part of the surface, and thus an inside point
    if (d <= tol) {
      isin[i] = true;
      continue;
    }

    // otherwise, we are going to shoot rays from the point in
    // arbitrary directions, and count the intersections with the
    // surface. An odd number of intersections means that the point is
    // inside. If a ray lies on an edge or facet, it will produce many
    // spurious intersections. That's why we use several rays, and
    // take the result of the majority

    // rays starting at the test point, in the ray directions
    Ray ray_query;
    unsigned int isin_vote = 0;
    for (unsigned int j = 0; j < direction.size(); ++j) {

      ray_query = Ray(xi, direction[j]);

      // if the number of intersections is odd, the point is inside the
      // surface
      isin_vote += tree.number_of_intersected_primitives(ray_query) % 2;
    }
    
    // the majority of rays decide whether the point is inside or outside
    isin[i] = isin_vote > (direction.size() / 2);
  }
  
}

#endif /* CGALINSURFACETRIANGULATION */
