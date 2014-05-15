/* CgalInSurfaceTriangulation.cpp
 *
 * CGAL_INSURFTRI  Find whether a point is inside or outside a closed
 * surface.
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
 *   TRI is a 3-column matrix. Each row represents the indices of the three
 *   vertices that form a triangle. TRI as a whole represents the closed
 *   surface.
 *
 *   X is a 3-column matrix. Each row represents the Cartesian coordinates
 *   of a vertex on the surface, indexed by TRI values.
 *
 *   XI is a 3-column matrix. Each row represents the Carterian coordinates
 *   of a point for which we want to find whether it's inside or outside the
 *   closed surface. Note that if you want to test all the voxels in an
 *   image, it is very slow and memory intensive to generate coordinates for
 *   each voxel. In that scenario, it is much better to use the cell array
 *   CI syntax shown below, and provide only values for the coordinate axes.
 *
 *   ISIN is a boolean vector with one element per point in XI. True means
 *   that the corresponding point is inside the closed surface (or on the
 *   surface), and false means that it's outside.
 *
 * ISIN = cgal_insurftri(TRI, X, CI)
 *
 *   CI is a cell array CI={XI, YI, ZI}, where XI, YI and ZI are row vectors
 *   that describe a rectangular grid. For example,
 *
 *     CI={linspace(-.25, .25, 5), ...
 *         linspace(-.25, .25, 4), ...
 *         linspace(-.25, .25, 3)};
 *
 *   describes a sampling grid of 4 rows x 5 columns x 3 slices (note that
 *   rows correspond to YI and columns to XI), of a domain 
 *   [-0.25, 0.25] x [-0.25, 0.25] x [-0.25, 0.25].
 *
 * ISIN = cgal_insurftri(..., DIRECTIONS, TOL)
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
 *          DIRECTIONS=rand(5, 3);
 *
 *   Warning! For the voting system to make sense, select an odd number of
 *   rays.
 *
 *   TOL is a scalar with the distance tolerance. Points at distance <= TOL
 *   are considered to be on the surface, and thus "inside". By default,
 *   TOL=1e-15.
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012-2013 University of Oxford
  * Version: 0.4.1
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

typedef MatlabImportFilter::MatlabInputPointer    MatlabInputPointer;

/*
 * pointIsIn(): auxiliary function to test whether a point is inside
 * or outside the surface
 */
bool pointIsIn(Point xi, Tree &tree, 
	       std::vector<Direction> &direction, double tol) {
  
  // minimum distance from the test point to the surface
  double d = tree.squared_distance(xi);
  
  // if the test point is close enough to the surface, we consider
  // it part of the surface, and thus an inside point
  if (d <= tol) {
    return true;
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
  return isin_vote > (direction.size() / 2);
}

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  enum InputIndexType {IN_TRI, IN_X, IN_XI, IN_DIRECTIONS, IN_TOL, InputIndexType_MAX};
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);

  // check that we have at least tri, x and xi
  matlabImport->CheckNumberOfArguments(3, InputIndexType_MAX);

  // register the inputs for this function at the import filter
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer;
  MatlabInputPointer inTRI =        matlabImport->RegisterInput(IN_TRI, "TRI");
  MatlabInputPointer inX =          matlabImport->RegisterInput(IN_X, "X");
  MatlabInputPointer inXI =         matlabImport->RegisterInput(IN_XI, "XI");
  MatlabInputPointer inDIRECTIONS = matlabImport->RegisterInput(IN_DIRECTIONS, "DIRECTIONS");
  MatlabInputPointer inTOL =        matlabImport->RegisterInput(IN_TOL, "TOL");

  // interface to deal with outputs to Matlab
  enum OutputIndexType {OUT_ISIN, OutputIndexType_MAX};
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->ConnectToMatlabFunctionOutput(nlhs, plhs);

  // check number of outputs the user is asking for
  matlabExport->CheckNumberOfArguments(0, OutputIndexType_MAX);

  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outISIN = matlabExport->RegisterOutput(OUT_ISIN, "ISIN");

  // if any of the inputs is empty, the output is empty too
  if (mxIsEmpty(prhs[IN_TRI]) || mxIsEmpty(prhs[IN_X]) || mxIsEmpty(prhs[IN_XI])) {
    matlabExport->CopyEmptyArrayToMatlab(outISIN);
    return;
  }

  // if user provides the ray directions, read them
  std::vector<Direction> directionDef; // ray query directions
  directionDef.push_back(Direction(1.0, 0.0, 0.0));    // default directions if 
  directionDef.push_back(Direction(-1.0, 1.0, 1.0));   // not provided by the
  directionDef.push_back(Direction(-1.0, -1.0, -1.0)); // user

  // dimensions of the directions matrix
  std::vector<Direction> direction
    = matlabImport->ReadVectorOfVectorsFromMatlab<K::RT, Direction>(inDIRECTIONS, directionDef);

  // distance tolerance value
  double tol = matlabImport->ReadScalarFromMatlab<double>(inTOL, 1e-15);

  // point coordinates with NaN values in case there's a problem reading them
  Point def(mxGetNaN(), mxGetNaN(), mxGetNaN());

  // get size of input matrix
  mwSize nrowsTri = mxGetM(prhs[IN_TRI]);
  mwSize nrowsXi = mxGetM(prhs[IN_XI]);
  mwSize ncolsTri = mxGetN(prhs[IN_TRI]);
  mwSize ncolsX = mxGetN(prhs[IN_X]);
  mwSize ncolsXi = mxGetN(prhs[IN_XI]);
  if ((ncolsTri != 3) || (ncolsX != 3) || (ncolsXi != 3)) {
    mexErrMsgIdAndTxt("Gerardus:CgalInSurfaceTriangulation:WrongInputFormat", 
		      "All input arguments must have 3 columns");
  }

  // read triangular mesh from function
  std::list<Triangle> triangles;
  mwIndex v0, v1, v2; // indices of the 3 vertices of each triangle
  Point x0, x1, x2; // coordinates of the 3 vertices of each triangle
  for (mwIndex i = 0; i < nrowsTri; ++i) {

    // exit if user pressed Ctrl+C
    ctrlcCheckPoint(__FILE__, __LINE__);

    // get indices of the 3 vertices of each triangle. These indices
    // follow Matlab's convention v0 = 1, 2, ..., n
    v0 = matlabImport->ReadScalarFromMatlab<mwIndex>(inTRI, i, 0, mxGetNaN());
    v1 = matlabImport->ReadScalarFromMatlab<mwIndex>(inTRI, i, 1, mxGetNaN());
    v2 = matlabImport->ReadScalarFromMatlab<mwIndex>(inTRI, i, 2, mxGetNaN());
    if (mxIsNaN(v0) || mxIsNaN(v1) || mxIsNaN(v2)) {
      mexErrMsgIdAndTxt("Gerardus:CgalInSurfaceTriangulation:WrongInputFormat", 
			"Parameter TRI: Vertex index is NaN");
    }
    
    // get coordinates of the 3 vertices (substracting 1 so that
    // indices follow the C++ convention 0, 1, ..., n-1)
    x0 = matlabImport->ReadRowVectorFromMatlab<void, Point>(inX, v0 - 1, def);
    x1 = matlabImport->ReadRowVectorFromMatlab<void, Point>(inX, v1 - 1, def);
    x2 = matlabImport->ReadRowVectorFromMatlab<void, Point>(inX, v2 - 1, def);

    // add triangle to the list of triangles in the surface
    triangles.push_back(Triangle(x0, x1, x2));
  }

  // construct AABB tree
  Tree tree(triangles.begin(),triangles.end());

  // construct internal data structure to accelerate distance queries
  if (!tree.accelerate_distance_queries()) {
    mexErrMsgTxt("Not enough memory to accelerate distance queries");
  }

  if (mxIsCell(inXI->pm)) { // xi is given by 3 vectors that describe a
                            // rectangular volume we want to test

    // check that the cell contains three vectors and get pointers to them
    if (mxGetN(prhs[IN_XI]) != 3) {
      mexErrMsgTxt("CI must be a cell array given as a row with 3 elements");
    }
    
    mxArray *pXi = mxGetCell(inXI->pm, 0);
    mxArray *pYi = mxGetCell(inXI->pm, 1);
    mxArray *pZi = mxGetCell(inXI->pm, 2);
    if (pXi == NULL || pYi == NULL || pZi == NULL) {
      mexErrMsgTxt("Cannot get pointer to vectors inside cell array CI");
    }
    
    // if any of the vectors is empty, we return an empty output
    if (mxIsEmpty(pXi) || mxIsEmpty(pYi) || mxIsEmpty(pZi)) {
      matlabExport->CopyEmptyArrayToMatlab(outISIN);
      return;
    }

    // check that the cell array contains vectors, instead of matrices
    if (mxGetM(pXi) != 1) {
      mexErrMsgTxt("XI contained in CI must be a row vector");
    }
    if (mxGetM(pYi) != 1) {
      mexErrMsgTxt("YI contained in CI must be a row vector");
    }
    if (mxGetM(pZi) != 1) {
      mexErrMsgTxt("ZI contained in CI must be a row vector");
    }

    // register the three cell components as inputs
    MatlabInputPointer inCIXI = matlabImport->RegisterInput(pXi, "CI{XI}");
    MatlabInputPointer inCIYI = matlabImport->RegisterInput(pYi, "CI{YI}");
    MatlabInputPointer inCIZI = matlabImport->RegisterInput(pZi, "CI{ZI}");

    // get length of each vector in CI
    size_t lenXi = std::max(mxGetM(pXi), mxGetN(pXi));
    size_t lenYi = std::max(mxGetM(pYi), mxGetN(pYi));
    size_t lenZi = std::max(mxGetM(pZi), mxGetN(pZi));
    
    // initialise output (note that rows correspond to Y coordinates,
    // and columns correspond to X coordinates)
    std::vector<mwSize> size;
    size.push_back(lenYi);
    size.push_back(lenXi);
    size.push_back(lenZi);
    bool *isin = matlabExport->AllocateNDArrayInMatlab<bool>(outISIN, size);

    // loop every point that is tested to see whether it's inside or
    // outside the surface
    Point xi; // test point coordinates
    for (mwIndex s = 0; s < lenZi; ++s) { // slice (slowest varying)

      // z-coordinate of the point to be tested
      double xi_z = matlabImport->ReadScalarFromMatlab<double>(inCIZI, 0, s, mxGetNaN());
      for (mwIndex c = 0; c < lenXi; ++c) { // column

    	// x-coordinate of the point to be tested
    	double xi_x = matlabImport->ReadScalarFromMatlab<double>(inCIXI, 0, c, mxGetNaN());
    	for (mwIndex r = 0; r < lenYi; ++r) { // row (fastest varying)
	  
    	  // exit if user pressed Ctrl+C
    	  ctrlcCheckPoint(__FILE__, __LINE__);
    
    	  // y-coordinate of the point to be tested
    	  double xi_y = matlabImport->ReadScalarFromMatlab<double>(inCIYI, 0, r, mxGetNaN());

    	  // test whether point is inside or outside the surface
    	  isin[sub2ind(size[0], size[1], size[2], r, c, s)] 
    	    = pointIsIn(Point(xi_x, xi_y, xi_z), tree, direction, tol);
    	} // r
      } // c
    } // s

  } else { // each row of xi is a point to test
    
    // initialise output
    bool *isin = matlabExport->AllocateColumnVectorInMatlab<bool>(outISIN, nrowsXi);
    
    // loop every point that is tested to see whether it's inside or
    // outside the surface
    Point xi; // test point coordinates
    for (mwIndex i = 0; i < nrowsXi; ++i) {
      
      // exit if user pressed Ctrl+C
      ctrlcCheckPoint(__FILE__, __LINE__);

      // get point coordinates to be tested
      xi = matlabImport->ReadRowVectorFromMatlab<void, Point>(inXI, i, def);

      // test whether point is inside or outside the surface
      isin[i] = pointIsIn(xi, tree, direction, tol);
    }

  } // ENDELSE: each row of xi is a point to test
  
}

#endif /* CGALINSURFACETRIANGULATION */
