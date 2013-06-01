/* CgalClosestFacetToPoint.cpp
 *
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2013 University of Oxford
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

#ifndef CGALCLOSESTFACETTOPOINT
#define CGALCLOSESTFACETTOPOINT

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
typedef K::FT                                     FT;
typedef CGAL::Point_3<K>                          Point;
typedef K::Segment_3                              Segment;
typedef CGAL::Triangle_3<K>                       Triangle;
typedef std::list<Triangle>::iterator             Iterator;
typedef CGAL::AABB_triangle_primitive<K,Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive>           AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits>     Tree;
typedef Tree::Object_and_primitive_id             Object_and_primitive_id;
typedef Tree::Point_and_primitive_id              Point_and_primitive_id;

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->SetMatlabArgumentsPointer(nrhs, prhs);

  // check that we have at least a filter name and input image
  matlabImport->CheckNumberOfArguments(3, 3);

  // interface to deal with outputs to Matlab
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->SetMatlabArgumentsPointer(nlhs, plhs);

  // check number of outputs the user is asking for
  matlabExport->CheckNumberOfArguments(0, 3);

  // if any of the inputs is empty, the output is empty too
  if (mxIsEmpty(prhs[0]) || mxIsEmpty(prhs[1]) || mxIsEmpty(prhs[2])) {
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    return;
  }

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

  // faces lookup table
  std::map<Primitive, int> vLookup;

  // read triangular mesh from function
  std::list<Triangle> triangles;
  mwIndex v0, v1, v2; // indices of the 3 vertices of each triangle
  Point x0, x1, x2; // coordinates of the 3 vertices of each triangle
  for (mwIndex i = 0; i < nrowsTri; ++i) {

    // exit if user pressed Ctrl+C
    ctrlcCheckPoint(__FILE__, __LINE__);

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
    x0 = matlabImport->GetRowVectorArgument<double, Point>(1, v0 - 1, "X", def);
    x1 = matlabImport->GetRowVectorArgument<double, Point>(1, v1 - 1, "X", def);
    x2 = matlabImport->GetRowVectorArgument<double, Point>(1, v2 - 1, "X", def);

    // add triangle to the list of triangles in the surface
    triangles.push_back(Triangle(x0, x1, x2));
  }

  // construct AABB tree
  Tree tree(triangles.begin(),triangles.end());

  std::cout << "tree size = " << tree.size() << std::endl;////////////
  tree.id();

  // construct internal data structure to accelerate distance queries
  if (!tree.accelerate_distance_queries()) {
    mexErrMsgTxt("Not enough memory to accelerate distance queries");
  }

  // initialise outputs
  plhs[0] = mxCreateNumericMatrix(nrowsXi, 1, mxDOUBLE_CLASS, mxREAL);
  if (plhs[0] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output 0");
  }
  plhs[1] = mxCreateNumericMatrix(nrowsXi, 3, mxDOUBLE_CLASS, mxREAL);
  if (plhs[1] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output 1");
  }
  plhs[2] = mxCreateNumericMatrix(nrowsXi, 1, mxDOUBLE_CLASS, mxREAL);
  if (plhs[2] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output 2");
  }
    
  // pointer to the outputs
  double *f = (double *)mxGetData(plhs[0]);
  if (f == NULL) {
    mexErrMsgTxt("Cannot get pointer to allocated output 0");
  }
  double *p = (double *)mxGetData(plhs[1]);
  if (p == NULL) {
    mexErrMsgTxt("Cannot get pointer to allocated output 1");
  }
  double *d = (double *)mxGetData(plhs[2]);
  if (d == NULL) {
    mexErrMsgTxt("Cannot get pointer to allocated output 2");
  }
  
  // loop every point to compute its distance to, intersection with
  // and closest facet of the surface
  Point xi; // test point coordinates
  for (mwIndex i = 0; i < nrowsXi; ++i) {
    
    // exit if user pressed Ctrl+C
    ctrlcCheckPoint(__FILE__, __LINE__);
    
    // get point coordinates to be tested
    xi = matlabImport->GetRowVectorArgument<double, Point>(2, i, "XI", def);
    
    // computes distance from query
    d[i] = sqrt(tree.squared_distance(xi));
    std::cout << "distance: " << d[i] << std::endl;

    // computes closest point and closest facet
    Point_and_primitive_id pp = tree.closest_point_and_primitive(xi);

    // closest point on the surface to the testing point
    p[i] = pp.first[0];
    p[i + nrowsXi] = pp.first[1];
    p[i + 2*nrowsXi] = pp.first[2];

    // closest facet
    Primitive fh = pp.second;
    // std::cout << "primitive = " << fh[0] << std::endl;///////////
    
  }
  
  
}

#endif /* CGALCLOSESTFACETTOPOINT */
