/* CgalTriFacetArea.cpp
 *
 * CGAL_TRIFACET_AREA  Area of the facets in a triangular mesh
 *
 * A = cgal_trifacet_area(TRI, X)
 *
 *   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
 *   triangular facet in the mesh.
 *
 *   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of the
 *   i-th node in the mesh.
 *
 *   A is a vector with the same number of rows as TRI. A(i) is the area of
 *   the triangle TRI(i).
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

#ifndef CGALTRIFACETAREA
#define CGALTRIFACETAREA

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>

/* Gerardus headers */
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double>            K;
typedef K::FT                                     FT;
typedef CGAL::Point_3<K>                          Point;
typedef K::Segment_3                              Segment;
typedef CGAL::Triangle_3<K>                       Triangle; // size 72 byte

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->SetMatlabArgumentsPointer(nrhs, prhs);

  // check that we have at least a filter name and input image
  matlabImport->CheckNumberOfArguments(2, 2);

  // interface to deal with outputs to Matlab
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->SetMatlabArgumentsPointer(nlhs, plhs);

  // check number of outputs the user is asking for
  matlabExport->CheckNumberOfArguments(0, 1);

  // default coordinates are NaN values, so that the user can spot
  // whether there was any problem reading them
  Point def(mxGetNaN(), mxGetNaN(), mxGetNaN());

  // if any of the inputs is empty, the output is empty too
  if (mxIsEmpty(prhs[0]) || mxIsEmpty(prhs[1])) {
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    return;
  }

  // get size of input matrix
  mwSize nrowsTri = mxGetM(prhs[0]);
  mwSize ncolsTri = mxGetN(prhs[0]);
  mwSize ncolsX = mxGetN(prhs[1]);
  if ((ncolsTri != 3) || (ncolsX != 3)) {
    mexErrMsgTxt("Both input arguments must have 3 columns");
  }

  // initialise output
  plhs[0] = mxCreateNumericMatrix(nrowsTri, 1, mxDOUBLE_CLASS, mxREAL);
  if (plhs[0] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output 0");
  }
    
  // pointer to the outputs
  double *area = (double *)mxGetData(plhs[0]);
  if (area == NULL) {
    mexErrMsgTxt("Cannot get pointer to allocated output 0");
  }
  
  // read triangular mesh from function
  Triangle tri;
  mwIndex v0, v1, v2; // indices of the 3 vertices of each triangle
  Point x0, x1, x2; // coordinates of the 3 vertices of each triangle

  for (mwIndex i = 0; i < nrowsTri; ++i) {

    // exit if user pressed Ctrl+C
    ctrlcCheckPoint(__FILE__, __LINE__);

    // get indices of the 3 vertices of each triangle. These indices
    // follow Matlab's convention v0 = 1, 2, ..., n
    v0 = matlabImport->GetScalarArgument<mwIndex>(0, i, 0, "TRI0", mxGetNaN());
    v1 = matlabImport->GetScalarArgument<mwIndex>(0, i, 1, "TRI1", mxGetNaN());
    v2 = matlabImport->GetScalarArgument<mwIndex>(0, i, 2, "TRI2", mxGetNaN());
    if (mxIsNaN(v0) || mxIsNaN(v1) || mxIsNaN(v2)) {
      mexErrMsgTxt("Parameter TRI: Vertex index is NaN");
    }
    
    // get coordinates of the 3 vertices (substracting 1 so that
    // indices follow the C++ convention 0, 1, ..., n-1)
    x0 = matlabImport->GetRowVectorArgument<double, Point>(1, v0 - 1, "X0", def);
    x1 = matlabImport->GetRowVectorArgument<double, Point>(1, v1 - 1, "X1", def);
    x2 = matlabImport->GetRowVectorArgument<double, Point>(1, v2 - 1, "X2", def);

    // create triangle from the vertices read at the input
    tri = Triangle(x0, x1, x2);

    // compute triangle area
    area[i] = std::sqrt(tri.squared_area());
    
    // // DEBUG
    // std::cout << x0 << "\t" << x1 << "\t" << x2 << "\t" << std::sqrt(tri.squared_area()) << std::endl;
  }


  
}

#endif /* CGALTRIFACETAREA */
