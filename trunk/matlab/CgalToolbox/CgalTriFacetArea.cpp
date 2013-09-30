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
  * Version: 0.3.1
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
  enum InputIndexType {IN_TRI, IN_X, InputIndexType_MAX};
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);

  // check that we have at least a filter name and input image
  matlabImport->CheckNumberOfArguments(2, InputIndexType_MAX);

  // register the inputs for this function at the import filter
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer;
  MatlabInputPointer inTRI = matlabImport->RegisterInput(IN_TRI, "TRI");
  MatlabInputPointer inX = matlabImport->RegisterInput(IN_X, "X");

  // interface to deal with outputs to Matlab
  enum OutputIndexType {OUT_A, OutputIndexType_MAX};
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->ConnectToMatlabFunctionOutput(nlhs, plhs);

  // check number of outputs the user is asking for
  matlabExport->CheckNumberOfArguments(0, OutputIndexType_MAX);

  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outA = matlabExport->RegisterOutput(OUT_A, "A");

  // default coordinates are NaN values, so that the user can spot
  // whether there was any problem reading them
  Point def(mxGetNaN(), mxGetNaN(), mxGetNaN());

  // if any of the inputs is empty, the output is empty too
  if (mxIsEmpty(prhs[IN_TRI]) || mxIsEmpty(prhs[IN_X])) {
    matlabExport->CopyEmptyArrayToMatlab(outA);
    return;
  }

  // get size of input matrix
  mwSize nrowsTri = mxGetM(prhs[IN_TRI]);
  mwSize ncolsTri = mxGetN(prhs[IN_TRI]);
  mwSize ncolsX = mxGetN(prhs[IN_X]);
  if ((ncolsTri != 3) || (ncolsX != 3)) {
    mexErrMsgTxt("Both input arguments must have 3 columns");
  }

  // initialise output
  double *area = matlabExport->AllocateColumnVectorInMatlab<double>(outA, nrowsTri);
  
  // read triangular mesh from function
  Triangle tri;
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
      mexErrMsgTxt("Input TRI: Vertex index is NaN");
    }
    
    // get coordinates of the 3 vertices (substracting 1 so that
    // indices follow the C++ convention 0, 1, ..., n-1)
    x0 = matlabImport->ReadRowVectorFromMatlab<void, Point>(inX, v0 - 1, def);
    x1 = matlabImport->ReadRowVectorFromMatlab<void, Point>(inX, v1 - 1, def);
    x2 = matlabImport->ReadRowVectorFromMatlab<void, Point>(inX, v2 - 1, def);

    // create triangle from the vertices read at the input
    tri = Triangle(x0, x1, x2);

    // compute triangle area
    area[i] = std::sqrt(tri.squared_area());
    
    // // DEBUG
    // std::cout << x0 << "\t" << x1 << "\t" << x2 << "\t" << std::sqrt(tri.squared_area()) << std::endl;
  }


  
}

#endif /* CGALTRIFACETAREA */
