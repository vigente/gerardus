/*
 * sparse_breakdown.cpp
 *
 * SPARSE_BREAKDOWN  Extract internal arrays from sparse matrix
 *
 * [IR, PR, JC] = SPARSE_BREAKDOWN(S)
 *
 *   S is a sparse matrix.
 *
 *   IR is the vector obtained with the C++ function mxGetIr().
 *   http://www.mathworks.com/access/helpdesk/help/techdoc/apiref/mxgetir.html
 *
 *   PR is the vector obtained with the C++ function mxGetPr().
 *   http://www.mathworks.com/access/helpdesk/help/techdoc/apiref/mxgetpr.html
 *
 *   JC is the vector obtained with the C++ function mxGetJc().
 *   http://www.mathworks.com/access/helpdesk/help/techdoc/apiref/mxgetjc.html
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.1.1
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


#include "mex.h"
#include "matrix.h"

#include <string.h>
#include <iostream>

// entry point for the MEX file
void mexFunction(int nlhs, // number of expected outputs
		 mxArray *plhs[], // array of pointers to outputs
		 int nrhs, // number of inputs
		 const mxArray *prhs[] // array of pointers to inputs
		 ) {

  // Syntax:
  //
  // [IR, PR, JC] = SPARSE_BREAKDOWN(S)

  // check arguments
  if (nrhs != 1) {
    mexErrMsgTxt("1 input argument required.");
  }
  if (nlhs > 3) {
    mexErrMsgTxt("Maximum of 3 output argument allowed.");
  }
  if (!mxIsSparse(prhs[0])) {
    mexErrMsgTxt("Input matrix must be sparse.");
  }

  // number of columns
  mwSize ncol = mxGetN(prhs[0]);

  // number of nonzero elements in the input sparse array
  mwSize nnz = *(mxGetJc(prhs[0]) + ncol);

  // create output arrays
  plhs[0] = mxCreateDoubleMatrix(nnz, 1, mxREAL);
  if (!plhs[0]) {
    mexErrMsgTxt("Not enough memory for output");
  }
  plhs[1] = mxCreateDoubleMatrix(nnz, 1, mxREAL);
  if (!plhs[1]) {
    mexErrMsgTxt("Not enough memory for output");
  }
  plhs[2] = mxCreateDoubleMatrix(ncol+1, 1, mxREAL);
  if (!plhs[2]) {
    mexErrMsgTxt("Not enough memory for output");
  }

  // pointers to outputs
  double *outIr = mxGetPr(plhs[0]);
  double *outPr = mxGetPr(plhs[1]);
  double *outJc = mxGetPr(plhs[2]);

  // pointers to arrays of input sparse matrix
  mwIndex *inIr = mxGetIr(prhs[0]);
  double *inPr = mxGetPr(prhs[0]);
  mwIndex *inJc = mxGetJc(prhs[0]);

  // populate the outputs
  memcpy(outPr, inPr, nnz*sizeof(double));
  for (mwIndex i = 0; i < nnz; ++i) {
    outIr[i] = (double)inIr[i];
  }
  for (mwIndex i = 0; i <= ncol; ++i) {
    outJc[i] = (double)inJc[i];
  }

}
