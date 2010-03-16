#include "mex.h"
#include "matrix.h"

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
