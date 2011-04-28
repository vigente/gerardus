/*
 * bwregiongrow_aux.cpp
 *
 * This is an auxiliary function to bwregiongrow.m, to implement part
 * of the algorithm in C++ so that we can avoid a costly loop in Matlab
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.1.0
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

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

/* mex headers */
#include <math.h>
#include <matrix.h>
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <limits>

/* run(): function in charge of processing. We cannot do this in
 *              the body of mexFunction() because we need to template
 *              the function so that we can operate with different
 *              types of Matlab matrices (double, uint8, etc)
 *
 */
template <class VoxelType>
void run(mxArray* &im,
	 const mxArray* &idxtodo,
	 const mxArray* &d,
	 const mxArray* &idict,
	 const mxArray* &TODO,
	 mxArray* &new_labels) {
  
  // get number of candidate voxels
  mwIndex N = mxGetNumberOfElements(idxtodo);

  // // DEBUG
  // std::cout << "Number of voxels = " << N << std::endl; 

  // create empty vector for output
  if (N == 0) {
    new_labels = mxCreateDoubleMatrix(0, 0, mxREAL);
    return;
  }

  // dealing with sparse array
  double *dpr = (double *)mxGetPr(d);
  mwIndex *dir   = mxGetIr(d);
  mwIndex *djc   = mxGetJc(d);

  // pointer to input arrays
  double *idxtodop = (double *)mxGetPr(idxtodo);
  double *idictp = (double *)mxGetPr(idict);
  VoxelType *imp = (VoxelType *)mxGetPr(im);
  VoxelType *TODOp = (VoxelType *)mxGetPr(TODO);

  // create output matrix for Matlab's result
  const mwSize dims[] = {1, N};
  new_labels = (mxArray *)mxCreateNumericArray(2, dims,
					       mxGetClassID(im),
					       mxREAL);
  if (new_labels == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output matrix");
  }
  VoxelType *new_labelsp =  (VoxelType *)mxGetPr(new_labels);

  // loop every candidate voxel
  for (mwIndex i = 0; i < N; ++i) {
    // get d matrix indices of voxels adjacent to the current candidate
    
    // vector pr contains all elements in the sparse matrix. We want
    // to know from which index to which index of pr are the d-matrix
    // values for this candidate
    //
    // the values are the distances. The corresponding values in dir
    // are the voxel indices
    mwIndex fromidx = djc[(mwIndex)idxtodop[i]-1];
    mwIndex toidx = djc[(mwIndex)idxtodop[i]] - 1;

    // initialize variables to identify the closest labelled voxel to
    // our candidate
    double dmin = std::numeric_limits<double>::max();
    new_labelsp[i] = 0;

    // loop adjacent voxels to candidate
    // // DEBUG
    // std::cout << "Candidate: " << idxtodop[i] << std::endl;
    for (mwIndex idx = fromidx; idx <= toidx; ++idx) {
      // // DEBUG
      // std::cout << "adjacentOf[" << idxtodop[i] << "] = " << dir[idx]+1
      // 		<< ", d = " << (double)dpr[idx] 
      // 		<< ", voxel = " << idictp[dir[idx]]
      // 		<< ", imlab = " 
      // 		<< (double)imp[(mwIndex)idictp[dir[idx]]-1]
      // 		<< std::endl;
      // we are looking for a voxel to copy its label, so ignore all
      // unlabelled voxels

      if ((double)imp[(mwIndex)idictp[dir[idx]]-1] != TODOp[0]) {
	// DEBUG
	// std::cout << "Tagged already" << std::endl;

	// if this neighbour is closer than what we had before, then
	// we want its label
	if (dpr[idx] < dmin) {
	  dmin = dpr[idx];
	  new_labelsp[i] = imp[(mwIndex)idictp[dir[idx]]-1];
	}
      }

      // NOTE!!! 
      // idictp starts at 1
      // idxtodop starts at 1
      // i starts at 0
      // dir[idx] starts at 0
      // dpr[idx] starts at 0

    }
  }

  // once the new labels are computed, add them to the actual
  // segmentation (it's important to wait until now, because
  // if we change the labels while still in the previous loop, we get
  // spill overs of some regions into other regions
  for (mwIndex i = 0; i < N; ++i) {
    imp[(mwIndex)idictp[(mwIndex)idxtodop[i]-1]-1] = new_labelsp[i];
  }  

  return;
}

// entry point for the mex function
//   prhs[0]: (in) im
//   prhs[1]: (in) idxtodo
//   prhs[2]: (in) d (double sparse)
//   prhs[3]: (in) idict
//   prhs[4]: (in) TODO
//   plhs[0]: (out) new_labels
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {
  // check number of input and output arguments
  if (nrhs != 5) {
    mexErrMsgTxt("Five input arguments required");
  }
  else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  // check arguments
  if (!mxIsDouble(prhs[1])) {
    mexErrMsgTxt("idxtodo must be of type double");
  }
  if (!mxIsSparse(prhs[2])) {
    mexErrMsgTxt("d is not a valid sparse matrix");
  }
  if (!mxIsDouble(prhs[2])) {
    mexErrMsgTxt("d must be of type double");
  }
  if (!mxIsDouble(prhs[3])) {
    mexErrMsgTxt("idict must be of type double");
  }
  // get input image class (double, uint8, etc)
  mxClassID  imClassId  = mxGetClassID(prhs[0]);
  if (imClassId != mxGetClassID(prhs[4])) {
    mexErrMsgTxt("im and TODO must be the same type");
  }

  // run function, templated according to the input matrix type
  switch(imClassId)  {
  case mxLOGICAL_CLASS:
    run<bool>(const_cast<mxArray*&>(prhs[0]), 
	      prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
    break;
  case mxDOUBLE_CLASS:
    run<double>(const_cast<mxArray*&>(prhs[0]), 
		prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
    break;
  case mxSINGLE_CLASS:
    run<float>(const_cast<mxArray*&>(prhs[0]), 
	       prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
    break;
  case mxINT8_CLASS:
    run<int8_T>(const_cast<mxArray*&>(prhs[0]), 
		prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
    break;
  case mxUINT8_CLASS:
    run<uint8_T>(const_cast<mxArray*&>(prhs[0]), 
		 prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
    break;
  case mxINT16_CLASS:
    run<int16_T>(const_cast<mxArray*&>(prhs[0]), 
		 prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
    break;
  case mxUINT16_CLASS:
    run<uint16_T>(const_cast<mxArray*&>(prhs[0]), 
		  prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
    break;
  case mxINT32_CLASS:
    run<int32_T>(const_cast<mxArray*&>(prhs[0]), 
		 prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
    break;
  // case mxUINT32_CLASS:
  // run<uint32_T>(const_cast<mxArray*&>(prhs[0]), 
  // 		prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
  //   break;
  case mxINT64_CLASS:
    run<int64_T>(const_cast<mxArray*&>(prhs[0]), 
		 prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
    break;
  // case mxUINT64_CLASS:
  //   run<uint64_T>(const_cast<mxArray*&>(prhs[0]), 
  // 		  prhs[1], prhs[2], prhs[3], prhs[4], plhs[0]);
  //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt("Input matrix has unknown type.");
    break;
  default:
    mexErrMsgTxt("Input matrix has invalid type.");
    break;
  }

  // exit successfully
  return;

}
