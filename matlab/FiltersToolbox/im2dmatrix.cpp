/*
 * im2dmatrix.cpp
 *
 * IM2DMATRIX  Sparse distance matrix between 2D or 3D image voxels, weighted
 * by voxel intensities
 *
 * A = im2dmatrix(IM)
 *
 *   IM is an image volume with dimensions (R, C, S).
 *
 *   A is a sparse matrix with dimensions (R*C*S, R*C*S), where element
 *   (i,j) is the mean intensity between voxels with linear indices i and j.
 *
 *   Voxels with an Inf intensity are skipped.
 *
 * ... = im2dmatrix(..., RES) [This option only available in the MEX version]
 *
 *   RES is a row vector with the voxel size of [row, column, slice] (2D) or
 *   [row, column, slice] (3D). By default, RES=[1.0 1.0 1.0].
 *
 *
 *   This function has a slow Matlab implementation (using loops), but a
 *   fast MEX version is provided with Gerardus too.
 *
 * See also: seg2dmat.
 */
/*
 * Author: Ramon Casero <rcasero@gmail.com>
 * Copyright © 2010-2011 University of Oxford
 * Version: 0.2.3
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "mex.h"
#include "matrix.h"

#include <iostream>
#include <vector>
#include <cmath>

#include "../GerardusCommon.h"

// entry point for the MEX file
void mexFunction(int nlhs, // number of expected outputs
		 mxArray *plhs[], // array of pointers to outputs
		 int nrhs, // number of inputs
		 const mxArray *prhs[] // array of pointers to inputs
		 ) {

  // variables
  mwSize R = 0, C = 0, S = 0; // number of rows, cols and slices of input image
  double *im; // pointer to the image volume
  double *out; // pointer to the output adjacency-distance sparse matrix

  // check arguments
  if ((nrhs < 1) || (nrhs > 2)) {
    mexErrMsgTxt("Wrong number of input arguments");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  // get image size
  const mwSize *sz = mxGetDimensions(prhs[0]);
  mwSize ndim = mxGetNumberOfDimensions(prhs[0]);
  if (ndim == 2) { // 2D image
    R = sz[0];
    C = sz[1];
    S = 1;
  } else if (ndim == 3) { // 3D image volume
    R = sz[0];
    C = sz[1];
    S = sz[2];
  } else {
    mexErrMsgTxt("Input argument has to be a 2D image or 3D image volume");
  }

  // defaults: voxel size (dR, dC, dS)
  std::vector<double> res; // voxel size
  if (nrhs < 2 || mxIsEmpty(prhs[1])) {
    res.push_back(1.0); // dR
    res.push_back(1.0); // dC
    res.push_back(1.0); // dS
  } else {
    if (!mxIsDouble(prhs[1]) || (mxGetM(prhs[1]) != 1) || (mxGetN(prhs[1]) != ndim)) {
      mexErrMsgTxt("Voxel size must be a row vector of class double with 1 element per image dimension");
    }
    double *pres = mxGetPr(prhs[1]);
    if (pres == NULL) {
      mexErrMsgTxt("Cannot get pointer to voxel size vector");
    }
    res.push_back(pres[0]); // dR
    res.push_back(pres[1]); // dC
    if (ndim == 2) { // 2D image
      res.push_back(1.0); // dS
    } else {
      res.push_back(pres[2]); // dS
    }
  }

  // pointer to input image, for convenience
  if (!mxIsDouble(prhs[0])) {
    mexErrMsgTxt("Input image array must be of type double");
  }
  im = mxGetPr(prhs[0]);

  // count number of non infinite voxels in the image
  mwSize nvox = 0;
  for (mwSize idx = 0; idx < R*C*S; ++idx) {
    if (!mxIsInf(im[idx])) {nvox++;}
  }

  // create sparse matrix for the output each voxel can be connected
  // to up to 26 voxels (imagine a (3,3,3)-cube of voxels, with our
  // voxel in the middle), so we allocate memory for that worse case
  // scenario
  plhs[0] = (mxArray *)mxCreateSparse(R*C*S, R*C*S, R*C*S*26, mxREAL);
  if (!plhs[0]) {mexErrMsgTxt("Not enough memory for output");}

  // pointer to the output matrix
  out = (double *)mxGetPr(plhs[0]);

  // pointer to the values in the sparse matrix (distance weighted by
  // mean voxel intensity)
  double *pr = mxGetPr(plhs[0]);
  if (!pr) {
    mexErrMsgTxt("Error loading vector pr from sparse matrix");
  }

  // pointer to vector jc in the sparse matrix (cumsum of number of
  // voxels connected to current voxel). If the sparse matrix has n
  // columns, then jc has n+1 elements. jc[i]=X means that in columns
  // 1 to i-1 there are a total of X non-zero elements. For example,
  // if column 0 has 3 nonzero entries and column 1 has 2, then jc=[0
  // 3 5]
  mwIndex *jc = mxGetJc(plhs[0]);
  if (!jc) {
    mexErrMsgTxt("Error loading vector jc from sparse matrix");
  }
  
  // pointer to the connected voxel indices
  //
  // in this case, indices follow the C++ convention
  //
  // ir is a vector that says where the non-zero entries are. For
  // example,
  //
  // s = 
  //   0     0     0
  //   1     0     2
  //   0     1     0
  //   0     0     0
  //   1     0     1
  //   0     0     1
  //   0     0     0
  //
  // ir = [1 4, 2, 1 4 5]
  mwIndex *ir = mxGetIr(plhs[0]);
  if (!ir) {
    mexErrMsgTxt("Error loading vector ir from sparse matrix");
  }

  // init jc so that we can add 1 every time we find a node connection
  for (mwSize idx = 0; idx < R*C*S+1; ++idx) {
    jc[idx] = 0;
  }
  
  // loop voxels searching for voxels connected to them
  //
  // here the indices follow the C convention, i.e. first index is 0,
  // and will be converted to Matlab's convention (first index is 1)
  // when necessary
  //
  // because of the way Jc is defined, we are going to store the
  // connectivity value for the first node in the second element of
  // jc[], for the 2nd node in the 3rd element, and so on
  mwSize outidx = 0; // index for ir, pr
  mwSize nedg = 0; // number of edges or connections between voxels
  mwSize idx = 0; // linear index for current voxel
  mwSize nnidx = 0; // linear index for adjacent voxels
  mwSize RC = R*C; // aux variable
  double ledge = 0.0; // length of edge linking two voxels
  for (mwSize s = 0; s < S; ++s) {
    for (mwSize c = 0; c < C; ++c) {
      for (mwSize r = 0; r < R; ++r) {

	// exit if user pressed Ctrl+C
	ctrlcCheckPoint(__FILE__, __LINE__);
	
	std::cout << "idx = " << idx << std::endl;////////////////

	// if current voxels is Inf, we don't include it in the
	// graph, so just skip to next iteration
	if (mxIsInf(im[idx])) {++idx; continue;}

	// init neighbour index
	nnidx = 0;

	// examine voxels surrounding the current voxel (up to 26
	// voxels); if a neighbour voxel is noninfinite, we say that
	// it is connected to the current voxel, and has to be added
	// to the sparse matrix
	for (mwSize nns = std::max((long int)0, (long int)s-1);
	     nns <= std::min(S-1, s+1); ++nns) {
	  for (mwSize nnc = std::max((long int)0, (long int)c-1);
	       nnc <= std::min(C-1, c+1); ++nnc) {
	    for (mwSize nnr = std::max((long int)0, (long int)r-1);
		 nnr <= std::min(R-1, r+1); ++nnr) {
	      
	      // don't connect current voxel to itself
	      if (nns == s && nnc == c && nnr == r) {continue;}
	      
	      // get linear index of neighbour voxel
	      nnidx = RC*nns + R*nnc + nnr;

	      // skip connected voxels that are Inf
	      if (mxIsInf(im[nnidx])) {continue;}

	      std::cout << "\tnnidx = " << nnidx << std::endl;////////////////

	      // length of edge linking two voxels (sqrt(dx^2 + dy^2 + dz^2))
	      ledge = abs(nnr - r) * res[0] * abs(nnr - r) * res[0];
	      ledge += abs(nnc - c) * res[1] * abs(nnc - c) * res[1];
	      ledge += abs(nns - s) * res[2] * abs(nns - s) * res[2];
	      ledge = sqrt(ledge);

	      std::cout << "\t\tledge = " << ledge << std::endl;////////////////

	      // the edge weight between the current voxel and the
	      // connected voxel is the mean between their node values
	      //
	      // this weight multiplies the connection length
	      pr[outidx] = (im[nnidx] + im[idx]) * 0.5 * ledge;

	      // instead of computing jc directly as the cumulative
	      // sum, for the moment we are going to keep track only
	      // of how many elements are in each column, i.e. how
	      // many neighbours the current voxel has (and note that
	      // jc[i+1] corresponds to column i)
	      //
	      // here we say that the voxel has another neighbour
	      jc[idx+1] += 1;

	      // row position of the neighbour (C++ convention)
	      ir[outidx] = nnidx;

	      // advance one position in ir, pr
	      ++outidx;

	    }
	  }
	}
	
	// this counter keeps track of the total number of nonzero
	// elements in the distance matrix
	nedg = outidx;

	// update current voxel index
	++idx;

      }
    }
  }
				
  // now vector jc contains the number of voxels connected to each
  // voxel (e.g. [0 4 1 0 2]). However, we need instead the
  // accumulated vector (e.g. [0 4 5 0 7])
  for (mwSize idx = 1; idx <= R*C*S; ++idx) {
    jc[idx] += jc[idx-1];
  }
  
  // correct the number of non-zero entries in the sparse matrix
  mxSetNzmax(plhs[0], nedg);
  
};
