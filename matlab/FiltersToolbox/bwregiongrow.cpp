/*
 * bwregiongrow.cpp
 *
 * BWREGIONGROW  Fast region grow labelling of binary image from multiple
 * seeds
 *
 * LAB = BWREGIONGROW(IM, TODO)
 *
 *   IM is a 2D matrix or 3D array with a multi-label segmentation. That is,
 *   IM contains a background and several objects, each object represented
 *   by all the connected voxels with the same label.
 *
 *   IM can have any Matlab numeric type (double, uint8, etc).
 *
 *   Numerical voxel values in IM are interpreted in the following way:
 *
 *     0:               background voxel, don't label
 *     TODO:            this voxel needs to be labelled using the region
 *                      grow algorithm
 *     Any other value: this voxel is a seed, and its value will be
 *                      propagated using the region grow algorithm
 *
 *   LAB has the same size as IM, and has the label values that partition
 *   the original binary image IM. Each partition has a different label.
 *   Partitions are computed with a region grow algorithm that expands the
 *   labels from the seeds.
 *
 *   At each iteration, partitions grow 1 voxel until the whole IM is
 *   labelled.
 *
 * LAB = BWREGIONGROW(..., RES)
 *
 *   RES is a 2-vector (in 2D) or 3-vector (in 3D) with the voxel size in
 *   each dimension. By default, it is assumed that RES=[1, 1, 1]. Voxel
 *   size is used to compute distances between voxels in the labelling
 *   process.
 *
 *   MAXITER is a scalar to tell the algorithm to stop after a number of
 *   region grow iterations. If MAXITER < 0, the algorithm iterates until
 *   all TODO voxels have been labelled. By default, MAXITER = -1.
*/

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.2.1
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
#include <vector>
#include <stdlib.h>

/*
 * sub2ind(): function that converts r, c, s indices to linear indices
 *            in a 3D array (same as Matlab's function sub2ind(),
 *            although in Matlab indices start at 1, and in C++, they
 *            start at 0)
 *
 */
mwIndex sub2ind(mwSize R, mwSize C, mwSize S,
		std::vector<mwIndex> rcs) {
  // check for out of range index
  if (
      (rcs[0] < 0) || (rcs[0] >= R) 
      || (rcs[1] < 0) || (rcs[1] >= C)
      || (rcs[2] < 0) || (rcs[2] >= S)
      ) {
    mexErrMsgTxt("Out of range index");
  }
  if ((R*C*S == 0) || (R < 0) || (C < 0) || (S < 0)) {
    mexErrMsgTxt("Size values cannot be 0 or negative");
  }

  // check that input vector has 3 elements
  if (rcs.size() != 3) {
    mexErrMsgTxt("Input vector must have 3 elements");
  }

  // convert r, c, s to linear index
  mwIndex idx = rcs[0] + rcs[1] * R + rcs[2] * R * C;

  // // DEBUG
  // std::cout << "idx = " << idx
  // 	    << std::endl;
  
  return idx;
}

/*
 * ind2sub(): function that converts linear indices in a 3D array to
 *            r, c, s indices (same as Matlab's function ind2sub(),
 *            although in Matlab indices start at 1, and in C++, they
 *            start at 0)
 *
 */
std::vector<mwIndex> ind2sub(mwSize R, mwSize C, mwSize S,
			     mwIndex idx) {
  // check for out of range index
  if (idx >= R*C*S || idx < 0) {
    mexErrMsgTxt("Out of range index");
  }
  if (R*C*S == 0 || R < 0 || C < 0 || S < 0) {
    mexErrMsgTxt("Size values cannot be 0 or negative");
  }

  // init output
  std::vector<mwIndex> rcs(3);
  
  // convert linear index to r, c, s 
  div_t divresult;

  divresult = div(idx, R*C);
  rcs[2] = divresult.quot; // slice value
  idx = divresult.rem;

  divresult = div(idx, R);
  rcs[1] = divresult.quot; // column value
  idx = divresult.rem;

  divresult = div(idx, 1);
  rcs[0] = divresult.quot; // row value

  // // DEBUG
  // std::cout << "rcs = " << rcs[0] << ", " 
  // 	    << rcs[1] << ", "
  // 	    << rcs[2]
  // 	    << std::endl;
  
  return rcs;
}

/*
 * getNeighbours(): function that gets a vector of indices that
 *                  correspond to the voxels adjacent to the current
 *                  one
 */
std::vector<mwIndex> getNeighbours(mwSize R, mwSize C, mwSize S,
				   mwIndex idx) {
  
  // output vector
  std::vector<mwIndex> neighbour;

  // convert linear index to image array indices
  std::vector<mwIndex> rcs = ind2sub(R, C, S, idx);
  
  std::vector<mwIndex> rcsn(3);
  if (S == 1) { // 2D image
    
    for (mwIndex c = rcs[1]-1; c <= rcs[1]+1; ++c) {
      for (mwIndex r = rcs[0]-1; r <= rcs[0]+1; ++r) {
	// input index cannot be an index of itself
	if ((r == rcs[0]) && (c == rcs[1])) {
	  continue;
	}

	// is the potential neighbour within the image bounds?
	if ((r >= 0) && (r < R) && (c >= 0) && (c < C)) {
	  // convert array indices to linear index
	  rcsn[0] = r;
	  rcsn[1] = c;
	  rcsn[2] = 0;
	  idx = sub2ind(R, C, 1, rcsn);

	  // add the neighbour to the output
	  neighbour.push_back(idx);
	}
      }
    }
    
  } else { // 3D image
    
    for (mwIndex s = rcs[2]-1; s <= rcs[2]+1; ++s) {
      for (mwIndex c = rcs[1]-1; c <= rcs[1]+1; ++c) {
	for (mwIndex r = rcs[0]-1; r <= rcs[0]+1; ++r) {
	  // input index cannot be an index of itself
	  if ((r == rcs[0]) && (c == rcs[1]) && (s == rcs[2])) {
	    continue;
	  }
	  
	  // is the potential neighbour within the image bounds?
	  if ((r >= 0) && (r < R) 
	      && (c >= 0) && (c < C)
	      && (s >= 0) && (s < S)) {
	    // convert array indices to linear index
	    rcsn[0] = r;
	    rcsn[1] = c;
	    rcsn[2] = s;
	    idx = sub2ind(R, C, S, rcsn);
	    
	    // add the neighbour to the output
	    neighbour.push_back(idx);
	  }
	}
      }
    }

  }

  return neighbour;
  
}

/* run(): function in charge of processing. We cannot do this in
 *              the body of mexFunction() because we need to template
 *              the function so that we can operate with different
 *              types of Matlab matrices (double, uint8, etc)
 *
 */
template <class VoxelType>
void run(mxArray* &im, const mxArray* _TODO,
	 const mxArray* _res, const mwIndex _maxiter=-1) {

  // local variables
  mwIndex maxiter = _maxiter;

  // full array pointers
  VoxelType *imp = (VoxelType *)mxGetPr(im);
  VoxelType *TODOp = (VoxelType *)mxGetPr(_TODO);

  // get image size
  const mwSize *dims = mxGetDimensions(im);

  mwSize S = 1; // number of slices in the image
  if (mxGetNumberOfDimensions(im) == 2) {
    S = 1;
  } else if (mxGetNumberOfDimensions(im) == 3) {
    S = dims[2];
  } else {
    mexErrMsgTxt("Input image must be 2D or 3D");
  }
  mwSize R = dims[0]; // number of rows in the image
  mwSize C = dims[1]; // number of columns in the image
  mwSize Nim = R*C*S; // number of voxels in the image

  // value of the todo label
  VoxelType TODO = TODOp[0];

  // get resolution
  std::vector<double> res;
  if (mxIsEmpty(_res)) { // make resolution = [1, 1] or [1, 1, 1]
    res.assign(mxGetNumberOfDimensions(im), 1);
  } else { // check and extract resolution
    if (!mxIsDouble(_res)) {
      mexErrMsgTxt("RES must be of type double");
    }
    double *resp = mxGetPr(_res);

    if ((mxGetM(_res) == 1) && (mxGetN(_res) == 2)) {
      res.push_back(resp[0]);
      res.push_back(resp[1]);
      res.push_back(0.0);
    } else if ((mxGetM(_res) == 1) && (mxGetN(_res) == 3)) {
      res.push_back(resp[0]);
      res.push_back(resp[1]);
      res.push_back(resp[2]);
    } else {
      mexErrMsgTxt("RES must be a row 2-vector or a row 3-vector ");
    }
  }

  // indices of current boundary voxels
  std::vector<mwIndex> boundary;

  // indices of neighbours to a boundary voxel
  std::vector<mwIndex> neighbour;

  // indices of new boundary voxels
  std::vector<mwIndex> newBoundary;

  // indices of labels for new boundary voxels
  std::vector<VoxelType> newLabel;

  // flags to say whether a voxel has already been added to the
  // boundary or not
  std::vector<bool> addedToNewBoundary(Nim, false);

  /*
   * Algorithm initialization
   *
   */

  // number of voxels that need to be labelled, not including the seed
  // voxels
  mwSize Ntodo = 0;

  // initialize the boundary with the seeds
  for (mwIndex i = 0; i < Nim; ++i) {

    // does this voxel need to be labelled?
    if (imp[i] == TODO) {
      Ntodo++;
    }

    // is this voxel a seed voxel?
    addedToNewBoundary[i] = (imp[i] != 0) && (imp[i] != TODO);
    if (addedToNewBoundary[i]) {
      // if so, add it to the list of boundary voxels
      boundary.push_back(i);
    }
  }

  // loop until all voxels have been labelled or until the maximum
  // number of iterations
  while (maxiter != 0) {

    // decrease maxiter counter if we have to stop after a maximum
    // number of iterations
    if (maxiter > 0) {
      maxiter--;
    }
    
    /*
     * expand the boundary by 1 voxel at every iteration
     *
     */
  
    // start with an empty new boundary
    newBoundary.clear();
    for (mwIndex i = 0; i < addedToNewBoundary.size(); ++i) {
      addedToNewBoundary[i] = false;
    }

    // loop every boundary voxel
    for (mwIndex i = 0; i < boundary.size(); ++i) {

      // get neighbours of current voxel
      neighbour = getNeighbours(R, C, S, boundary[i]);
     
      // loop every neighbour
      for (mwIndex j = 0; j < neighbour.size(); ++j) {

	// ignore any neighbour that either belongs to the background
	// or that is already labelled or that has already been added
	// to the new boundary
	if ((imp[neighbour[j]] == TODO)
	    && !addedToNewBoundary[neighbour[j]]) {
	  // add this neighbour to the new boundary
	  newBoundary.push_back(neighbour[j]);
	  addedToNewBoundary[neighbour[j]] = true;
	}
	
      }
    }

    // if the algorithm has labelled all voxels already, terminate
    if (newBoundary.size() == 0) {
      return;
    }

    /*
     * find labels for the voxels in the new boundary. For each new
     * voxel, we need to find the closest voxel in the old boundary,
     * and copy its label
     *
     */
    
    // create a label vector of the same size as the new boundary. As
    // we are going to overwrite every label value, we don't need to
    // clear the vector first
    newLabel.resize(newBoundary.size());

    // squared distance between new boundary voxel and its closest
    // neighbour
    double d2min;

    // array indices of boundary voxel and neighbour
    std::vector<mwIndex> rcs(3), rcsn(3);

    // loop each new boundary voxel
    for (mwIndex i = 0; i < newBoundary.size(); ++i) {
      
      // initialize minimum distance to infinity
      d2min = std::numeric_limits<double>::max();

      // linear index to array indices
      rcs = ind2sub(R, C, S, newBoundary[i]);

      // if this voxel has a label already, we skip it
      if (imp[newBoundary[i]] != TODO) {
	continue;
      }

      // find neighbours
      neighbour = getNeighbours(R, C, S, newBoundary[i]);

      // loop every neighbour
      for (mwIndex j = 0; j < neighbour.size(); ++j) {
	// if the neighbour belongs to the background or doesn't have
	// a label, skip it
	if ((imp[neighbour[j]] != 0) 
	    && (imp[neighbour[j]] != TODO)) {

	  // linear index to array indices
	  rcsn = ind2sub(R, C, S, neighbour[j]);

	  // compute squared distance between boundary voxel and this
	  // neighbour, reusing rcsn
	  rcsn[0] -= rcs[0];
	  rcsn[1] -= rcs[1];
	  rcsn[2] -= rcs[2];
	  rcsn[0] *= res[0];
	  rcsn[1] *= res[1];
	  rcsn[2] *= res[2];
	  double d2 = rcsn[0]*rcsn[0] + rcsn[1]*rcsn[1] 
	    + rcsn[2]*rcsn[2];

	  // if this neighbour is closer than what we had before, it
	  // becomes the closest neighbour and we copy its label
	  if (d2 < d2min) {
	    d2min = d2;
	    newLabel[i] = imp[neighbour[j]];
	  }

	}
      }
      
    }

    /*
     * transfer the new labels to the image. It's important that we do
     * it now, and not in the previous step of the algorithm,
     * otherwise we can get labels spilling over to other regions
     *
     */
    
    // loop each new label
    for (mwIndex i = 0; i < newLabel.size(); ++i) {
      imp[newBoundary[i]] = newLabel[i];
    }

    // make the just labelled boundary the starting boundary for the
    // next step of the algorithm
    boundary = newBoundary;

  }
  
  return;
}

// entry point for the mex function
//   prhs[0]: (in) im: input image
//   prhs[1]: (in) TODO: label of the voxels that need to be labelled
//   prhs[2]: (in) res: 3-vector with resolution values
//   prhs[3]: (in) maxiter: maximum number of iterations
//   plhs[0]: (out) lab: labelled input image
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {
  // check number of input and output arguments
  if ((nrhs < 2) || (nrhs > 4)) {
    mexErrMsgTxt("Two or four input arguments required");
  }
  else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  // get input image class (double, uint8, etc)
  mxClassID  imClassId  = mxGetClassID(prhs[0]);
  if (imClassId != mxGetClassID(prhs[1])) {
    mexErrMsgTxt("im and TODO must be the same type");
  }

  // duplicate input image. This is going to be the output data, and
  // where the processing will happen
  plhs[0] = mxDuplicateArray(prhs[0]);

  // if res was not provided, or is empty, create an empty array to
  // pass it to run()
  mxArray *res = NULL;
  if (nrhs < 3 || mxIsEmpty(prhs[2])) {
    res = mxCreateDoubleMatrix(0, 0, mxREAL);
  } else {
    res = const_cast<mxArray *>(prhs[2]);
  }

  // defaults
  mwIndex maxiter = -1;
  if (nrhs < 4 || mxIsEmpty(prhs[3])) {
    maxiter = -1;
  } else {
    if (!mxIsDouble(prhs[3])) {
      mexErrMsgTxt("MAXITER must be a double scalar");
    }
    maxiter = (mwIndex)mxGetPr(prhs[3])[0];
  }
  
  // run function, templated according to the input matrix type
  switch(imClassId)  {
  case mxLOGICAL_CLASS:
    run<bool>(plhs[0], prhs[1], 
	      const_cast<const mxArray *>(res), maxiter);
    break;
  case mxDOUBLE_CLASS:
    run<double>(plhs[0], prhs[1], 
		const_cast<const mxArray *>(res), maxiter);
    break;
  case mxSINGLE_CLASS:
    run<float>(plhs[0], prhs[1], 
	       const_cast<const mxArray *>(res), maxiter);
    break;
  case mxINT8_CLASS:
    run<int8_T>(plhs[0], prhs[1], 
		const_cast<const mxArray *>(res), maxiter);
    break;
  case mxUINT8_CLASS:
    run<uint8_T>(plhs[0], prhs[1], 
		 const_cast<const mxArray *>(res), maxiter);
    break;
  case mxINT16_CLASS:
    run<int16_T>(plhs[0], prhs[1], 
		 const_cast<const mxArray *>(res), maxiter);
    break;
  case mxUINT16_CLASS:
    run<uint16_T>(plhs[0], prhs[1], 
		  const_cast<const mxArray *>(res), maxiter);
    break;
  case mxINT32_CLASS:
    run<int32_T>(plhs[0], prhs[1], 
		 const_cast<const mxArray *>(res), maxiter);
    break;
  // case mxUINT32_CLASS:
  // run<uint32_T>(plhs[0], prhs[1], const_cast<const mxArray *>(res), maxiter);
  //   break;
  case mxINT64_CLASS:
    run<int64_T>(plhs[0], prhs[1], 
		 const_cast<const mxArray *>(res), maxiter);
    break;
  // case mxUINT64_CLASS:
  //   run<uint64_T>(plhs[0], prhs[1], const_cast<const mxArray *>(res), maxiter);
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
