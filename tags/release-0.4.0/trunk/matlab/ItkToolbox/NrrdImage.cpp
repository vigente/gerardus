/*
 * NrrdImage.cpp
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.1.2
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

#ifndef NRRDIMAGE_CPP
#define NRRDIMAGE_CPP

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <math.h>
#include <matrix.h>
#include <vector>

/* Gerardus headers */
#include "NrrdImage.hpp"

// The constructor works as a parser. It figures out whether the input
// data is a normal array with the image, or an NRRD struct. In the
// latter case, it reads the struct fields to extract the spacing and
// offset of the image
NrrdImage::NrrdImage(const mxArray * nrrd) {

  // initialize memory for member variables
  spacing.resize(Dimension);
  min.resize(Dimension);
  size.resize(Dimension);

  // check whether the input image is an SCI NRRD struct instead of
  // only the array. The SCI NRRD struct has information about the
  // scaling, offset, etc. of the image
  if (mxIsStruct(nrrd)) { // input is a struct
    // read struct fields...

    // ... data field (where the image is contained)
    data = mxGetField(nrrd, 0, "data");
    if (data == NULL) {
      mexErrMsgTxt("NRRD format error: Missing or invalid data field");
    }

    // ... axis field
    mxArray *axis = mxGetField(nrrd, 0, "axis");
    if (axis == NULL) {
      mexErrMsgTxt("NRRD format error: Missing or invalid axis field");
    }
    if (!mxIsStruct(axis)) {
      mexErrMsgTxt("NRRD format error: axis field is not a struct");
    }
    if (mxGetM(axis) != 3) {
      mexErrMsgTxt("NRRD format error: axis field must be a 3x1 struct array");
    }

    // ... spacing field (image resolution)
    // ... min field (image origin)
    mxArray *spacingMx, *minMx;
    double *spacingMxp, *minMxp;
    for (mwIndex i = 0; i < 3; ++i) {
      spacingMx = mxGetField(axis, i, "spacing");
      if (spacingMx == NULL) {
	mexErrMsgTxt("NRRD format error: Missing or invalid axis.spacing field");
      }
      minMx = mxGetField(axis, i, "min");
      if (minMx == NULL) {
	mexErrMsgTxt("NRRD format error: Missing or invalid axis.min field");
      }
      spacingMxp = (double *)mxGetData(spacingMx);
      minMxp = (double *)mxGetData(minMx);
      if (mxIsFinite(spacingMxp[0]) && !mxIsNaN(spacingMxp[0])) {
	spacing[i] = spacingMxp[0];
      } else {
	spacing[i] = 1.0;
      }
      if (mxIsFinite(minMxp[0]) && !mxIsNaN(minMxp[0])) {
	min[i] = minMxp[0];
      } else {
	min[i] = 0.0;
      }
    }

  } else { // input is not a struct, but just an image array
    data = const_cast<mxArray *>(nrrd);
    spacing[0] = 1.0;
    spacing[1] = 1.0;
    spacing[2] = 1.0;
    min[0] = 0.0;
    min[1] = 0.0;
    min[2] = 0.0;
  }


  // get number of dimensions in the input image
  if (!mxIsNumeric(data) 
      && !mxIsLogical(data)) {
    mexErrMsgTxt("IM must be a 2D or 3D numeric or boolean matrix");
  }
  ndim = mxGetNumberOfDimensions(data);
  
  // get size of input arguments
  dims = const_cast<mwSize *>(mxGetDimensions(data));
  if (dims == NULL) {
    mexErrMsgTxt("Invalid input image dimensions array");
  }
  size[0] = dims[0];
  if (ndim > 3) {
    mexErrMsgTxt("Input segmentation mask must be 2D or 3D");
  } else if (ndim == 2) {
    size[1] = dims[1];
    size[2] = 1;
  } else if (ndim == 3) {
    size[1] = dims[1];
    size[2] = dims[2];
  } else {
    mexErrMsgTxt("Assertion fail: number of dimensions is " + ndim);
  }

}

// compute the maximum distance between any two voxels in this image
// (in voxel units). This is the length of the largest diagonal in the
// cube
double NrrdImage::maxVoxDistance() {
    return std::sqrt((double)(size[0]-1)*(size[0]-1)
		   + (double)(size[1]-1)*(size[1]-1)
		   + (double)(size[2]-1)*(size[2]-1));
}

// compute the number of voxels in the image volume
mwSize NrrdImage::numEl() {
  return size[0]*size[1]*size[2];
}

#endif /* NRRDIMAGE_CPP */
