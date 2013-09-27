/*
 * MatlabImageHeader.h
 *
 * Class to provide an interface to read the metadata from a an input
 * image. This allows to create a header class that can be used to
 * deal with the image more easily.
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012-2013 University of Oxford
  * Version: 0.2.1
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

#ifndef MATLABIMAGEHEADER_H
#define MATLABIMAGEHEADER_H

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <vector>
#include <string>

class MatlabImageHeader {

 public:

  mxArray *data;               // pointer to the image voxels
  mxClassID type;              // pixel type
  std::vector<mwSize> size;    // number of voxels (row, col, slice)
  std::vector<double> spacing; // voxel size (row, col, slice)
  std::vector<double> origin;  // coordinates of centre of first voxel (row, col, slice)

  MatlabImageHeader(const mxArray *arg, std::string paramName);

  // get number of dimensions of the image
  size_t GetNumberOfDimensions() {
    return this->size.size();
  }
  
};

// constructor of the auxiliary class that preprocesses a Matlab
// argument that corresponds to an image
MatlabImageHeader::MatlabImageHeader(const mxArray *arg, std::string paramName) {

  if (arg == NULL) {
    mexErrMsgTxt(("Parameter " + paramName + " is a NULL pointer").c_str());
  }

  // get pointer to image data
  if (mxIsEmpty(arg)) { // input is empty

    this->data = NULL;
    this->type = mxUNKNOWN_CLASS;
    this->size.resize(0);
    this->spacing.resize(0);
    this->origin.resize(0);
    return;

  } else if (mxIsStruct(arg)) { // input is a SCI MAT struct
    // read struct fields...

    // ... data field (where the image is contained)
    this->data = mxGetField(arg, 0, "data");
    if (this->data == NULL) {
      mexErrMsgTxt(("Parameter " + paramName + 
		    ": Struct format error: Missing or invalid data field").c_str());
    }
  } else { // input is not a struct, but just an image array

    this->data = const_cast<mxArray *>(arg);

  }

  // get image type
  this->type = mxGetClassID(data);

  mwSize ndim; // number of dimensions in the image
  mwSize *dims; // dimensions array

  // get number of dimensions in the input image
  if (!mxIsNumeric(data) 
      && !mxIsLogical(data)) {
    mexErrMsgTxt(("Parameter " + paramName + 
		  " must be a numeric or boolean array").c_str());
  }
  ndim = mxGetNumberOfDimensions(data);
  
  // get image size
  dims = const_cast<mwSize *>(mxGetDimensions(data));
  if (dims == NULL) {
    mexErrMsgTxt(("Parameter " + paramName + 
		  ": Matlab cannot obtain vector of dimensions from the data").c_str());
  }
  for (unsigned int i = 0; i < ndim; ++i) {
    this->size.push_back(dims[i]);
  }

  // get rest of image info (spacing, origin)
  if (mxIsStruct(arg)) { // input is a SCI MAT struct{
    // ... axis field
    mxArray *axis = mxGetField(arg, 0, "axis");
    if (axis == NULL) {
      mexErrMsgTxt(("Parameter " + paramName + 
  		    ": Struct format error: Missing or invalid axis field").c_str());
    }
    if (!mxIsStruct(axis)) {
      mexErrMsgTxt(("Parameter " + paramName + 
  		    ": Struct format error: axis field is not a struct").c_str());
    }
    if (mxGetM(axis) < ndim) {
      mexErrMsgTxt(("Parameter " + paramName + 
  		    ": Struct format error: not enough elements in the axis field vector").c_str());
    }

    // ... spacing field (image resolution)
    // ... min field (note: "min"    = left corner of first voxel,
    //                      "origin" = centre of first voxel)
    mxArray *spacingMx, *minMx;
    double *spacingMxp, *minMxp;
    for (unsigned int i = 0; i < ndim; ++i) {
      spacingMx = mxGetField(axis, i, "spacing");
      if (spacingMx == NULL) {
  	mexErrMsgTxt(("Parameter " + paramName + 
  		      ": Struct format error: Missing or invalid axis.spacing field").c_str());
      }
      minMx = mxGetField(axis, i, "min");
      if (minMx == NULL) {
  	mexErrMsgTxt(("Parameter " + paramName + 
  		      ": Struct format error: Missing or invalid axis.min field").c_str());
      }
      spacingMxp = (double *)mxGetData(spacingMx);
      minMxp = (double *)mxGetData(minMx);
      if (mxIsFinite(spacingMxp[0]) && !mxIsNaN(spacingMxp[0])) {
  	spacing.push_back(spacingMxp[0]);
      } else {
  	spacing.push_back(1.0);
      }
      if (mxIsFinite(minMxp[0]) && !mxIsNaN(minMxp[0])) {
  	origin.push_back(minMxp[0] + spacing[i] / 2.0);
      } else {
  	origin.push_back(-spacing[i] / 2.0);
      }
    }

  } else { // input is not a struct, but just an image array

    for (unsigned int i = 0; i < ndim; ++i) {
      spacing.push_back(1.0);
      origin.push_back(0.0);
    }

  }

}

#endif /* MATLABIMAGEHEADER_H */
