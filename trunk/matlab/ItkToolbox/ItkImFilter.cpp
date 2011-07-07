/* ItkImFilter.cpp
 *
 * ITK_IMFILTER: Run ITK filter on a 2D or 3D image
 *
 * This MEX function is a multiple-purpose wrapper to be able to run
 * all ITK filters that inherit from itk::ImageToImageFilter on a
 * Matlab 2D image or 3D image volume.
 *
 * B = ITK_IMFILTER(TYPE, A)
 *
 *   TYPE is a string with the filter we want to run. Currently, only
 *   the following options are implemented:
 *
 *     'skel':    (BinaryThinningImageFilter3D) Skeletonize a
 *                binary mask
 *
 *                B has the same size and class as A
 *
 *     'dandist': (DanielssonDistanceMapImageFilter) Compute unsigned
 *                distance map for a binary mask. Distance values are
 *                given in voxel coordinates
 *
 *                B has the same size as A. B has a type large enough
 *                to store the maximum distance in the image. The
 *                largest available type is double. If this is not
 *                enough, a warning message is displayed, and double
 *                is used as the output type
 *
 *     'maudist': (SignedMaurerDistanceMapImageFilter) Compute signed
 *                distance map for a binary mask. Distance values are
 *                given in real world coordinates, if the input image
 *                is given as an NRRD struct, or in voxel units, if
 *                the input image is a normal array. The output type
 *                is always double.
 *
 *   A is a 2D matrix or 3D volume with the image or
 *   segmentation. Currently, A can be of any of the following
 *   Matlab classes:
 *
 *     boolean
 *     double
 *     single
 *     int8
 *     uint8
 *     int16
 *     uint16
 *     int32
 *     int64
 *
 *   A can also be a SCI NRRD struct, A = nrrd, with the following fields:
 *
 *     nrrd.data: 2D or 3D array with the image or segmentation, as above
 *     nrrd.axis: 3x1 struct array with fields:
 *       nnrd.axis.size:    number of voxels in the image
 *       nnrd.axis.spacing: voxel size, image resolution
 *       nnrd.axis.min:     real world coordinates of image origin
 *       nnrd.axis.max:     ignored
 *       nnrd.axis.center:  ignored
 *       nnrd.axis.label:   ignored
 *       nnrd.axis.unit:    ignored
 *
 *   An SCI NRRD struct is the output of Matlab's function
 *   scinrrd_load(), also available from Gerardus.
 *
 *   B has the same size and class as the image in A, and contains the
 *   filtered image or segmentation mask.
 *
 * This function must be compiled before it can be used from Matlab.
 * If Gerardus' root directory is e.g. ~/gerardus, type from a
 * linux shell
 *
 *    $ cd ~/gerardus/matlab
 *    $ mkdir bin
 *    $ cd bin
 *    $ cmake ..
 *    $ make install
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.3.10
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

#ifndef ITKIMFILTER_CPP
#define ITKIMFILTER_CPP

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <cmath>
#include <matrix.h>
#include <vector>

/* ITK headers */
#include "itkImage.h"
#include "itkBinaryThinningImageFilter3D.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

/* Gerardus headers */
#include "NrrdImage.hpp"
#include "BaseFilter.hpp"
#include "DanielssonFilter.hpp"
#include "SignedMaurerFilter.hpp"
#include "ThinningFilter.hpp"

/*
 * Argument Parsers
 *
 * parseInputTypeToTemplate()
 * parseOutputTypeToTemplate< InVoxelType >()
 * parseFilterTypeToTemplate< InVoxelType, OutVoxelType >()
 *
 * These functions are used to be able to map between the input/output
 * data types that are only know at run-time, and the input/output
 * data templates that ITK requires and must be know at compilation
 * time.
 *
 * To avoid a nesting nightmare:
 *
 * switch FilterType {
 *   switch InputDataType {
 *     switch OutputDataType {
 *     }
 *   }
 * }
 *
 * we split the conversion in 3 steps. The first function,
 * parseInputTypeToTemplate() reads the input data type and the
 * filter type, and instantiates
 * parseOutputTypeToTemplate<InVoxelType>() for each input data
 * type.
 *
 * This way, we have effectively mapped the input type from a run-time
 * variable to a set of compilation time templates.
 *
 * Then parseOutputTypeToTemplate<InVoxelType>() decides on the
 * output data type depending on the filter and the input, and
 * instantiates one
 * parseFilterTypeToTemplate<InVoxelType, OutVoxelType>() for each
 * output data type.
 *
 * This way, we have instantiated all combinations of input/output
 * data types, without having to code them explicitly.
 *
 * Finally, parseFilterTypeToTemplate<InVoxelType, OutVoxelType>()
 * checks that the requested filter is implemented and maps the filter
 * variable to all filter templates.
 *
 * Note that the actual filtering is done using derived filter classes
 * from BaseFilter, that is instantiated from
 * parseFilterTypeToTemplate.
 */

// parseFilterTypeToTemplate<InVoxelType, OutVoxelType>()
template <class InVoxelType, class OutVoxelType>
void parseFilterTypeToTemplate(NrrdImage nrrd,
			       int nargout,
			       mxArray** argOut,
			       char *filterName) {

  // image type definitions
  typedef double TScalarType; // data type for scalars
  typedef itk::Image< InVoxelType, Dimension > 
    InImageType;
  typedef itk::Image< OutVoxelType, Dimension > 
    OutImageType;

  // pointer to the filter object (we are using polymorphism)
  BaseFilter<InVoxelType, OutVoxelType> *filter = NULL;

  // convert run-time filter string to template
  if (!strcmp(filterName, "skel")) {
    
    filter = new ThinningFilter<InVoxelType, 
      OutVoxelType>(nrrd, nargout, argOut);

  }  else if (!strcmp(filterName, "dandist")) {
    
    filter = new DanielssonFilter<InVoxelType, 
      OutVoxelType>(nrrd, nargout, argOut);

  }  else if (!strcmp(filterName, "maudist")) {

    filter = new SignedMaurerFilter<InVoxelType, 
      OutVoxelType>(nrrd, nargout, argOut);

  } else {
    mexErrMsgTxt("Filter type not implemented");
  }
  
  // set up and run filter
  filter->CopyMatlabInputsToItkImages();
  filter->FilterSetup();
  filter->RunFilter();
  filter->CopyAllFilterOutputsToMatlab();
}

// parseOutputTypeToTemplate<InVoxelType>()
template <class InVoxelType>
void parseOutputTypeToTemplate(NrrdImage nrrd,
			       int nargout,
			       mxArray** argOut,
			       char *filter) {

  // make it easier to remember the different cases for the output
  // voxel type
  enum OutVoxelType {
    SAME, BOOL, UINT8, UINT16, SINGLE, DOUBLE
  };

  // establish output voxel type according to the filter
  OutVoxelType outVoxelType = DOUBLE;
  if (!strcmp(filter, "skel")) {

    outVoxelType = SAME;

  } else if (!strcmp(filter, "dandist")) {
    // find how many bits we need to represent the maximum distance
    // that two voxels can have between them (in voxel units)
      mwSize nbit = (mwSize)ceil(log(nrrd.maxVoxDistance()) / log(2.0));

    // select an output voxel size enough to save the maximum distance
    // value
    if (nbit <= 2) {
      outVoxelType = BOOL;
    } else if (nbit <= 8) {
      outVoxelType = UINT8;
    } else if (nbit <= 16) {
      outVoxelType = UINT16;
    } else if (nbit <= 128) {
      outVoxelType = SINGLE;
    } else {
      outVoxelType = DOUBLE;
    }
    
  } else if (!strcmp(filter, "maudist")) {

    outVoxelType = DOUBLE;

  } else {
    mexErrMsgTxt("Filter type not implemented");
  }

  switch(outVoxelType) {
  case SAME:
    parseFilterTypeToTemplate<InVoxelType, 
			      InVoxelType>(nrrd, nargout, argOut, filter);
    break;
  case BOOL:
    parseFilterTypeToTemplate<InVoxelType, 
			      bool>(nrrd, nargout, argOut, filter);
    break;
  case UINT8:
    parseFilterTypeToTemplate<InVoxelType, 
			      uint8_T>(nrrd, nargout, argOut, filter);
    break;
  case UINT16:
    parseFilterTypeToTemplate<InVoxelType, 
			      uint16_T>(nrrd, nargout, argOut, filter);
    break;
  case SINGLE:
    parseFilterTypeToTemplate<InVoxelType, 
			      float>(nrrd, nargout, argOut, filter);
    break;
  case DOUBLE:
    parseFilterTypeToTemplate<InVoxelType, 
			      double>(nrrd, nargout, argOut, filter);
    break;
  default:
    mexErrMsgTxt("Invalid output type.");
    break;
  }
}

// parseInputTypeToTemplate()
void parseInputTypeToTemplate(NrrdImage nrrd,
			      int nargout,
			      mxArray** argOut,
			      char *filter) {
  
  // input image type
  mxClassID inputVoxelClassId = mxGetClassID(nrrd.getData());

  switch(inputVoxelClassId)  { // swith input image type
  case mxLOGICAL_CLASS:
    parseOutputTypeToTemplate<bool>(nrrd, nargout, argOut, filter);
    break;
  case mxDOUBLE_CLASS:
    parseOutputTypeToTemplate<double>(nrrd, nargout, argOut, filter);
    break;
  case mxSINGLE_CLASS:
    parseOutputTypeToTemplate<float>(nrrd, nargout, argOut, filter);
    break;
  case mxINT8_CLASS:
    parseOutputTypeToTemplate<int8_T>(nrrd, nargout, argOut, filter);
    break;
  case mxUINT8_CLASS:
    parseOutputTypeToTemplate<uint8_T>(nrrd, nargout, argOut, filter);
    break;
  case mxINT16_CLASS:
    parseOutputTypeToTemplate<int16_T>(nrrd, nargout, argOut, filter);
    break;
  case mxUINT16_CLASS:
    parseOutputTypeToTemplate<uint16_T>(nrrd, nargout, argOut, filter);
    break;
  case mxINT32_CLASS:
    parseOutputTypeToTemplate<int32_T>(nrrd, nargout, argOut, filter);
    break;
  // case mxUINT32_CLASS:
  //   break;
  case mxINT64_CLASS:
    parseOutputTypeToTemplate<int64_T>(nrrd, nargout, argOut, filter);
    break;
  // case mxUINT64_CLASS:
  //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt("Input matrix has unknown type.");
    break;
  default:
    mexErrMsgTxt("Input matrix has invalid type.");
    break;
  }
}

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {
  // check number of input and output arguments
  if (nrhs != 2) {
    mexErrMsgTxt("Two input arguments required");
  }

  // read image and its parameters, whether it's in NRRD format, or
  // just a 2D or 3D array
  NrrdImage nrrd(prhs[1]);

  // get type of filter
  char *filter = mxArrayToString(prhs[0]);
  if (filter == NULL) {
    mexErrMsgTxt("Invalid FILTER string");
  }

  // run filter (this function starts a cascade of functions designed
  // to translate the run-time type variables like inputVoxelClassId
  // to templates, so that we don't need to nest lots of "switch" or
  // "if" statements)
  parseInputTypeToTemplate(nrrd,
			   nlhs,
			   plhs,
			   filter);

  // exit successfully
  return;

}

#endif /* ITKIMFILTER_CPP */
