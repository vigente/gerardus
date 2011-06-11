/*
 * ArgumentParsers.cpp
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

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.2.2
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

#ifndef ARGUMENTPARSERS_CPP
#define ARGUMENTPARSERS_CPP

/* C++ headers */
#include <vector>

/* mex headers */
#include <mex.h>

/* ITK headers */
#include "itkImage.h"

/* Gerardus headers */
#include "NrrdImage.hpp"
#include "BaseFilter.hpp"
#include "DanielssonFilter.hpp"
#include "SignedMaurerFilter.hpp"
#include "ThinningFilter.hpp"

// parseFilterTypeToTemplate<InVoxelType, OutVoxelType>()
template <class InVoxelType, class OutVoxelType>
void parseFilterTypeToTemplate(char *filterName,
			       NrrdImage nrrd,
			       int nargout,
			       mxArray** argOut) {

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
void parseOutputTypeToTemplate(char *filter,
			       NrrdImage nrrd,
			       int nargout,
			       mxArray** argOut) {

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
      InVoxelType>(filter, nrrd, nargout, argOut);
    break;
  case BOOL:
    parseFilterTypeToTemplate<InVoxelType, 
      bool>(filter, nrrd, nargout, argOut);
    break;
  case UINT8:
    parseFilterTypeToTemplate<InVoxelType, 
      uint8_T>(filter, nrrd, nargout, argOut);
    break;
  case UINT16:
    parseFilterTypeToTemplate<InVoxelType, 
      uint16_T>(filter, nrrd, nargout, argOut);
    break;
  case SINGLE:
    parseFilterTypeToTemplate<InVoxelType, 
      float>(filter, nrrd, nargout, argOut);
    break;
  case DOUBLE:
    parseFilterTypeToTemplate<InVoxelType, 
      double>(filter, nrrd, nargout, argOut);
    break;
  default:
    mexErrMsgTxt("Invalid output type.");
    break;
  }
}

// parseInputTypeToTemplate()
void parseInputTypeToTemplate(mxClassID inputVoxelClassId, 
			      char *filter,
			      NrrdImage nrrd,
			      int nargout,
			      mxArray** argOut) {
  
  switch(inputVoxelClassId)  { // swith input image type
  case mxLOGICAL_CLASS:
    parseOutputTypeToTemplate<bool>(filter, nrrd, nargout, argOut);
    break;
  case mxDOUBLE_CLASS:
    parseOutputTypeToTemplate<double>(filter, nrrd, nargout, argOut);
    break;
  case mxSINGLE_CLASS:
    parseOutputTypeToTemplate<float>(filter, nrrd, nargout, argOut);
    break;
  case mxINT8_CLASS:
    parseOutputTypeToTemplate<int8_T>(filter, nrrd, nargout, argOut);
    break;
  case mxUINT8_CLASS:
    parseOutputTypeToTemplate<uint8_T>(filter, nrrd, nargout, argOut);
    break;
  case mxINT16_CLASS:
    parseOutputTypeToTemplate<int16_T>(filter, nrrd, nargout, argOut);
    break;
  case mxUINT16_CLASS:
    parseOutputTypeToTemplate<uint16_T>(filter, nrrd, nargout, argOut);
    break;
  case mxINT32_CLASS:
    parseOutputTypeToTemplate<int32_T>(filter, nrrd, nargout, argOut);
    break;
  // case mxUINT32_CLASS:
  //   break;
  case mxINT64_CLASS:
    parseOutputTypeToTemplate<int64_T>(filter, nrrd, nargout, argOut);
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

#endif /* ARGUMENTPARSERS_CPP */
