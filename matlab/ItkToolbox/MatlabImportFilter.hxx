/*
 * MatlabImportFilter.hxx
 *
 * Class to provide an interface to import data from Matlab mxArrays
 * into ITK
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012 University of Oxford
  * Version: 0.1.0
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

#ifndef MATLABIMPORTFILTER_HXX
#define MATLABIMPORTFILTER_HXX

/* ITK headers */
#include "itkImportImageFilter.h"

/* Gerardus headers */
#include "MatlabImportFilter.h"

// constructor
MatlabImportFilter::MatlabImportFilter() {
  this->args = NULL;
  this->numArgs = 0;
}

// destructor
// this class is only an interface to handle pointer, so we must not
// attempt to delete the argument list provided by Matlab
MatlabImportFilter::~MatlabImportFilter() {}

// function to check that number of input arguments is within
// certain limits
void MatlabImportFilter::CheckNumberOfArguments(unsigned int min, unsigned int max) {
  if (this->numArgs < min) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (this->numArgs > max) {
    mexErrMsgTxt("Too many input arguments");
  }
}

// function to get the value of input arguments that are strings
std::string MatlabImportFilter::GetStringArgument(unsigned int idx,
						       std::string paramName,
						       std::string def) {
  
  // if user didn't provide a value, or provided an empty array, return the default
  if (idx >= this->numArgs || mxIsEmpty(this->args[idx])) {
    return def;
  }

  // if user provided a value, check that it's a string
  if (!mxIsChar(this->args[idx])) {
    mexErrMsgTxt(("Parameter " + paramName + " must be a string").c_str());
  }
  
  // import string
  return mxArrayToString(this->args[idx]);

}

// function to get the value of input parameters that are numeric
// scalars from the array of input arguments
template <class ParamType>
ParamType MatlabImportFilter::GetScalarArgument(unsigned int idx, 
						std::string paramName,
						ParamType def) {
  
  // if user didn't provide a value, or provided an empty array, return the default
  if (idx >= this->numArgs || mxIsEmpty(this->args[idx])) {
    return def;
  }
  
  // check for null pointer
  if (this->args[idx] == NULL) {
    mexErrMsgTxt(("Parameter " + paramName + " provided, but NULL pointer.").c_str());
  }

  // if user provided a value, check that it's a scalar, whether in
  // numeric or logical form
  if ((!mxIsNumeric(this->args[idx])
       && !mxIsLogical(this->args[idx]))
      || mxGetNumberOfElements(this->args[idx]) != 1) {
    mexErrMsgTxt(("Parameter " + paramName + " must be a scalar").c_str());
  }
  
  // output
  ParamType value = 0;
  
  // input image type
  mxClassID inputVoxelClassId = mxGetClassID(this->args[idx]);
  
  // macro to make the code in the switch statement cleaner
#define GETVALUE(Tx)						\
  {								\
    Tx *valuep = (Tx *)mxGetData(this->args[idx]);		\
    value = (ParamType)valuep[0];				\
  }
  
  // cast the class type provided by Matlab to the type requested by
  // the user
  switch(inputVoxelClassId)  {
  case mxLOGICAL_CLASS:
    GETVALUE(mxLogical);
    break;
  case mxDOUBLE_CLASS:
    GETVALUE(double);
    break;
  case mxSINGLE_CLASS:
    GETVALUE(float);
    break;
  case mxINT8_CLASS:
    GETVALUE(int8_T);
    break;
  case mxUINT8_CLASS:
    GETVALUE(uint8_T);
    break;
  case mxINT16_CLASS:
    GETVALUE(int16_T);
    break;
  case mxUINT16_CLASS:
    GETVALUE(uint16_T);
    break;
  case mxINT32_CLASS:
    GETVALUE(int32_T);
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    GETVALUE(int64_T);
    break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt(("Parameter " + paramName + " has unknown type.").c_str());
    break;
  default:
    mexErrMsgTxt(("Parameter " + paramName + " has invalid type.").c_str());
    break;
  }
  
#undef GETVALUE
  
  return value;
}

// function to get the an input argument that is a vector of scalars
// from the array of input arguments
template <class ParamType, class ParamValueType>
ParamType MatlabImportFilter::GetVectorArgument(unsigned int idx, 
						std::string paramName,
						ParamType def) {

  // if user didn't provide a value, or provided an empty array,
  // return default
  if (idx >= this->numArgs || mxIsEmpty(this->args[idx])) {
    return def;
  }

  // check for null pointer
  if (this->args[idx] == NULL) {
    mexErrMsgTxt(("Parameter " + paramName + " provided, but NULL pointer.").c_str());
  }

  // check that we have a row or column vector, numeric or boolean
  if ((mxGetM(this->args[idx]) != 1) && (mxGetN(this->args[idx]) != 1)) {
    mexErrMsgTxt(("Parameter " + paramName + " must be a vector.").c_str());
  }
  if (!mxIsNumeric(this->args[idx]) && !mxIsLogical(this->args[idx])) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must be a numeric or logical vector.").c_str());
  }

  // vector to be returned
  ParamType param;

  // input vector length
  size_t len = std::max(mxGetM(this->args[idx]), mxGetN(this->args[idx]));
  if (len != param.GetSizeDimension()) {
    mexErrMsgTxt(("Parameter " + paramName 
    		  + " has the wrong length.").c_str());
  }

  // input image type
  mxClassID inputVoxelClassId = mxGetClassID(this->args[idx]);
  
  // macro to make the code in the switch statement cleaner
#define GETVALUE(Tx)					\
  {							\
    Tx *valuep = (Tx *)mxGetData(this->args[idx]);	\
    for (size_t i = 0; i < len; ++i) {			\
      param[i] = (ParamValueType)valuep[i];		\
    }							\
  }
  
  // cast the class type provided by Matlab to the type requested by
  // the user
  switch(inputVoxelClassId)  { 
  case mxLOGICAL_CLASS:
    GETVALUE(mxLogical);
    break;
  case mxDOUBLE_CLASS:
    GETVALUE(double);
    break;
  case mxSINGLE_CLASS:
    GETVALUE(float);
    break;
  case mxINT8_CLASS:
    GETVALUE(int8_T);
    break;
  case mxUINT8_CLASS:
    GETVALUE(uint8_T);
    break;
  case mxINT16_CLASS:
    GETVALUE(int16_T);
    break;
  case mxUINT16_CLASS:
    GETVALUE(uint16_T);
    break;
  case mxINT32_CLASS:
    GETVALUE(int32_T);
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    GETVALUE(int64_T);
    break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt(("Parameter " + paramName + " has unknown type.").c_str());
    break;
  default:
    mexErrMsgTxt(("Parameter " + paramName + " has invalid type.").c_str());
    break;
  }

#undef GETVALUE

  return param;
}

// function to get an input argument that is an image
template <class TPixel, unsigned int VImageDimension>
typename itk::Image<TPixel, VImageDimension>::Pointer
MatlabImportFilter::GetImageArgument(unsigned int idx, std::string paramName) {
  
  // note that:
  //
  // 1) in ITK we have X,Y,Z indices, while in Matlab we have R,C,S
  // (row, column, slice)
  //
  // 2) matrices in ITK are read by columns, while in Matlab
  // they are read by rows 
  //
  // So imagine we have this (2, 3) matrix in Matlab
  //
  //   a b   |
  //   c d   | y-axis (resolution 1.0)
  //   e f   |
  //   ---
  //   x-axis (resolution 0.5)
  //
  //   size = [3 2 1]
  //
  // The C-style array is going to be (reading by rows)
  //
  //   im = [a c e b d f]
  //
  // ITK is going to read by colums, thinking that the size is
  //
  //   size = [sx sy sz] = [3 2 1]
  //
  //   a c e   |
  //   b d f   | y-axis (resolution 0.5)
  //   -----
  //   x-axis (resolution 1.0)
  //
  // Note that the matrix has been transposed, but this is not a
  // problem, because the resolution values have been "transposed"
  // too
  //
  // Having the matrix transposed may make us feel a bit uneasy, but
  // it has the advantage that Matlab and ITK can use the same C-style
  // array, without having to rearrange its elements

  // instantiate the filter that will act as an interface between the
  // Matlab image array and the ITK filter
  typedef itk::ImportImageFilter<TPixel, VImageDimension> ImportImageFilterType;
  typename ImportImageFilterType::RegionType region;
  typename ImportImageFilterType::SizeType size;
  typename ImportImageFilterType::IndexType start;
  typename ImportImageFilterType::SpacingType spacing;
  typename ImportImageFilterType::OriginType origin;
  typename ImportImageFilterType::Pointer importFilter = ImportImageFilterType::New();

  // create a header for the image
  MatlabImageHeader imageHeader(this->args[idx], paramName);

  // convert image header parameters to a format that can be passed to
  // the import filter
  for (size_t i = 0; i < imageHeader.GetNumberOfDimensions(); ++i) {
    start[i] = 0;
    size[i] = imageHeader.size[i];
    spacing[i] = imageHeader.spacing[i];
    origin[i] = imageHeader.origin[i];
  }
  
  // set image metainformation
  region.SetIndex(start);
  region.SetSize(size);
  importFilter->SetRegion(region);
  importFilter->SetSpacing(spacing);
  importFilter->SetOrigin(origin);

  // get pointer to input segmentation mask
  const TPixel *im = (TPixel *)mxGetData(imageHeader.data);

  // pass pointer to Matlab image to the import filter, and tell it to
  // NOT attempt to delete the memory when it's destructor is
  // called. This is important, because the input image still has to
  // live in Matlab's memory after running the filter
  const bool importFilterWillOwnTheBuffer = false;
  importFilter->SetImportPointer(const_cast<TPixel *>(im),
				 mxGetNumberOfElements(this->args[idx]),
				 importFilterWillOwnTheBuffer);

  // actually import the image
  importFilter->Update();

  // succesful exit
  return importFilter->GetOutput();

}

#endif /* MATLABIMPORTFILTER_HXX */
