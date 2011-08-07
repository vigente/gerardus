/* 
 * BaseFilter.hpp
 *
 * BaseFilter<InVoxelType, OutVoxelType>: This is where
 * the code to actually run the filter on the image lives.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.4.0
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

#ifndef BASEFILTER_HPP
#define BASEFILTER_HPP

/* ITK headers */
#include "itkImageToImageFilter.h"

/* Gerardus headers */
#include "GerardusCommon.hpp"
#include "NrrdImage.hpp"

/* 
 * BaseFilter
 */
template <class InVoxelType, class OutVoxelType>
class BaseFilter {
private:
  typedef double TScalarType; // data type for scalars
  
protected:
  typedef typename itk::Image< InVoxelType, Dimension > InImageType;
  typedef typename itk::Image< OutVoxelType, Dimension > OutImageType;
  typedef typename itk::ImageToImageFilter<InImageType, OutImageType> BaseFilterType;
  typename InImageType::Pointer image;
  NrrdImage nrrd;
  int nargout;
  mxArray** argOut;
  mwSize nparam;
  const mxArray** argParam;
  typename BaseFilterType::Pointer filter;

  // function to get the value of input parameters that are numeric
  // scalars from the array of input arguments
  template <class ParamType>
  ParamType getScalarParamValue(std::string paramName,
				mwIndex idx, ParamType def);

public:

  // constructor when the filter takes user-provided parameters
  // (e.g. dilation radius) apart from the input arguments with the
  // filter type and image
  BaseFilter(const NrrdImage &_nrrd, int _nargout, mxArray** _argOut,
	     const int _nargin, const mxArray** _argIn)
    : nrrd(_nrrd), nargout(_nargout), argOut(_argOut) {
    // number of parameters
    nparam = (mwSize)(_nargin - 2);
    if (nparam < 0) {
      mexErrMsgTxt("Assertion error: Number of parameters cannot be negative");
    } else if (nparam == 0) {
      argParam = NULL;
    } else {
      argParam = &_argIn[2];
    }
  }
  
  // constructor when the filter takes no user-provided parameters
  BaseFilter(const NrrdImage &_nrrd, int _nargout, mxArray** _argOut)
    : nrrd(_nrrd), nargout(_nargout), argOut(_argOut) {
    nparam = 0;
    argParam = NULL;
  }

  // constructor when a derived filter needs to be excluded from
  // instantiation
  BaseFilter() {;}

  // functions to create the ITK images, filter it and return a Matlab result
  virtual void CopyMatlabInputsToItkImages();
  virtual void FilterSetup();
  virtual void RunFilter();
  virtual void CopyAllFilterOutputsToMatlab();
  void CopyFilterImageOutputToMatlab();
};

/*
 * Templated member functions. They need to go here in the header
 * file, because if they go in the .cpp, then they need instance
 * specialization for every <InVoxelType, OutVoxelType, ParamType>
 *
 */

// function to get the value of input parameters that are numeric
// scalars from the array of input arguments
template <class InVoxelType, class OutVoxelType>
template <class ParamType>
ParamType BaseFilter<InVoxelType, 
		     OutVoxelType>::getScalarParamValue(std::string paramName,
				mwIndex idx, ParamType def) {
  
  // if user didn't provide a value, or provided an empty array, return the default
  if (idx >= this->nparam || mxIsEmpty(this->argParam[idx])) {
    return def;
  }

  // if user provided a value, check that it's a scalar, whether in
  // numeric or logical form
  if ((!mxIsNumeric(this->argParam[idx])
       && !mxIsLogical(this->argParam[idx]))
      || mxGetNumberOfElements(this->argParam[idx]) != 1) {
    mexErrMsgTxt((std::string("Parameter ") + paramName 
		  + std::string(" must be a scalar")).c_str());
  }
  
  // output
  ParamType value = 0;
  
  // input image type
  mxClassID inputVoxelClassId = mxGetClassID(this->argParam[idx]);

  // macro to make the code in the switch statement cleaner
#define GETVALUE(Tx, x, Ty, y)						\
  if (x) {								\
    Tx *valuep = (Tx *)mxGetData(x);					\
    y = (Ty)valuep[0];							\
  } else {								\
    mexErrMsgTxt((std::string("Parameter ") + paramName			\
		  + std::string(" couldn't be read")).c_str());		\
}
  
    // cast the class type provided by Matlab to the type requested by
    // the user
  switch(inputVoxelClassId)  { 
  case mxLOGICAL_CLASS:
    GETVALUE(mxLogical, this->argParam[idx], ParamType, value);
    break;
  case mxDOUBLE_CLASS:
    GETVALUE(double, this->argParam[idx], ParamType, value);
    break;
  case mxSINGLE_CLASS:
    GETVALUE(float, this->argParam[idx], ParamType, value);
    break;
  case mxINT8_CLASS:
    GETVALUE(int8_T, this->argParam[idx], ParamType, value);
    break;
  case mxUINT8_CLASS:
    GETVALUE(uint8_T, this->argParam[idx], ParamType, value);
    break;
  case mxINT16_CLASS:
    GETVALUE(int16_T, this->argParam[idx], ParamType, value);
    break;
  case mxUINT16_CLASS:
    GETVALUE(uint16_T, this->argParam[idx], ParamType, value);
    break;
  case mxINT32_CLASS:
    GETVALUE(int32_T, this->argParam[idx], ParamType, value);
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    GETVALUE(int64_T, this->argParam[idx], ParamType, value);
    break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt("Input parameter has unknown type.");
    break;
  default:
    mexErrMsgTxt("Input parameter has invalid type.");
      break;
  }

#undef GETVALUE
  
  return value;
}

#endif /* BASEFILTER_HPP */
