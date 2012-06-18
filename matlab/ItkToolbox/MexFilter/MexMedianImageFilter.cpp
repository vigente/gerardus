/*
 * MexMedianImageFilter.cpp
 *
 * Code that is specific to itk::MedianImageFilter.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012 University of Oxford
  * Version: 0.1.1
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

#ifndef MEXMEDIANIMAGEFILTER_CPP
#define MEXMEDIANIMAGEFILTER_CPP

/* C++ headers */

/* ITK headers */
#include "itkBoxImageFilter.h"

/* Gerardus headers */
#include "MexMedianImageFilter.hpp"

/*
 * strings that the user can type to invoke this filter in itk_imfilter()
 */
const std::string MexMedianImageFilter<std::string, 
		  std::string>::longname = "MedianImageFilter";
const std::string MexMedianImageFilter<std::string, 
		  std::string>::shortname = "median";

/* 
 * constructor (here we instantiate the filter and process the
 * user-provided input parameters, if any)
 */
template <class InVoxelType, class OutVoxelType>
MexMedianImageFilter<InVoxelType, OutVoxelType>::MexMedianImageFilter(
                                const NrrdImage &_nrrd, 
				int _numArgOut, mxArray** _argOut,
				const int _numArgIn, const mxArray** _argIn) :
  MexBaseFilter<InVoxelType, OutVoxelType>(_nrrd, _numArgOut, _argOut, _numArgIn, _argIn) {

  // instantiate filter in this derived class, but on the base class
  // pointer, thanks to polimorphism. This way, we can run methods on
  // the derived class from the base class
  this->filter = DerivedImageToImageFilterType::New();

  // get a pointer to the filter in this derived class. We cannot use
  // this->filter if we want to access methods that are only in the
  // derived class, because this->filter points to the filter in the
  // base class
  derivedFilter = 
    dynamic_cast<DerivedImageToImageFilterType *>(this->filter.GetPointer());

  // check number of user-provided parameters (user-provided
  // parameters are the extra input arguments apart from the filter
  // type and input image)
  if (this->numParam < 0) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (this->numParam > 3) {
    mexErrMsgTxt("Too many input arguments");
  }
  if (this->numParam > 0 && this->argParam == NULL) {
    mexErrMsgTxt("Assertion fail: There is at least one parameter, but pointer to parameter array is NULL");
  }
  
  // get user-provided parameters: 
  //    parameter name
  //    index (0 = first parameter)
  //    default value
  // Example:
  this->radiusR = this->template 
    GetScalarParamValue<mwSize>("RadiusR", 0, 0);
  this->radiusC = this->template 
    GetScalarParamValue<mwSize>("RadiusC", 1, 0);
  this->radiusS = this->template 
    GetScalarParamValue<mwSize>("RadiusS", 2, 0);

}

/* 
 * overriding of MexBaseFilter virtual methods, if needed
 */

// // add code here if you need to pass user-provided parameters to the
// // filter, or perform any other kind of filter setup
template <class InVoxelType, class OutVoxelType>
void MexMedianImageFilter<InVoxelType, 
			    OutVoxelType>::FilterAdvancedSetup() {

  // set half size of the filter's box
  typedef typename itk::BoxImageFilter<
    itk::Image<InVoxelType, Dimension>,
    itk::Image<OutVoxelType, Dimension> > BoxFilterType;
  typename BoxFilterType::RadiusType radius;
  radius[0] = this->radiusR;
  radius[1] = this->radiusC;
  radius[2] = this->radiusS;
  derivedFilter->SetRadius(radius);

}

// // this method overrides MexBaseFilter::CheckNumberOfOutputs(). Change
// // the code here if this filter provides more than one output
// template <class InVoxelType, class OutVoxelType>
// void MexMedianImageFilter<InVoxelType, 
// 			    OutVoxelType>::CheckNumberOfOutputs() {
//
//   // prevent the user from asking for too many output arguments
//   if (this->numArgOut > 1) {
//     mexErrMsgTxt("Too many output arguments");
//   }
//
// }

// // this method overrides
// // MexBaseFilter::ExportOtherFilterOutputsToMatlab(). Add code here
// // if you want to extract extra outputs from the filter, apart from
// // the main filtered image output
// template <class InVoxelType, class OutVoxelType>
// void MexMedianImageFilter<InVoxelType, 
// 			    OutVoxelType>::ExportOtherFilterOutputsToMatlab() {
//
// }

/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#define FILTERINST(T1, T2)				\
  template class MexMedianImageFilter<T1, T2>;

FILTERINST(mxLogical, mxLogical)
// FILTERINST(mxLogical, uint8_T)
// FILTERINST(mxLogical, int8_T)
// FILTERINST(mxLogical, uint16_T)
// FILTERINST(mxLogical, int16_T)
// FILTERINST(mxLogical, int32_T)
// FILTERINST(mxLogical, int64_T)
// FILTERINST(mxLogical, float)
// FILTERINST(mxLogical, double)

// FILTERINST(uint8_T, mxLogical)
FILTERINST(uint8_T, uint8_T)
// FILTERINST(uint8_T, int8_T)
// FILTERINST(uint8_T, uint16_T)
// FILTERINST(uint8_T, int16_T)
// FILTERINST(uint8_T, int32_T)
// FILTERINST(uint8_T, int64_T)
// FILTERINST(uint8_T, float)
// FILTERINST(uint8_T, double)

// FILTERINST(int8_T, mxLogical)
// FILTERINST(int8_T, uint8_T)
FILTERINST(int8_T, int8_T)
// FILTERINST(int8_T, uint16_T)
// FILTERINST(int8_T, int16_T)
// FILTERINST(int8_T, int32_T)
// FILTERINST(int8_T, int64_T)
// FILTERINST(int8_T, float)
// FILTERINST(int8_T, double)

// FILTERINST(uint16_T, mxLogical)
// FILTERINST(uint16_T, uint8_T)
// FILTERINST(uint16_T, int8_T)
FILTERINST(uint16_T, uint16_T)
// FILTERINST(uint16_T, int16_T)
// FILTERINST(uint16_T, int32_T)
// FILTERINST(uint16_T, int64_T)
// FILTERINST(uint16_T, float)
// FILTERINST(uint16_T, double)

// FILTERINST(int16_T, mxLogical)
// FILTERINST(int16_T, uint8_T)
// FILTERINST(int16_T, int8_T)
// FILTERINST(int16_T, uint16_T)
FILTERINST(int16_T, int16_T)
// FILTERINST(int16_T, int32_T)
// FILTERINST(int16_T, int64_T)
// FILTERINST(int16_T, float)
// FILTERINST(int16_T, double)

// FILTERINST(int32_T, mxLogical)
// FILTERINST(int32_T, uint8_T)
// FILTERINST(int32_T, int8_T)
// FILTERINST(int32_T, uint16_T)
// FILTERINST(int32_T, int16_T)
FILTERINST(int32_T, int32_T)
// FILTERINST(int32_T, int64_T)
// FILTERINST(int32_T, float)
// FILTERINST(int32_T, double)

// FILTERINST(int64_T, mxLogical)
// FILTERINST(int64_T, uint8_T)
// FILTERINST(int64_T, int8_T)
// FILTERINST(int64_T, uint16_T)
// FILTERINST(int64_T, int16_T)
// FILTERINST(int64_T, int32_T)
FILTERINST(int64_T, int64_T)
// FILTERINST(int64_T, float)
// FILTERINST(int64_T, double)

// FILTERINST(float, mxLogical)
// FILTERINST(float, uint8_T)
// FILTERINST(float, int8_T)
// FILTERINST(float, uint16_T)
// FILTERINST(float, int16_T)
// FILTERINST(float, int32_T)
// FILTERINST(float, int64_T)
FILTERINST(float, float)
// FILTERINST(float, double)

// FILTERINST(double, mxLogical)
// FILTERINST(double, uint8_T)
// FILTERINST(double, int8_T)
// FILTERINST(double, uint16_T)
// FILTERINST(double, int16_T)
// FILTERINST(double, int32_T)
// FILTERINST(double, int64_T)
// FILTERINST(double, float)
FILTERINST(double, double)

#undef FILTERINST

#endif /* MEXMEDIANIMAGEFILTER_CPP */
