/*
 * MexBinaryDilateImageFilter.cpp
 *
 * Code that is specific to itk::BinaryDilateImageFilter. Support for
 * radius and foreground value arguments. Structuring element is a
 * ball.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.0.2
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

#ifndef MEXBINARYDILATEIMAGEFILTER_CPP
#define MEXBINARYDILATEIMAGEFILTER_CPP

/* C++ headers */
#include <limits>

/* Gerardus headers */
#include "MexBinaryDilateImageFilter.hpp"

/* 
 * MexBinaryDilateImageFilter::FilterSetup
 */
template <class InVoxelType, class OutVoxelType>
void MexBinaryDilateImageFilter<InVoxelType, 
			OutVoxelType>::FilterSetup() {
  
  // check input parameters
  if (this->nparam > 2) {
    mexErrMsgTxt("Too many input arguments");
  }
  if (this->nparam > 0 && this->argParam == NULL) {
    mexErrMsgTxt("Assertion fail: There is at least one parameter, but pointer to parameter array is NULL");
  }
  
  // list of user-provided parameters for this filter (compulsory or
  // optional)
  unsigned long radius;   // (comp) radius of the ball in voxels
  InVoxelType foreground; // (opt) voxels with this value will be
                          // dilated. Default, maximum value of the
                          // pixel type

  // get parameters
  if (this->nparam < 1) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (this->nparam > 2) {
    mexErrMsgTxt("Too many input arguments");
  }
  radius  = this->template 
    getScalarParamValue<unsigned long>("RADIUS", 0, 0);
  foreground = this->template 
    getScalarParamValue<InVoxelType>("FOREGROUND", 1, std::numeric_limits<InVoxelType>::max());
  
  // create a local pointer to the filter so that we can use
  // methods that are not part of the MexBaseFilter
  typename FilterType::Pointer localFilter = 
    dynamic_cast<typename MexBinaryDilateImageFilter<InVoxelType,
			     OutVoxelType>::FilterType *>(this->filter.GetPointer());
  
  // instantiate structuring element
  StructuringElementType structuringElement;
  structuringElement.SetRadius(radius);
  structuringElement.CreateStructuringElement();
  localFilter->SetKernel(structuringElement);
  
  // pass parameters to filter
  localFilter->SetForegroundValue(foreground);

  // pass image to filter
  localFilter->SetInput(this->image);
  
}

/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#define FILTERINST(T1, T2)				\
  template class MexBinaryDilateImageFilter<T1, T2>;

FILTERINST(mxLogical, mxLogical)
FILTERINST(uint8_T, uint8_T)
FILTERINST(int8_T, int8_T)
FILTERINST(uint16_T, uint16_T)
FILTERINST(int16_T, int16_T)
FILTERINST(int32_T, int32_T)
FILTERINST(int64_T, int64_T)
FILTERINST(float, float)
FILTERINST(double, double)

#undef FILTERINST

#endif /* MEXBINARYDILATEIMAGEFILTER_CPP */
