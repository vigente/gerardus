/*
 * MexBinaryErodeImageFilter.cpp
 *
 * Code that is specific to itk::BinaryErodeImageFilter. Support for
 * radius and foreground value arguments. Structuring element is a
 * ball.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
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

#ifndef MEXBINARYERODEIMAGEFILTER_CPP
#define MEXBINARYERODEIMAGEFILTER_CPP

/* C++ headers */
#include <limits>

/* Gerardus headers */
#include "MexBinaryErodeImageFilter.hpp"

/*
 * strings that the user can type to invoke this filter in itk_imfilter()
 */
const std::string MexBinaryErodeImageFilter<std::string, 
		  std::string>::longname = "BinaryErodeImageFilter";
const std::string MexBinaryErodeImageFilter<std::string, 
		  std::string>::shortname = "bwerode";

/* 
 * constructor
 */
template <class InVoxelType, class OutVoxelType>
MexBinaryErodeImageFilter<InVoxelType, OutVoxelType>::MexBinaryErodeImageFilter(
                                const NrrdImage &_nrrd, int _numArgOut, mxArray** _argOut,
				const int _numArgIn, const mxArray** _argIn)  :
  MexBaseFilter<InVoxelType, OutVoxelType>(_nrrd, _numArgOut, _argOut,
					   _numArgIn, _argIn) {

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

  // check number of input parameters
  if (this->numParam < 1) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (this->numParam > 2) {
    mexErrMsgTxt("Too many input arguments");
  }
  if (this->numParam > 0 && this->argParam == NULL) {
    mexErrMsgTxt("Assertion fail: There is at least one parameter, but pointer to parameter array is NULL");
  }
  
  // get user-provided parameters
  this->radius  = this->template 
    GetScalarParamValue<unsigned long>("RADIUS", 0, 0);
  this->foreground = this->template 
    GetScalarParamValue<InVoxelType>("FOREGROUND", 1, std::numeric_limits<InVoxelType>::max());

}

/* 
 * if this particular filter needs to redifine one or more BaseFilter
 * virtual methods, the corresponding definitions go here
 */

template <class InVoxelType, class OutVoxelType>
void MexBinaryErodeImageFilter<InVoxelType, OutVoxelType>::FilterAdvancedSetup() {
  
  // instantiate structuring element
  StructuringElementType structuringElement;
  structuringElement.SetRadius(this->radius);
  structuringElement.CreateStructuringElement();
  derivedFilter->SetKernel(structuringElement);
  
  // pass other parameters to filter
  derivedFilter->SetForegroundValue(this->foreground);

}


/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#define FILTERINST(T1, T2)				\
  template class MexBinaryErodeImageFilter<T1, T2>;

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

#endif /* MEXBINARYERODEIMAGEFILTER_CPP */
