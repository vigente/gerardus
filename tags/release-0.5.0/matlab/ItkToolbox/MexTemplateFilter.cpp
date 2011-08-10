/*
 * MexTemplateImageFilter.cpp
 *
 * Code that is specific to itk::TemplateImageFilter.
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

#ifndef MEXTEMPLATEIMAGEFILTER_CPP
#define MEXTEMPLATEIMAGEFILTER_CPP

/* C++ headers */

/* Gerardus headers */
#include "MexTemplateImageFilter.hpp"

/*
 * strings that the user can type to invoke this filter in itk_imfilter()
 */
const std::string MexTemplateImageFilter<std::string, 
		  std::string>::longname = "TemplateImageFilter";
const std::string MexTemplateImageFilter<std::string, 
		  std::string>::shortname = "template";

/* 
 * if this particular filter needs to redifine one or more BaseFilter
 * virtual methods, the corresponding definitions go here
 */

// template <class InVoxelType, class OutVoxelType>
// void MexTemplateImageFilter<InVoxelType, OutVoxelType>::FilterSetup() {
  
//   // pass image to filter
//   localFilter->SetInput(this->image);
  
// }

// template <class InVoxelType, class OutVoxelType>
// void MexTemplateFilter<InVoxelType, OutVoxelType>::RunFilter() {
  
//   // run filter
//   filter->Update();
  
// }

// template <class InVoxelType, class OutVoxelType>
// void MexTemplateFilter<InVoxelType, OutVoxelType>::CopyAllFilterOutputsToMatlab() {
  
//   // by default, we assume that all filters produce at least 1 main
//   // output
//   this->CopyFilterImageOutputToMatlab();

//   // prevent the user from asking for too many output arguments
//   if (nargout > 1) {
//     mexErrMsgTxt("Too many output arguments");
//   }

// }

/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#error FILTERINST types cannot be automatically determined by add_itk_imfilter_template.sh
#define FILTERINST(T1, T2)				\
  template class MexTemplateImageFilter<T1, T2>;

FILTERINST(mxLogical, mxLogical)
FILTERINST(mxLogical, uint8_T)
FILTERINST(mxLogical, int8_T)
FILTERINST(mxLogical, uint16_T)
FILTERINST(mxLogical, int16_T)
FILTERINST(mxLogical, int32_T)
FILTERINST(mxLogical, int64_T)
FILTERINST(mxLogical, float)
FILTERINST(mxLogical, double)

FILTERINST(uint8_T, mxLogical)
FILTERINST(uint8_T, uint8_T)
FILTERINST(uint8_T, int8_T)
FILTERINST(uint8_T, uint16_T)
FILTERINST(uint8_T, int16_T)
FILTERINST(uint8_T, int32_T)
FILTERINST(uint8_T, int64_T)
FILTERINST(uint8_T, float)
FILTERINST(uint8_T, double)

FILTERINST(int8_T, mxLogical)
FILTERINST(int8_T, uint8_T)
FILTERINST(int8_T, int8_T)
FILTERINST(int8_T, uint16_T)
FILTERINST(int8_T, int16_T)
FILTERINST(int8_T, int32_T)
FILTERINST(int8_T, int64_T)
FILTERINST(int8_T, float)
FILTERINST(int8_T, double)

FILTERINST(uint16_T, mxLogical)
FILTERINST(uint16_T, uint8_T)
FILTERINST(uint16_T, int8_T)
FILTERINST(uint16_T, uint16_T)
FILTERINST(uint16_T, int16_T)
FILTERINST(uint16_T, int32_T)
FILTERINST(uint16_T, int64_T)
FILTERINST(uint16_T, float)
FILTERINST(uint16_T, double)

FILTERINST(int16_T, mxLogical)
FILTERINST(int16_T, uint8_T)
FILTERINST(int16_T, int8_T)
FILTERINST(int16_T, uint16_T)
FILTERINST(int16_T, int16_T)
FILTERINST(int16_T, int32_T)
FILTERINST(int16_T, int64_T)
FILTERINST(int16_T, float)
FILTERINST(int16_T, double)

FILTERINST(int32_T, mxLogical)
FILTERINST(int32_T, uint8_T)
FILTERINST(int32_T, int8_T)
FILTERINST(int32_T, uint16_T)
FILTERINST(int32_T, int16_T)
FILTERINST(int32_T, int32_T)
FILTERINST(int32_T, int64_T)
FILTERINST(int32_T, float)
FILTERINST(int32_T, double)

FILTERINST(int64_T, mxLogical)
FILTERINST(int64_T, uint8_T)
FILTERINST(int64_T, int8_T)
FILTERINST(int64_T, uint16_T)
FILTERINST(int64_T, int16_T)
FILTERINST(int64_T, int32_T)
FILTERINST(int64_T, int64_T)
FILTERINST(int64_T, float)
FILTERINST(int64_T, double)

FILTERINST(float, mxLogical)
FILTERINST(float, uint8_T)
FILTERINST(float, int8_T)
FILTERINST(float, uint16_T)
FILTERINST(float, int16_T)
FILTERINST(float, int32_T)
FILTERINST(float, int64_T)
FILTERINST(float, float)
FILTERINST(float, double)

FILTERINST(double, mxLogical)
FILTERINST(double, uint8_T)
FILTERINST(double, int8_T)
FILTERINST(double, uint16_T)
FILTERINST(double, int16_T)
FILTERINST(double, int32_T)
FILTERINST(double, int64_T)
FILTERINST(double, float)
FILTERINST(double, double)

#undef FILTERINST

#endif /* MEXTEMPLATEIMAGEFILTER_CPP */
