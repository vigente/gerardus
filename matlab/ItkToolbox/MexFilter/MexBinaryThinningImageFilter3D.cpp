/*
 * MexBinaryThinningImageFilter3D.cpp
 *
 * Code that is specific to itk::BinaryThinningImageFilter3D
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

#ifndef MEXBINARYTHINNINGIMAGEFILTER3D_CPP
#define MEXBINARYTHINNINGIMAGEFILTER3D_CPP

#include "MexBinaryThinningImageFilter3D.hpp"

/*
 * strings that the user can use to invoke this filter in itk_imfilter()
 */
const std::string MexBinaryThinningImageFilter3D<std::string, 
		  std::string>::longname = "BinaryThinningImageFilter3D";
const std::string MexBinaryThinningImageFilter3D<std::string, 
		  std::string>::shortname = "skel";

/* 
 * constructor (here we instantiate the filter and process the
 * user-provided input parameters, if any)
 */
template <class InVoxelType, class OutVoxelType>
MexBinaryThinningImageFilter3D<InVoxelType, OutVoxelType>::MexBinaryThinningImageFilter3D(
                                const NrrdImage &_nrrd, 
				int _nargout, mxArray** _argOut,
				const int _nargin, const mxArray** _argIn) :
  MexBaseFilter<InVoxelType, OutVoxelType>(_nrrd, _nargout, _argOut, _nargin, _argIn) {

  // instantiate filter
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
  if (this->nparam < 0) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (this->nparam > 0) {
    mexErrMsgTxt("Too many input arguments");
  }
  if (this->nparam > 0 && this->argParam == NULL) {
    mexErrMsgTxt("Assertion fail: There is at least one parameter, but pointer to parameter array is NULL");
  }
  
  // get user-provided parameters
  // Example: 
  // this->foreground = this->template 
  //   GetScalarParamValue<InVoxelType>("FOREGROUND", 1, std::numeric_limits<InVoxelType>::max());

}

/* 
 * MexBinaryThinningImageFilter3D : MexBaseFilter
 */

/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#define FILTERINST(T1, T2)					\
  template class MexBinaryThinningImageFilter3D<T1, T2>;

FILTERINST(uint8_T, uint8_T)
FILTERINST(int8_T, int8_T)
FILTERINST(uint16_T, uint16_T)
FILTERINST(int16_T, int16_T)
FILTERINST(int32_T, int32_T)
FILTERINST(int64_T, int64_T)
FILTERINST(float, float)
FILTERINST(double, double)

#undef FILTERINST

#endif /* MEXBINARYTHINNINGIMAGEFILTER3D_CPP */
