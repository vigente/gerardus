/*
 * MexSignedMaurerDistanceMapImageFilter.cpp
 *
 * Code that is specific to itk::SignedMaurerDistanceMapImageFilter
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.3.1
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

#ifndef MEXSIGNEDMAURERDISTANCEMAPIMAGEFILTER_CPP
#define MEXSIGNEDMAURERDISTANCEMAPIMAGEFILTER_CPP

#include "MexSignedMaurerDistanceMapImageFilter.hpp"

/*
 * strings that the user can use to invoke this filter in itk_imfilter()
 */
const std::string MexSignedMaurerDistanceMapImageFilter<std::string, 
		  std::string>::longname = "SignedMaurerDistanceMapImageFilter";
const std::string MexSignedMaurerDistanceMapImageFilter<std::string, 
		  std::string>::shortname = "maudist";

/* 
 * constructor (here we instantiate the filter and process the
 * user-provided input parameters, if any)
 */
template <class InVoxelType, class OutVoxelType>
MexSignedMaurerDistanceMapImageFilter<InVoxelType, OutVoxelType>::MexSignedMaurerDistanceMapImageFilter(
                                const NrrdImage &_nrrd, 
				int _numArgOut, mxArray** _argOut,
				const int _numArgIn, const mxArray** _argIn) :
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

  // check number of user-provided parameters (user-provided
  // parameters are the extra input arguments apart from the filter
  // type and input image)
  if (this->numParam < 0) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (this->numParam > 0) {
    mexErrMsgTxt("Too many input arguments");
  }
  if (this->numParam > 0 && this->argParam == NULL) {
    mexErrMsgTxt("Assertion fail: There is at least one parameter, but pointer to parameter array is NULL");
  }
  
}

/* 
 * overriding of MexBaseFilter virtual methods, if needed
 */

template <class InVoxelType, class OutVoxelType>
void MexSignedMaurerDistanceMapImageFilter<InVoxelType, 
					   OutVoxelType>::FilterAdvancedSetup() {
  
  // compute distances using real world coordinates, instead of voxel
  // indices
  derivedFilter->SetUseImageSpacing(true);

  // give output as actual distances, instead of squared distances
  derivedFilter->SquaredDistanceOff();
}

/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#define FILTERINST(T1, T2)					\
  template class MexSignedMaurerDistanceMapImageFilter<T1, T2>;

FILTERINST(mxLogical, double)
FILTERINST(uint8_T, double)
FILTERINST(int8_T, double)
FILTERINST(uint16_T, double)
FILTERINST(int16_T, double)
FILTERINST(int32_T, double)
FILTERINST(int64_T, double)
FILTERINST(float, double)
FILTERINST(double, double)

#undef FILTERINST

#endif /* MEXSIGNEDMAURERDISTANCEMAPIMAGEFILTER_CPP */
