/*
 * DanielssonFilter.cpp
 *
 * Code that is specific to the DanielssonDistanceMapImageFilter
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.3.0
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

#ifndef DANIELSSONFILTER_CPP
#define DANIELSSONFILTER_CPP

/* ITK headers */
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"

/* Gerardus headers */
#include "GerardusCommon.hpp"
#include "DanielssonFilter.hpp"

template <class InVoxelType, class OutVoxelType>
void DanielssonFilter<InVoxelType,
		      OutVoxelType>::CopyAllFilterOutputsToMatlab() {
  
  // by default, we assume that all filters produce at least 1 main
  // output
  this->CopyFilterImageOutputToMatlab();

  if (this->nargout <= 2) {
    CopyFilterNearestOutputToMatlab();
  } else {
    mexErrMsgTxt("Too many output arguments");
  }

}

/*
 * DanielssonFilter::CopyFilterNearestOutputToMatlab()
 *
 * Pass to Matlab an array of the same size as the image. Each element
 * has the linear index of the closest object voxel to each image
 * voxel. The distance between both voxels is the distance returned in
 * the distance map. This is a reformatting of the output provided by
 * itk::DanielssonDistanceMapImageFilter::GetVectorDistanceMap()
 */
template <class InVoxelType, class OutVoxelType>
void DanielssonFilter<InVoxelType, 
		      OutVoxelType>::CopyFilterNearestOutputToMatlab() {

  typedef double OutType;

  // if the input image is empty, create empty segmentation mask for
  // output, and we don't need to do any further processing
  if (this->nrrd.getR() == 0 || this->nrrd.getC() == 0) {
    this->argOut[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
    return;
  }

  // output class ID
  // convert output data type to output class ID
  mxClassID outputVoxelClassId = mxUNKNOWN_CLASS;
  if (TypeIsBool< OutType >::value) {
    outputVoxelClassId = mxLOGICAL_CLASS;
  } else if (TypeIsUint8< OutType >::value) {
    outputVoxelClassId = mxUINT8_CLASS;
  } else if (TypeIsUint16< OutType >::value) {
    outputVoxelClassId = mxUINT16_CLASS;
  } else if (TypeIsFloat< OutType >::value) {
    outputVoxelClassId = mxSINGLE_CLASS;
  } else if (TypeIsDouble< OutType >::value) {
    outputVoxelClassId = mxDOUBLE_CLASS;
  } else {
    mexErrMsgTxt("Assertion fail: Unrecognised output voxel type");
  }

  // create output matrix for Matlab's result
  this->argOut[1] = (mxArray *)mxCreateNumericArray( this->nrrd.getNdim(), 
						     this->nrrd.getDims(),
						     outputVoxelClassId,
						     mxREAL);
  if (this->argOut[1] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output matrix");
  }
  OutType *imOutp =  (OutType *)mxGetPr(this->argOut[1]);
  
  // populate output image
  typedef typename DanielssonFilter<InVoxelType, 
    OutVoxelType>::FilterType::VectorImageType OffsetImageType;

  // the filter member variable is declared in BaseFilter as a general
  // ImageToImageFilter, but we want to use some methods that belong
  // only to the derived filter class
  // DanielssonDistanceMapImageFilter. In order to do this, we need
  // to declare a local filter variable that is of type
  // DanielssonDistanceMapImageFilter, and dynamic cast it to filter
  // in the BaseFilter class
  typename FilterType::Pointer localFilter = 
    dynamic_cast<typename DanielssonFilter<InVoxelType, 
    OutVoxelType>::FilterType *>(this->filter.GetPointer());

  typedef itk::ImageRegionConstIterator<OffsetImageType> 
    OutConstIteratorType;

  OutConstIteratorType citer(localFilter->GetVectorDistanceMap(),
	     localFilter->GetVectorDistanceMap()->GetLargestPossibleRegion());

  itk::Offset<Dimension> idx3;
  mwIndex i = 0; // voxel linear index
  for (citer.GoToBegin(), i = 0; !citer.IsAtEnd(); ++citer, ++i) {

    // current voxel in the image
    // image linear index => r, c, s indices
    idx3 = ind2sub_itkOffset(this->nrrd.getR(), this->nrrd.getC(), 
			     this->nrrd.getS(), i);

    // compute coordinates of the closest object voxel to the current
    // voxel
    idx3 += citer.Get();

    // convert image r, c, s indices => linear index
    // Note: we need to add 1 to the index to follow Matlab's convention
    imOutp[i] = sub2ind(this->nrrd.getR(), this->nrrd.getC(), 
			this->nrrd.getS(), idx3) + 1;
    
  }

}

/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#define FILTERINST(T1, T2)						\
  template class DanielssonFilter<T1, T2>;				\

FILTERINST(bool, bool);
FILTERINST(bool, uint8_T)
FILTERINST(bool, uint16_T)
FILTERINST(bool, float)
FILTERINST(bool, double)

FILTERINST(uint8_T, bool);
FILTERINST(uint8_T, uint8_T)
FILTERINST(uint8_T, uint16_T)
FILTERINST(uint8_T, float)
FILTERINST(uint8_T, double)

FILTERINST(int8_T, bool);
FILTERINST(int8_T, uint8_T)
FILTERINST(int8_T, uint16_T)
FILTERINST(int8_T, float)
FILTERINST(int8_T, double)

FILTERINST(uint16_T, bool);
FILTERINST(uint16_T, uint8_T)
FILTERINST(uint16_T, uint16_T)
FILTERINST(uint16_T, float)
FILTERINST(uint16_T, double)

FILTERINST(int16_T, bool);
FILTERINST(int16_T, uint8_T)
FILTERINST(int16_T, uint16_T)
FILTERINST(int16_T, float)
FILTERINST(int16_T, double)

FILTERINST(int32_T, bool);
FILTERINST(int32_T, uint8_T)
FILTERINST(int32_T, uint16_T)
FILTERINST(int32_T, float)
FILTERINST(int32_T, double)

FILTERINST(int64_T, bool);
FILTERINST(int64_T, uint8_T)
FILTERINST(int64_T, uint16_T)
FILTERINST(int64_T, float)
FILTERINST(int64_T, double)

FILTERINST(float, bool);
FILTERINST(float, uint8_T)
FILTERINST(float, uint16_T)
FILTERINST(float, float)
FILTERINST(float, double)

FILTERINST(double, bool);
FILTERINST(double, uint8_T)
FILTERINST(double, uint16_T)
FILTERINST(double, float)
FILTERINST(double, double)

#undef FILTERINST

#endif /* DANIELSSONFILTER_CPP */
