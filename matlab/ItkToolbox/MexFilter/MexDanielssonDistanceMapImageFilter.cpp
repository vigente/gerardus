/*
 * MexDanielssonDistanceMapImageFilter.cpp
 *
 * Code that is specific to itk::DanielssonDistanceMapImageFilter
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.5.1
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

#ifndef MEXDANIELSSONDISTANCEMAPIMAGEFILTER_CPP
#define MEXDANIELSSONDISTANCEMAPIMAGEFILTER_CPP

/* ITK headers */
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"

/* Gerardus headers */
#include "GerardusCommon.hpp"
#include "MexDanielssonDistanceMapImageFilter.hpp"

/*
 * strings that the user can use to invoke this filter in itk_imfilter()
 */
const std::string MexDanielssonDistanceMapImageFilter<std::string, 
		  std::string>::longname = "DanielssonDistanceMapImageFilter";
const std::string MexDanielssonDistanceMapImageFilter<std::string, 
		  std::string>::shortname = "dandist";

/* 
 * constructor (here we instantiate the filter and process the
 * user-provided input parameters, if any)
 */
template <class InVoxelType, class OutVoxelType>
MexDanielssonDistanceMapImageFilter<InVoxelType, OutVoxelType>::MexDanielssonDistanceMapImageFilter(
                                const NrrdImage &_nrrd, 
				int _nargout, mxArray** _argOut,
				const int _nargin, const mxArray** _argIn) :
  MexBaseFilter<InVoxelType, OutVoxelType>(_nrrd, _nargout, _argOut,
					   _nargin, _argIn) {

  // instantiate filter in this derived class, but on the base class
  // pointer, thanks to polimorphism. This way, we can run methods on
  // the derived class from the base class
  this->filter = DerivedImageToImageFilterType::New();

  // std::cout << "this->filter 0 = " 
  // 	    << this->filter->GetOutput(0) << std::endl;//////////
  // std::cout << "this->filter 1 = " 
  // 	    << this->filter->GetOutput(1) << std::endl;//////////
  // std::cout << "this->filter 2 = " 
  // 	    << this->filter->GetOutput(2) << std::endl;//////////

  // std::cout << "this->filter NumberOfOutputs = " 
  // 	    << this->filter->GetNumberOfOutputs() << std::endl;//////////

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
 * Definition of methods from BaseFilter that this filter needs to
 * override
 */

template <class InVoxelType, class OutVoxelType>
void MexDanielssonDistanceMapImageFilter<InVoxelType, 
					 OutVoxelType>::CheckNumberOfOutputs() {
  
  // prevent the user from asking for too many output arguments
  if (this->nargout > 2) {
    mexErrMsgTxt("Too many output arguments");
  }

}

template <class InVoxelType, class OutVoxelType>
void MexDanielssonDistanceMapImageFilter<InVoxelType,
		      OutVoxelType>::ExportOtherFilterOutputsToMatlab() {

  // convert the 3-vector format to index of nearest segmented
  // voxel. This way, we can give the output as a matrix of the same
  // size as the input
  if (this->nargout > 1) {
    CopyFilterNearestOutputToMatlab();

    // link the filter's second output to ITK's second output
    // typedef typename MexBaseFilter<InVoxelType, OutVoxelType>::InImageType::OffsetType OffsetType;
    // typedef typename OffsetType::OffsetValueType OffsetValueType;
    // std::cout << "OffsetType = " << print_T<OffsetType>() << std::endl;
    // std::cout << "OffsetValueType = " << print_T<typename OffsetType::OffsetValueType>() << std::endl;

    // this->template MallocMatlabOutputBuffer<OffsetValueType, Dimension>(2);
    // std::cout << "before 0 = " 
    // 	      << derivedFilter->GetOutput(0) << std::endl;//////////
    // std::cout << "before 1 = " 
    // 	      << derivedFilter->GetOutput(1) << std::endl;//////////
    // std::cout << "before 2 = " 
    // 	      << derivedFilter->GetOutput(2) << std::endl;//////////

    // derivedFilter->itk::ProcessObject::GetOutput(2);

    // std::cout << "ProcessObject 2 = " 
    // 	      << derivedFilter->ProcessObject::GetOutput(2) << std::endl;//////////

  //   std::cout << "derivedFilter GetVectorDistanceMap = " 
  // 	      << derivedFilter->GetVectorDistanceMap() << std::endl;//////////

  //   std::cout << "derivedFilter NumberOfOutputs = " 
  // 	      << derivedFilter->GetNumberOfOutputs() << std::endl;//////////

  //   this->template foo<OffsetValueType, Dimension>(2, 
  // 	  dynamic_cast<typename MexBaseFilter<InVoxelType, OutVoxelType>::ImageToImageFilterType *>(derivedFilter.GetPointer()));
  }

}

/*
 * MexDanielssonDistanceMapImageFilter::CopyFilterNearestOutputToMatlab()
 *
 * Pass to Matlab an array of the same size as the image. Each element
 * has the linear index of the closest object voxel to each image
 * voxel. The distance between both voxels is the distance returned in
 * the distance map. This is a reformatting of the output provided by
 * itk::DanielssonDistanceMapImageFilter::GetVectorDistanceMap()
 */
template <class InVoxelType, class OutVoxelType>
void MexDanielssonDistanceMapImageFilter<InVoxelType, 
					 OutVoxelType>::CopyFilterNearestOutputToMatlab() {

  typedef double OutputType;

  // allocate memory for the output buffer
  this->template MallocMatlabOutputBuffer<OutputType>(1);

  OutputType *imOutp =  (OutputType *)mxGetData(this->argOut[1]);
  
  // populate output image
  typedef typename MexDanielssonDistanceMapImageFilter<InVoxelType, 
    OutVoxelType>::DerivedImageToImageFilterType::VectorImageType OffsetImageType;

  typedef itk::ImageRegionConstIterator<OffsetImageType> 
    OutConstIteratorType;

  OutConstIteratorType citer(derivedFilter->GetVectorDistanceMap(),
	     derivedFilter->GetVectorDistanceMap()->GetLargestPossibleRegion());

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

#define FILTERINST(T1, T2)					\
  template class MexDanielssonDistanceMapImageFilter<T1, T2>;

FILTERINST(mxLogical, mxLogical);
FILTERINST(mxLogical, uint8_T)
FILTERINST(mxLogical, uint16_T)
FILTERINST(mxLogical, float)
FILTERINST(mxLogical, double)

FILTERINST(uint8_T, mxLogical);
FILTERINST(uint8_T, uint8_T)
FILTERINST(uint8_T, uint16_T)
FILTERINST(uint8_T, float)
FILTERINST(uint8_T, double)

FILTERINST(int8_T, mxLogical);
FILTERINST(int8_T, uint8_T)
FILTERINST(int8_T, uint16_T)
FILTERINST(int8_T, float)
FILTERINST(int8_T, double)

FILTERINST(uint16_T, mxLogical);
FILTERINST(uint16_T, uint8_T)
FILTERINST(uint16_T, uint16_T)
FILTERINST(uint16_T, float)
FILTERINST(uint16_T, double)

FILTERINST(int16_T, mxLogical);
FILTERINST(int16_T, uint8_T)
FILTERINST(int16_T, uint16_T)
FILTERINST(int16_T, float)
FILTERINST(int16_T, double)

FILTERINST(int32_T, mxLogical);
FILTERINST(int32_T, uint8_T)
FILTERINST(int32_T, uint16_T)
FILTERINST(int32_T, float)
FILTERINST(int32_T, double)

FILTERINST(int64_T, mxLogical);
FILTERINST(int64_T, uint8_T)
FILTERINST(int64_T, uint16_T)
FILTERINST(int64_T, float)
FILTERINST(int64_T, double)

FILTERINST(float, mxLogical);
FILTERINST(float, uint8_T)
FILTERINST(float, uint16_T)
FILTERINST(float, float)
FILTERINST(float, double)

FILTERINST(double, mxLogical);
FILTERINST(double, uint8_T)
FILTERINST(double, uint16_T)
FILTERINST(double, float)
FILTERINST(double, double)

#undef FILTERINST

#endif /* MEXDANIELSSONDISTANCEMAPIMAGEFILTER_CPP */
