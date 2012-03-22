/* 
 * MexBaseFilter.cpp
 *
 * MexBaseFilter<InVoxelType, OutVoxelType>: This is where the code to
 * actually run the filter on the image lives.
 *
 * The reason is that template explicit specialization is only
 * possible in classes, not in functions. We need explicit
 * specialization to prevent the compiler from compiling certain
 * input/output image data types for some filters that don't accept
 * them.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.5.2
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

#ifndef MEXBASEFILTER_CPP
#define MEXBASEFILTER_CPP

/* C++ headers */
#include <iostream>
#include <vector>

/* mex headers */
#include <mex.h>

/* ITK headers */

/* Gerardus headers */
#include "GerardusCommon.hpp"
#include "NrrdImage.hpp"
#include "MexBaseFilter.hpp"

/*
 * MexBaseFilter<InVoxelType, OutVoxelType>: This is where
 * the code to actually run the filter on the image lives.
 */

// BaseFilter cannot be invoked by the user, but defining these
// static strings is necessary when we use EXCLUDEFILTER with
// derived filters
const std::string MexBaseFilter<std::string, std::string>::longname = "BaseFilter";
const std::string MexBaseFilter<std::string, std::string>::shortname = "base";

// constructor when the filter takes user-provided parameters
// (e.g. dilation radius) apart from the input arguments with the
// filter type and image
template <class InVoxelType, class OutVoxelType>
MexBaseFilter<InVoxelType, OutVoxelType>::MexBaseFilter(const NrrdImage &_nrrd, 
							int _nargout, mxArray** _argOut,
							const int _nargin, const mxArray** _argIn)
  : nrrd(_nrrd), nargout(_nargout), argOut(_argOut) {

  // check that we have at least filter type and input image as input
  // arguments
  if (_nargin < 2) {
    mexErrMsgTxt("Not enough input arguments");
  }

  // number of extra parameters. Every specific filter has different
  // requirements in terms of this number
  nparam = (mwSize)(_nargin - 2);
  if (nparam < 0) {
    mexErrMsgTxt("Assertion error: Number of parameters cannot be negative");
  } else if (nparam == 0) {
    argParam = NULL;
  } else {
    argParam = &_argIn[2];
  }
}

// destructor to deallocate any memory that the Matlab garbage
// collector won't free automatically. The destructor has to be
// virtual so that derived classes can add their own local memory
// deallocation steps
template <class InVoxelType, class OutVoxelType>
MexBaseFilter<InVoxelType, OutVoxelType>::~MexBaseFilter() {
  
  // assigning NULL to a smart pointer is the correct way of asking
  // for the object to be destroyed. This will deallocate the memory
  // buffers that the filter uses but have not been mummified by us
  this->filter = NULL;
  
}

// check numer of outputs requested by the user. By default, the
// function provides 0 or 1 output (the filtered image), but this
// method can be overriden in child filters with more outputs
template <class InVoxelType, class OutVoxelType>
void MexBaseFilter<InVoxelType, OutVoxelType>::CheckNumberOfOutputs() {
  
  // prevent the user from asking for too many output arguments
  if (nargout > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

}

// put the pointer to the Matlab image buffer into an
// itk::ImportImageFilter. We can feed the ImportImageFilter directly
// to the filter. This way, we don't need to copy the buffer
// to an itk::Image, which saves time, and avoids duplicating the
// image in memory
template <class InVoxelType, class OutVoxelType>
void MexBaseFilter<InVoxelType, OutVoxelType>::GraftMatlabInputBufferIntoItkImportFilter() {
  
  // note that:
  //
  // 1) in ITK we have X,Y,Z indices, while in Matlab we have R,C,S
  //
  // 2) matrices in ITK are read by columns, while in Matlab
  // they are read by rows 
  //
  // So imagine we have this (2, 3) matrix in Matlab, in the NRRD
  //
  //   a b   |
  //   c d   | y-axis (resolution 1.0)
  //   e f   |
  //   ---
  //   x-axis (resolution 0.5)
  //
  //   [nrrd.axis.size] = [3 2 1]
  //
  // The C-style array is going to be (reading by rows)
  //
  //   im = [a c e b d f]
  //
  // ITK is going to read by colums, thinking that the size is
  //
  //   image.size = [sx sy sz] = [3 2 1]
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

  // get pointer to input segmentation mask
  const InVoxelType *im = (InVoxelType *)mxGetData(nrrd.getData());

  // init the filter that will act as an interface between the Matlab
  // image array and the ITK filter
  importFilter = ImportFilterType::New();

  // create ITK image to hold the segmentation mask
  typename ImportFilterType::RegionType region;
  typename ImportFilterType::SizeType size;
  typename ImportFilterType::IndexType start;
  typename ImportFilterType::SpacingType spacing;
  typename InImageType::PointType origin;

  // get image parameters for each dimension
  for (mwIndex i = 0; i < Dimension; ++i) {
    // the region of interest is the whole image
    start[i] = 0;
    size[i] = nrrd.getSize()[i];
    spacing[CAST2MWSIZE(i)] = nrrd.getSpacing()[i];
    origin[CAST2MWSIZE(i)] = nrrd.getMin()[i] 
      + (nrrd.getSpacing()[i] / 2.0); // note that in NRRD, "min" is the
                                      // edge of the voxel, while in
                                      // ITK, "origin" is the centre
                                      // of the voxel

  }
  region.SetIndex(start);
  region.SetSize(size);

  // pass input region parameters to the import filter
  importFilter->SetRegion(region);
  importFilter->SetSpacing(spacing);
  importFilter->SetOrigin(origin);

  // pass pointer to Matlab image to the import filter, and tell it to
  // NOT attempt to delete the memory when it's destructor is
  // called. This is important, because the input image still has to
  // live in Matlab's memory after running the filter
  const bool importImageFilterWillOwnTheBuffer = false;
  importFilter->SetImportPointer(const_cast<InVoxelType *>(im),
				 mxGetNumberOfElements(nrrd.getData()),
				 importImageFilterWillOwnTheBuffer);

  importFilter->Update();
  
}

// filter setup code common to all filters: pass image to the
// filter, allocate memory for the Matlab output, and graft the
// Matlab output into the filter output
template <class InVoxelType, class OutVoxelType>
void MexBaseFilter<InVoxelType, OutVoxelType>::FilterBasicSetup() {

  // pass input image to filter
  this->filter->SetInput(this->importFilter->GetOutput());

  // link the filter main output (the filtered image) to the first
  // Matlab output buffer
  this->template MallocMatlabOutputBuffer<OutVoxelType>(0);
  this->template GraftMatlabOutputBufferIntoItkFilterOutput<OutVoxelType>(0);

}

// by default, this method doesn't do anything, but can be overriden
// when a filter needs some extra setput steps (e.g. passing parameters)
template <class InVoxelType, class OutVoxelType>
void MexBaseFilter<InVoxelType, OutVoxelType>::FilterAdvancedSetup() {

}

// filter the image
template <class InVoxelType, class OutVoxelType>
void MexBaseFilter<InVoxelType, OutVoxelType>::RunFilter() {
  
  // run filter
  filter->Update();
  
}

// prevent the C++ destructor from deleting the data of a filter
// output when program returns. This is necessary in order to use
// that output from Matlab
template <class InVoxelType, class OutVoxelType>
void MexBaseFilter<InVoxelType,
		   OutVoxelType>::MummifyFilterOutput(unsigned int idx) {

  // mummify filter output buffer for Matlab
  typename OutImageType::PixelContainer * container;
  container = this->filter->GetOutput(idx)->GetPixelContainer();
  container->SetContainerManageMemory(false);
  
}

// by default, this method doesn't do anything, but can be overriden
// when a child filter provides other outputs apart from the
// filtered image
template <class InVoxelType, class OutVoxelType>
void MexBaseFilter<InVoxelType, OutVoxelType>::ExportOtherFilterOutputsToMatlab() {
  
}

/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#define FILTERINST(T1, T2)			\
  template class MexBaseFilter<T1, T2>;

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

#endif /* MEXBASEFILTER_CPP */
