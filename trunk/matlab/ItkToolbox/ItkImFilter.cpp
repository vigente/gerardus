/* ItkImFilter.cpp
 *
 * ITK_IMFILTER: Run ITK filter on a 2D or 3D image
 *
 * This MEX function is a multiple-purpose wrapper to be able to run
 * all ITK filters that inherit from itk::ImageToImageFilter on a
 * Matlab 2D image or 3D image volume.
 *
 * B = ITK_IMFILTER(TYPE, A)
 *
 *   TYPE is a string with the filter we want to run. Currently, only
 *   the following options are implemented:
 *
 *     'skel': (BinaryThinningImageFilter3D) Skeletonize a
 *             segmentation mask
 *
 *   A is a 2D matrix or 3D volume with the image or
 *   segmentation. Currently, A can be of any of the following
 *   Matlab classes:
 *
 *     boolean
 *     double
 *     single
 *     int8
 *     uint8
 *     int16
 *     uint16
 *     int32
 *     int64
 *
 *   B has the same size and class as A, and contains the filtered
 *   image or segmentation mask.
 *
 * This function must be compiled before it can be used from Matlab.
 * If Gerardus' root directory is e.g. ~/gerardus, type from a
 * linux shell
 *
 *    $ cd ~/gerardus/matlab
 *    $ mkdir bin
 *    $ cd bin
 *    $ cmake ..
 *    $ make install
 *
 * If cmake throws an error because it cannot find Matlab, then edit
 * gerardus/matlab/CMakeLists.txt, and where it says
 *
 *    SET(MATLAB_ROOT "/usr/local/matlab/R2010b/")
 *
 *  change to your own Matlab root path.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.1.0
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

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

/* mex headers */
#include <math.h>
#include <matrix.h>
#include <mex.h>

/* C++ headers */
#include <iostream>

/* ITK headers */
#include "itkImage.h"
#include "itkBinaryThinningImageFilter3D.h"

#ifndef ITKIMFILTER_CPP
#define ITKIMFILTER_CPP

/*
 * runFilter(): function in charge of processing. We cannot do this in
 *              the body of mexFunction() because we need to template
 *              the function so that we can operate with different
 *              types of Matlab matrices (double, uint8, etc)
 *
 * filter: string with the name of the filter we want to run
 * imIn:   input image, read-only
 * imOut:  output image, taken by reference so that we can allocate
 *         memory and populate it from here
 */
template <class VoxelType>
void runFilter(char *filterType, const mxArray *imIn, mxArray* &imOut) {

  // get number of dimensions in the input segmentation
  if (!mxIsNumeric(imIn) && !mxIsLogical(imIn)) {
    mexErrMsgTxt("IM must be a 2D or 3D numeric or boolean matrix");
  }
  mwSize ndim = mxGetNumberOfDimensions(imIn);

  // get size of input arguments
  mwSize Rin = 0; // rows
  mwSize Cin = 0; // cols
  mwSize Sin = 0; // slices
  const mwSize *dims = mxGetDimensions(imIn);
  Rin = dims[0];
  if (ndim > 3) {
    mexErrMsgTxt("Input segmentation mask must be 2D or 3D");
  } else if (ndim == 2) {
    Cin = dims[1];
    Sin = 1;
  } else if (ndim == 3) {
    Cin = dims[1];
    Sin = dims[2];
  } else {
    mexErrMsgTxt("Assertion fail: number of dimensions is " + ndim);
  }

  // create empty segmentation mask for output
  if (Rin == 0 || Cin == 0) {
    imOut = mxCreateDoubleMatrix(0, 0, mxREAL);
    return;
  }

  // get pointer to input segmentation mask
  const VoxelType *im = (VoxelType *)mxGetPr(imIn);
  
  // image type definitions
  static const unsigned int Dimension = 3; // volume data dimension
                                           // (3D volume)
  typedef double TScalarType; // data type for scalars
  typedef itk::Image< VoxelType, Dimension > 
    ImageType;
  typedef itk::ImageRegionIterator< ImageType > 
    IteratorType;
  typedef itk::ImageRegionConstIterator< ImageType > 
    ConstIteratorType;

  // create ITK image to hold the segmentation mask
  typename ImageType::Pointer image = ImageType::New();
  typename ImageType::RegionType region;
  typename ImageType::IndexType start;
  typename ImageType::SizeType size; 

  // note that:
  //
  // 1) in ITK we have X,Y,Z indices, while in Matlab we have R,C,S
  //
  // 2) matrices in ITK are read by columns, while in Matlab
  // they are read by rows, 
  //
  // So both things cancel each other.
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  size[0] = Rin;
  size[1] = Cin;
  size[2] = Sin;
  region.SetIndex(start);
  region.SetSize(size);
  image->SetRegions(region);
  image->Allocate();
  image->Update();

  // loop through every voxel in the image and copy it to the ITK image
  IteratorType iter(image, image->GetLargestPossibleRegion());
  mwIndex i = 0; // index for the input image
  for (iter.GoToBegin(), i=0; !iter.IsAtEnd(); ++iter, ++i) {
    iter.Set(im[i]);
    
    // Note: If you need to debug this loop, take into account that
    //   std::cout << im[i] << std::endl;
    // doesn't work as expected if im is e.g. a 
    // (uint8_T *) pointer. Instead, you need to do
    //   std::cout << (double)im[i] << std::endl;

  }

  // filter image
  typedef itk::ImageToImageFilter< ImageType, ImageType > 
    FilterType;
  typename FilterType::Pointer filter;
  if (!strcmp(filterType, "skel")) {
    // skeletonize image
    typedef itk::BinaryThinningImageFilter3D< ImageType, ImageType > 
      ThinningFilterType;
    filter = ThinningFilterType::New();
    filter->SetInput(image);
    filter->Update();
  } else {
    mexErrMsgTxt("Filter not implemented");
  }
  
  // create output matrix for Matlab's result
  imOut = (mxArray *)mxCreateNumericArray(ndim, dims,
					    mxGetClassID(imIn),
					    mxREAL);
  if (imOut == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output matrix");
  }
  VoxelType *imOutp =  (VoxelType *)mxGetPr(imOut);

  // populate output image
  ConstIteratorType citer(filter->GetOutput(), 
			  filter->GetOutput()->GetLargestPossibleRegion());
  for (citer.GoToBegin(), i=0; i < Rin*Cin*Sin; ++citer, ++i) {
    imOutp[i] = (VoxelType)citer.Get();
  }

  return;
}

// entry point for the mex function
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {
  // check number of input and output arguments
  if (nrhs != 2) {
    mexErrMsgTxt("Two input arguments required");
  }
  else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  // get type of filter
  char *filter = mxArrayToString(prhs[0]);
  if (filter == NULL) {
    mexErrMsgTxt("Invalid FILTER string");
  }

  // run filter, templated according to the input matrix type
  mxClassID  ExternalVoxelClassId  = mxGetClassID(prhs[1]);
  switch(ExternalVoxelClassId)  {
  case mxLOGICAL_CLASS:
    runFilter<bool>(filter, prhs[1], plhs[0]);
    break;
  case mxDOUBLE_CLASS:
    runFilter<double>(filter, prhs[1], plhs[0]);
    break;
  case mxSINGLE_CLASS:
    runFilter<float>(filter, prhs[1], plhs[0]);
    break;
  case mxINT8_CLASS:
    runFilter<int8_T>(filter, prhs[1], plhs[0]);
    break;
  case mxUINT8_CLASS:
    runFilter<uint8_T>(filter, prhs[1], plhs[0]);
    break;
  case mxINT16_CLASS:
    runFilter<int16_T>(filter, prhs[1], plhs[0]);
    break;
  case mxUINT16_CLASS:
    runFilter<uint16_T>(filter, prhs[1], plhs[0]);
    break;
  case mxINT32_CLASS:
    runFilter<int32_T>(filter, prhs[1], plhs[0]);
    break;
  // case mxUINT32_CLASS:
  //   runFilter<uint32_T>(filter, prhs[1], plhs[0]);
  //   break;
  case mxINT64_CLASS:
    runFilter<int64_T>(filter, prhs[1], plhs[0]);
    break;
  // case mxUINT64_CLASS:
  //   runFilter<uint64_T>(filter, prhs[1], plhs[0]);
  //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt("Input matrix has unknown type.");
    break;
  default:
    mexErrMsgTxt("Input matrix has invalid type.");
    break;
  }

  // exit successfully
  return;

}

#endif /* ITKIMFILTER_CPP */
