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
 *     'skel':    (BinaryThinningImageFilter3D) Skeletonize a
 *                binary mask
 *
 *                B has the same size and class as A
 *
 *     'dandist': (DanielssonDistanceMapImageFilter) Compute unsigned
 *                distance map for a binary mask. Distance values are
 *                given in voxel coordinates
 *
 *                B has the same size as A. B has a type large enough
 *                to store the maximum distance in the image. The
 *                largest available type is double. If this is not
 *                enough, a warning message is displayed, and double
 *                is used as the output type
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
 * Cmake has equivalents for Windows and MacOS, but I have not tried
 * them.
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
  * Version: 0.2.0
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
#include "itkDanielssonDistanceMapImageFilter.h"

#ifndef ITKIMFILTER_CPP
#define ITKIMFILTER_CPP

/*
 * TypeIsUint16(): Block of functions to allow testing of template
 *                 types
 */
template< class T >
struct TypeIsBool
{ static const bool value = false; };

template<>
struct TypeIsBool< bool >
{ static const bool value = true; };

template< class T >
struct TypeIsUint8
{ static const bool value = false; };

template<>
struct TypeIsUint8< uint8_T >
{ static const bool value = true; };

template< class T >
struct TypeIsUint16
{ static const bool value = false; };

template<>
struct TypeIsUint16< uint16_T >
{ static const bool value = true; };

template< class T >
struct TypeIsFloat
{ static const bool value = false; };

template<>
struct TypeIsFloat< float >
{ static const bool value = true; };

template< class T >
struct TypeIsDouble
{ static const bool value = false; };

template<>
struct TypeIsDouble< double >
{ static const bool value = true; };


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

/* 
 * runFilter<VoxelType> is for filters that produce an output image
 * with the same type as the input image (e.g. thinning algorithm)
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

  // select filter
  typedef itk::ImageToImageFilter< ImageType, ImageType > 
    FilterType;
  typename FilterType::Pointer filter;
  if (!strcmp(filterType, "skel")) {
    // select skeletonize filter
    typedef itk::BinaryThinningImageFilter3D< ImageType, ImageType > 
      ThinningFilterType;
    filter = ThinningFilterType::New();
  } else {
    mexErrMsgTxt("Filter not implemented");
  }

  // run filter on input image
  filter->SetInput(image);
  filter->Update();
  
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

/* 
 * runFilter<InVoxelType, OutVoxelType> is for filters that produce an
 * output image with a different type than the input image
 * (e.g. Danielsson distance).
 *
 * The reason why we have to split this into two definitions is that
 * the thinning filter will not compile if it is templated with
 * different types for the input and output
 */
template <class InVoxelType, class OutVoxelType>
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
  const InVoxelType *im = (InVoxelType *)mxGetPr(imIn);
  
  // image type definitions
  static const unsigned int Dimension = 3; // volume data dimension
                                           // (3D volume)
  typedef double TScalarType; // data type for scalars
  typedef itk::Image< InVoxelType, Dimension > 
    InImageType;
  typedef itk::Image< OutVoxelType, Dimension > 
    OutImageType;
  typedef itk::ImageRegionIterator< InImageType > 
    InIteratorType;
  typedef itk::ImageRegionConstIterator< OutImageType > 
    OutConstIteratorType;

  // create ITK image to hold the segmentation mask
  typename InImageType::Pointer image = InImageType::New();
  typename InImageType::RegionType region;
  typename InImageType::IndexType start;
  typename InImageType::SizeType size; 

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
  InIteratorType iter(image, image->GetLargestPossibleRegion());
  mwIndex i = 0; // index for the input image
  for (iter.GoToBegin(), i=0; !iter.IsAtEnd(); ++iter, ++i) {
    iter.Set(im[i]);
    
    // Note: If you need to debug this loop, take into account that
    //   std::cout << im[i] << std::endl;
    // doesn't work as expected if im is e.g. a 
    // (uint8_T *) pointer. Instead, you need to do
    //   std::cout << (double)im[i] << std::endl;

  }

  // select filter
  typedef itk::ImageToImageFilter< InImageType, OutImageType > 
    FilterType;
  typename FilterType::Pointer filter;
  if (!strcmp(filterType, "dandist")) {
    
    // select Danielsson unsigned distance map
    typedef itk::DanielssonDistanceMapImageFilter< InImageType, 
      OutImageType > DanielssonDistanceFilterType;
    filter = DanielssonDistanceFilterType::New();

    // only needed if we want the Voronoi regions of the different
    // nonzero pixels
    // dynamic_cast<DanielssonDistanceFilterType *>(
    // 		 filter.GetPointer())->InputIsBinaryOn();
    
  } else {
    mexErrMsgTxt("Filter not implemented");
  }

  // run filter on input image
  filter->SetInput(image);
  filter->Update();
  
  // convert output data type to output class ID
  mxClassID outputVoxelClassId;
  if (TypeIsBool< OutVoxelType >::value) {
    outputVoxelClassId = mxLOGICAL_CLASS;
  } else if (TypeIsUint8< OutVoxelType >::value) {
    outputVoxelClassId = mxUINT8_CLASS;
  } else if (TypeIsUint16< OutVoxelType >::value) {
    outputVoxelClassId = mxUINT16_CLASS;
  } else if (TypeIsFloat< OutVoxelType >::value) {
    outputVoxelClassId = mxSINGLE_CLASS;
  } else if (TypeIsDouble< OutVoxelType >::value) {
    outputVoxelClassId = mxDOUBLE_CLASS;
  } else {
    mexErrMsgTxt("Assertion fail: Unrecognised output voxel type");
  }

  // create output matrix for Matlab's result
  imOut = (mxArray *)mxCreateNumericArray(ndim, dims,
					    outputVoxelClassId,
					    mxREAL);
  if (imOut == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output matrix");
  }
  OutVoxelType *imOutp =  (OutVoxelType *)mxGetPr(imOut);

  // populate output image
  OutConstIteratorType citer(filter->GetOutput(), 
			  filter->GetOutput()->GetLargestPossibleRegion());
  for (citer.GoToBegin(), i=0; i < Rin*Cin*Sin; ++citer, ++i) {
    imOutp[i] = (OutVoxelType)citer.Get();
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

  // get number of dimensions in the input segmentation
  if (!mxIsNumeric(prhs[1]) && !mxIsLogical(prhs[1])) {
    mexErrMsgTxt("IM must be a 2D or 3D numeric or boolean matrix");
  }
  mwSize ndim = mxGetNumberOfDimensions(prhs[1]);

  // get size of input arguments
  mwSize Rin = 0; // rows
  mwSize Cin = 0; // cols
  mwSize Sin = 0; // slices
  const mwSize *dims = mxGetDimensions(prhs[1]);
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

  // get type of filter
  char *filter = mxArrayToString(prhs[0]);
  if (filter == NULL) {
    mexErrMsgTxt("Invalid FILTER string");
  }

  // establish output voxel type according to the filter
  enum OutVoxelType {SAME, BOOL, UINT8, UINT16, SINGLE, DOUBLE};
  OutVoxelType outVoxelType = DOUBLE;
  if (!strcmp(filter, "skel")) {
    outVoxelType = SAME;
  } else if (!strcmp(filter, "dandist")) {
    // compute the maximum distance in voxel units that the image can have
    mwSize dmax = std::sqrt(Rin*Rin + Cin*Cin + Sin*Sin);

    // find how many bits we need to represent that distance
    mwSize nbit = (mwSize)ceil(log2(dmax));

    // select an output voxel size enough to save the maximum distance
    // value
    if (nbit <= 2) {
      outVoxelType = BOOL;
    } else if (nbit <= 8) {
      outVoxelType = UINT8;
    } else if (nbit <= 16) {
      outVoxelType = UINT16;
    } else if (nbit <= 128) {
      outVoxelType = SINGLE;
    } else {
      outVoxelType = DOUBLE;
    }
    
  } else {
    mexErrMsgTxt("Invalid FILTER string");
  }

  // run filter, templated according to the input and output image types
  mxClassID inputVoxelClassId = mxGetClassID(prhs[1]);
  switch(inputVoxelClassId)  { // swith input image type
  case mxLOGICAL_CLASS:
    if (outVoxelType == SAME) { // if block with output image type
      runFilter<bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == BOOL) {
      runFilter<bool, bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT8) {
      runFilter<bool, uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT16) {
      runFilter<bool, uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == SINGLE) {
      runFilter<bool, float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == DOUBLE) {
      runFilter<bool, double>(filter, prhs[1], plhs[0]);
    }
    break;
  case mxDOUBLE_CLASS:
    if (outVoxelType == SAME) {
      runFilter<double>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == BOOL) {
      runFilter<double, bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT8) {
      runFilter<double, uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT16) {
      runFilter<double, uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == SINGLE) {
      runFilter<double, float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == DOUBLE) {
      runFilter<double, double>(filter, prhs[1], plhs[0]);
    }
    break;
  case mxSINGLE_CLASS:
    if (outVoxelType == SAME) {
      runFilter<float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == BOOL) {
      runFilter<float, bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT8) {
      runFilter<float, uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT16) {
      runFilter<float, uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == SINGLE) {
      runFilter<float, float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == DOUBLE) {
      runFilter<float, double>(filter, prhs[1], plhs[0]);
    }
    break;
  case mxINT8_CLASS:
    if (outVoxelType == SAME) {
      runFilter<int8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == BOOL) {
      runFilter<int8_T, bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT8) {
      runFilter<int8_T, uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT16) {
      runFilter<int8_T, uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == SINGLE) {
      runFilter<int8_T, float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == DOUBLE) {
      runFilter<int8_T, double>(filter, prhs[1], plhs[0]);
    }
    break;
  case mxUINT8_CLASS:
    if (outVoxelType == SAME) {
      runFilter<uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == BOOL) {
      runFilter<uint8_T, bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT8) {
      runFilter<uint8_T, uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT16) {
      runFilter<uint8_T, uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == SINGLE) {
      runFilter<uint8_T, float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == DOUBLE) {
      runFilter<uint8_T, double>(filter, prhs[1], plhs[0]);
    }
    break;
  case mxINT16_CLASS:
    if (outVoxelType == SAME) {
      runFilter<int16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == BOOL) {
      runFilter<int16_T, bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT8) {
      runFilter<int16_T, uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT16) {
      runFilter<int16_T, uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == SINGLE) {
      runFilter<int16_T, float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == DOUBLE) {
      runFilter<int16_T, double>(filter, prhs[1], plhs[0]);
    }
    break;
  case mxUINT16_CLASS:
    if (outVoxelType == SAME) {
      runFilter<uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == BOOL) {
      runFilter<uint16_T, bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT8) {
      runFilter<uint16_T, uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT16) {
      runFilter<uint16_T, uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == SINGLE) {
      runFilter<uint16_T, float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == DOUBLE) {
      runFilter<uint16_T, double>(filter, prhs[1], plhs[0]);
    }
    break;
  case mxINT32_CLASS:
    if (outVoxelType == SAME) {
      runFilter<int32_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == BOOL) {
      runFilter<int32_T, bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT8) {
      runFilter<int32_T, uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT16) {
      runFilter<int32_T, uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == SINGLE) {
      runFilter<int32_T, float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == DOUBLE) {
      runFilter<int32_T, double>(filter, prhs[1], plhs[0]);
    }
    break;
  // case mxUINT32_CLASS:
  //   break;
  case mxINT64_CLASS:
    if (outVoxelType == SAME) {
      runFilter<int64_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == BOOL) {
      runFilter<int64_T, bool>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT8) {
      runFilter<int64_T, uint8_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == UINT16) {
      runFilter<int64_T, uint16_T>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == SINGLE) {
      runFilter<int64_T, float>(filter, prhs[1], plhs[0]);
    } else if (outVoxelType == DOUBLE) {
      runFilter<int64_T, double>(filter, prhs[1], plhs[0]);
    }
    break;
  // case mxUINT64_CLASS:
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
