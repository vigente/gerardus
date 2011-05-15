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
 *     'maudist': (SignedMaurerDistanceMapImageFilter) Compute signed
 *                distance map for a binary mask. Distance values are
 *                given in real world coordinates, if the input image
 *                is given as an NRRD struct, or in voxel units, if
 *                the input image is a normal array. The output type
 *                is always double.
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
 *   A can also be a SCI NRRD struct, A = nrrd, with the following fields:
 *
 *     nrrd.data: 2D or 3D array with the image or segmentation, as above
 *     nrrd.axis: 3x1 struct array with fields:
 *       nnrd.axis.size:    number of voxels in the image
 *       nnrd.axis.spacing: voxel size, image resolution
 *       nnrd.axis.min:     real world coordinates of image origin
 *       nnrd.axis.max:     ignored
 *       nnrd.axis.center:  ignored
 *       nnrd.axis.label:   ignored
 *       nnrd.axis.unit:    ignored
 *
 *   An SCI NRRD struct is the output of Matlab's function
 *   scinrrd_load(), also available from Gerardus.
 *
 *   B has the same size and class as the image in A, and contains the
 *   filtered image or segmentation mask.
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
  * Version: 0.3.4
  * $Rev: 376 $
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
#include "itkSignedMaurerDistanceMapImageFilter.h"

#ifndef ITKIMFILTER_CPP
#define ITKIMFILTER_CPP

/*
 * Global variables and declarations
 */
static const unsigned int Dimension = 3; // image dimension (we assume
					 // a 3D volume also for 2D images)


/*
 * Block of functions to allow testing of template types
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
 * NrrdImage: class to parse an image that follows the SCI NRRD format
 * obtained from saving an image in Seg3D to a Matlab file
 */
class NrrdImage {
private:
  
  mxArray *data; // pointer to the image data in Matlab format
  mwSize ndim; // number of elements in the dimensions array dims
  mwSize *dims; // dimensions array
  // size[0] -> rows
  // size[1] -> columns
  // size[2] -> slices
  std::vector<mwSize> size; // number of voxels in each dimension
  std::vector<double> spacing; // voxel size in each dimension
  std::vector<double> min; // real world coordinates of the image origin

public:
  NrrdImage(const mxArray * nrrd);
  mxArray * getData() {return data;}
  std::vector<mwSize> getSize() {return size;}
  std::vector<double> getSpacing() {return spacing;}
  std::vector<double> getMin() {return min;}
  mwSize getR() {return size[0];}
  mwSize getC() {return size[1];}
  mwSize getS() {return size[2];}
  mwSize getDr() {return spacing[0];}
  mwSize getDc() {return spacing[1];}
  mwSize getDs() {return spacing[2];}
  mwSize getDx() {return spacing[1];}
  mwSize getDy() {return spacing[0];}
  mwSize getDz() {return spacing[2];}
  mwSize getMinR() {return min[0];}
  mwSize getMinC() {return min[1];}
  mwSize getMinS() {return min[2];}
  mwSize getMinX() {return min[1];}
  mwSize getMinY() {return min[0];}
  mwSize getMinZ() {return min[2];}
  mwSize getNdim() {return ndim;}
  mwSize *getDims() {return dims;}
  mwSize maxVoxDistance();
  mwSize numEl();
};

// The constructor works as a parser. It figures out whether the input
// data is a normal array with the image, or an NRRD struct. In the
// latter case, it reads the struct fields to extract the spacing and
// offset of the image
NrrdImage::NrrdImage(const mxArray * nrrd) {

  // initialize memory for member variables
  spacing.resize(Dimension);
  min.resize(Dimension);
  size.resize(Dimension);

  // check whether the input image is an SCI NRRD struct instead of
  // only the array. The SCI NRRD struct has information about the
  // scaling, offset, etc. of the image
  if (mxIsStruct(nrrd)) { // input is a struct
    // read struct fields...

    // ... data field (where the image is contained)
    data = mxGetField(nrrd, 0, "data");
    if (data == NULL) {
      mexErrMsgTxt("NRRD format error: Missing or invalid data field");
    }

    // ... axis field
    mxArray *axis = mxGetField(nrrd, 0, "axis");
    if (axis == NULL) {
      mexErrMsgTxt("NRRD format error: Missing or invalid axis field");
    }
    if (!mxIsStruct(axis)) {
      mexErrMsgTxt("NRRD format error: axis field is not a struct");
    }
    if (mxGetM(axis) != 3) {
      mexErrMsgTxt("NRRD format error: axis field must be a 3x1 struct array");
    }

    // ... spacing field (image resolution)
    // ... min field (image origin)
    mxArray *spacingMx, *minMx;
    double *spacingMxp, *minMxp;
    for (mwIndex i = 0; i < 3; ++i) {
      spacingMx = mxGetField(axis, i, "spacing");
      if (spacingMx == NULL) {
	mexErrMsgTxt("NRRD format error: Missing or invalid axis.spacing field");
      }
      minMx = mxGetField(axis, i, "min");
      if (minMx == NULL) {
	mexErrMsgTxt("NRRD format error: Missing or invalid axis.min field");
      }
      spacingMxp = (double *)mxGetData(spacingMx);
      minMxp = (double *)mxGetData(minMx);
      if (mxIsFinite(spacingMxp[0]) && !mxIsNaN(spacingMxp[0])) {
	spacing[i] = spacingMxp[0];
      } else {
	spacing[i] = 1.0;
      }
      if (mxIsFinite(minMxp[0]) && !mxIsNaN(minMxp[0])) {
	min[i] = minMxp[0];
      } else {
	min[i] = 0.0;
      }
    }

  } else { // input is not a struct, but just an image array
    data = const_cast<mxArray *>(nrrd);
    spacing[0] = 1.0;
    spacing[1] = 1.0;
    spacing[2] = 1.0;
    min[0] = 0.0;
    min[1] = 0.0;
    min[2] = 0.0;
  }


  // get number of dimensions in the input image
  if (!mxIsNumeric(data) 
      && !mxIsLogical(data)) {
    mexErrMsgTxt("IM must be a 2D or 3D numeric or boolean matrix");
  }
  ndim = mxGetNumberOfDimensions(data);
  
  // get size of input arguments
  dims = const_cast<mwSize *>(mxGetDimensions(data));
  if (dims == NULL) {
    mexErrMsgTxt("Invalid input image dimensions array");
  }
  size[0] = dims[0];
  if (ndim > 3) {
    mexErrMsgTxt("Input segmentation mask must be 2D or 3D");
  } else if (ndim == 2) {
    size[1] = dims[1];
    size[2] = 1;
  } else if (ndim == 3) {
    size[1] = dims[1];
    size[2] = dims[2];
  } else {
    mexErrMsgTxt("Assertion fail: number of dimensions is " + ndim);
  }

}

// compute the maximum distance between any two voxels in this image
// (in voxel units). This is the length of the largest diagonal in the
// cube
mwSize NrrdImage::maxVoxDistance() {
  return std::sqrt((size[0]-1)*(size[0]-1)
		   + (size[1]-1)*(size[1]-1)
		   + (size[2]-1)*(size[2]-1));
}

// compute the number of voxels in the image volume
mwSize NrrdImage::numEl() {
  return size[0]*size[1]*size[2];
}

/* 
 * FilterFactory<InVoxelType, OutVoxelType, FilterType>: This is where
 * the code to actually run the filter on the image lives.
 *
 * Instead of having a function (e.g. runFilter), we have the code in
 * the constructor of class FilterFactory.
 *
 * The reason is that template explicit specialization is only
 * possible in classes, not in functions. We need explicit
 * specialization to prevent the compiler from compiling certain
 * input/output image data types for some filters that don't accept
 * them.
 */

template <class InVoxelType, class OutVoxelType, class FilterType>
class FilterFactory {
public:
  FilterFactory(char *filterType, NrrdImage &nrrd, mxArray* &imOut);
};

/*
 * FilterParamFactory: class to pass parameters specific to one filter
 * but not the others
 */

// default: any filter without an explicit specialization
template <class InVoxelType, class OutVoxelType, 
	  class FilterType>
class FilterParamFactory {
public:
  FilterParamFactory(typename FilterType::Pointer filter) {
    // by default, we assume that filters do not need parameters. If a
    // filter needs some specific parameters, or setting any flags, we
    // need to declare a explicit specialization of this class, and
    // put the corresponding code there
    //
    // Hence, this constructor is empty
    ;
  }
};

// SignedMaurerDistanceMapImageFilter: Specific parameters
template <class InVoxelType, class OutVoxelType>
class FilterParamFactory<InVoxelType, OutVoxelType,
			 itk::SignedMaurerDistanceMapImageFilter< 
			   itk::Image<InVoxelType, Dimension>,
			   itk::Image<OutVoxelType, Dimension> > 
			 > {
public:
  FilterParamFactory(typename itk::SignedMaurerDistanceMapImageFilter< 
		     itk::Image<InVoxelType, Dimension>,
		     itk::Image<OutVoxelType, Dimension> >::Pointer filter) {
    // we want the output in real world coordinates by default. If the
    // user wants voxel units, then provide a plain image at the
    // input, or make the spacing in the NRRD struct = [1.0, 1.0, 1.0]
    filter->SetUseImageSpacing(true);
  
    // we want actual Euclidean distances, not squared ones
    filter->SquaredDistanceOff();
  }
};

// avoid compiling combinations of input/output data types that are
// not available for some filters.  Read the header help of
// FilterFactoryExclusions.hpp for more info
#import "FilterFactoryExclusions.hpp"

// Constructor: where the actual filtering code lives
template <class InVoxelType, class OutVoxelType, class FilterType>
FilterFactory<InVoxelType, OutVoxelType, FilterType>::FilterFactory
(char *filterName, NrrdImage &nrrd, mxArray* &imOut) {
  
  // if the input image is empty, create empty segmentation mask for
  // output. We don't need to do any processing
  if (nrrd.getR() == 0 || nrrd.getC() == 0) {
    imOut = mxCreateDoubleMatrix(0, 0, mxREAL);
    return;
  }
  
  // get pointer to input segmentation mask
  const InVoxelType *im = (InVoxelType *)mxGetPr(nrrd.getData());
  
  // image type definitions
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
  typename InImageType::SpacingType spacing;
  
  // note that:
  //
  // 1) in ITK we have X,Y,Z indices, while in Matlab we have R,C,S
  //
  // 2) matrices in ITK are read by columns, while in Matlab
  // they are read by rows, 
  //
  // So both things cancel each other.
  for (mwIndex i = 0; i < Dimension; ++i) {
    start[i] = nrrd.getMin()[i];
    size[i] = nrrd.getSize()[i];
    spacing[i] = nrrd.getSpacing()[i];
  }
  region.SetIndex(start);
  region.SetSize(size);
  image->SetRegions(region);
  image->SetSpacing(spacing);
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
  
  // instantiate filter
  typename FilterType::Pointer filter = FilterType::New();
  
  // pass any parameters specific to this filter. We have to do it
  // with explicit specialization of class FilterParamFactory to
  // prevent the compiler from compiling e.g. filter->SetAlpha1(0.3)
  // for filters that don't declare that member function
  FilterParamFactory<InVoxelType, OutVoxelType, 
    FilterType> filterParam(filter);

  // run filter on input image
  filter->SetInput(image);
  filter->Update();
  
  // convert output data type to output class ID
  mxClassID outputVoxelClassId = mxUNKNOWN_CLASS;
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
  imOut = (mxArray *)mxCreateNumericArray(nrrd.getNdim(), nrrd.getDims(),
					  outputVoxelClassId,
					  mxREAL);
  if (imOut == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output matrix");
  }
  OutVoxelType *imOutp =  (OutVoxelType *)mxGetPr(imOut);
  
  // populate output image
  OutConstIteratorType citer(filter->GetOutput(), 
			     filter->GetOutput()->GetLargestPossibleRegion());
  for (citer.GoToBegin(), i=0; i < nrrd.numEl(); ++citer, ++i) {
    imOutp[i] = (OutVoxelType)citer.Get();
  }
  

}

/*
 * runTimeInputTypeToTemplate()
 * runTimeOutputTypeToTemplate< InVoxelType >()
 * runTimeFilterTypeToTemplate< InVoxelType, OutVoxelType >()
 *
 * These functions are used to be able to map between the input/output
 * data types that are only know at run-time, and the input/output
 * data templates that ITK requires and must be know at compilation
 * time.
 *
 * To avoid a nesting nightmare:
 *
 * switch FilterType {
 *   switch InputDataType {
 *     switch OutputDataType {
 *     }
 *   }
 * }
 *
 * we split the conversion in 3 steps. The first function,
 * runTimeInputTypeToTemplate() reads the input data type and the
 * filter type, and instantiates
 * runTimeOutputTypeToTemplate<InVoxelType>() for each input data
 * type.
 *
 * This way, we have effectively mapped the input type from a run-time
 * variable to a set of compilation time templates.
 *
 * Then runTimeOutputTypeToTemplate<InVoxelType>() decides on the
 * output data type depending on the filter and the input, and
 * instantiates one
 * runTimeFilterTypeToTemplate<InVoxelType, OutVoxelType>() for each
 * output data type.
 *
 * This way, we have instantiated all combinations of input/output
 * data types, without having to code them explicitly.
 *
 * Finally, runTimeFilterTypeToTemplate<InVoxelType, OutVoxelType>()
 * checks that the requested filter is implemented and maps the filter
 * variable to all filter templates.
 *
 * Note that the actual filtering is done in the constructor of
 * FilterFactory, that is instantiated from
 * runTimeFilterTypeToTemplate.
 */

// runTimeFilterTypeToTemplate<InVoxelType, OutVoxelType>()
template <class InVoxelType, class OutVoxelType>
void runTimeFilterTypeToTemplate(char *filter,
				 NrrdImage nrrd,
				 mxArray* &imOut) {

  // image type definitions
  static const unsigned int Dimension = 3; // volume data dimension
                                           // (3D volume)
  typedef double TScalarType; // data type for scalars
  typedef itk::Image< InVoxelType, Dimension > 
    InImageType;
  typedef itk::Image< OutVoxelType, Dimension > 
    OutImageType;


  // convert run-time filter string to template
  if (!strcmp(filter, "skel")) {
    FilterFactory<InVoxelType, OutVoxelType, 
      itk::BinaryThinningImageFilter3D< InImageType, OutImageType >
      > filterFactory(filter, nrrd, imOut);
  }  else if (!strcmp(filter, "dandist")) {
    FilterFactory<InVoxelType, OutVoxelType, 
      itk::DanielssonDistanceMapImageFilter< InImageType, OutImageType >
      > filterFactory(filter, nrrd, imOut);
  }  else if (!strcmp(filter, "maudist")) {
    FilterFactory<InVoxelType, OutVoxelType, 
      itk::SignedMaurerDistanceMapImageFilter< InImageType, OutImageType >
      > filterFactory(filter, nrrd, imOut);
  } else {
    mexErrMsgTxt("Filter type not implemented");
  }
  
}

// runTimeOutputTypeToTemplate<InVoxelType>()
template <class InVoxelType>
void runTimeOutputTypeToTemplate(char *filter,
				 NrrdImage nrrd,
				 mxArray* &imOut) {

  // make it easier to remember the different cases for the output
  // voxel type
  enum OutVoxelType {
    SAME, BOOL, UINT8, UINT16, SINGLE, DOUBLE
  };

  // establish output voxel type according to the filter
  OutVoxelType outVoxelType = DOUBLE;
  if (!strcmp(filter, "skel")) {

    outVoxelType = SAME;

  } else if (!strcmp(filter, "dandist")) {
    // find how many bits we need to represent the maximum distance
    // that two voxels can have between them (in voxel units)
    mwSize nbit = (mwSize)ceil(log2(nrrd.maxVoxDistance()));

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
    
  } else if (!strcmp(filter, "maudist")) {

    outVoxelType = DOUBLE;

  } else {
    mexErrMsgTxt("Filter type not implemented");
  }

  switch(outVoxelType) {
  case SAME:
    runTimeFilterTypeToTemplate<InVoxelType, 
      InVoxelType>(filter, nrrd, imOut);
    break;
  case BOOL:
    runTimeFilterTypeToTemplate<InVoxelType, 
      bool>(filter, nrrd, imOut);
    break;
  case UINT8:
    runTimeFilterTypeToTemplate<InVoxelType, 
      uint8_T>(filter, nrrd, imOut);
    break;
  case UINT16:
    runTimeFilterTypeToTemplate<InVoxelType, 
      uint16_T>(filter, nrrd, imOut);
    break;
  case SINGLE:
    runTimeFilterTypeToTemplate<InVoxelType, 
      float>(filter, nrrd, imOut);
    break;
  case DOUBLE:
    runTimeFilterTypeToTemplate<InVoxelType, 
      double>(filter, nrrd, imOut);
    break;
  default:
    mexErrMsgTxt("Invalid output type.");
    break;
  }
}

// runTimeInputTypeToTemplate()
void runTimeInputTypeToTemplate(mxClassID inputVoxelClassId, 
				char *filter,
				NrrdImage nrrd,
				mxArray* &imOut) {
  
  switch(inputVoxelClassId)  { // swith input image type
  case mxLOGICAL_CLASS:
    runTimeOutputTypeToTemplate<bool>(filter, nrrd, imOut);
    break;
  case mxDOUBLE_CLASS:
    runTimeOutputTypeToTemplate<double>(filter, nrrd, imOut);
    break;
  case mxSINGLE_CLASS:
    runTimeOutputTypeToTemplate<float>(filter, nrrd, imOut);
    break;
  case mxINT8_CLASS:
    runTimeOutputTypeToTemplate<int8_T>(filter, nrrd, imOut);
    break;
  case mxUINT8_CLASS:
    runTimeOutputTypeToTemplate<uint8_T>(filter, nrrd, imOut);
    break;
  case mxINT16_CLASS:
    runTimeOutputTypeToTemplate<int16_T>(filter, nrrd, imOut);
    break;
  case mxUINT16_CLASS:
    runTimeOutputTypeToTemplate<uint16_T>(filter, nrrd, imOut);
    break;
  case mxINT32_CLASS:
    runTimeOutputTypeToTemplate<int32_T>(filter, nrrd, imOut);
    break;
  // case mxUINT32_CLASS:
  //   break;
  case mxINT64_CLASS:
    runTimeOutputTypeToTemplate<int64_T>(filter, nrrd, imOut);
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
}

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {
  // check number of input and output arguments
  if (nrhs != 2) {
    mexErrMsgTxt("Two input arguments required");
  }
  else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  // read image and its parameters, whether it's in NRRD format, or
  // just a 2D or 3D array
  NrrdImage nrrd(prhs[1]);

  // get type of filter
  char *filter = mxArrayToString(prhs[0]);
  if (filter == NULL) {
    mexErrMsgTxt("Invalid FILTER string");
  }

  // input image type
  mxClassID inputVoxelClassId = mxGetClassID(nrrd.getData());

  // run filter (this function starts a cascade of functions design to
  // translate the run-time type variables like inputVoxelClassId to
  // templates, so that we don't need to nest lots of switch of if
  // clauses
  runTimeInputTypeToTemplate(inputVoxelClassId, 
			     filter,
			     nrrd,
			     plhs[0]);

  // exit successfully
  return;

}

#endif /* ITKIMFILTER_CPP */
