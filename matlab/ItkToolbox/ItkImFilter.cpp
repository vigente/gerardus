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

#ifndef ITKIMFILTER_CPP
#define ITKIMFILTER_CPP

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <cmath>
#include <matrix.h>
#include <vector>

/* ITK headers */
#include "itkImage.h"
#include "itkBinaryThinningImageFilter3D.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

/* Gerardus headers */
#include "NrrdImage.hpp"

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
 * Global variables and declarations
 */
static const unsigned int Dimension = 3; // image dimension (we assume
					 // a 3D volume also for 2D images)

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

/*
 * FilterOutputFactory<FilterType>: Create outputs that are specific
 * to each filter.
 */

template <class FilterType>
class FilterOutputFactory {
public:
  FilterOutputFactory(typename FilterType::Pointer filter, 
		      mxArray** &argOut) {;}
};

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
  FilterFactory(char *filterType, NrrdImage &nrrd, mxArray** &argOut);
};

// avoid compiling combinations of input/output data types that are
// not available for some filters.  Read the header help of
// FilterFactoryExclusions.hpp for more info
#import "FilterFactoryExclusions.hpp"

// Constructor: where the actual filtering code lives
template <class InVoxelType, class OutVoxelType, class FilterType>
FilterFactory<InVoxelType, OutVoxelType, FilterType>::FilterFactory
(char *filterName, NrrdImage &nrrd, mxArray** &argOut) {
  
  // if the input image is empty, create empty segmentation mask for
  // output. We don't need to do any processing
  if (nrrd.getR() == 0 || nrrd.getC() == 0) {
    argOut[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
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
  argOut[0] = (mxArray *)mxCreateNumericArray(nrrd.getNdim(), nrrd.getDims(),
					  outputVoxelClassId,
					  mxREAL);
  if (argOut[0] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output matrix");
  }
  OutVoxelType *imOutp =  (OutVoxelType *)mxGetPr(argOut[0]);
  
  // populate output image
  OutConstIteratorType citer(filter->GetOutput(), 
			     filter->GetOutput()->GetLargestPossibleRegion());
  for (citer.GoToBegin(), i=0; i < nrrd.numEl(); ++citer, ++i) {
    imOutp[i] = (OutVoxelType)citer.Get();
  }
  
  // create outputs that are specific to each filter
  // FilterOutputFactory<FilterType> filterOutput(filter, argOut);
  FilterOutputFactory<FilterType> filterOutput(filter, argOut);

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
				 mxArray** &argOut) {

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
      > filterFactory(filter, nrrd, argOut);
  }  else if (!strcmp(filter, "dandist")) {
    FilterFactory<InVoxelType, OutVoxelType, 
      itk::DanielssonDistanceMapImageFilter< InImageType, OutImageType >
      > filterFactory(filter, nrrd, argOut);
  }  else if (!strcmp(filter, "maudist")) {
    FilterFactory<InVoxelType, OutVoxelType, 
      itk::SignedMaurerDistanceMapImageFilter< InImageType, OutImageType >
      > filterFactory(filter, nrrd, argOut);
  } else {
    mexErrMsgTxt("Filter type not implemented");
  }
  
}

// runTimeOutputTypeToTemplate<InVoxelType>()
template <class InVoxelType>
void runTimeOutputTypeToTemplate(char *filter,
				 NrrdImage nrrd,
				 mxArray** &argOut) {

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
      InVoxelType>(filter, nrrd, argOut);
    break;
  case BOOL:
    runTimeFilterTypeToTemplate<InVoxelType, 
      bool>(filter, nrrd, argOut);
    break;
  case UINT8:
    runTimeFilterTypeToTemplate<InVoxelType, 
      uint8_T>(filter, nrrd, argOut);
    break;
  case UINT16:
    runTimeFilterTypeToTemplate<InVoxelType, 
      uint16_T>(filter, nrrd, argOut);
    break;
  case SINGLE:
    runTimeFilterTypeToTemplate<InVoxelType, 
      float>(filter, nrrd, argOut);
    break;
  case DOUBLE:
    runTimeFilterTypeToTemplate<InVoxelType, 
      double>(filter, nrrd, argOut);
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
				mxArray** &argOut) {
  
  switch(inputVoxelClassId)  { // swith input image type
  case mxLOGICAL_CLASS:
    runTimeOutputTypeToTemplate<bool>(filter, nrrd, argOut);
    break;
  case mxDOUBLE_CLASS:
    runTimeOutputTypeToTemplate<double>(filter, nrrd, argOut);
    break;
  case mxSINGLE_CLASS:
    runTimeOutputTypeToTemplate<float>(filter, nrrd, argOut);
    break;
  case mxINT8_CLASS:
    runTimeOutputTypeToTemplate<int8_T>(filter, nrrd, argOut);
    break;
  case mxUINT8_CLASS:
    runTimeOutputTypeToTemplate<uint8_T>(filter, nrrd, argOut);
    break;
  case mxINT16_CLASS:
    runTimeOutputTypeToTemplate<int16_T>(filter, nrrd, argOut);
    break;
  case mxUINT16_CLASS:
    runTimeOutputTypeToTemplate<uint16_T>(filter, nrrd, argOut);
    break;
  case mxINT32_CLASS:
    runTimeOutputTypeToTemplate<int32_T>(filter, nrrd, argOut);
    break;
  // case mxUINT32_CLASS:
  //   break;
  case mxINT64_CLASS:
    runTimeOutputTypeToTemplate<int64_T>(filter, nrrd, argOut);
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

  // run filter (this function starts a cascade of functions designed
  // to translate the run-time type variables like inputVoxelClassId
  // to templates, so that we don't need to nest lots of "switch" or
  // "if" statements)
  runTimeInputTypeToTemplate(inputVoxelClassId, 
			     filter,
			     nrrd,
			     plhs);

  // exit successfully
  return;

}

#endif /* ITKIMFILTER_CPP */
