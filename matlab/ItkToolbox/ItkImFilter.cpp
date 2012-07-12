/* ItkImFilter.cpp
 *
 * ITK_IMFILTER: Run ITK filter on a 2D, 3D, 4D or 5D image
 *
 * This MEX function is a multiple-purpose wrapper to be able to run
 * all ITK filters that inherit from itk::ImageToImageFilter on a
 * Matlab 2D image or 3D, 4D or 5D image volume.
 *
 * B = itk_imfilter(TYPE, A, [FILTER PARAMETERS])
 *
 *   TYPE is a string with the filter we want to run. See below for a whole
 *   list of options.
 *
 *   A is a 2D matrix or 3D, 4D or 5D volume with the image or
 *   segmentation. Currently, A can be of any of the following Matlab
 *   classes:
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
 *   A can also be a SCI MAT struct, A = scimat, with the following fields:
 *
 *     scimat.data: 2D or 3D array with the image or segmentation, as above
 *     scimat.axis: 3x1 struct array with fields:
 *       scimat.axis.size:    number of voxels in the image
 *       scimat.axis.spacing: voxel size, image resolution
 *       scimat.axis.min:     real world coordinates of "left" edge of first voxel
 *       scimat.axis.max:     ignored
 *       scimat.axis.center:  ignored
 *       scimat.axis.label:   ignored
 *       scimat.axis.unit:    ignored
 *
 *   (An SCI MAT struct is the output of Matlab's function scinrrd_load(),
 *   also available from Gerardus.)
 *
 *   [FILTER PARAMETERS] is an optional list of parameters, specific for
 *   each filter. See below for details.
 *
 *   B has the same size as the image in A, and contains the filtered image
 *   or segmentation mask. It's type depends on the type of A and the filter
 *   used, and it's computed automatically.
 *
 *
 * Supported filters:
 * -------------------------------------------------------------------------
 *
 * B = itk_imfilter('skel', A)
 *
 *   (itk::BinaryThinningImageFilter3D) Skeletonize a binary mask.
 *
 *   B has the same size and class as A
 *
 * [B, V, W] = itk_imfilter('dandist', A).
 *
 *   (itk::DanielssonDistanceMapImageFilter) Compute unsigned distance map
 *   for a binary mask. Distance values are given in voxel coordinates.
 *
 *   B has the same size as A and type double. Each element in B
 *   contains an approximation to the Euclidean distance of that voxel
 *   to the closest foreground voxel, in index units.
 *
 *   V has the same size and type as A. V is a Voronoi partition of A,
 *   using the same indices.
 *
 *   W has size (3,R,C,S) if A has size (R,C,S), and type int64. Each
 *   3-vector W(:,i,j,k) is a vector pointing to the closest
 *   foreground voxel from A(i,j,k).
 *
 * B = itk_imfilter('maudist', A)
 *
 *   (itk::SignedMaurerDistanceMapImageFilter) Compute signed distance
 *   map for a binary mask. Distance values are given in real world
 *   coordinates, if the input image is given as an SCI MAT struct, or
 *   in voxel units, if the input image is a normal array.
 *
 *   The output type is always double.
 *
 * B = itk_imfilter('bwdilate', A, RADIUS, FOREGROUND)
 * B = itk_imfilter('bwerode', A, RADIUS, FOREGROUND)
 *
 *   (itk::BinaryDilateImageFilter). Binary dilation. The structuring
 *   element is a ball.
 *   (itk::BinaryErodeImageFilter). Binary erosion. The structuring
 *   element is a ball.
 *
 *   RADIUS is a scalar with the radius of the ball in voxel units. If a
 *   non-integer number is provided, then floor(RADIUS) is used. By default,
 *   RADIUS = 0 and no dilation is performed.
 *
 *   FOREGROUND is a scalar. Voxels with that value will be the only ones
 *   dilated. By default, FOREGROUND is the maximum value allowed for the
 *   type, e.g. FOREGROUND=255 if the image is uint8. This is the default in
 *   ITK, so we respect it even if we would rather have 1 as the default.
 *
 * B = itk_imfilter('advess', A, SIGMAMIN, SIGMAMAX, NUMSIGMASTEPS, ISSIGMASTEPLOG,
 *                  NUMITERATIONS, WSTRENGTH, SENSITIVITY, TIMESTEP, EPSILON)
 *
 *   (itk::AnisotropicDiffusionVesselEnhancementImageFilter)
 *   Anisotropic difussion vessel enhancement.
 *
 *   Enquobahrie A., Ibanez L., Bullitt E., Aylward S. "Vessel
 *   Enhancing Diffusion Filter", Insight Journal,
 *   2007. http://hdl.handle.net/1926/558.
 *
 *   B has the same size and class as A.
 *
 *   Note: A should have a signed type (e.g. int16, single). Images
 *   with unsigned types (e.g. uint16) will cause intermediate results
 *   that should be negative to be truncated to 0, and the result will
 *   not be meaningful. The best compromise between accuracy and
 *   saving memory seems to be type single (= float).
 *
 *   Note: While it is possible to run the filter on a SCI MAT struct,
 *   results seem better if run directly on the image. The filter
 *   doesn't seem to be spacing invariant.
 *
 *   SIGMAMIN, SIGMAMAX are scalars with the limits of the multiscale
 *   scheme, in the same units as the image. They should be set to
 *   roughly the diameters of the smallest and largest vessels in the
 *   image. By default, SIGMAMIN=0.2, SIGMAMAX=2.0.
 *
 *   NUMSIGMASTEPS is a scalar with the number of scales for the
 *   analysis. The scales change exponentially, not linearly. Casual
 *   testing suggests that the final result does not depend heavily on
 *   this parameter. By default, NUMSIGMASTEPS=10.
 *
 *   ISSIGMASTEPLOG is a boolean that determines whether the
 *   intermediate scales between SIGMAMIN to SIGMAMAX are distributed
 *   logarithmically (true) or linearly (false). The latter seems to
 *   work better for small ranges.
 *
 *   NUMITERATIONS is a scalar with the number of times the multiscale
 *   anisotropic difussion method is run. In practice, a higher number
 *   of iterations means more blurring along the vessels, which is
 *   usually desirable. The result will depend heavily on the number
 *   of iterations chosen. By default, NUMITERATIONS=10.
 *
 *   WSTRENGTH is a scalar that indicates the strength of anisotropic
 *   diffusion. Casual testing suggests that the result doesn't depend
 *   much on this value. By default, WSTRENGTH=25.0.
 *
 *   SENSITIVITY is a scalar that indicates the sensitivity to the
 *   vesselness response. Casual testing suggests that the result
 *   doesn't depend much on this value. By default, SENSITIVITY=5.0.
 *
 *   TIMESTEP is a scalar with the time step size in the difussion
 *   process. It needs to be small enough to avoid divergence, but
 *   otherwise casual testing suggests that the result doesn't depend
 *   much on this value. For 3D images, TIMESTEP < 0.0625. By default,
 *   TIMESTEP=0.001.
 *
 *   EPSILON is a scalar. It's a small number to ensure the positive
 *   definiteness of the diffusion tensor. By default, EPSILON=0.01.
 *
 * B = itk_imfilter('hesves', A, SIGMAMIN, SIGMAMAX, NUMSIGMASTEPS, ISSIGMASTEPLOG)
 *
 *   (itk::MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter)
 *   Vesselness measure from a multiscale scheme based on
 *   eigenanalysis of the Hessian.
 *
 *   Enquobahrie A., Ibanez L., Bullitt E., Aylward S. "Vessel
 *   Enhancing Diffusion Filter", Insight Journal,
 *   2007. http://hdl.handle.net/1926/558.
 *
 *   B has the same size as A, but is always of type double.
 *
 *   Input arguments are the same as the four first input arguments of
 *   filter "advess" above.
 *
 * B = itk_imfilter('median', A, RADIUS)
 *
 *   (itk::MedianImageFilter)
 *   Median of a rectangular neighbourhood.
 *
 *   B has the same size and class as A.
 *
 *   RADIUS is a vector of scalars with the half-size of the filter's
 *   box in each dimension. E.g. RADIUS=[2, 3, 4] means that the
 *   median is computed in a rectangular neighbourhood of [5, 7, 9]
 *   voxels.
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 1.0.0
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
#include "itkMedianImageFilter.h"
#include "itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.h"
#include "itkAnisotropicDiffusionVesselEnhancementImageFilter.h"
#include "itkBinaryThinningImageFilter3D.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"

/* Gerardus headers */
#include "MatlabImageHeader.h"
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

// list of supported filters. It has to be an enum so that we can pass
// it as a template constant parameter
enum SupportedFilter {
  nMedianImageFilter,
  nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter,
  nAnisotropicDiffusionVesselEnhancementImageFilter,
  nBinaryThinningImageFilter3D,
  nDanielssonDistanceMapImageFilter,
  nSignedMaurerDistanceMapImageFilter,
  nBinaryDilateImageFilter,
  nBinaryErodeImageFilter
};

// FilterWrapper():
//
// This block contains one FilterWrapper partial specialisation per
// filter. In this class, the filter acquires the inputs from Matlab,
// parameters are set, and the outputs are grafted onto Matlab.
//
// The reason to use an encapsulating class like this FilterWrapper is
// because some filters do not accept certains dimensions or input
// types. Thus, we can use partial template specialisation to give a
// runtime error in those cases and avoid instantiating the ITK
// filter, which would give a compilation error
template <class TPixelIn, class TPixelOut, unsigned int VImageDimension,
	  unsigned int FilterEnum>
class FilterWrapper {
public:
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    mexErrMsgTxt("Unsupported filter type");
  }
};

// MedianImageFilter
template <class TPixelIn, class TPixelOut, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, TPixelOut, VImageDimension,
		    nMedianImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // instantiate the filter
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::MedianImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 5);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->
		     GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));
    
    // set half size of the filter's box
    typedef typename itk::BoxImageFilter<
      itk::Image<TPixelIn, VImageDimension>,
      itk::Image<TPixelOut, VImageDimension> > BoxFilterType;
    typename BoxFilterType::RadiusType radius;
    radius.Fill(0);
    filter->SetRadius(matlabImport->
		      GetVectorArgument<typename BoxFilterType::RadiusType, 
					typename BoxFilterType::RadiusValueType>
		      (2, "RADIUS", radius));
  
    // connect ITK filter outputs to Matlab outputs
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->GraftItkImageOntoMatlab<TPixelOut, VImageDimension>
	(filter->GetOutputs()[0], im.size, 0, "0");
    }

    // run filter
    filter->Update();

  }
};

// MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter
template <class TPixelIn, class TPixelOut, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, TPixelOut, VImageDimension,
		    nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // instantiate the filter
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef typename itk::Image<TPixelOut, VImageDimension> OutImageType;
    typedef itk::MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter
      <InImageType, OutImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 6);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->
		     GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));
    
    // filter parameters
    filter->SetSigmaMin(matlabImport->template
			GetScalarArgument<double>(2, "SIGMAMIN", 0.2));
    filter->SetSigmaMax(matlabImport->template
			GetScalarArgument<double>(3, "SIGMAMAX", 2.0));
    filter->SetNumberOfSigmaSteps(matlabImport->template
			GetScalarArgument<int>(4, "NUMSIGMASTEPS", 10));
    filter->SetIsSigmaStepLog(matlabImport->template
			GetScalarArgument<bool>(5, "ISSIGMASTEPLOG", true));

    // connect ITK filter outputs to Matlab outputs
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->GraftItkImageOntoMatlab<TPixelOut, VImageDimension>
	(filter->GetOutputs()[0], im.size, 0, "0");
    }

    // run filter
    filter->Update();

  }
};

template <class TPixelIn, class TPixelOut>
class FilterWrapper<TPixelIn, TPixelOut, 2,
		    nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter only accepts 3D input images");
  }
};

template <class TPixelIn, class TPixelOut>
class FilterWrapper<TPixelIn, TPixelOut, 4,
		    nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter only accepts 3D input images");
  }
};

// AnisotropicDiffusionVesselEnhancementImageFilter
template <class TPixelIn, class TPixelOut, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, TPixelOut, VImageDimension,
		    nAnisotropicDiffusionVesselEnhancementImageFilter> {
public:

  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // instantiate the filter
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::AnisotropicDiffusionVesselEnhancementImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 11);
    matlabExport->CheckNumberOfArguments(0, 1);

    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));

    filter->SetSigmaMin(matlabImport->
		       GetScalarArgument<double>(2,  "SIGMAMIN", 0.2));
    filter->SetSigmaMax(matlabImport->
		       GetScalarArgument<double>(3,  "SIGMAMAX", 2.0));
    filter->SetNumberOfSigmaSteps(matlabImport->
		       GetScalarArgument<int>   (4,  "NUMSIGMASTEPS", 10));
    filter->SetIsSigmaStepLog(matlabImport->
		       GetScalarArgument<bool>  (5,  "ISSIGMASTEPLOG", true));
    filter->SetNumberOfIterations(matlabImport->
		       GetScalarArgument<int>   (6,  "NUMITERATIONS", 1));
    filter->SetWStrength(matlabImport->
		       GetScalarArgument<double>(7,  "WSTRENGTH", 25.0));
    filter->SetSensitivity(matlabImport->
		       GetScalarArgument<double>(8,  "SENSITIVITY", 5.0));
    filter->SetTimeStep(matlabImport->
		       GetScalarArgument<double>(9,  "TIMESTEP", 1e-3));
    filter->SetEpsilon(matlabImport->
		       GetScalarArgument<double>(10, "EPSILON", 1e-2));
    
    // connect ITK filter outputs to Matlab outputs
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->GraftItkImageOntoMatlab<TPixelOut, VImageDimension>
	(filter->GetOutputs()[0], im.size, 0, "0");
    }

    // run filter
    filter->Update();

  }
};

template <class TPixelIn, class TPixelOut>
class FilterWrapper<TPixelIn, TPixelOut, 2,
		    nAnisotropicDiffusionVesselEnhancementImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("AnisotropicDiffusionVesselEnhancementImageFilter only accepts 3D input images");
  }
};

template <class TPixelIn, class TPixelOut>
class FilterWrapper<TPixelIn, TPixelOut, 4,
	      nAnisotropicDiffusionVesselEnhancementImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("AnisotropicDiffusionVesselEnhancementImageFilter only accepts 3D input images");
  }
};

// BinaryThinningImageFilter3D
template <class TPixelIn, class TPixelOut, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, TPixelOut, VImageDimension,
		    nBinaryThinningImageFilter3D> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // instantiate the filter
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::BinaryThinningImageFilter3D<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 2);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->
		     GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));

    // connect ITK filter outputs to Matlab outputs
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->GraftItkImageOntoMatlab<TPixelOut, VImageDimension>
	(filter->GetOutputs()[0], im.size, 0, "0");
    }

    // run filter
    filter->Update();

  }
};

template <class TPixelIn, class TPixelOut>
class FilterWrapper<TPixelIn, TPixelOut, 2,
		    nBinaryThinningImageFilter3D> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("BinaryThinningImageFilter3D only accepts 3D input images");
  }
};

template <class TPixelIn, class TPixelOut>
class FilterWrapper<TPixelIn, TPixelOut, 4,
		    nBinaryThinningImageFilter3D> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("BinaryThinningImageFilter3D only accepts 3D input images");
  }
};

// DanielssonDistanceMapImageFilter
template <class TPixelIn, class TPixelOut, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, TPixelOut, VImageDimension,
		    nDanielssonDistanceMapImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // instantiate the filter
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef typename itk::Image<TPixelOut, VImageDimension> OutImageType;
    typedef itk::DanielssonDistanceMapImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 2);
    matlabExport->CheckNumberOfArguments(0, 3);
    
    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->
		     GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));

    // connect ITK filter outputs to Matlab outputs

    // distance map
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->GraftItkImageOntoMatlab<TPixelOut, VImageDimension>
	(filter->GetOutputs()[0], im.size, 0, "0");
    }
    // Voronoi map
    if (matlabExport->GetNumberOfArguments() >= 2) {
      matlabExport->GraftItkImageOntoMatlab<TPixelIn, VImageDimension>
	(filter->GetOutputs()[1], im.size, 1, "1");
    }
    // vectors pointing to closest foreground voxel
    if (matlabExport->GetNumberOfArguments() >= 3) {
      matlabExport->GraftItkImageOntoMatlab<typename InImageType::OffsetType::OffsetValueType,
					    VImageDimension,
					    typename InImageType::OffsetType::OffsetType>
	(filter->GetOutputs()[2], im.size, 2, "2");
    }

    // run filter
    filter->Update();

  }
};

// SignedMaurerDistanceMapImageFilter
template <class TPixelIn, class TPixelOut, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, TPixelOut, VImageDimension,
		    nSignedMaurerDistanceMapImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // instantiate the filter
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef typename itk::Image<TPixelOut, VImageDimension> OutImageType;
    typedef itk::SignedMaurerDistanceMapImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 2);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // compute distances using real world coordinates, instead of voxel
    // indices
    filter->SetUseImageSpacing(true);
    
    // give output as actual distances
    filter->SquaredDistanceOff();
    
    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->
		     GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));

    // connect ITK filter outputs to Matlab outputs

    // distance map
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->GraftItkImageOntoMatlab<TPixelOut, VImageDimension>
	(filter->GetOutputs()[0], im.size, 0, "0");
    }

    // run filter
    filter->Update();

  }
};

template <class TPixelOut, unsigned int VImageDimension>
class FilterWrapper<mxLogical, TPixelOut, VImageDimension,
		    nSignedMaurerDistanceMapImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("SignedMaurerDistanceMapImageFilter does not accept input image with type boolean");
  }
};

// BinaryDilateImageFilter
template <class TPixelIn, class TPixelOut, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, TPixelOut, VImageDimension,
		    nBinaryDilateImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // instantiate the filter
    typedef itk::BinaryBallStructuringElement<TPixelIn, VImageDimension>
      StructuringElementType;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::BinaryDilateImageFilter<InImageType, OutImageType, StructuringElementType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 4);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->
		     GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));
    
    // instantiate structuring element
    // (comp) radius of the ball in voxels
    StructuringElementType structuringElement;
    structuringElement.SetRadius(matlabImport->
				 GetScalarArgument<unsigned long>(2, "RADIUS", 0));
    structuringElement.CreateStructuringElement();
    filter->SetKernel(structuringElement);

    // pass other parameters to filter
    // (opt) voxels with this value will be dilated. Default, maximum
    // value of the pixel type (this is the ITK default, so we
    // reproduce it here, even if it "1" would be more convenient)
    filter->SetForegroundValue(matlabImport->template
			       GetScalarArgument<TPixelIn>(3, "FOREGROUND", std::numeric_limits<TPixelIn>::max()));

    // connect ITK filter outputs to Matlab outputs
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->GraftItkImageOntoMatlab<TPixelOut, VImageDimension>
	(filter->GetOutputs()[0], im.size, 0, "0");
    }

    // run filter
    filter->Update();

  }
};

// BinaryErodeImageFilter
template <class TPixelIn, class TPixelOut, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, TPixelOut, VImageDimension,
		    nBinaryErodeImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // instantiate the filter
    typedef itk::BinaryBallStructuringElement<TPixelIn, VImageDimension>
      StructuringElementType;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::BinaryErodeImageFilter<InImageType, OutImageType, StructuringElementType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 4);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->
		     GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));
    
    // instantiate structuring element
    // (comp) radius of the ball in voxels
    StructuringElementType structuringElement;
    structuringElement.SetRadius(matlabImport->
				 GetScalarArgument<unsigned long>(2, "RADIUS", 0));
    structuringElement.CreateStructuringElement();
    filter->SetKernel(structuringElement);

    // pass other parameters to filter
    // (opt) voxels with this value will be dilated. Default, maximum
    // value of the pixel type (this is the ITK default, so we
    // reproduce it here, even if it "1" would be more convenient)
    filter->SetForegroundValue(matlabImport->template
			       GetScalarArgument<TPixelIn>(3, "FOREGROUND", std::numeric_limits<TPixelIn>::max()));

    // connect ITK filter outputs to Matlab outputs
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->GraftItkImageOntoMatlab<TPixelOut, VImageDimension>
	(filter->GetOutputs()[0], im.size, 0, "0");
    }

    // run filter
    filter->Update();

  }
};

/*
 * Argument Parsers
 *
 * These functions are used to be able to map between the input/output
 * data types that are only know at run-time, and the input/output
 * data templates that ITK requires and must be know at compilation
 * time.
 */

// parseOutputImageTypeToTemplate()
template <class TPixelIn, unsigned int VImageDimension>
void parseOutputImageTypeToTemplate(MatlabImportFilter::Pointer matlabImport,
				    MatlabExportFilter::Pointer matlabExport,
				    MatlabImageHeader &im) {

  // input image type
  typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;

  // name of the filter
  std::string filterName = matlabImport->GetStringArgument(0, "0", "");

  // select the output type corresponding to each filter
  if (filterName == "median" 
	     || filterName == "MedianImageFilter") {

    FilterWrapper<TPixelIn, TPixelIn, VImageDimension, 
		  nMedianImageFilter> 
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "advess" 
      || filterName == "AnisotropicDiffusionVesselEnhancementImageFilter") {

    FilterWrapper<TPixelIn, TPixelIn, VImageDimension, 
		  nAnisotropicDiffusionVesselEnhancementImageFilter> 
      wrapper(matlabImport, matlabExport, im);

  } else if (filterName == "bwdilate" 
	     || filterName == "BinaryDilateImageFilter") {

    FilterWrapper<TPixelIn, TPixelIn, VImageDimension, 
		  nBinaryDilateImageFilter> 
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "bwerode" 
	     || filterName == "BinaryErodeImageFilter") {

    FilterWrapper<TPixelIn, TPixelIn, VImageDimension, 
		  nBinaryErodeImageFilter> 
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "skel" 
	     || filterName == "BinaryThinningImageFilter3D") {

    FilterWrapper<TPixelIn, TPixelIn, VImageDimension, 
		  nBinaryThinningImageFilter3D>
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "dandist" 
	     || filterName == "DanielssonDistanceMapImageFilter") {

    FilterWrapper<TPixelIn, double, VImageDimension, 
		  nDanielssonDistanceMapImageFilter>
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "hesves" 
	     || filterName == "MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter") {

    FilterWrapper<TPixelIn, double, VImageDimension, 
		  nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter>
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "maudist" 
	     || filterName == "SignedMaurerDistanceMapImageFilter") {

    FilterWrapper<TPixelIn, int, VImageDimension, 
		  nSignedMaurerDistanceMapImageFilter>
      filterWrapper(matlabImport, matlabExport, im);

  } else {
    mexErrMsgTxt("Invalid filter type");    
  }  

}
  
// parseInputImageTypeToTemplate()
template <unsigned int VImageDimension>
void parseInputImageTypeToTemplate(MatlabImportFilter::Pointer matlabImport,
				   MatlabExportFilter::Pointer matlabExport,
				   MatlabImageHeader &im) {
  
  // input image type
  switch(im.type)  {
  case mxLOGICAL_CLASS:
    parseOutputImageTypeToTemplate<mxLogical, VImageDimension>(matlabImport, matlabExport, im);
    break;
  case mxDOUBLE_CLASS:
    parseOutputImageTypeToTemplate<double, VImageDimension>(matlabImport, matlabExport, im);
    break;
  case mxSINGLE_CLASS:
    parseOutputImageTypeToTemplate<float, VImageDimension>(matlabImport, matlabExport, im);
    break;
  case mxINT8_CLASS:
    parseOutputImageTypeToTemplate<int8_T, VImageDimension>(matlabImport, matlabExport, im);
    break;
  case mxUINT8_CLASS:
    parseOutputImageTypeToTemplate<uint8_T, VImageDimension>(matlabImport, matlabExport, im);
    break;
  case mxINT16_CLASS:
    parseOutputImageTypeToTemplate<int16_T, VImageDimension>(matlabImport, matlabExport, im);
    break;
  case mxUINT16_CLASS:
    parseOutputImageTypeToTemplate<uint16_T, VImageDimension>(matlabImport, matlabExport, im);
    break;
  case mxINT32_CLASS:
    parseOutputImageTypeToTemplate<int32_T, VImageDimension>(matlabImport, matlabExport, im);
    break;
  // case mxUINT32_CLASS:
  //   break;
  case mxINT64_CLASS:
    parseOutputImageTypeToTemplate<int64_T, VImageDimension>(matlabImport, matlabExport, im);
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

void parseInputImageDimensionToTemplate(MatlabImportFilter::Pointer matlabImport,
					MatlabExportFilter::Pointer matlabExport) {
  
  // the 2nd input argument is the input image. It can be given as an
  // array, or a SCI MAT struct, so it's necessary to pre-process the
  // pointer to do checks and extract the meta information
  MatlabImageHeader im(matlabImport->GetArg(1), "1");

  switch (im.GetNumberOfDimensions()) {
  case 2:
    parseInputImageTypeToTemplate<2>(matlabImport, matlabExport, im);
    break;
  case 3:
    parseInputImageTypeToTemplate<3>(matlabImport, matlabExport, im);
    break;
  case 4:
    parseInputImageTypeToTemplate<4>(matlabImport, matlabExport, im);
    break;
  // case 5:
  //   parseInputImageTypeToTemplate<5>(matlabImport, matlabExport, im);
  //   break;
  default:
    mexErrMsgTxt("Input image can only have 2 to 5 dimensions");
    break;
  }

}


/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->SetMatlabArgumentsPointer(nrhs, prhs);

  // check that we have at least a filter name and input image
  matlabImport->CheckNumberOfArguments(2, UINT_MAX);

  // interface to deal with output arguments from Matlab
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->SetMatlabArgumentsPointer(nlhs, plhs);

  // run filter (this function starts a cascade of functions designed
  // to translate the run-time type variables like inputVoxelClassId
  // to templates, so that we don't need to nest lots of "switch" or
  // "if" statements)
  parseInputImageDimensionToTemplate(matlabImport, matlabExport);

  // exit successfully
  return;

}

#endif /* ITKIMFILTER_CPP */
