/* ItkImFilter.cpp
 *
 * ITK_IMFILTER: Run ITK filter on a 2D, 3D or 4D image
 *
 * This MEX function is a multiple-purpose wrapper to be able to run
 * all ITK filters that inherit from itk::ImageToImageFilter on a
 * Matlab 2D image or 3D or 4D image volume.
 *
 * B = ITK_IMFILTER(TYPE, A, [FILTER PARAMETERS])
 *
 *   TYPE is a string with the filter we want to run. See below for a whole
 *   list of options.
 *
 *   A is a 2D matrix, or 3D or 4D volume with the image or
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
 *   A can also be a SCI MAT struct, A = scimat, with the following fields:
 *
 *     scimat.data: 2D or 3D array with the image or segmentation, as above
 *     scimat.axis: 3x1 struct array with fields:
 *       scimat.axis.size:    number of voxels in the image
 *       scimat.axis.spacing: voxel size, image resolution
 *       scimat.axis.min:     real world coordinates of image origin
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
 * B = ITK_IMFILTER('skel', A)
 *
 *   (itk::BinaryThinningImageFilter3D)
 *   Skeletonize a binary mask
 *
 *   A is a segmentation.
 *
 *   B has the same size and class as A
 *
 * -------------------------------------------------------------------------
 *
 * [B, V, W] = itk_imfilter('dandist', A).
 * [B, V, W] = itk_imfilter('signdandist', A).
 *
 *   (itk::DanielssonDistanceMapImageFilter)
 *   (itk::SignedDanielssonDistanceMapImageFilter)
 *   Compute unsigned/signed distance map for a binary mask. Distance values are
 *   given in voxel coordinates.
 *
 *   A is a segmentation.
 *
 *   B has the same size as A and type float. Each element in B
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
 * -------------------------------------------------------------------------
 *
 * B = ITK_IMFILTER('maudist', A)
 *
 *   (itk::SignedMaurerDistanceMapImageFilter)
 *   Compute signed distance map for a binary mask. Distance values are
 *   given in real world coordinates, if the input image is given as a SCI
 *   MAT struct, or in voxel units, if the input image is a normal array. 
 *
 *   A is a segmentation.
 *
 *   B has the same size as A and type float.
 *
 * -------------------------------------------------------------------------
 *
 * B = ITK_IMFILTER('appsigndist', A)
 *
 *   (itk::ApproximateSignedDistanceMapImageFilter) 
 *   Compute signed distance map for a binary mask. Distance values
 *   are given in real world coordinates, if the input image is given
 *   as a SCIMAT struct, or in voxel units, if the input image is a
 *   plain array. The distances computed by this filter are Chamfer
 *   distances, which are only an approximation to Euclidian
 *   distances, and are not as exact approximations as those
 *   calculated by the DanielssonDistanceMapImageFilter. On the other
 *   hand, this filter is faster.
 *
 *   A is a segmentation.
 *
 *   B has the same size as A and type float.
 *
 * -------------------------------------------------------------------------
 *
 * B = ITK_IMFILTER('bwdilate', A, RADIUS, FOREGROUND)
 * B = ITK_IMFILTER('bwerode', A, RADIUS, FOREGROUND)
 *
 *   (itk::BinaryDilateImageFilter). 
 *   Binary dilation. The structuring element is a ball.
 *   (itk::BinaryErodeImageFilter).
 *   Binary erosion. The structuring element is a ball.
 *
 *   A is a segmentation.
 *
 *   RADIUS is a scalar with the radius of the ball in voxel units. If a
 *   non-integer number is provided, then floor(RADIUS) is used. By default,
 *   RADIUS = 0 and no dilation is performed.
 *
 *   FOREGROUND is a scalar. Voxels with that value will be the only ones
 *   dilated. By default, FOREGROUND is the maximum value allowed for the
 *   type, e.g. FOREGROUND=255 if the image is uint8. This is the default in
 *   ITK, so we respect it.
 *
 * -------------------------------------------------------------------------
 *
 * B = ITK_IMFILTER('advess', A, SIGMAMIN, SIGMAMAX, NUMSIGMASTEPS, NUMITERATIONS,
 *                  WSTRENGTH, SENSITIVITY, TIMESTEP, EPSILON)
 *
 *   (itk::AnisotropicDiffusionVesselEnhancementImageFilter)
 *   Anisotropic difussion vessel enhancement.
 *
 *   Enquobahrie A., Ibanez L., Bullitt E., Aylward S. "Vessel
 *   Enhancing Diffusion Filter", Insight Journal,
 *   2007. http://hdl.handle.net/1926/558.
 *
 *   A is a grayscale image.
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
 *   results seem better if run directly on the image. The
 *   filter doesn't seem to be spacing invariant.
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
 * -------------------------------------------------------------------------
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
 *   A is an image.
 *
 *   B has the same size as A and type double.
 *
 *   Input arguments are the same as the four first input arguments of
 *   filter "advess" above.
 *
 * -------------------------------------------------------------------------
 *
 * B = itk_imfilter('median', A, RADIUS)
 *
 *   (itk::MedianImageFilter)
 *   Median of a rectangular neighbourhood.
 *
 *   A is an image.
 *
 *   B has the same size and class as A.
 *
 *   RADIUS is a vector of scalars with the half-size of the filter's
 *   box in each dimension. E.g. RADIUS=[2, 3, 4] means that the
 *   median is computed in a rectangular neighbourhood of [5, 7, 9]
 *   voxels.
 *
 * -------------------------------------------------------------------------
 *
 * B = itk_imfilter('mrf', A, MU)
 *
 *   (itk::MRFImageFilter)
 *   Markov Random Field segmentation.
 *
 *   This filter can be used to improve a previous segmentation. A Markov
 *   Random Field (MRF) filter imposes the constraint that neighbouring
 *   voxels are more likely to have the same label. For
 *   example, gmth_seg() can be used to obtain a previous rough 2-label
 *   segmentation of an object over a background, and then the MRF filter
 *   applied to the computed Gaussian-mixture model mean values to smooth
 *   the segmentation.
 *
 *   A is an image.
 *
 *   B is a segmentation with the same size as A and type uint8.
 *
 *   MU is a row vector with the mean intensity values (centroids) of each
 *   label. If MU has n elements, then B will have n different labels. 
 *
 * B = itk_imfilter(..., WEIGHTS, SMOOTH, NITER, TOL)
 *
 *   WEIGHTS is an array with the same dimension as A. This array defines
 *   the neighbourhood of each pixel, and the relative importance each
 *   neighbouring pixel contributes to the labelling of the central pixel.
 *   By default, WEIGHTS is an array of all 1.0 (except for the central
 *   element, that is 0.0) and size 3x3, 3x3x3 or 3x3x3x3, depending on the
 *   image dimension.
 *
 *   SMOOTH is a scalar that represents the tradeoff between fidelity to the
 *   observed image and the smoothness of the segmented image. Typical
 *   smoothing factors have values from 1 to 5. This factor will multiply
 *   the weights that define the influence of neighbours on the
 *   classification of a given pixel.  The higher the value, the more
 *   uniform will be the regions resulting from the classification
 *   refinement. By default, SMOOTH=1e-7 and almost no smoothing is applied.
 *
 *   NITER is a scalar with the number of iterations the filter will run. By
 *   default, NITER=100.
 *
 *   TOL is a scalar with the error tolerance that will be used as a
 *   criterion for convergence. By default, TOL=1e-7.
 *
 * -------------------------------------------------------------------------
 *
 * B = itk_imfilter('voteholefill', A)
 *
 *   (itk::VotingBinaryIterativeHoleFillingImageFilter)
 *   Fills in holes and cavities by iteratively applying a voting operation.
 *
 *   A is a binary image.
 *
 *   B is a binary image of the same size and type as A.
 *
 * B = itk_imfilter(..., RADIUS, THR, BACKGROUND, FOREGROUND)
 *
 *   RADIUS is an array with the same dimension as A. RADIUS gives the
 *   radius of the box around the current voxel in each dimension. Each
 *   voxel within the box counts as a vote for whether the current
 *   background voxel should be flipped to foreground. By default RADIUS is
 *   1 in all dimensions, i.e. a box of side = 3.
 *
 *   THR is the majority threshold, i.e. the number of pixels over 50% that
 *   will decide whether a background pixel will become foreground or not.
 *
 *   BACKGROUND, FOREGROUND are the voxel values for background and
 *   foreground voxels, respectively. By default, BACKGROUND=0,
 *   FOREGROUND=1.
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011-2012 University of Oxford
  * Version: 1.3.0
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

/* ITK dependencies */
#include "itkImage.h"
#include "itkComposeImageFilter.h"
#include "itkFixedArray.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkMinimumDecisionRule.h"

/* ITK filter headers */
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.h"
#include "itkAnisotropicDiffusionVesselEnhancementImageFilter.h"
#include "itkBinaryThinningImageFilter3D.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkMRFImageFilter.h"

/* Gerardus headers */
#include "GerardusCommon.h"
#include "MatlabImageHeader.h"
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

// list of supported filters. It has to be an enum so that we can pass
// it as a template constant parameter
enum SupportedFilter {
  nVotingBinaryIterativeHoleFillingImageFilter,
  nApproximateSignedDistanceMapImageFilter,
  nMedianImageFilter,
  nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter,
  nAnisotropicDiffusionVesselEnhancementImageFilter,
  nBinaryThinningImageFilter3D,
  nSignedDanielssonDistanceMapImageFilter,
  nDanielssonDistanceMapImageFilter,
  nSignedMaurerDistanceMapImageFilter,
  nBinaryDilateImageFilter,
  nBinaryErodeImageFilter,
  nMRFImageFilter
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
template <class TPixelIn, unsigned int VImageDimension,
	  unsigned int FilterEnum>
class FilterWrapper {
public:
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    mexErrMsgTxt("Unsupported filter type");
  }
};

// VotingBinaryIterativeHoleFillingImageFilter
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nVotingBinaryIterativeHoleFillingImageFilter> {

public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 5);
    matlabExport->CheckNumberOfArguments(0, 1);

    // instantiate the filter
    typedef TPixelIn TPixelOut;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter<InImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->
		     GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));

    // default parameters
    typename InImageType::SizeType radiusDef;
    radiusDef.Fill(1);

    // filter parameters
    filter->SetRadius(matlabImport->template
		      GetRowVectorArgument<typename InImageType::SizeValueType,
					   typename InImageType::SizeType>(2, "RADIUS", radiusDef));
    filter->SetMaximumNumberOfIterations(matlabImport->template
					 GetScalarArgument<unsigned int>(3, "MAXITER", 1));
    filter->SetMajorityThreshold(matlabImport->template
				 GetScalarArgument<unsigned int>(4, "THR", 2));
    filter->SetBackgroundValue(matlabImport->template
			       GetScalarArgument<TPixelIn>(5, "BACKGROUND", 0));
    filter->SetForegroundValue(matlabImport->template
			       GetScalarArgument<TPixelIn>(6, "FOREGROUND", 1));

    // run filter
    filter->Update();

    // copy ITK filter outputs to Matlab outputs
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->CopyItkImageToMatlab<TPixelOut, VImageDimension>
	(filter->GetOutputs()[0], im.size, 0, "0");
    }

  }
};

// ApproximateSignedDistanceMapImageFilter
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nApproximateSignedDistanceMapImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 2);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // instantiate the filter
    typedef float TPixelOut;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef typename itk::Image<TPixelOut, VImageDimension> OutImageType;
    typedef itk::ApproximateSignedDistanceMapImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    // expect segmented object of 1s over background of 0s
    filter->SetInsideValue(1);
    filter->SetOutsideValue(0);
    
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

// MedianImageFilter
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension, 
		    nMedianImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 5);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // instantiate the filter
    typedef TPixelIn TPixelOut;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::MedianImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
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
		      GetRowVectorArgument<typename BoxFilterType::RadiusValueType, 
					   typename BoxFilterType::RadiusType>
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
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 6);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // instantiate the filter
    typedef double TPixelOut;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef typename itk::Image<TPixelOut, VImageDimension> OutImageType;
    typedef itk::MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter
      <InImageType, OutImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
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

template <class TPixelIn>
class FilterWrapper<TPixelIn, 2,
		    nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter only accepts 3D input images");
  }
};

template <class TPixelIn>
class FilterWrapper<TPixelIn, 4,
		    nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter only accepts 3D input images");
  }
};

// AnisotropicDiffusionVesselEnhancementImageFilter
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nAnisotropicDiffusionVesselEnhancementImageFilter> {
public:

  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 11);
    matlabExport->CheckNumberOfArguments(0, 1);

    // instantiate the filter
    typedef TPixelIn TPixelOut;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::AnisotropicDiffusionVesselEnhancementImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
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

template <class TPixelIn>
class FilterWrapper<TPixelIn, 2,
		    nAnisotropicDiffusionVesselEnhancementImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("AnisotropicDiffusionVesselEnhancementImageFilter only accepts 3D input images");
  }
};

template <class TPixelIn>
class FilterWrapper<TPixelIn, 4,
	      nAnisotropicDiffusionVesselEnhancementImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("AnisotropicDiffusionVesselEnhancementImageFilter only accepts 3D input images");
  }
};

// BinaryThinningImageFilter3D
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nBinaryThinningImageFilter3D> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 2);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // instantiate the filter
    typedef TPixelIn TPixelOut;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::BinaryThinningImageFilter3D<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
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

template <class TPixelIn>
class FilterWrapper<TPixelIn, 2,
		    nBinaryThinningImageFilter3D> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("BinaryThinningImageFilter3D only accepts 3D input images");
  }
};

template <class TPixelIn>
class FilterWrapper<TPixelIn, 4,
		    nBinaryThinningImageFilter3D> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("BinaryThinningImageFilter3D only accepts 3D input images");
  }
};

// SignedDanielssonDistanceMapImageFilter
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nSignedDanielssonDistanceMapImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 2);
    matlabExport->CheckNumberOfArguments(0, 3);
    
    // instantiate the filter
    typedef float TPixelOut;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef typename itk::Image<TPixelOut, VImageDimension> OutImageType;
    typedef itk::SignedDanielssonDistanceMapImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
    // connect Matlab inputs to ITK filter
    filter->SetInput(matlabImport->
		     GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));

    // connect ITK filter outputs to Matlab outputs

    // distance map
    if (matlabExport->GetNumberOfArguments() >= 1) {
      matlabExport->GraftItkImageOntoMatlab<TPixelOut, VImageDimension>
	// (filter->GetOutputs()[0], im.size, 0, "0");
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

// DanielssonDistanceMapImageFilter
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nDanielssonDistanceMapImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 2);
    matlabExport->CheckNumberOfArguments(0, 3);
    
    // instantiate the filter
    typedef double TPixelOut;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef typename itk::Image<TPixelOut, VImageDimension> OutImageType;
    typedef itk::DanielssonDistanceMapImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
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
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nSignedMaurerDistanceMapImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 2);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // instantiate the filter
    typedef float TPixelOut;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef typename itk::Image<TPixelOut, VImageDimension> OutImageType;
    typedef itk::SignedMaurerDistanceMapImageFilter<InImageType, OutImageType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
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

template <unsigned int VImageDimension>
class FilterWrapper<mxLogical, VImageDimension,
		    nSignedMaurerDistanceMapImageFilter> {
public:
  FilterWrapper(MatlabImportFilter::Pointer, MatlabExportFilter::Pointer,
		MatlabImageHeader &) {
    mexErrMsgTxt("SignedMaurerDistanceMapImageFilter does not accept input image with type boolean");
  }
};

// BinaryDilateImageFilter
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nBinaryDilateImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 4);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // instantiate the filter
    typedef TPixelIn TPixelOut;
    typedef itk::BinaryBallStructuringElement<TPixelIn, VImageDimension>
      StructuringElementType;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::BinaryDilateImageFilter<InImageType, OutImageType, StructuringElementType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
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
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nBinaryErodeImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(2, 4);
    matlabExport->CheckNumberOfArguments(0, 1);
    
    // instantiate the filter
    typedef TPixelIn TPixelOut;
    typedef itk::BinaryBallStructuringElement<TPixelIn, VImageDimension>
      StructuringElementType;
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    typedef InImageType OutImageType;
    typedef itk::BinaryErodeImageFilter<InImageType, OutImageType, StructuringElementType>
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    
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

// MRFImageFilter
template <class TPixelIn, unsigned int VImageDimension>
class FilterWrapper<TPixelIn, VImageDimension,
		    nMRFImageFilter> {
public:
  
  FilterWrapper(MatlabImportFilter::Pointer matlabImport,
		MatlabExportFilter::Pointer matlabExport,
		MatlabImageHeader &im) {
    
    // check number of input and output arguments
    matlabImport->CheckNumberOfArguments(3, 7);
    matlabExport->CheckNumberOfArguments(0, 1);

    /* type definitions */

    // input image
    typedef typename itk::Image<TPixelIn, VImageDimension> InImageType;
    
    // segmentation masks
    typedef unsigned char LabelPixelType;
    typedef itk::Image<LabelPixelType, VImageDimension> LabelImageType;

    // output pixel type
    typedef LabelPixelType TPixelOut;

    // dummy compose filter to convert the scalar image into a 1-vector image
    typedef itk::FixedArray<TPixelIn, 1> ArrayPixelType;
    typedef itk::Image<ArrayPixelType, VImageDimension> ArrayImageType;
    typedef itk::ComposeImageFilter<
      InImageType, ArrayImageType> ScalarToArrayFilterType;

    // filter
    typedef itk::MRFImageFilter<ArrayImageType, LabelImageType>
      FilterType;

    // classifier
    typedef itk::ImageClassifierBase<ArrayImageType, LabelImageType> SupervisedClassifierType;

    // decision rule
    typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;

    // membership function
    typedef itk::Statistics::DistanceToCentroidMembershipFunction<ArrayPixelType>
      MembershipFunctionType;
    typedef typename MembershipFunctionType::Pointer MembershipFunctionPointer;

    /* filter actions */    

    // instantiate the filter
    typename FilterType::Pointer filter = FilterType::New();

    /*    
     * get input arguments (grouped here for clarity)
     */
    
    // from the ITK guide: "Since the Markov Random Field algorithm is
    // defined in general for images whose pixels have multiple
    // components, that is, images of vector type, we must adapt our
    // scalar image in order to satisfy the interface expected by the
    // \code{MRFImageFilter}. We do this by using the
    // \doxygen{ComposeImageFilter}. With this filter we will present
    // our scalar image as a vector image whose vector pixels contain
    // a single component"
    typename ScalarToArrayFilterType::Pointer
      scalarToArrayFilter = ScalarToArrayFilterType::New();
    scalarToArrayFilter->SetInput(matlabImport->
    				  GetImageArgument<TPixelIn, VImageDimension>(1, "IM"));

    // vector of centroids
    std::vector<TPixelIn> centroid = matlabImport->template
      GetRowVectorArgument<TPixelIn, std::vector<TPixelIn> >(2, "MU", std::vector<TPixelIn>(0));
    unsigned int numberOfClasses = centroid.size();

    // by default, the neighbourhood is a hypercube with 1 voxel to
    // either side of the centre, i.e. a hypercube with side 3. All
    // elements of the default hypercube are 1.0, except for the
    // central pixel, that is 0.0
    mwSize neighLength = (mwSize)std::pow(3.0, (double)VImageDimension);
    std::vector<double> weights(neighLength, 1.0);
    weights[(neighLength-1)/2] = 0.0;
    typename InImageType::SizeType neighHalfSize;
    neighHalfSize.Fill(1);

    // read neighbourhood weights provided by the user, but as a vector
    weights = matlabImport->template
      GetArrayArgumentAsVector<std::vector<double> >(3, "WEIGHT", weights);
    
    // get size of neighbourhood weights array as provided by the
    // user. We get the half-size, as required by this filter (size =
    // 2 * halfsize + 1)
    neighHalfSize = matlabImport->template
      GetArrayHalfSize<typename InImageType::SizeValueType, 
		       typename InImageType::SizeType>(3, "WEIGHT", neighHalfSize);

    double smoothingFactor = matlabImport->template
      GetScalarArgument<double>(4, "SMOOTH", 1e-7);
    unsigned int maximumNumberOfIterations = matlabImport->template
      GetScalarArgument<unsigned int>(5, "NITER", 100);
    double errorTolerance = matlabImport->template
      GetScalarArgument<double>(6, "TOL", 1e-7);

    // ITK guide: "number of classes to be used during the
    // classification, the maximum number of iterations to be run in
    // this filter and the error tolerance that will be used as a
    // criterion for convergence"
    //
    // ITK guide: "the smoothing factor represents the tradeoff
    // between fidelity to the observed image and the smoothness of
    // the segmented image. Typical smoothing factors have values
    // between 1~5. This factor will multiply the weights that define
    // the influence of neighbors on the classification of a given
    // pixel.  The higher the value, the more uniform will be the
    // regions resulting from the classification refinement"
    filter->SetNumberOfClasses(numberOfClasses);
    filter->SetSmoothingFactor(smoothingFactor);
    filter->SetMaximumNumberOfIterations(maximumNumberOfIterations);
    filter->SetErrorTolerance(errorTolerance);

    // ITK guide: "Given that the MRF filter need to continually
    // relabel the pixels, it needs access to a set of membership
    // functions that will measure to what degree every pixel belongs
    // to a particular class.  The classification is performed by the
    // \doxygen{ImageClassifierBase} class, that is instantiated using
    // the type of the input vector image and the type of the labeled
    // image
    typename SupervisedClassifierType::Pointer classifier 
      = SupervisedClassifierType::New();

    // ITK guide: "The classifier needs a decision rule to be set by
    // the user. Note that we must use \code{GetPointer()} in the call
    // of the \code{SetDecisionRule()} method because we are passing a
    // SmartPointer, and smart pointers cannot perform polymorphism,
    // we must then extract the raw pointer that is associated to the
    // smart pointer. This extraction is done with the GetPointer()
    // method"
    //
    // MinimumDecisionRule returns the class label with the smallest
    // discriminant score
    typename DecisionRuleType::Pointer classifierDecisionRule 
      = DecisionRuleType::New();
    classifier->SetDecisionRule(classifierDecisionRule.GetPointer());

    // ITK guide: "we now instantiate the membership functions. In
    // this case we use the
    // \subdoxygen{Statistics}{DistanceToCentroidMembershipFunction}
    // class templated over the pixel type of the vector image, that
    // in our example happens to be a vector of dimension 1"
    double meanDistance = 0.0;
    typename MembershipFunctionType::CentroidType centroidAux(1);
    for(unsigned int i=0; i < numberOfClasses; i++) {
      MembershipFunctionPointer membershipFunction =
    	MembershipFunctionType::New();
      
      centroidAux[0] = centroid[i];
      
      membershipFunction->SetCentroid(centroidAux);
      
      classifier->AddMembershipFunction(membershipFunction);
      meanDistance += static_cast<double>(centroid[i]);
    }
    meanDistance /= numberOfClasses;
    
    // ITK guide: "and we set the neighborhood radius that will define
    // the size of the clique to be used in the computation of the
    // neighbors' influence in the classification of any given
    // pixel. Note that despite the fact that we call this a radius,
    // it is actually the half size of an hypercube. That is, the
    // actual region of influence will not be circular but rather an
    // N-Dimensional box. For example, a neighborhood radius of 2 in a
    // 3D image will result in a clique of size 5x5x5 pixels, and a
    // radius of 1 will result in a clique of size 3x3x3 pixels."
    filter->SetNeighborhoodRadius(neighHalfSize);

    // ITK guide: "We now scale weights so that the smoothing function
    // and the image fidelity functions have comparable value. This is
    // necessary since the label image and the input image can have
    // different dynamic ranges. The fidelity function is usually
    // computed using a distance function, such as the
    // \doxygen{DistanceToCentroidMembershipFunction} or one of the
    // other membership functions. They tend to have values in the
    // order of the means specified."
    double totalWeight = 0;
    for(std::vector<double>::const_iterator wcIt = weights.begin();
	wcIt != weights.end(); ++wcIt ) {
      totalWeight += *wcIt;
    }
    for(std::vector<double>::iterator wIt = weights.begin();
	wIt != weights.end(); wIt++) {
      *wIt = static_cast<double> ((*wIt) * meanDistance / (2 * totalWeight));
    }

    filter->SetMRFNeighborhoodWeight(weights);
    
    // ITK guide: "Finally, the classifier class is connected to the Markof Random Fields filter."
    filter->SetClassifier(classifier);

    // connect Matlab inputs to ITK filter
    filter->SetInput(scalarToArrayFilter->GetOutput());
    
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
  if (filterName == "appsigndist" 
  	     || filterName == "ApproximateSignedDistanceMapImageFilter") {

    FilterWrapper<TPixelIn, VImageDimension, 
		  nApproximateSignedDistanceMapImageFilter> 
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "median" 
  	     || filterName == "MedianImageFilter") {

    FilterWrapper<TPixelIn, VImageDimension, 
  		  nMedianImageFilter> 
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "advess" 
      || filterName == "AnisotropicDiffusionVesselEnhancementImageFilter") {

    FilterWrapper<TPixelIn, VImageDimension, 
  		  nAnisotropicDiffusionVesselEnhancementImageFilter> 
      wrapper(matlabImport, matlabExport, im);

  } else if (filterName == "bwdilate" 
  	     || filterName == "BinaryDilateImageFilter") {

    FilterWrapper<TPixelIn, VImageDimension, 
  		  nBinaryDilateImageFilter> 
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "bwerode" 
  	     || filterName == "BinaryErodeImageFilter") {

    FilterWrapper<TPixelIn, VImageDimension, 
  		  nBinaryErodeImageFilter> 
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "skel" 
  	     || filterName == "BinaryThinningImageFilter3D") {

    FilterWrapper<TPixelIn, VImageDimension, 
  		  nBinaryThinningImageFilter3D>
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "signdandist" 
  	     || filterName == "SignedDanielssonDistanceMapImageFilter") {

    FilterWrapper<TPixelIn, VImageDimension, 
  		  nSignedDanielssonDistanceMapImageFilter>
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "dandist" 
  	     || filterName == "DanielssonDistanceMapImageFilter") {

    FilterWrapper<TPixelIn, VImageDimension, 
  		  nDanielssonDistanceMapImageFilter>
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "hesves" 
  	     || filterName == "MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter") {

    FilterWrapper<TPixelIn, VImageDimension, 
  		  nMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter>
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "maudist" 
  	     || filterName == "SignedMaurerDistanceMapImageFilter") {

    FilterWrapper<TPixelIn, VImageDimension, 
  		  nSignedMaurerDistanceMapImageFilter>
      filterWrapper(matlabImport, matlabExport, im);

  } else if (filterName == "mrf" 
      || filterName == "MRFImageFilter") {
    
    FilterWrapper<TPixelIn, VImageDimension, 
  		  nMRFImageFilter> 
      filterWrapper(matlabImport, matlabExport, im);
    
  } else if (filterName == "voteholefill" 
      || filterName == "VotingBinaryIterativeHoleFillingImageFilter") {
    
    FilterWrapper<TPixelIn, VImageDimension, 
  		  nVotingBinaryIterativeHoleFillingImageFilter> 
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
  default:
    mexErrMsgTxt("Input image can only have 2 to 4 dimensions");
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
