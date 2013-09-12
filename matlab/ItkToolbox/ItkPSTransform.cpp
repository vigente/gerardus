/* ItkPSTransform.cpp
 *
 * ITK_PSTRANSFORM: Spatial transformation (i.e. warp) of a set of points,
 * defined from a known landmark correspondence
 *
 * This MEX-function provides a Matlab interface to run ITK's Kernel
 * and B-spline transforms:
 *
 *   itk::ElasticBodySplineKernelTransform
 *   itk::ElasticBodyReciprocalSplineKernelTransform
 *   itk::ThinPlateSplineKernelTransform
 *   itk::ThinPlateR2LogRSplineKernelTransform
 *   itk::VolumeSplineKernelTransform
 *   itk::BSplineScatteredDataPointSetToImageFilter
 *
 *
 * YI = itk_pstransform(TRANSFORM, X, Y, XI)
 *
 *   X, Y are 2-column (2D) or 3-column (3D) matrices. Each row has
 *   the coordinates of a point. The warp is defined so that
 *   X(i,:)->Y(i,:).
 *
 *   XI is a matrix with the same number of columns as X, Y. Each row
 *   has the coordinates of a point to be warped.
 *
 *   YI has the same dimensions as XI. YI contains the coordinates of
 *   the warped points.
 *
 *   TRANSFORM is a string that allows to select the type of warp (no
 *   defaults):
 *
 * YI = itk_pstransform('elastic', X, Y, XI)
 * YI = itk_pstransform('elasticr', X, Y, XI)
 * YI = itk_pstransform('tps', X, Y, XI)
 * YI = itk_pstransform('tpsr2', X, Y, XI)
 * YI = itk_pstransform('volume', X, Y, XI)
 *
 *   'elastic':  itk::ElasticBodySplineKernelTransform
 *   'elasticr': itk::ElasticBodyReciprocalSplineKernelTransform
 *   'tps':      itk::ThinPlateSplineKernelTransform
 *   'tpsr2':    itk::ThinPlateR2LogRSplineKernelTransform
 *   'volume':   itk::VolumeSplineKernelTransform
 *
 *   Note that 'tpsr2' produces the same result as our Matlab implementation
 *   pts_tps_map(), as it implements the classic kernel proposed by
 *   Bookstein, r^2 ln(r^2).
 *
 * YI = itk_pstransform('bspline', X, Y, XI, ORDER, LEVELS)
 *
 *   'bspline':  itk::BSplineScatteredDataPointSetToImageFilter
 *
 *   By Nicholas J. Tustison, James C. Gee in the Insight Journal
 *   paper: http://hdl.handle.net/1926/140
 *
 *   ORDER is an integer with the B-spline order. By default, ORDER=3,
 *   and the B-spline is cubic.
 *
 *   LEVELS is an integer with the number of multi-resolution levels
 *   in the algorithm. A higher number of levels will make the spline
 *   more flexible and match the landmarks better. By default, LEVELS=5.
 *
 * See also: pts_tps_map, pts_tps_weights.
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011-2013 University of Oxford
  * Version: 0.5.0
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

#ifndef ITKPSTRANSFORM_CPP
#define ITKPSTRANSFORM_CPP

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>

/* ITK headers */
#include "itkImage.h"
#include "itkKernelTransform.h"
#include "itkElasticBodySplineKernelTransform.h"
#include "itkElasticBodyReciprocalSplineKernelTransform.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "itkVolumeSplineKernelTransform.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#if ITK_VERSION_MAJOR>=4
#include "itkBSplineControlPointImageFunction.h"
#endif

/* Gerardus headers */
#include "GerardusCommon.h"
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

/* Inputs/outputs interfaces */
enum InputIndexType {IN_TRANSFORM, IN_X, IN_Y, IN_XI, 
		     IN_ORDER, IN_LEVELS, InputIndexType_MAX}; // IN_ORDER, IN_LEVELS only for B-spline
enum OutputIndexType {OUT_YI, OutputIndexType_MAX};

/* Functions */

// runBSplineTransform<TScalarType, Dimension>()
template <class TScalarType, unsigned int Dimension>
void runBSplineTransform(MatlabImportFilter::Pointer matlabImport,
			 MatlabExportFilter::Pointer matlabExport) {

  // retrieve pointers to the inputs that we are going to need here
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer; 
  MatlabInputPointer inX         = matlabImport->GetRegisteredInput("X");
  MatlabInputPointer inY         = matlabImport->GetRegisteredInput("Y");
  MatlabInputPointer inXI        = matlabImport->GetRegisteredInput("XI");
  MatlabInputPointer inORDER     = matlabImport->GetRegisteredInput("ORDER");
  MatlabInputPointer inLEVELS    = matlabImport->GetRegisteredInput("LEVELS");

  // register the output for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outYI = matlabExport->RegisterOutput(OUT_YI, "YI");

  // spline order (input argument): default or user-provided
  unsigned int splineOrder = matlabImport->ReadScalarFromMatlab<unsigned int>(inORDER, 3);

  // number of levels (input argument): default or user-provided
  unsigned int numOfLevels = matlabImport->ReadScalarFromMatlab<unsigned int>(inLEVELS, 5);

  // get size of input arguments
  mwSize Mx = mxGetM(inX->pm); // number of source points
  mwSize Mxi = mxGetM(inXI->pm); // number of points to be warped

  // pointers to input matrices
  TScalarType *x 
    = (TScalarType *)mxGetData(inX->pm); // source points
  TScalarType *y 
    = (TScalarType *)mxGetData(inY->pm); // target points
  TScalarType *xi 
    = (TScalarType *)mxGetData(inXI->pm); // points to be warped
  if (x == NULL) {
    mexErrMsgTxt("Cannot get a pointer to input X");
  }
  if (y == NULL) {
    mexErrMsgTxt("Cannot get a pointer to input Y");
  }
  if (xi == NULL) {
    mexErrMsgTxt("Cannot get a pointer to input XI");
  }

  // type definitions for the BSPline transform
  typedef itk::Vector<TScalarType, Dimension> DataType;
  typedef itk::PointSet<DataType, Dimension> PointSetType;
  typedef itk::Image<DataType, Dimension> ImageType;
  typedef typename 
    itk::BSplineScatteredDataPointSetToImageFilter<PointSetType, 
						   ImageType> TransformType;

  // variables to store the input points
  typename PointSetType::Pointer pointSet = PointSetType::New();
  typename PointSetType::PointType xParam;
  DataType v; // v = y-x, i.e. displacement vector between source and
              // target landmark

  // init variables to contain the limits of a bounding box that
  // contains all the points
  typename ImageType::PointType orig, term;
  for (mwSize col=0; col < (mwSize)Dimension; ++col) {
    orig[CAST2MWSIZE(col)] = std::numeric_limits<TScalarType>::max();
    term[CAST2MWSIZE(col)] = std::numeric_limits<TScalarType>::min();
  }

  // find bounding box limits
  for (mwSize row=0; row < Mx; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      orig[CAST2MWSIZE(col)] = std::min((TScalarType)orig[CAST2MWSIZE(col)], x[Mx * col + row]);
      term[CAST2MWSIZE(col)] = std::max((TScalarType)term[CAST2MWSIZE(col)], x[Mx * col + row]);
      orig[CAST2MWSIZE(col)] = std::min((TScalarType)orig[CAST2MWSIZE(col)], y[Mx * col + row]);
      term[CAST2MWSIZE(col)] = std::max((TScalarType)term[CAST2MWSIZE(col)], y[Mx * col + row]);
    }
  }
  for (mwSize row=0; row < Mxi; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      orig[CAST2MWSIZE(col)] = std::min((TScalarType)orig[CAST2MWSIZE(col)], xi[Mxi * col + row]);
      term[CAST2MWSIZE(col)] = std::max((TScalarType)term[CAST2MWSIZE(col)], xi[Mxi * col + row]);
    }
  }

  // compute length of each size of the bounding box
  DataType len = term - orig;
  TScalarType lenmax = std::numeric_limits<TScalarType>::min();
  for (mwSize col=0; col < (mwSize)Dimension; ++col) {
    lenmax = std::max(lenmax, len[CAST2MWSIZE(col)]);
  }

  // duplicate the input x and y matrices to PointSet format so that
  // we can pass it to the ITK function
  //
  // we also translate and scale all points so that the bounding box
  // fits within the domain [0, 1] x [0, 1] x [0,1]. We need to do
  // this because the BSpline function requires the parametric domain
  // to be within [0, 1] x [0, 1] x [0,1]
  for (mwSize row=0; row < Mx; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      v[CAST2MWSIZE(col)] = (y[Mx * col + row] - x[Mx * col + row]) / lenmax;
      xParam[CAST2MWSIZE(col)] = (x[Mx * col + row] - orig[CAST2MWSIZE(col)]) / lenmax;
    }
    pointSet->SetPoint(row, xParam);
    pointSet->SetPointData(row, v);
  }

  // instantiate and set-up transform
  typename TransformType::Pointer transform = TransformType::New();
  transform->SetGenerateOutputImage(false);
  transform->SetInput(pointSet);
  transform->SetSplineOrder(splineOrder);
  typename TransformType::ArrayType ncps ;
  ncps.Fill(splineOrder + 1);
  transform->SetNumberOfControlPoints(ncps);
  transform->SetNumberOfLevels(numOfLevels);

  // note that closedim, spacing, sz and orig are all refered to the
  // parametric domain, i.e. the domain of x and xi
  typename TransformType::ArrayType closedim;
  typename ImageType::SpacingType spacing;
  typename ImageType::SizeType sz;

  // the parametric domain is not periodic in any dimension
  closedim.Fill(0);

  // as we are not creating the image, we don't need to provide a
  // sensible number of voxels. But size has to be at least 2 voxels
  // to avoid a run-time error
  sz.Fill(2);

  // because the parameterization is in [0, 1] x [0, 1] x [0,1], and
  // we have only size = 2 voxels in every dimension, the spacing will
  // be 1.0 / (2 - 1) = 1.0
  spacing.Fill(1.0);

  // because of the reparameterization, the origin we have to pass to
  // the transform is not the origin of the real points, but the
  // origin of the [0, 1] x [0, 1] x [0,1] bounding box
  typename ImageType::PointType origZero;
  origZero.Fill(0.0);
  
  transform->SetCloseDimension(closedim);
  transform->SetSize(sz);
  transform->SetSpacing(spacing);
  transform->SetOrigin(origZero);

  // run transform
  transform->Update();

  // create output vector and pointer to populate it
  mwSize ndimxi = mxGetNumberOfDimensions(inXI->pm); 
  const mwSize *dimsxi = mxGetDimensions(inXI->pm);
  std::vector<mwSize> size;
  for (mwIndex i = 0; i < ndimxi; ++i) {
    size.push_back(dimsxi[i]);
  }
  TScalarType *yi 
    = matlabExport->AllocateNDArrayInMatlab<TScalarType>(outYI, size);

  // from ITK v4.x, we need to instantiate a function to evaluate
  // points of the B-spline, as the Evaluate() method has been removed
  // from the TransformType
#if ITK_VERSION_MAJOR>=4
  // Note: in the following, we have to use TCoordRep=double, because
  // ITK gives a compilation error of an abstract class not having
  // been implemented. Otherwise, we would use
  // TCoordRep=TScalar=float, as in the rest of this program
  typedef typename 
    itk::BSplineControlPointImageFunction<ImageType, double> EvalFunctionType;
  typename EvalFunctionType::Pointer function = EvalFunctionType::New();

  function->SetSplineOrder(splineOrder);
  function->SetOrigin(origZero);
  function->SetSpacing(spacing);
  function->SetSize(sz);
  function->SetInputImage(transform->GetPhiLattice());
#endif

  // sample the warp field
  DataType vi; // warp field sample
  typename PointSetType::PointType xiParam; // sampling coordinates
  for (mwSize row=0; row < Mxi; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      xiParam[CAST2MWSIZE(col)] = (xi[Mxi * col + row] - orig[CAST2MWSIZE(col)]) / lenmax;
    }
#if ITK_VERSION_MAJOR<4
    transform->Evaluate(xiParam, vi);
#else
    vi = function->Evaluate(xiParam);
#endif
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      yi[Mxi * col + row] = xi[Mxi * col + row] + vi[CAST2MWSIZE(col)] * lenmax;
    }
  }

  // exit function
  return;

}

// runKernelTransform<TScalarType, Dimension, TransformType>()
template <class TScalarType, unsigned int Dimension, class TransformType>
void runKernelTransform(MatlabImportFilter::Pointer matlabImport,
			MatlabExportFilter::Pointer matlabExport) {

  // check number of input arguments (the kernel transform syntax
  // accepts up to 4 arguments only. Thus, we cannot use InputIndexType_MAX)
  matlabImport->CheckNumberOfArguments(4, 4);

  // retrieve pointers to the inputs that we are going to need here
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer; 
  MatlabInputPointer inX         = matlabImport->GetRegisteredInput("X");
  MatlabInputPointer inY         = matlabImport->GetRegisteredInput("Y");
  MatlabInputPointer inXI        = matlabImport->GetRegisteredInput("XI");

  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outYI = matlabExport->RegisterOutput(OUT_YI, "YI");

  // get size of input arguments
  mwSize Mx = mxGetM(inX->pm); // number of source points
  mwSize Mxi = mxGetM(inXI->pm); // number of points to be warped
  mwSize ndimxi; // number of dimension of points to be warped
  const mwSize *dimsxi; // dimensions vector of array of points to be warped

  // create pointers to input matrices
  TScalarType *x 
    = (TScalarType *)mxGetData(inX->pm); // source points
  TScalarType *y 
    = (TScalarType *)mxGetData(inY->pm); // target points
  TScalarType *xi 
    = (TScalarType *)mxGetData(inXI->pm); // points to be warped
  if (x == NULL) {
    mexErrMsgTxt("Cannot get a pointer to input X");
  }
  if (y == NULL) {
    mexErrMsgTxt("Cannot get a pointer to input Y");
  }
  if (xi == NULL) {
    mexErrMsgTxt("Cannot get a pointer to input XI");
  }

  // type definitions and variables to store points for the kernel transform
  typedef typename TransformType::PointSetType PointSetType;
  typename PointSetType::Pointer fixedPointSet = PointSetType::New();
  typename PointSetType::Pointer movingPointSet = PointSetType::New();
  typename PointSetType::Pointer toWarpPointSet = PointSetType::New();
  typedef typename PointSetType::PointsContainer PointsContainer;
  typename PointsContainer::Pointer fixedPointContainer = PointsContainer::New();
  typename PointsContainer::Pointer movingPointContainer = PointsContainer::New();
  typedef typename PointSetType::PointType PointType;
  PointType fixedPoint;
  PointType movingPoint;
  PointType toWarpPoint;
  PointType warpedPoint;

  // duplicate the input x and y matrices to PointSet format so that
  // we can pass it to the ITK function
  mwSize pointId=0;
  for (mwSize row=0; row < Mx; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      fixedPoint[CAST2MWSIZE(col)] = y[Mx * col + row];
      movingPoint[CAST2MWSIZE(col)] = x[Mx * col + row];
    }
    fixedPointContainer->InsertElement(pointId, fixedPoint);
    movingPointContainer->InsertElement(pointId, movingPoint);
    ++pointId;
  }
  fixedPointSet->SetPoints(fixedPointContainer);
  movingPointSet->SetPoints(movingPointContainer);

  // compute the transform
  typename TransformType::Pointer transform;
  transform = TransformType::New();
  
  transform->SetSourceLandmarks(movingPointSet);
  transform->SetTargetLandmarks(fixedPointSet);
  transform->ComputeWMatrix();
  
  // create output vector and pointer to populate it
  ndimxi = mxGetNumberOfDimensions(inXI->pm);
  dimsxi = mxGetDimensions(inXI->pm);
  std::vector<mwSize> size;
  for (mwIndex i = 0; i < ndimxi; ++i) {
    size.push_back(dimsxi[i]);
  }

  TScalarType *yi 
    = matlabExport->AllocateNDArrayInMatlab<TScalarType>(outYI, size);

  // transform points
  for (mwSize row=0; row < Mxi; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      toWarpPoint[CAST2MWSIZE(col)] = xi[Mxi * col + row];
    }
    warpedPoint = transform->TransformPoint(toWarpPoint);
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      yi[Mxi * col + row] = warpedPoint[CAST2MWSIZE(col)];
    }
  }
  
  // exit function
  return;
  
}

// parseTransformType<TScalarType, Dimension>()
template <class TScalarType, unsigned int Dimension>
void parseTransformType(MatlabImportFilter::Pointer matlabImport,
			MatlabExportFilter::Pointer matlabExport) {


  // retrieve pointers to the inputs that we are going to need here
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer; 
  MatlabInputPointer inTRANSFORM = matlabImport->GetRegisteredInput("TRANSFORM");

  // get type of transform
  char *transform = mxArrayToString(inTRANSFORM->pm);
  if (transform == NULL) {
    mexErrMsgTxt("Cannot read TRANSFORM string");
  }
  
  // kernel transform types
  typedef itk::ElasticBodySplineKernelTransform<TScalarType, Dimension> 
    ElasticTransformType;
  typedef itk::ElasticBodyReciprocalSplineKernelTransform<TScalarType, Dimension> 
    ElasticReciprocalTransformType;
  typedef itk::ThinPlateSplineKernelTransform<TScalarType, Dimension> 
    TpsTransformType;
  typedef itk::ThinPlateR2LogRSplineKernelTransform<TScalarType, Dimension> 
    TpsR2LogRTransformType;
  typedef itk::VolumeSplineKernelTransform<TScalarType, Dimension> 
    VolumeTransformType;

  // select transform function
  if (!strcmp(transform, "elastic")) {
    runKernelTransform<TScalarType, Dimension, 
		       ElasticTransformType>(matlabImport, matlabExport);
  } else if (!strcmp(transform, "elasticr")) {
    runKernelTransform<TScalarType, Dimension, 
		       ElasticReciprocalTransformType>(matlabImport, matlabExport);
  } else if (!strcmp(transform, "tps")) {
    runKernelTransform<TScalarType, Dimension, 
		       TpsTransformType>(matlabImport, matlabExport);
  } else if (!strcmp(transform, "tpsr2")) {
    runKernelTransform<TScalarType, Dimension, 
		       TpsR2LogRTransformType>(matlabImport, matlabExport);
  } else if (!strcmp(transform, "volume")) {
    runKernelTransform<TScalarType, Dimension, 
		       VolumeTransformType>(matlabImport, matlabExport);
  } else if (!strcmp(transform, "bspline")) {
    runBSplineTransform<TScalarType, Dimension>(matlabImport, matlabExport);
  } else if (!strcmp(transform, "")) {
    std::cout << 
      "Implemented transform types: elastic, elasticr, tps, tpsr2, volume, bspline" 
	      << std::endl;
    return;
  } else {
    mexErrMsgTxt("Transform not implemented");
  }

  // exit function
  return;
  
}

// parseDimensionToTemplate<TScalarType>()
template <class TScalarType>
void parseDimensionToTemplate(MatlabImportFilter::Pointer matlabImport,
			      MatlabExportFilter::Pointer matlabExport) {

  // retrieve pointers to the inputs that we are going to need here
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer; 
  MatlabInputPointer inX         = matlabImport->GetRegisteredInput("X");

  // dimension of points
  mwSize Nx = mxGetN(inX->pm);

  // parse the dimension value
  switch (Nx) {
  case 2:
    parseTransformType<TScalarType, 2>(matlabImport, matlabExport);
    break;
  case 3:
    parseTransformType<TScalarType, 3>(matlabImport, matlabExport);
    break;
  default:
    mexErrMsgTxt("Input points can only have dimensions 2 or 3");
    break;
  }

  // exit function
  return;

}

// parseInputTypeToTemplate()
void parseInputTypeToTemplate(MatlabImportFilter::Pointer matlabImport,
			      MatlabExportFilter::Pointer matlabExport) {

  // retrieve pointers to the inputs that we are going to need here
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer; 
  MatlabInputPointer inX         = matlabImport->GetRegisteredInput("X");
  MatlabInputPointer inY         = matlabImport->GetRegisteredInput("Y");
  MatlabInputPointer inXI        = matlabImport->GetRegisteredInput("XI");

  // point coordinate type
  mxClassID pointCoordClassId = mxGetClassID(inX->pm);

  // check that all point coordinates have the same type (it simplifies
  // things with templates)
  if ((pointCoordClassId != mxGetClassID(inY->pm))
      | (pointCoordClassId != mxGetClassID(inXI->pm))) {
    mexErrMsgTxt("Input arguments X, Y and XI must have the same type");
  }
  
  // swith input point type
  switch(pointCoordClassId) {
  case mxDOUBLE_CLASS:
    parseDimensionToTemplate<double>(matlabImport, matlabExport);
    break;
  case mxSINGLE_CLASS:
    parseDimensionToTemplate<float>(matlabImport, matlabExport);
    break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt("Point coordinates have unknown type");
    break;
  default:
    mexErrMsgTxt("Point coordinates can only be of type single or double");
    break;
  }

  // exit function
  return;
  
}

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);

  // register all possible inputs for this function at the import filter
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer;
  MatlabInputPointer inTRANSFORM = matlabImport->RegisterInput(IN_TRANSFORM, "TRANSFORM");
  MatlabInputPointer inX         = matlabImport->RegisterInput(IN_X, "X");
  MatlabInputPointer inY         = matlabImport->RegisterInput(IN_Y, "Y");
  MatlabInputPointer inXI        = matlabImport->RegisterInput(IN_XI, "XI");
  MatlabInputPointer inORDER     = matlabImport->RegisterInput(IN_ORDER, "ORDER");
  MatlabInputPointer inLEVELS    = matlabImport->RegisterInput(IN_LEVELS, "LEVELS");

  // interface to deal with output arguments from Matlab
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->ConnectToMatlabFunctionOutput(nlhs, plhs);
    
  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outYI = matlabExport->RegisterOutput(OUT_YI, "YI");

  // check number of input and output arguments
  matlabImport->CheckNumberOfArguments(4, InputIndexType_MAX);
  matlabExport->CheckNumberOfArguments(0, OutputIndexType_MAX);
    
  // if there are no points to warp, return empty array
  if (mxIsEmpty(inXI->pm)) {
    matlabExport->CopyEmptyArrayToMatlab(outYI);
    return;
  }

  // check size of input arguments
  mwSize Mx = mxGetM(inX->pm); // number of source points
  mwSize My = mxGetM(inY->pm); // number of target points
  mwSize Dimension = mxGetN(inX->pm); // dimension of source points
  mwSize dimy = mxGetN(inY->pm); // dimension of target points
  mwSize dimxi = mxGetN(inXI->pm); // dimension of points to be warped
  
  // the landmark arrays must have the same number of points
  // (degenerate case, both are empty)
  if (Mx != My) {
    mexErrMsgTxt("X and Y must have the same number of points (i.e. rows).");
  }

  // if there are no landmarks, we apply no transformation to the
  // points to warp
  if (mxIsEmpty(inX->pm)) {
    *outYI->ppm = mxDuplicateArray(inXI->pm);
    return;
  }

  // if there are landmarks and points to warp, all must have the same dimension
  if (Dimension != dimy || Dimension != dimxi) {
    mexErrMsgTxt("X, Y and XI must all have the same dimension (i.e. number of columns).");
  }

  // run filter (this function starts a cascade of functions designed
  // to translate the run-time type variables like inputVoxelClassId
  // to templates, so that we don't need to nest lots of "switch" or
  // "if" statements)
  parseInputTypeToTemplate(matlabImport, matlabExport);

  // exit successfully
  return;

}

#endif /* ITKPSTRANSFORM_CPP */
