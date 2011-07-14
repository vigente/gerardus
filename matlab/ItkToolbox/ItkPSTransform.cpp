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
 * YI = ITK_PSTRANSFORM(TRANSFORM, X, Y, XI)
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
 * YI = ITK_PSTRANSFORM('elastic', X, Y, XI)
 * YI = ITK_PSTRANSFORM('elasticr', X, Y, XI)
 * YI = ITK_PSTRANSFORM('tps', X, Y, XI)
 * YI = ITK_PSTRANSFORM('tpsr2', X, Y, XI)
 * YI = ITK_PSTRANSFORM('volume', X, Y, XI)
 *
 *   'elastic':  itk::ElasticBodySplineKernelTransform
 *   'elasticr': itk::ElasticBodyReciprocalSplineKernelTransform
 *   'tps':      itk::ThinPlateSplineKernelTransform
 *   'tpsr2':    itk::ThinPlateR2LogRSplineKernelTransform
 *   'volume':   itk::VolumeSplineKernelTransform
 *
 * YI = ITK_PSTRANSFORM('bspline', X, Y, XI, ORDER, LEVELS)
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
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.2.3
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

// this definition is necessary for ITK v3.20.0 to avoid an error when trying to
// compile itk::FixedArray::operator[](unsigned __int64) for Windows 64 bit, but
// maybe we can remove it when ITK v4.0.0 is released
#ifdef _WIN64
#define CAST2MWSIZE(x) static_cast<unsigned long long>(x)
#else
#define CAST2MWSIZE(x) x
#endif

/* Gerardus headers */

/* Functions */

// runBSplineTransform<TScalarType, Dimension>()
template <class TScalarType, unsigned int Dimension>
void runBSplineTransform(int nArgIn, const mxArray** argIn,
			 mxArray** argOut) {

  // check number of input arguments
  if (nArgIn > 6) {
    mexErrMsgTxt("Too many input arguments");
  }

  // spline order (input argument): default or user-provided
  unsigned int splineOrder = 3;
  if ((nArgIn > 4) && !mxIsEmpty(argIn[4])) {
    if (!mxIsDouble(argIn[4]) 
	|| mxIsComplex(argIn[4])
	|| (mxGetM(argIn[4]) != 1) 
	|| (mxGetN(argIn[4]) != 1)) {
      mexErrMsgTxt("ORDER must be an integer scalar of type double >= 0");
    }
    double *pSplineOrder = mxGetPr(argIn[4]);
    if ((pSplineOrder[0] < 0) 
	|| (pSplineOrder[0] != std::floor(pSplineOrder[0] + 0.5))) {
      mexErrMsgTxt("ORDER must be an integer scalar of type double >= 0");
    }
    splineOrder = (unsigned int)pSplineOrder[0];
  }  

  // number of levels (input argument): default or user-provided
  unsigned int numOfLevels = 5;
  if ((nArgIn > 5) && !mxIsEmpty(argIn[5])) {
    if (!mxIsDouble(argIn[5]) 
	|| mxIsComplex(argIn[5])
	|| (mxGetM(argIn[5]) != 1) 
	|| (mxGetN(argIn[5]) != 1)) {
      mexErrMsgTxt("LEVELS must be an integer scalar of type double >= 1");
    }
    double *pNumOfLevels = mxGetPr(argIn[5]);
    if ((pNumOfLevels[0] < 1) 
	|| (pNumOfLevels[0] != std::floor(pNumOfLevels[0] + 0.5))) {
      mexErrMsgTxt("LEVELS must be an integer scalar of type double >= 1");
    }
    numOfLevels = (unsigned int)pNumOfLevels[0];
  }  

  // get size of input arguments
  mwSize Mx = mxGetM(argIn[1]); // number of source points
  mwSize Mxi = mxGetM(argIn[3]); // number of points to be warped

  // create pointers to input matrices
  TScalarType *x = (TScalarType *)mxGetData(argIn[1]); // source points
  TScalarType *y = (TScalarType *)mxGetData(argIn[2]); // target points
  TScalarType *xi = (TScalarType *)mxGetData(argIn[3]); // points to be warped

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
    orig[col] = std::numeric_limits<TScalarType>::max();
    term[col] = std::numeric_limits<TScalarType>::min();
  }

  // find bounding box limits
  for (mwSize row=0; row < Mx; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      orig[col] = std::min((TScalarType)orig[col], x[Mx * col + row]);
      term[col] = std::max((TScalarType)term[col], x[Mx * col + row]);
      orig[col] = std::min((TScalarType)orig[col], y[Mx * col + row]);
      term[col] = std::max((TScalarType)term[col], y[Mx * col + row]);
    }
  }
  for (mwSize row=0; row < Mxi; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      orig[col] = std::min((TScalarType)orig[col], xi[Mxi * col + row]);
      term[col] = std::max((TScalarType)term[col], xi[Mxi * col + row]);
    }
  }

  // compute length of each size of the bounding box
  DataType len = term - orig;
  TScalarType lenmax = std::numeric_limits<TScalarType>::min();
  for (mwSize col=0; col < (mwSize)Dimension; ++col) {
    lenmax = std::max(lenmax, len[col]);
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
      v[col] = (y[Mx * col + row] - x[Mx * col + row]) / lenmax;
      xParam[col] = (x[Mx * col + row] - orig[col]) / lenmax;
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

  // number of dimension of points to be warped
  mwSize ndimxi = mxGetNumberOfDimensions(argIn[3]); 
  // dimensions vector of array of points to be warped
  const mwSize *dimsxi = mxGetDimensions(argIn[3]);

  // create output vector and pointer to populate it
  argOut[0] = (mxArray *)mxCreateNumericArray(ndimxi, 
  					      dimsxi, 
  					      mxGetClassID(argIn[1]), 
  					      mxREAL);
  TScalarType *yi = (TScalarType *)mxGetData(argOut[0]);

  // sample the warp field
  DataType vi; // warp field sample
  typename PointSetType::PointType xiParam; // sampling coordinates
  for (mwSize row=0; row < Mxi; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      xiParam[col] = (xi[Mxi * col + row] - orig[col]) / lenmax;
    }
    transform->Evaluate(xiParam, vi);
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      yi[Mxi * col + row] = xi[Mxi * col + row] + vi[col] * lenmax;
    }
  }

  // exit function
  return;

}

// runKernelTransform<TScalarType, Dimension, TransformType>()
template <class TScalarType, unsigned int Dimension, class TransformType>
void runKernelTransform(int nArgIn, const mxArray** argIn,
			mxArray** argOut) {

  // check number of input arguments
  if (nArgIn > 4) {
    mexErrMsgTxt("Too many input arguments");
  }

  // get size of input arguments
  mwSize Mx = mxGetM(argIn[1]); // number of source points
  mwSize Mxi = mxGetM(argIn[3]); // number of points to be warped
  mwSize ndimxi; // number of dimension of points to be warped
  const mwSize *dimsxi; // dimensions vector of array of points to be warped

  // create pointers to input matrices
  TScalarType *x = (TScalarType *)mxGetData(argIn[1]); // source points
  TScalarType *y = (TScalarType *)mxGetData(argIn[2]); // target points
  TScalarType *xi = (TScalarType *)mxGetData(argIn[3]); // points to be warped

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
  ndimxi = mxGetNumberOfDimensions(argIn[3]);
  dimsxi = mxGetDimensions(argIn[3]);
  argOut[0] = (mxArray *)mxCreateNumericArray(ndimxi, 
  					      dimsxi, 
  					      mxGetClassID(argIn[1]), 
  					      mxREAL);
  TScalarType *yi = (TScalarType *)mxGetData(argOut[0]);

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
void parseTransformType(int nArgIn, const mxArray** argIn,
			mxArray** argOut) {
  
  // get type of transform
  char *transform = mxArrayToString(argIn[0]);
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
		       ElasticTransformType>(nArgIn, argIn, argOut);
  } else if (!strcmp(transform, "elasticr")) {
    runKernelTransform<TScalarType, Dimension, 
		       ElasticReciprocalTransformType>(nArgIn, argIn, argOut);
  } else if (!strcmp(transform, "tps")) {
    runKernelTransform<TScalarType, Dimension, 
		       TpsTransformType>(nArgIn, argIn, argOut);
  } else if (!strcmp(transform, "tpsr2")) {
    runKernelTransform<TScalarType, Dimension, 
		       TpsR2LogRTransformType>(nArgIn, argIn, argOut);
  } else if (!strcmp(transform, "volume")) {
    runKernelTransform<TScalarType, Dimension, 
		       VolumeTransformType>(nArgIn, argIn, argOut);
  } else if (!strcmp(transform, "bspline")) {
    runBSplineTransform<TScalarType, Dimension>(nArgIn, argIn, argOut);
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
void parseDimensionToTemplate(int nArgIn, const mxArray** argIn,
			      mxArray** argOut) {

  // dimension of points
  mwSize Nx = mxGetN(argIn[1]);

  // parse the dimension value
  switch (Nx) {
  case 2:
    parseTransformType<TScalarType, 2>(nArgIn, argIn, argOut);
    break;
  case 3:
    parseTransformType<TScalarType, 3>(nArgIn, argIn, argOut);
    break;
  default:
    mexErrMsgTxt("Input points can only have dimensions 2 or 3");
    break;
  }

  // exit function
  return;

}

// parseInputTypeToTemplate()
void parseInputTypeToTemplate(int nArgIn, const mxArray** argIn,
			      mxArray** argOut) {

  // point coordinate type
  mxClassID pointCoordClassId = mxGetClassID(argIn[1]);

  // check that all point coordinates have the same type (it simplifies
  // things with templates)
  if ((pointCoordClassId != mxGetClassID(argIn[2]))
      | (pointCoordClassId != mxGetClassID(argIn[3]))) {
    mexErrMsgTxt("Input arguments X, Y and XI must have the same type");
  }
  
  // swith input image type
  switch(pointCoordClassId) {
  case mxDOUBLE_CLASS:
    parseDimensionToTemplate<double>(nArgIn, argIn, argOut);
    break;
  case mxSINGLE_CLASS:
    parseDimensionToTemplate<float>(nArgIn, argIn, argOut);
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

  // check number of input and output arguments
  if (nrhs < 4) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  // if there are no points to warp, return empty array
  if (mxIsEmpty(prhs[3])) {
    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    return;
  }

  // check size of input arguments
  mwSize Mx = mxGetM(prhs[1]); // number of source points
  mwSize My = mxGetM(prhs[2]); // number of target points
  mwSize Dimension = mxGetN(prhs[1]); // dimension of source points
  mwSize dimy = mxGetN(prhs[2]); // dimension of target points
  mwSize dimxi = mxGetN(prhs[3]); // dimension of points to be warped
  
  // the landmark arrays must have the same number of points
  // (degenerate case, both are empty)
  if (Mx != My) mexErrMsgTxt("X and Y must have the same number of points (i.e. rows).");

  // if there are no landmarks, we apply no transformation to the
  // points to warp
  if (mxIsEmpty(prhs[1])) {
    plhs[0] = mxDuplicateArray(prhs[3]);
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
  parseInputTypeToTemplate(nrhs, prhs, plhs);

  // exit successfully
  return;

}

#endif /* ITKPSTRANSFORM_CPP */
