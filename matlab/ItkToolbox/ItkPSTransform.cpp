/* ItkPSTransform.cpp
 *
 * ITK_PSTRANSFORM: Spatial transformation or warp on a point set
 * defined from a known landmark correspondence
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.1.1
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
// #include <cmath>
// #include <matrix.h>
// #include <vector>

/* ITK headers */
#include "itkElasticBodySplineKernelTransform.h"
#include "itkElasticBodyReciprocalSplineKernelTransform.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "itkVolumeSplineKernelTransform.h"

/* Gerardus headers */

/* Functions */

// runKernelTransform<TScalarType, Dimension, TransformType>()
template <class TScalarType, unsigned int Dimension, class TransformType>
void runKernelTransform(const mxArray** argIn,
	       mxArray** argOut) {

  // get size of input arguments
  mwSize Mx = mxGetM(argIn[1]); // number of source points
  mwSize Mxi = mxGetM(argIn[3]); // number of points to be warped
  mwSize ndimxi; // number of dimension of points to be warped
  const mwSize *dimsxi; // dimensions vector of array of points to be warped

  // create pointers to input matrices
  TScalarType *x = (TScalarType *)mxGetPr(argIn[1]); // source points
  TScalarType *y = (TScalarType *)mxGetPr(argIn[2]); // target points
  TScalarType *xi = (TScalarType *)mxGetPr(argIn[3]); // points to be warped

  // duplicate the input x and y matrices to PointSet format so that
  // we can pass it to the ITK function

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

  mwSize pointId=0;
  for (mwSize row=0; row < Mx; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      fixedPoint[col] = y[Mx * col + row];
      movingPoint[col] = x[Mx * col + row];
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
  TScalarType *yi = (TScalarType *)mxGetPr(argOut[0]);

  // transform points
  for (mwSize row=0; row < Mxi; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      toWarpPoint[col] = xi[Mxi * col + row];
    }
    warpedPoint = transform->TransformPoint(toWarpPoint);
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      yi[Mxi * col + row] = warpedPoint[col];
    }
  }
  
  // exit function
  return;
  
}

// parseTransformType<TScalarType, Dimension>()
template <class TScalarType, unsigned int Dimension>
void parseTransformType(const mxArray** argIn,
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
		       ElasticTransformType>(argIn, argOut);
  } else if (!strcmp(transform, "elasticr")) {
    runKernelTransform<TScalarType, Dimension, 
		       ElasticReciprocalTransformType>(argIn, argOut);
  } else if (!strcmp(transform, "tps")) {
    runKernelTransform<TScalarType, Dimension, 
		       TpsTransformType>(argIn, argOut);
  } else if (!strcmp(transform, "tpsr2")) {
    runKernelTransform<TScalarType, Dimension, 
		       TpsR2LogRTransformType>(argIn, argOut);
  } else if (!strcmp(transform, "volume")) {
    runKernelTransform<TScalarType, Dimension, 
		       VolumeTransformType>(argIn, argOut);
  } else if (!strcmp(transform, "bspline")) {
    mexErrMsgTxt("BSpline transform not implemented yet");
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
void parseDimensionToTemplate(const mxArray** argIn,
			      mxArray** argOut) {

  // dimension of points
  mwSize Nx = mxGetN(argIn[1]);

  // parse the dimension value
  switch (Nx) {
  case 2:
    parseTransformType<TScalarType, 2>(argIn, argOut);
    break;
  case 3:
    parseTransformType<TScalarType, 3>(argIn, argOut);
    break;
  default:
    mexErrMsgTxt("Input points can only have dimensions 2 or 3");
    break;
  }

  // exit function
  return;

}

// parseInputTypeToTemplate()
void parseInputTypeToTemplate(const mxArray** argIn,
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
    parseDimensionToTemplate<double>(argIn, argOut);
    break;
  case mxSINGLE_CLASS:
    parseDimensionToTemplate<float>(argIn, argOut);
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
  if (nrhs != 4) {
    mexErrMsgTxt("Four input arguments required");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  // check size of input arguments
  mwSize Mx = mxGetM(prhs[1]); // number of source points
  mwSize My = mxGetM(prhs[2]); // number of target points
  mwSize Dimension = mxGetN(prhs[1]); // dimension of source points
  mwSize dimy = mxGetN(prhs[2]); // dimension of target points
  mwSize dimxi = mxGetN(prhs[3]); // dimension of points to be warped
  
  if (Mx != My) mexErrMsgTxt("X and Y must have the same number of points (i.e. rows).");
  if (Dimension != 3 || Dimension != dimy || Dimension != dimxi) {
    mexErrMsgTxt("X, Y and XI must all have the same dimension (i.e. number of columns).");
  }

  // run filter (this function starts a cascade of functions designed
  // to translate the run-time type variables like inputVoxelClassId
  // to templates, so that we don't need to nest lots of "switch" or
  // "if" statements)
  parseInputTypeToTemplate(prhs, plhs);

  // exit successfully
  return;

}

#endif /* ITKPSTRANSFORM_CPP */
