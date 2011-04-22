/*
 * ItkKernelTransform.cpp
 *
 * ITK_KERNEL_TRANSFORM  ITK warps between 3D point sets with a known
 * correspondence
 *
 * This MEX-function provides a Matlab interface to run ITK's Kernel
 * Transforms:
 *
 *   itk::ElasticBodySplineKernelTransform
 *   itk::ElasticBodyReciprocalSplineKernelTransform
 *   itk::ThinPlateSplineKernelTransform
 *   itk::ThinPlateR2LogRSplineKernelTransform
 *   itk::VolumeSplineKernelTransform
 *
 * YI = ITK_KERNEL_TRANSFORM(X, Y, XI, TYPE)
 *
 *   X, Y are 3-column matrices with N rows. Each row has the
 *   coordinates of a point. The warp is defined so that
 *   X(i,:)->Y(i,:).
 *
 *   XI is a 3-column matrix with M rows. Each row has the coordinates
 *   of a point to be warped.
 *
 *   YI has the same dimensions as XI. YI contains the coordinates of
 *   the warped points.
 *
 *   TYPE is a string that allows to select the type of warp (no defaults):
 *
 *   'elastic':  itk::ElasticBodySplineKernelTransform
 *   'elasticr': itk::ElasticBodyReciprocalSplineKernelTransform
 *   'tps':      itk::ThinPlateSplineKernelTransform
 *   'tpsr2':    itk::ThinPlateR2LogRSplineKernelTransform
 *   'volume':   itk::VolumeSplineKernelTransform
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
 * change to your own Matlab root path.
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
#include "itkTransformBase.h"
#include "itkElasticBodySplineKernelTransform.h"
#include "itkElasticBodyReciprocalSplineKernelTransform.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "itkVolumeSplineKernelTransform.h"
#include "itkPointSet.h"

#ifndef KERNELTRANSFORM_CPP
#define KERNELTRANSFORM_CPP

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // check number of input and output arguments
  if (nrhs != 4) {
    mexErrMsgTxt("Four input arguments required");
  }
  else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  // check size of input arguments
  mwSize Mx = mxGetM(prhs[0]); // number of source points
  mwSize My = mxGetM(prhs[1]); // number of target points
  mwSize Mxi = mxGetM(prhs[2]); // number of points to be warped
  mwSize dimx = mxGetN(prhs[0]); // dimension of source points
  mwSize dimy = mxGetN(prhs[1]); // dimension of target points
  mwSize dimxi = mxGetN(prhs[2]); // dimension of points to be warped
  
  if (Mx != My) mexErrMsgTxt("X and Y must have the same number of points (i.e. rows).");
  if (dimx != dimy || dimx != dimxi) mexErrMsgTxt("X, Y and XI must have all dimension 3 (i.e. 3 columns).");

  // create output vector and pointer to populate it
  plhs[0] = mxCreateDoubleMatrix(Mxi, dimxi, mxREAL);
  double *yi = mxGetPr(plhs[0]); // warped points

  // create pointers to input matrices
  double *x = mxGetPr(prhs[0]); // source points
  double *y = mxGetPr(prhs[1]); // target points
  double *xi = mxGetPr(prhs[2]); // points to be warped

  // read type of transform (TPS, elastic, etc.)
  char *type;
  type = mxArrayToString(prhs[3]);
  if (type == NULL) {
    mexErrMsgTxt("Cannot read transform TYPE string");
  }

  // duplicate the input x and y matrices to PointSet format so that
  // we can pass it to the ITK function
  const unsigned int Dimension = 3;
  typedef double TScalarType; // data type for scalars (e.g. point
                              // coordinates)
  typedef itk::KernelTransform< TScalarType, Dimension > 
    KernelTransformType;
  typedef itk::ElasticBodySplineKernelTransform< TScalarType, Dimension > 
    ElasticTransformType;
  typedef itk::ElasticBodyReciprocalSplineKernelTransform< TScalarType, Dimension > 
    ElasticReciprocalTransformType;
  typedef itk::ThinPlateSplineKernelTransform< TScalarType, Dimension > 
    TpsTransformType;
  typedef itk::ThinPlateR2LogRSplineKernelTransform< TScalarType, Dimension > 
    TpsR2LogRTransformType;
  typedef itk::VolumeSplineKernelTransform< TScalarType, Dimension > 
    VolumeTransformType;

  typedef KernelTransformType::PointSetType PointSetType;

  PointSetType::Pointer fixedPointSet = PointSetType::New();
  PointSetType::Pointer movingPointSet = PointSetType::New();
  PointSetType::Pointer toWarpPointSet = PointSetType::New();
  typedef PointSetType::PointsContainer PointsContainer;
  PointsContainer::Pointer fixedPointContainer = PointsContainer::New();
  PointsContainer::Pointer movingPointContainer = PointsContainer::New();
  typedef PointSetType::PointType PointType;
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

  // select warp function
  KernelTransformType::Pointer transform;
  if (!strcmp(type, "elastic")) {
    transform = ElasticTransformType::New();
  } else if (!strcmp(type, "elasticr")) {
    transform = ElasticReciprocalTransformType::New();
  } else if (!strcmp(type, "tps")) {
    transform = TpsTransformType::New();
  } else if (!strcmp(type, "tpsr2")) {
    transform = TpsR2LogRTransformType::New();
  } else if (!strcmp(type, "volume")) {
    transform = VolumeTransformType::New();
  } else if (!strcmp(type, "")) {
    std::cout << "Implemented transform types: elastic, elasticr, tps, tpsr2, volume" << std::endl;
    return;
  } else {
    mexErrMsgTxt("Transform TYPE not implemented");
  }
  transform->SetSourceLandmarks(movingPointSet);
  transform->SetTargetLandmarks(fixedPointSet);
  transform->ComputeWMatrix();

  // warp points
  for (mwSize row=0; row < Mxi; ++row) {
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      toWarpPoint[col] = xi[Mxi * col + row];
    }
    warpedPoint = transform->TransformPoint(toWarpPoint);
    for (mwSize col=0; col < (mwSize)Dimension; ++col) {
      yi[Mxi * col + row] = warpedPoint[col];
    }
  }

  return;

}

#endif /* KERNELTRANSFORM_CPP */
