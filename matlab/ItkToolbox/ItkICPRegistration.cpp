/**
 * ItkICPRegistration.cpp
 *
 * ITK_ICP_REGISTRATION  Iterative Closest Point registration
 *
 * This function is a derived work of IterativeClosestPoint3.cxx
 * https://github.com/Kitware/ITK/blob/master/Examples/Registration/IterativeClosestPoint3.cxx
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2013 University of Oxford
  * Version: 0.0.5
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

#ifndef ITKICPREGISTRATION
#define ITKICPREGISTRATION

#define DEBUG

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>

/* Gerardus headers */
#include "MatlabImportFilter.h"
#include "MatlabExportFilter.h"

/* ITK headers */
#include "itkTranslationTransform.h"
#include "itkEuclideanDistancePointMetric.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkPointSet.h"
#include "itkPointSetToPointSetRegistrationMethod.h"
//#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkPointSetToImageFilter.h"

// type definitions
static const unsigned int Dimension = 3;
typedef double CoordinateType;
typedef itk::PointSet<CoordinateType, Dimension> PointSetType;
typedef PointSetType::PointType PointType;

/*
 * mexFunction(): entry point for the mex function
 */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // interface to deal with input arguments from Matlab
  enum InputIndexType {IN_X, IN_Y, IN_TRANSFORM, 
		       IN_NITER, IN_GRADTOL, IN_VALTOL, IN_EPSFUN, InputIndexType_MAX};
  MatlabImportFilter::Pointer matlabImport = MatlabImportFilter::New();
  matlabImport->ConnectToMatlabFunctionInput(nrhs, prhs);

  // check the number of input arguments
  matlabImport->CheckNumberOfArguments(2, InputIndexType_MAX);

  // register the inputs for this function at the import filter
  typedef MatlabImportFilter::MatlabInputPointer MatlabInputPointer;
  MatlabInputPointer inX = matlabImport->RegisterInput(IN_X, "X");
  MatlabInputPointer inY = matlabImport->RegisterInput(IN_Y, "Y");
  MatlabInputPointer inTRANSFORM = matlabImport->RegisterInput(IN_TRANSFORM, "TRANSFORM");
  MatlabInputPointer inNITER = matlabImport->RegisterInput(IN_NITER, "NITER");
  MatlabInputPointer inGRADTOL = matlabImport->RegisterInput(IN_GRADTOL, "GRADTOL");
  MatlabInputPointer inVALTOL = matlabImport->RegisterInput(IN_VALTOL, "VALTOL");
  MatlabInputPointer inEPSFUN = matlabImport->RegisterInput(IN_EPSFUN, "EPSFUN");


  // interface to deal with outputs to Matlab
  enum OutputIndexType {OUT_YY, OUT_T, OutputIndexType_MAX};
  MatlabExportFilter::Pointer matlabExport = MatlabExportFilter::New();
  matlabExport->ConnectToMatlabFunctionOutput(nlhs, plhs);
  
  // check that the number of outputs the user is asking for is valid
  matlabExport->CheckNumberOfArguments(0, OutputIndexType_MAX);

  // register the outputs for this function at the export filter
  typedef MatlabExportFilter::MatlabOutputPointer MatlabOutputPointer;
  MatlabOutputPointer outYY = matlabExport->RegisterOutput(OUT_YY, "Y2");
  MatlabOutputPointer outT  = matlabExport->RegisterOutput(OUT_T, "T");

  // if any input point set is empty, the outputs are empty too
  if (mxIsEmpty(prhs[IN_X]) || mxIsEmpty(prhs[IN_Y])) {
    matlabExport->CopyEmptyArrayToMatlab(outYY);
    matlabExport->CopyEmptyArrayToMatlab(outT);
    return;
  }

  // get size of input matrix with the points
  mwSize nrowsX = mxGetM(prhs[IN_X]);
  mwSize ncolsX = mxGetN(prhs[IN_X]);
  mwSize nrowsY = mxGetM(prhs[IN_Y]);
  mwSize ncolsY = mxGetN(prhs[IN_Y]);
  if (ncolsX != Dimension || ncolsY != Dimension) {
    mexErrMsgTxt("X and Y must have 3 columns");
  }

  // if there's some problem reading the point, default is NaN
  PointType def;
  def.Fill(mxGetNaN());

  // read point sets
  PointSetType::Pointer fixedPointSet = PointSetType::New();
  // fixedPointSet->GetPoints()->CastToSTLContainer().resize(nrowsX);
  PointSetType::Pointer xDef = PointSetType::New();
  // fixedPointSet->GetPoints()->CastToSTLContainer()
  //   = matlabImport->ReadVectorOfVectorsFromMatlab<PointType::ValueType, PointType>
  //   (IN_X, "X", xDef->CastToSTLContainer());
  //@@
  // matlabImport->ReadMatrixFromMatlabIntoVectorOfVectors<PointType::ValueType, std::vector<PointType> >
  //   (IN_X, "X", 
  //    fixedPointSet->GetPoints()->CastToSTLContainer(), fixedPointSet->GetNumberOfPoints(), Dimension);

  PointSetType::Pointer movingPointSet = PointSetType::New();
  movingPointSet->GetPoints()->CastToSTLContainer().resize(nrowsY);
  //@@
  // matlabImport->ReadMatrixFromMatlabIntoVectorOfVectors<PointType::ValueType, std::vector<PointType> >
  //   (IN_Y, "Y", 
  //    movingPointSet->GetPoints()->CastToSTLContainer(), movingPointSet->GetNumberOfPoints(), Dimension);

#ifdef DEBUG
  // debug
  std::cout << "Number of points in X = " << fixedPointSet->GetNumberOfPoints() << std::endl;
  std::cout << "Number of points in Y = " << movingPointSet->GetNumberOfPoints() << std::endl;
#endif
  
  /*
   * registration method
   */

  typedef itk::PointSetToPointSetRegistrationMethod<PointSetType, 
						    PointSetType> RegistrationType;

  RegistrationType::Pointer registration = RegistrationType::New();

  registration->SetFixedPointSet(fixedPointSet);
  registration->SetMovingPointSet(movingPointSet);

  /*
   * metric
   */

  // generic metric to use for common methods
  typedef itk::PointSetToPointSetMetric<PointSetType, PointSetType> GenericPointSetToPointSetMetricType;
  GenericPointSetToPointSetMetricType::Pointer metric;

  // specific metrics that can be chosen by the user
  static enum MetricEnum {
    METRIC_Unknown,
    METRIC_EuclideanDistancePoint
  } metricLabel = METRIC_Unknown;

  typedef itk::EuclideanDistancePointMetric<PointSetType, PointSetType> 
    EuclideanDistancePointMetricType;

  EuclideanDistancePointMetricType::Pointer metricEuclideanDistancePoint;

  // read metric chosen by the user
  metricLabel = METRIC_EuclideanDistancePoint;

  switch (metricLabel) {
  case METRIC_EuclideanDistancePoint:
    metricEuclideanDistancePoint = EuclideanDistancePointMetricType::New();
    registration->SetMetric(metricEuclideanDistancePoint);
    metric = dynamic_cast<GenericPointSetToPointSetMetricType*>(metricEuclideanDistancePoint.GetPointer());
    break;
  default:
    mexErrMsgTxt("Metric not implemented");
    break;
  }

  /*
   * transform
   */

  // generic transform to use for common methods
  typedef itk::Transform<CoordinateType, Dimension> GenericTransformType;
  GenericTransformType::Pointer transform;

  // specific transforms that can be chosen by the user
  static enum TransformEnum {
    TRANSFORM_Unknown,
    TRANSFORM_Translation
  } transformLabel = TRANSFORM_Unknown;

  typedef itk::TranslationTransform<CoordinateType, Dimension> TranslationTransformType;

  TranslationTransformType::Pointer transformTranslation;

  // read transform chosen by the user
  transformLabel = TRANSFORM_Translation;

  switch (transformLabel) {
  case TRANSFORM_Translation:
    transformTranslation = TranslationTransformType::New();
    transformTranslation->SetIdentity();
    registration->SetTransform(transformTranslation);
    registration->SetInitialTransformParameters(transformTranslation->GetParameters());
    transform = dynamic_cast<GenericTransformType *>(transformTranslation.GetPointer());
    break;
  default:
    mexErrMsgTxt("Transform not implemented");
    break;
  }

  /*
   * optimizer
   */

  // generic optimizer to use for common methods
  typedef itk::Optimizer GenericOptimizerType;
  GenericOptimizerType::Pointer optimizer;
  GenericOptimizerType::ScalesType scales(transform->GetNumberOfParameters());

  // specific optimizers that can be chosen by the user
  static enum OptimizerEnum {
    OPTIMIZER_Unknown,
    OPTIMIZER_LevenbergMarquardt
  } optimizerLabel = OPTIMIZER_Unknown;

  typedef itk::LevenbergMarquardtOptimizer LevenbergMarquardtOptimizerType;

  LevenbergMarquardtOptimizerType::Pointer optimizerLevenbergMarquardt;

  // read optimizer chosen by the user
  optimizerLabel = OPTIMIZER_LevenbergMarquardt;

  switch (optimizerLabel) {
  case OPTIMIZER_LevenbergMarquardt:
    // register the fields from the optimizer struct
    //    matlabImport->RegisterInputArgumentFromMatlab(mxGetField());


    optimizerLevenbergMarquardt = LevenbergMarquardtOptimizerType::New();
    optimizerLevenbergMarquardt->SetUseCostFunctionGradient(false);



    // optimizerLevenbergMarquardt->SetNumberOfIterations
    //   (matlabImport->ReadScalarFromMatlab<double>()
    //    );


    optimizerLevenbergMarquardt->SetValueTolerance(1e-5);
    optimizerLevenbergMarquardt->SetGradientTolerance(1e-5);
    optimizerLevenbergMarquardt->SetEpsilonFunction(1e-6);
    scales.Fill(1.0);
    optimizerLevenbergMarquardt->SetScales(scales);
    registration->SetOptimizer(optimizerLevenbergMarquardt);
    optimizer = dynamic_cast<GenericOptimizerType *>(optimizerLevenbergMarquardt.GetPointer());
    break;
  default:
    mexErrMsgTxt("Optimizer not implemented");
    break;
  }

  /*
   * run registration
   */

  try 
    {
      registration->Update();
    }
  catch( itk::ExceptionObject & e )
    {
      mexErrMsgTxt(e.GetDescription());
    }

  /*
   * export results
   */
  
  std::cout << "Solution = " << transform->GetParameters() << std::endl;
  std::cout << "Stop condition = " << optimizer->GetStopConditionDescription() << std::endl;

  // warp the moving points according to the solution
  PointSetType::Pointer warpedMovingPointSet = movingPointSet;
  for (mwIndex i = 0; i < movingPointSet->GetNumberOfPoints(); ++i) {
    warpedMovingPointSet->SetPoint(i, transform->TransformPoint(movingPointSet->GetPoint(i)));
  }

  // copy to Matlab the moving points after warping
  matlabExport->CopyVectorOfVectorsToMatlab<CoordinateType, 
  					    const std::vector<PointType> 
					    >
    (outYY,
     warpedMovingPointSet->GetPoints()->CastToSTLConstContainer(), 
     warpedMovingPointSet->GetNumberOfPoints(), Dimension);

  // registration parameters
  matlabExport->CopyVectorOfScalarsToMatlab<CoordinateType, 
					    GenericTransformType::ParametersType>
    (outT, 
     registration->GetTransform()->GetParameters(), 
     transform->GetNumberOfParameters());



}

#endif /* ITKICPREGISTRATION */
