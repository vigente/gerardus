/*
 * MexMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.hpp
 *
 * Code that is specific to itk::MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012 University of Oxford
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

#ifndef MEXMULTISCALEHESSIANSMOOTHED3DTOVESSELNESSMEASUREIMAGEFILTER_HPP
#define MEXMULTISCALEHESSIANSMOOTHED3DTOVESSELNESSMEASUREIMAGEFILTER_HPP

/* mex headers */
#include <mex.h>

/* ITK headers */
#include "itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.h"

/* Gerardus headers */
#include "MexBaseFilter.hpp"

/* 
 * MexMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter : MexBaseFilter
 */
template <class InVoxelType, class OutVoxelType>
class MexMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter : 
  public MexBaseFilter<InVoxelType, OutVoxelType> {

private:
  
  typedef itk::MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter< 
    itk::Image<InVoxelType, Dimension>,
    itk::Image<OutVoxelType, Dimension> > DerivedImageToImageFilterType;

  // this->filter is the base filter, and doesn't provide access to
  // methods that are exclusive to this derived filter. derivedFilter
  // is a pointer that provides said access
  typename DerivedImageToImageFilterType::Pointer derivedFilter;

protected:

  // user-provided input parameters
  double sigmaMin;
  double sigmaMax;
  int    numSigmaSteps;
  bool   isSigmaStepLog;

public:

  // constructor (to instantiate the filter and process the
  // user-provided input parameters, if any)
  MexMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter(const NrrdImage &_nrrd, 
			 int _nargout, mxArray** _argOut,
			 const int _nargin, const mxArray** _argIn);

  // if this particular filter needs to redefine one or more MexBaseFilter
  // virtual methods, the corresponding declarations go here
  // void CheckNumberOfOutputs();
  void FilterAdvancedSetup();
  // void ExportOtherFilterOutputsToMatlab();

};

// input/output voxel type is never going to be a string, so we are
// going to use this specialization just to have a container to put
// the text strings with the two names the user can type to select
// this filter
template <>
class MexMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter<std::string, std::string> {
public:
  
  static const std::string longname;
  static const std::string shortname;
  
};

#endif /* MEXMULTISCALEHESSIANSMOOTHED3DTOVESSELNESSMEASUREIMAGEFILTER_HPP */