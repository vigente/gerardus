/*
 * MexTemplateImageFilter.hpp
 *
 * Code that is specific to itk::TemplateImageFilter.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.1.0
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

#ifndef MEXTEMPLATEIMAGEFILTER_HPP
#define MEXTEMPLATEIMAGEFILTER_HPP

/* mex headers */
#include <mex.h>

/* ITK headers */
#include "itkImage.h"
#include "itkTemplateImageFilter.h"

/* Gerardus headers */
#include "MexBaseFilter.hpp"

/* 
 * MexTemplateImageFilter : MexBaseFilter
 */
template <class InVoxelType, class OutVoxelType>
class MexTemplateImageFilter : 
  public MexBaseFilter<InVoxelType, OutVoxelType> {

private:
  
  typedef itk::TemplateImageFilter< 
    itk::Image<InVoxelType, Dimension>,
    itk::Image<OutVoxelType, Dimension> > FilterType;

protected:

public:

  // constructor for filters that take user-defined parameters
  MexTemplateImageFilter(const NrrdImage &_nrrd, int _nargout, mxArray** _argOut,
			     const int _nargin, const mxArray** _argIn) :
    MexBaseFilter<InVoxelType, OutVoxelType>(_nrrd, _nargout, _argOut,
					     _nargin, _argIn) {

  // constructor for filters without parameters

  // MexTemplateImageFilter(const NrrdImage &_nrrd, int _nargout, mxArray** _argOut) :
  //   MexBaseFilter<InVoxelType, OutVoxelType>(_nrrd, _nargout, _argOut) {

    // instantiate filter
    this->filter = FilterType::New();
  }

  // if this particular filter needs to redifine one or more BaseFilter
  // virtual methods, the corresponding declarations go here

  // void FilterSetup();
  // void RunFilter();
  // void CopyAllFilterOutputsToMatlab();

};

// input/output voxel type is never going to be a string, so we are
// going to use this specialization just to have a container to put
// the text strings with the two names the user can type to select
// this filter
template <>
class MexTemplateImageFilter< std::string, std::string > {
public:
  
  static const std::string longname;
  static const std::string shortname;
  
};

/*
 * Filter exclusions: input/output data type combinations that are not
 * allowed or not going to be used for this filter
 */

// #define EXCLUDEFILTER(T1, T2)						\
//   template <>								\
//   class MexTemplateImageFilter< T1, T2 > :				\
//     public MexBaseFilter<T1, T2> {					\
//   public:								\
//     MexTemplateImageFilter(const NrrdImage &, int, mxArray**,		\
// 		       const int, const mxArray **) {;}			\
//     MexTemplateImageFilter(const NrrdImage &, int, mxArray**) {;}	\
//     void CopyMatlabInputToItkImage() {;}				\
//     void FilterSetup() {;}						\
//     void RunFilter() {;}						\
//     void CopyAllFilterOutputsToMatlab() {;}				\
//     void CopyFilterImageOutputToMatlab() {;}				\
//   };

// EXCLUDEFILTER(mxLogical, mxLogical)
// EXCLUDEFILTER(mxLogical, uint8_T)
// EXCLUDEFILTER(mxLogical, int8_T)
// EXCLUDEFILTER(mxLogical, uint16_T)
// EXCLUDEFILTER(mxLogical, int16_T)
// EXCLUDEFILTER(mxLogical, int32_T)
// EXCLUDEFILTER(mxLogical, int64_T)
// EXCLUDEFILTER(mxLogical, float)
// EXCLUDEFILTER(mxLogical, double)

// EXCLUDEFILTER(uint8_T, mxLogical)
// EXCLUDEFILTER(uint8_T, uint8_T)
// EXCLUDEFILTER(uint8_T, int8_T)
// EXCLUDEFILTER(uint8_T, uint16_T)
// EXCLUDEFILTER(uint8_T, int16_T)
// EXCLUDEFILTER(uint8_T, int32_T)
// EXCLUDEFILTER(uint8_T, int64_T)
// EXCLUDEFILTER(uint8_T, float)
// EXCLUDEFILTER(uint8_T, double)

// EXCLUDEFILTER(int8_T, mxLogical)
// EXCLUDEFILTER(int8_T, uint8_T)
// EXCLUDEFILTER(int8_T, int8_T)
// EXCLUDEFILTER(int8_T, uint16_T)
// EXCLUDEFILTER(int8_T, int16_T)
// EXCLUDEFILTER(int8_T, int32_T)
// EXCLUDEFILTER(int8_T, int64_T)
// EXCLUDEFILTER(int8_T, float)
// EXCLUDEFILTER(int8_T, double)

// EXCLUDEFILTER(uint16_T, mxLogical)
// EXCLUDEFILTER(uint16_T, uint8_T)
// EXCLUDEFILTER(uint16_T, int8_T)
// EXCLUDEFILTER(uint16_T, uint16_T)
// EXCLUDEFILTER(uint16_T, int16_T)
// EXCLUDEFILTER(uint16_T, int32_T)
// EXCLUDEFILTER(uint16_T, int64_T)
// EXCLUDEFILTER(uint16_T, float)
// EXCLUDEFILTER(uint16_T, double)

// EXCLUDEFILTER(int16_T, mxLogical)
// EXCLUDEFILTER(int16_T, uint8_T)
// EXCLUDEFILTER(int16_T, int8_T)
// EXCLUDEFILTER(int16_T, uint16_T)
// EXCLUDEFILTER(int16_T, int16_T)
// EXCLUDEFILTER(int16_T, int32_T)
// EXCLUDEFILTER(int16_T, int64_T)
// EXCLUDEFILTER(int16_T, float)
// EXCLUDEFILTER(int16_T, double)

// EXCLUDEFILTER(int32_T, mxLogical)
// EXCLUDEFILTER(int32_T, uint8_T)
// EXCLUDEFILTER(int32_T, int8_T)
// EXCLUDEFILTER(int32_T, uint16_T)
// EXCLUDEFILTER(int32_T, int16_T)
// EXCLUDEFILTER(int32_T, int32_T)
// EXCLUDEFILTER(int32_T, int64_T)
// EXCLUDEFILTER(int32_T, float)
// EXCLUDEFILTER(int32_T, double)

// EXCLUDEFILTER(int64_T, mxLogical)
// EXCLUDEFILTER(int64_T, uint8_T)
// EXCLUDEFILTER(int64_T, int8_T)
// EXCLUDEFILTER(int64_T, uint16_T)
// EXCLUDEFILTER(int64_T, int16_T)
// EXCLUDEFILTER(int64_T, int32_T)
// EXCLUDEFILTER(int64_T, int64_T)
// EXCLUDEFILTER(int64_T, float)
// EXCLUDEFILTER(int64_T, double)

// EXCLUDEFILTER(float, mxLogical)
// EXCLUDEFILTER(float, uint8_T)
// EXCLUDEFILTER(float, int8_T)
// EXCLUDEFILTER(float, uint16_T)
// EXCLUDEFILTER(float, int16_T)
// EXCLUDEFILTER(float, int32_T)
// EXCLUDEFILTER(float, int64_T)
// EXCLUDEFILTER(float, float)
// EXCLUDEFILTER(float, double)

// EXCLUDEFILTER(double, mxLogical)
// EXCLUDEFILTER(double, uint8_T)
// EXCLUDEFILTER(double, int8_T)
// EXCLUDEFILTER(double, uint16_T)
// EXCLUDEFILTER(double, int16_T)
// EXCLUDEFILTER(double, int32_T)
// EXCLUDEFILTER(double, int64_T)
// EXCLUDEFILTER(double, float)
// EXCLUDEFILTER(double, double)

// #undef EXCLUDEFILTER

#endif /* MEXTEMPLATEIMAGEFILTER_HPP */
