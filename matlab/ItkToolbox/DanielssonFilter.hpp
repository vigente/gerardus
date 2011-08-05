/*
 * DanielssonFilter.hpp
 *
 * Code that is specific to the DanielssonDistanceMapImageFilter
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.3.1
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

#ifndef DANIELSSONFILTER_HPP
#define DANIELSSONFILTER_HPP

/* mex headers */
#include <mex.h>

/* ITK headers */
#include "itkImage.h"
#include "itkDanielssonDistanceMapImageFilter.h"

/* Gerardus headers */
#include "BaseFilter.hpp"

/* 
 * DanielssonFilter : BaseFilter
 */
template <class InVoxelType, class OutVoxelType>
class DanielssonFilter : 
  public BaseFilter<InVoxelType, OutVoxelType> {
private:
  typedef itk::DanielssonDistanceMapImageFilter< 
  itk::Image<InVoxelType, Dimension>,
  itk::Image<OutVoxelType, Dimension> > FilterType;

protected:

public:
  DanielssonFilter(const NrrdImage &nrrd, int _nargout, mxArray** argOut) :
    BaseFilter<InVoxelType, OutVoxelType>(nrrd, _nargout, argOut) {
    // instantiate filter
    this->filter = FilterType::New();
  }
  void CopyAllFilterOutputsToMatlab();
  void CopyFilterNearestOutputToMatlab();
};

/*
 * Filter exclusions: input/output data type combinations that are not
 * allowed for this filter
 */

#define EXCLUDEFILTER(T1, T2)						\
  template <>								\
  class DanielssonFilter< T1, T2 > :					\
    public BaseFilter<T1, T2> {						\
  public:								\
    DanielssonFilter(const NrrdImage &, int, mxArray**) {;}		\
    void CopyMatlabInputsToFilter() {;}					\
    void FilterSetup() {;}						\
    void RunFilter() {;}						\
    void CopyAllFilterOutputsToMatlab() {;}				\
  };

EXCLUDEFILTER(bool, int8_T);
EXCLUDEFILTER(bool, int16_T)
EXCLUDEFILTER(bool, int32_T)
EXCLUDEFILTER(bool, int64_T)

EXCLUDEFILTER(uint8_T, int8_T)
EXCLUDEFILTER(uint8_T, int16_T)
EXCLUDEFILTER(uint8_T, int32_T)
EXCLUDEFILTER(uint8_T, int64_T)

EXCLUDEFILTER(int8_T, int8_T)
EXCLUDEFILTER(int8_T, int16_T)
EXCLUDEFILTER(int8_T, int32_T)
EXCLUDEFILTER(int8_T, int64_T)

EXCLUDEFILTER(uint16_T, int8_T)
EXCLUDEFILTER(uint16_T, int16_T)
EXCLUDEFILTER(uint16_T, int32_T)
EXCLUDEFILTER(uint16_T, int64_T)

EXCLUDEFILTER(int16_T, int8_T)
EXCLUDEFILTER(int16_T, int16_T)
EXCLUDEFILTER(int16_T, int32_T)
EXCLUDEFILTER(int16_T, int64_T)

EXCLUDEFILTER(int32_T, int8_T)
EXCLUDEFILTER(int32_T, int16_T)
EXCLUDEFILTER(int32_T, int32_T)
EXCLUDEFILTER(int32_T, int64_T)

EXCLUDEFILTER(int64_T, int8_T)
EXCLUDEFILTER(int64_T, int16_T)
EXCLUDEFILTER(int64_T, int32_T)
EXCLUDEFILTER(int64_T, int64_T)

EXCLUDEFILTER(float, int8_T)
EXCLUDEFILTER(float, int16_T)
EXCLUDEFILTER(float, int32_T)
EXCLUDEFILTER(float, int64_T)

EXCLUDEFILTER(double, int8_T)
EXCLUDEFILTER(double, int16_T)
EXCLUDEFILTER(double, int32_T)
EXCLUDEFILTER(double, int64_T)

#undef EXCLUDEFILTER

#endif /* DANIELSSONFILTER_HPP */
