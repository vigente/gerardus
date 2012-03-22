/*
 * MexBinaryErodeImageFilter.hpp
 *
 * Code that is specific to itk::BinaryErodeImageFilter. Support for
 * radius and foreground value arguments. Structuring element is a
 * ball.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.2.0
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

#ifndef MEXBINARYERODEIMAGEFILTER_HPP
#define MEXBINARYERODEIMAGEFILTER_HPP

/* mex headers */
#include <mex.h>

/* ITK headers */
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

/* Gerardus headers */
#include "MexBaseFilter.hpp"

/* 
 * MexBinaryErodeImageFilter : MexBaseFilter
 */
template <class InVoxelType, class OutVoxelType>
class MexBinaryErodeImageFilter : 
  public MexBaseFilter<InVoxelType, OutVoxelType> {

private:
  
  typedef itk::BinaryBallStructuringElement<InVoxelType, Dimension >
  StructuringElementType;
  
  typedef itk::BinaryErodeImageFilter< 
    itk::Image<InVoxelType, Dimension>,
    itk::Image<OutVoxelType, Dimension>, StructuringElementType > FilterType;

protected:

  // user-provided input arguments
  unsigned long radius;   // (comp) radius of the ball in voxels
  InVoxelType foreground; // (opt) voxels with this value will be
                          // dilated. Default, maximum value of the
                          // pixel type

public:

  // constructor for filters that take user-defined parameters
  MexBinaryErodeImageFilter(const NrrdImage &_nrrd, 
			    int _nargout, mxArray** _argOut,
			    const int _nargin, const mxArray** _argIn);

  // if this particular filter needs to redefine one or more BaseFilter
  // virtual methods, the corresponding declarations go here
  void FilterAdvancedSetup();

};

// input/output voxel type is never going to be a string, so we are
// going to use this specialization just to have a container to put
// the text strings with the two names the user can type to select
// this filter
template <>
class MexBinaryErodeImageFilter< std::string, std::string > {
public:
  
  static const std::string longname;
  static const std::string shortname;
  
};

#endif /* MEXBINARYERODEIMAGEFILTER_HPP */
