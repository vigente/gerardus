/*
 * ArgumentParsers.hpp
 *
 * parseInputTypeToTemplate()
 * parseOutputTypeToTemplate< InVoxelType >()
 * parseFilterTypeToTemplate< InVoxelType, OutVoxelType >()
 *
 * These functions are used to be able to map between the input/output
 * data types that are only know at run-time, and the input/output
 * data templates that ITK requires and must be know at compilation
 * time.
 *
 * To avoid a nesting nightmare:
 *
 * switch FilterType {
 *   switch InputDataType {
 *     switch OutputDataType {
 *     }
 *   }
 * }
 *
 * we split the conversion in 3 steps. The first function,
 * parseInputTypeToTemplate() reads the input data type and the
 * filter type, and instantiates
 * parseOutputTypeToTemplate<InVoxelType>() for each input data
 * type.
 *
 * This way, we have effectively mapped the input type from a run-time
 * variable to a set of compilation time templates.
 *
 * Then parseOutputTypeToTemplate<InVoxelType>() decides on the
 * output data type depending on the filter and the input, and
 * instantiates one
 * parseFilterTypeToTemplate<InVoxelType, OutVoxelType>() for each
 * output data type.
 *
 * This way, we have instantiated all combinations of input/output
 * data types, without having to code them explicitly.
 *
 * Finally, parseFilterTypeToTemplate<InVoxelType, OutVoxelType>()
 * checks that the requested filter is implemented and maps the filter
 * variable to all filter templates.
 *
 * Note that the actual filtering is done in the constructor of
 * BaseFilter, that is instantiated from
 * parseFilterTypeToTemplate.
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

#ifndef ARGUMENTPARSERS_HPP
#define ARGUMENTPARSERS_HPP

/* mex headers */
#include <mex.h>

/* Gerardus headers */
#import "NrrdImage.hpp"

// parseFilterTypeToTemplate<InVoxelType, OutVoxelType>()
template <class InVoxelType, class OutVoxelType>
void parseFilterTypeToTemplate(char *filterName,
			       NrrdImage nrrd,
			       int nargout,
			       mxArray** argOut);

// parseOutputTypeToTemplate<InVoxelType>()
template <class InVoxelType>
void parseOutputTypeToTemplate(char *filter,
			       NrrdImage nrrd,
			       int nargout,
			       mxArray** argOut);
// parseInputTypeToTemplate()
void parseInputTypeToTemplate(mxClassID inputVoxelClassId, 
			      char *filter,
			      NrrdImage nrrd,
			      int nargout,
			      mxArray** argOut);

#endif /* ARGUMENTPARSERS_HPP */
