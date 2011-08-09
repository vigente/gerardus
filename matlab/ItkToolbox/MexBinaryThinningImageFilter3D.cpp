/*
 * MexBinaryThinningImageFilter3D.cpp
 *
 * Code that is specific to itk::BinaryThinningImageFilter3D
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.1.4
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

#ifndef MEXBINARYTHINNINGIMAGEFILTER3D_CPP
#define MEXBINARYTHINNINGIMAGEFILTER3D_CPP

#include "MexBinaryThinningImageFilter3D.hpp"

/*
 * strings that the user can use to invoke this filter in itk_imfilter()
 */
const std::string MexBinaryThinningImageFilter3D<std::string, 
		  std::string>::longname = "BinaryThinningImageFilter3D";
const std::string MexBinaryThinningImageFilter3D<std::string, 
		  std::string>::shortname = "skel";

/* 
 * MexBinaryThinningImageFilter3D : MexBaseFilter
 */

/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#define FILTERINST(T1, T2)					\
  template class MexBinaryThinningImageFilter3D<T1, T2>;

FILTERINST(uint8_T, uint8_T)
FILTERINST(int8_T, int8_T)
FILTERINST(uint16_T, uint16_T)
FILTERINST(int16_T, int16_T)
FILTERINST(int32_T, int32_T)
FILTERINST(int64_T, int64_T)
FILTERINST(float, float)
FILTERINST(double, double)

#undef FILTERINST

#endif /* MEXBINARYTHINNINGIMAGEFILTER3D_CPP */
