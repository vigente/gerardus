/*
 * DanielssonFilter.cpp
 *
 * Code that is specific to the DanielssonDistanceMapImageFilter
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

#ifndef DANIELSSONFILTER_CPP
#define DANIELSSONFILTER_CPP

#import "DanielssonFilter.hpp"

/*
 * Instantiate filter with all the input/output combinations that it
 * accepts. This is necessary for the linker. The alternative is to
 * have all the code in the header files, but this makes compilation
 * slower and maybe the executable larger
 */

#define FILTERINST(T1, T2)						\
  template class DanielssonFilter<T1, T2>;				\

FILTERINST(bool, bool);
FILTERINST(bool, uint8_T)
FILTERINST(bool, uint16_T)
FILTERINST(bool, float)
FILTERINST(bool, double)

FILTERINST(uint8_T, bool);
FILTERINST(uint8_T, uint8_T)
FILTERINST(uint8_T, uint16_T)
FILTERINST(uint8_T, float)
FILTERINST(uint8_T, double)

FILTERINST(int8_T, bool);
FILTERINST(int8_T, uint8_T)
FILTERINST(int8_T, uint16_T)
FILTERINST(int8_T, float)
FILTERINST(int8_T, double)

FILTERINST(uint16_T, bool);
FILTERINST(uint16_T, uint8_T)
FILTERINST(uint16_T, uint16_T)
FILTERINST(uint16_T, float)
FILTERINST(uint16_T, double)

FILTERINST(int16_T, bool);
FILTERINST(int16_T, uint8_T)
FILTERINST(int16_T, uint16_T)
FILTERINST(int16_T, float)
FILTERINST(int16_T, double)

FILTERINST(int32_T, bool);
FILTERINST(int32_T, uint8_T)
FILTERINST(int32_T, uint16_T)
FILTERINST(int32_T, float)
FILTERINST(int32_T, double)

FILTERINST(int64_T, bool);
FILTERINST(int64_T, uint8_T)
FILTERINST(int64_T, uint16_T)
FILTERINST(int64_T, float)
FILTERINST(int64_T, double)

FILTERINST(float, bool);
FILTERINST(float, uint8_T)
FILTERINST(float, uint16_T)
FILTERINST(float, float)
FILTERINST(float, double)

FILTERINST(double, bool);
FILTERINST(double, uint8_T)
FILTERINST(double, uint16_T)
FILTERINST(double, float)
FILTERINST(double, double)

#undef FILTERINST

#endif /* DANIELSSONFILTER_CPP */
