/*
 * GerardusCommon.hpp
 *
 * Miscellaneous functions of general use.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.2.1
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

#ifndef GERARDUSCOMMON_HPP
#define GERARDUSCOMMON_HPP

/* mex headers */
#include <mex.h>

/* C++ headers */
//#import <vector>

/* ITK headers */
#include "itkOffset.h"

/* Global constants */
static const unsigned int Dimension = 3; // volume data dimension
                                         // (3D volume)

/* 
 * CAST2MWSIZE(): macro to cast to mwSize type. This definition is
 *                necessary for ITK v3.20.0 to avoid an error when
 *                trying to compile
 *                itk::FixedArray::operator[](unsigned __int64) for
 *                Windows 64 bit, but * maybe we can remove it when
 *                ITK v4.0.0 is released
 */
#ifdef _WIN64
#define CAST2MWSIZE(x) static_cast<unsigned long>(x)
#else
#define CAST2MWSIZE(x) static_cast<mwSize>(x)
#endif

/*
 * sub2ind(): function that converts r, c, s indices to linear indices
 *            in a 3D array (same as Matlab's function sub2ind(),
 *            although in Matlab indices start at 1, and in C++, they
 *            start at 0)
 *
 */
mwIndex sub2ind(mwSize R, mwSize C, mwSize S, std::vector<mwIndex> rcs);
mwIndex sub2ind(mwSize R, mwSize C, mwSize S, itk::Offset<Dimension> rcs);

/*
 * ind2sub(): function that converts linear indices in a 3D array to
 *            r, c, s indices (same as Matlab's function ind2sub(),
 *            although in Matlab indices start at 1, and in C++, they
 *            start at 0)
 *
 */
std::vector<mwIndex> ind2sub(mwSize R, mwSize C, mwSize S, mwIndex idx);
itk::Offset<Dimension> ind2sub_itkOffset(mwSize R, mwSize C, mwSize S, 
					 mwIndex idx);

#endif /* GERARDUSCOMMON_HPP */
