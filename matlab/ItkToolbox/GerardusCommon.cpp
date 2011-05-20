/*
 * GerardusCommon.cpp
 *
 * Miscellaneous functions of general use.
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

#ifndef GERARDUSCOMMON_CPP
#define GERARDUSCOMMON_CPP

/* mex headers */
#include <mex.h>

/* C++ headers */
#import <vector>

/* ITK headers */
#include "itkOffset.h"

/* Gerardus headers */
#import "GerardusCommon.hpp"

/*
 * sub2ind(): function that converts r, c, s indices to linear indices
 *            in a 3D array (same as Matlab's function sub2ind(),
 *            although in Matlab indices start at 1, and in C++, they
 *            start at 0)
 *
 */
mwIndex sub2ind(mwSize R, mwSize C, mwSize S,
		std::vector<mwIndex> rcs) {
  // check for out of range index
  if (
      (rcs[0] < 0) || (rcs[0] >= R) 
      || (rcs[1] < 0) || (rcs[1] >= C)
      || (rcs[2] < 0) || (rcs[2] >= S)
      ) {
    mexErrMsgTxt("Out of range index");
  }
  if ((R*C*S == 0) || (R < 0) || (C < 0) || (S < 0)) {
    mexErrMsgTxt("Size values cannot be 0 or negative");
  }

  // check that input vector has 3 elements
  if (rcs.size() != 3) {
    mexErrMsgTxt("Input vector must have 3 elements");
  }

  // convert r, c, s to linear index
  mwIndex idx = rcs[0] + rcs[1] * R + rcs[2] * R * C;

  // // DEBUG
  // std::cout << "idx = " << idx
  // 	    << std::endl;
  
  return idx;
}

mwIndex sub2ind(mwSize R, mwSize C, mwSize S,
		itk::Offset<Dimension> rcs) {
  // check for out of range index
  if (
      ((mwSize)rcs[0] < 0) || ((mwSize)rcs[0] >= R) 
      || ((mwSize)rcs[1] < 0) || ((mwSize)rcs[1] >= C)
      || ((mwSize)rcs[2] < 0) || ((mwSize)rcs[2] >= S)
      ) {
    mexErrMsgTxt("Out of range index");
  }
  if ((R*C*S == 0) || (R < 0) || (C < 0) || (S < 0)) {
    mexErrMsgTxt("Size values cannot be 0 or negative");
  }

  // convert r, c, s to linear index
  mwIndex idx = rcs[0] + rcs[1] * R + rcs[2] * R * C;

  // // DEBUG
  // std::cout << "idx = " << idx
  // 	    << std::endl;
  
  return idx;
}

/*
 * ind2sub(): function that converts linear indices in a 3D array to
 *            r, c, s indices (same as Matlab's function ind2sub(),
 *            although in Matlab indices start at 1, and in C++, they
 *            start at 0)
 *
 */
std::vector<mwIndex> ind2sub(mwSize R, mwSize C, mwSize S,
			     mwIndex idx) {
  // check for out of range index
  if (idx >= R*C*S || idx < 0) {
    mexErrMsgTxt("Out of range index");
  }
  if (R*C*S == 0 || R < 0 || C < 0 || S < 0) {
    mexErrMsgTxt("Size values cannot be 0 or negative");
  }

  // init output
  std::vector<mwIndex> rcs(3);
  
  // convert linear index to r, c, s 
  div_t divresult;

  divresult = div(idx, R*C);
  rcs[2] = divresult.quot; // slice value
  idx = divresult.rem;

  divresult = div(idx, R);
  rcs[1] = divresult.quot; // column value
  idx = divresult.rem;

  divresult = div(idx, 1);
  rcs[0] = divresult.quot; // row value

  // // DEBUG
  // std::cout << "rcs = " << rcs[0] << ", " 
  // 	    << rcs[1] << ", "
  // 	    << rcs[2]
  // 	    << std::endl;
  
  return rcs;
}

#endif /* GERARDUSCOMMON_CPP */
