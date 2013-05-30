/*
 * VectorWrapper.hxx
 *
 * Code definitions of the declarations in VectorWrapper.h
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012 University of Oxford
  * Version: 0.3.0
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

#ifndef VECTORWRAPPER_HXX
#define VECTORWRAPPER_HXX

/* ITK headers */
#include "itkSize.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>

/* Boost headers */
#include <boost/lexical_cast.hpp>

/* Gerardus headers */
#include "VectorWrapper.h"

/*
 * By default, VectorWrapper assumes that we want to put Matlab's row data into an
 * std::vector<type>
 *
 * VectorSize: this template is ignored
 */

// read a row from a Matlab matrix
template<class VectorValueType, class VectorType, class MatlabValueType, mwSize VectorSize>
VectorType VectorWrapper<VectorValueType, VectorType, MatlabValueType, VectorSize>::ReadRowVector
  (const mxArray *pm, mwIndex row, std::string paramName) {

  // check that the pointer is valid
  if (pm == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": pointer to Matlab input argument is NULL").c_str());
  }
  
  // matrix dimensions
  mwSize nrows = mxGetM(pm);
  mwSize ncols = mxGetN(pm);
  
  // check that row index is within range
  if (row < 0 || row >= nrows) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": row index out of bounds").c_str());
  }
  
  // init vector with NaN values
  VectorType v(ncols, mxGetNaN());
  
  // get pointer to the data in the mxArray
  MatlabValueType *valuep = (MatlabValueType *)mxGetData(pm);
  
  if (valuep == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": pointer to content of Matlab input argument is NULL").c_str());
  }
  
  // read the row from Matlab into the vector. Note that valuep is
  // of type MatlabType, so we need to cast it to the VectorValueType
  for (size_t col = 0; col < ncols; ++col) {
    v[col] = (VectorValueType)valuep[col * nrows + row];
  }
  
  // return output vector
  return v;

};

// read a whole array into a vector
template<class VectorValueType, class VectorType, class MatlabValueType, mwSize VectorSize>
VectorType VectorWrapper<VectorValueType, VectorType, MatlabValueType, VectorSize>::ReadArrayAsVector
  (const mxArray *pm, std::string paramName) {
  
  // check that the pointer is valid
  if (pm == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": pointer to Matlab input argument is NULL").c_str());
  }
  
  // matrix dimensions
  mwSize numel = mxGetNumberOfElements(pm);
  
  // get pointer to the data in the mxArray
  MatlabValueType *valuep = (MatlabValueType *)mxGetData(pm);
  
  if (valuep == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": pointer to content of Matlab input argument is NULL").c_str());
  }
  
  // create output vector
  VectorType v;
  v.assign(valuep, valuep + numel);
  
  // return vector
  return v;
    
};



/*
 * Partial specialisation if we want to put Matlab's row data into an
 * itk::Size<Dimension>::SizeType vector-like class
 */

// ReadRowVector
//
// partial specialization
template<class MatlabValueType, mwSize VectorSize>
typename itk::Size<VectorSize>::SizeType
VectorWrapper<typename itk::Size<VectorSize>::SizeValueType, 
	      typename itk::Size<VectorSize>::SizeType, MatlabValueType, VectorSize>
::ReadRowVector(const mxArray *pm, mwIndex row, std::string paramName) {

  // check that the pointer is valid
  if (pm == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": pointer to Matlab input argument is NULL").c_str());
  }
  
  // matrix dimensions
  mwSize nrows = mxGetM(pm);
  mwSize ncols = mxGetN(pm);

  // check that the row has the right number of elements
  if (ncols != VectorSize) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must have " 
		  + boost::lexical_cast<std::string>(VectorSize) + " columns").c_str());
  }
  
  // check that row index is within range
  if (row < 0 || row >= nrows) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": row index out of bounds").c_str());
  }

  // instantiate output vector
  typename itk::Size<VectorSize>::SizeType v;
  
  // get pointer to the data in the mxArray
  MatlabValueType *valuep = (MatlabValueType *)mxGetData(pm);
  
  if (valuep == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": pointer to content of Matlab input argument is NULL").c_str());
  }
  
  // read the row from Matlab into the vector
  for (size_t col = 0; col < VectorSize; ++col) {
    v[col] = (MatlabValueType)valuep[col * nrows + row];
  }

  // return output vector
  return v;
  
}

// ReadSize
//
// partial specialization
template<class MatlabValueType, mwSize VectorSize>
typename itk::Size<VectorSize>::SizeType
VectorWrapper<typename itk::Size<VectorSize>::SizeValueType, 
	      typename itk::Size<VectorSize>::SizeType, MatlabValueType, VectorSize>
::ReadSize(const mxArray *pm, std::string paramName) {

  // check for null pointer
  if (pm == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " provided, but pointer is NULL").c_str());
  }
  
  // get number of dimensions
  mwSize ndim = mxGetNumberOfDimensions(pm);

  if (ndim != VectorSize) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": Cannot read parameter size. Output vector has wrong length.").c_str());
  }
  
  // get dimensions array
  const mwSize *dims = mxGetDimensions(pm);

  // init output
  typename itk::Size<VectorSize>::SizeType size;

  // copy dimensions to output vector
  for (mwIndex i = 0; i < VectorSize; ++i) {
    size[i] = (typename itk::Size<VectorSize>::SizeValueType)dims[i];
  }

  // return output
  return size;

}

// ReadHalfSize
//
// partial specialization
template<class MatlabValueType, mwSize VectorSize>
typename itk::Size<VectorSize>::SizeType
VectorWrapper<typename itk::Size<VectorSize>::SizeValueType, 
	      typename itk::Size<VectorSize>::SizeType, MatlabValueType, VectorSize>
::ReadHalfSize(const mxArray *pm, std::string paramName) {

  // check for null pointer
  if (pm == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " provided, but pointer is NULL").c_str());
  }
  
  // get number of dimensions
  mwSize ndim = mxGetNumberOfDimensions(pm);

  if (ndim != VectorSize) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": Cannot read parameter size. Output vector has wrong length.").c_str());
  }
  
  // get dimensions array
  const mwSize *dims = mxGetDimensions(pm);

  // init output
  typename itk::Size<VectorSize>::SizeType halfsize;

  // copy dimensions to output vector
  for (mwIndex i = 0; i < VectorSize; ++i) {
    if (dims[i] % 2) {
      halfsize[i] = (typename itk::Size<VectorSize>::SizeValueType)((dims[i] - 1)/2);
    } else {
      mexErrMsgTxt(("Parameter " + paramName 
		    + ": All values must be odd, in order to define a box around the central pixel").c_str());
    }
  }

  // return output
  return halfsize;

}

/*
 * Partial specialisation if we want to put Matlab's row data into a
 * CGAL::Point_3<CGAL::Simple_cartesian<type> >
 * CGAL::Direction_3<CGAL::Simple_cartesian<type> > 
 * vector-like class
 */

// ReadCgalRowVector
//
// auxiliary function so that we don't need to rewrite this code in
// every partial specialization
template<class VectorValueType, class VectorType, class MatlabValueType>
VectorType ReadCgalRowVector(const mxArray *pm, mwIndex row, std::string paramName) {
    
    // check that the pointer is valid
    if (pm == NULL) {
      mexErrMsgTxt(("Parameter " + paramName 
		    + ": pointer to Matlab input argument is NULL").c_str());
    }
    
    // matrix dimensions
    mwSize nrows = mxGetM(pm);
    mwSize ncols = mxGetN(pm);
    
    // check that the row has the right number of elements
    if (ncols != 3) {
      mexErrMsgTxt(("Parameter " + paramName 
		    + " must have 3 columns").c_str());
    }
    
    // check that row index is within range
    if (row < 0 || row >= nrows) {
      mexErrMsgTxt(("Parameter " + paramName 
		    + ": row index out of bounds").c_str());
    }
    
    // get pointer to the data in the mxArray
    MatlabValueType *valuep = (MatlabValueType *)mxGetData(pm);
    
    if (valuep == NULL) {
      mexErrMsgTxt(("Parameter " + paramName 
		    + ": pointer to content of Matlab input argument is NULL").c_str());
    }
    
    // instantiate and read output vector. Both have to be done at the
    // same time, because the only way to populate a CGAL::Point_3 is through
    // its constructor
    VectorType v( (VectorValueType)valuep[row],
		  (VectorValueType)valuep[nrows + row],
		  (VectorValueType)valuep[2 * nrows + row]
		  );
    
    // return output vector
    return v;

}

#endif /* VECTORWRAPPER_HXX */
