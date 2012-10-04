/*
 * MatlabImportFilter.hxx
 *
 * Class to provide an interface to import data from Matlab mxArrays
 * into ITK
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012 University of Oxford
  * Version: 0.3.2
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

#ifndef MATLABIMPORTFILTER_HXX
#define MATLABIMPORTFILTER_HXX

/* ITK headers */
#include "itkImportImageFilter.h"

/* Gerardus headers */
#include "GerardusCommon.h"
#include "MatlabImageHeader.h"
#include "MatlabImportFilter.h"

// constructor
MatlabImportFilter::MatlabImportFilter() {
  this->args.clear();
}

// destructor
// this class is only an interface to handle pointer, so we must not
// attempt to delete the argument list provided by Matlab
MatlabImportFilter::~MatlabImportFilter() {}

// function to check that number of input arguments is within
// certain limits
void MatlabImportFilter::CheckNumberOfArguments(unsigned int min, unsigned int max) {
  if (this->args.size() < min) {
    mexErrMsgTxt("Not enough input arguments");
  }
  if (this->args.size() > max) {
    mexErrMsgTxt("Too many input arguments");
  }
}

// function to get the value of input arguments that are strings
std::string MatlabImportFilter::GetStringArgument(unsigned int idx,
						       std::string paramName,
						       std::string def) {
  
  // if user didn't provide a value, or provided an empty array, return the default
  if (idx >= this->args.size() || mxIsEmpty(this->args[idx])) {
    return def;
  }

  // if user provided a value, check that it's a string
  if (!mxIsChar(this->args[idx])) {
    mexErrMsgTxt(("Parameter " + paramName + " must be a string").c_str());
  }
  
  // import string
  return mxArrayToString(this->args[idx]);

}

// function to get the value of input parameters that are numeric
// scalars from the array of input arguments
template <class ParamType>
ParamType MatlabImportFilter::GetScalarArgument(unsigned int idx, 
						std::string paramName,
						ParamType def) {
  return MatlabImportFilter::GetScalarArgument<ParamType>(idx, 0, 0, paramName, def);
}

// function to get one scalar value from an input argument that is a matrix
//
// idx: argument index
// row: matrix row index of the scalar
// col: matrix column index of the scalar
// def: value returned by default if argument is empty or not provided
template <class ParamType>
ParamType MatlabImportFilter::GetScalarArgument(unsigned int idx,
						mwIndex row,
						mwIndex col,
						std::string paramName,
						ParamType def) {

  // if user didn't provide a value, or provided an empty array, return the default
  if (idx >= this->args.size() || mxIsEmpty(this->args[idx])) {
    return def;
  }
  
  // check for null pointer
  if (this->args[idx] == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " provided, but NULL pointer.").c_str());
  }

  // if user provided a parameter, check that it's a scalar, whether
  // in numeric or logical form
  if (!mxIsNumeric(this->args[idx])
       && !mxIsLogical(this->args[idx])) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must be of scalar or logical type").c_str());
  }
  
  // get size of input matrix
  mwSize nrows = mxGetM(this->args[idx]);
  mwSize ncols = mxGetN(this->args[idx]);

  // check that requested row and column are within the matrix range
  if (row < 0 || row >= nrows) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": row index out of bounds").c_str());
  }
  if (col < 0 || col >= ncols) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": column index out of bounds").c_str());
  }

  // output
  ParamType value = 0;
  
  // input image type
  mxClassID inputVoxelClassId = mxGetClassID(this->args[idx]);
  
  // macro to make the code in the switch statement cleaner
#define GETVALUE(Tx)						\
  {								\
    Tx *valuep = (Tx *)mxGetData(this->args[idx]);		\
    value = (ParamType)valuep[col * nrows + row];		\
  }
  
  // cast the class type provided by Matlab to the type requested by
  // the user
  switch(inputVoxelClassId)  {
  case mxLOGICAL_CLASS:
    GETVALUE(mxLogical);
    break;
  case mxDOUBLE_CLASS:
    GETVALUE(double);
    break;
  case mxSINGLE_CLASS:
    GETVALUE(float);
    break;
  case mxINT8_CLASS:
    GETVALUE(int8_T);
    break;
  case mxUINT8_CLASS:
    GETVALUE(uint8_T);
    break;
  case mxINT16_CLASS:
    GETVALUE(int16_T);
    break;
  case mxUINT16_CLASS:
    GETVALUE(uint16_T);
    break;
  case mxINT32_CLASS:
    GETVALUE(int32_T);
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    GETVALUE(int64_T);
    break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt(("Parameter " + paramName + " has unknown type.").c_str());
    break;
  default:
    mexErrMsgTxt(("Parameter " + paramName + " has invalid type.").c_str());
    break;
  }
  
#undef GETVALUE
  
  return value;
}

// function to get an input argument as a vector of scalars. The
// argument itself can be a row vector, or a 2D matrix. In the latter
// case, the user has to select one of the rows of the matrix
//
// this function is designed to deal with all types that are
// conceptually a vector, even if they are not a std::vector. Read the
// help of the VectorWrapper class defined in GerardusCommon.h for
// more details
template <class ParamType, class ParamValueType>
ParamType MatlabImportFilter::GetRowVectorArgument(unsigned int idx, 
						mwIndex row,
						std::string paramName,
						ParamType def) {

  // if user didn't provide a value, or provided an empty array,
  // return default
  if (idx >= this->args.size() || mxIsEmpty(this->args[idx])) {
    return def;
  }
  
  // check for null pointer
  if (this->args[idx] == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " provided, but NULL pointer.").c_str());
  }

  // check that we have a 2D matrix, numeric or boolean
  if (mxGetNumberOfDimensions(this->args[idx]) > 2) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must be a 2D matrix.").c_str());
  }
  if (!mxIsNumeric(this->args[idx]) && !mxIsLogical(this->args[idx])) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must be a numeric or logical matrix.").c_str());
  }

  // matrix dimensions
  mwSize nrows = mxGetM(this->args[idx]);
  mwSize ncols = mxGetN(this->args[idx]);

  // vector to be returned
  ParamType param;

  // check that the vector can have exactly the length of one row of
  // the input matrix
  VectorWrapper<ParamType, ParamValueType> paramWrap(param);
  paramWrap.Resize(ncols);

  // check that row index is within range
  if (row < 0 || row >= nrows) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": row index out of bounds.").c_str());
  }

  // input image type
  mxClassID inputVoxelClassId = mxGetClassID(this->args[idx]);
  
  // macro to make the code in the switch statement cleaner
#define GETVALUE(Tx)						\
  {								\
    Tx *valuep = (Tx *)mxGetData(this->args[idx]);		\
    for (size_t col = 0; col < ncols; ++col) {			\
      param[col] = (ParamValueType)valuep[col * nrows + row];	\
    }								\
  }
  
  // cast the class type provided by Matlab to the type requested by
  // the user
  switch(inputVoxelClassId)  { 
  case mxLOGICAL_CLASS:
    GETVALUE(mxLogical);
    break;
  case mxDOUBLE_CLASS:
    GETVALUE(double);
    break;
  case mxSINGLE_CLASS:
    GETVALUE(float);
    break;
  case mxINT8_CLASS:
    GETVALUE(int8_T);
    break;
  case mxUINT8_CLASS:
    GETVALUE(uint8_T);
    break;
  case mxINT16_CLASS:
    GETVALUE(int16_T);
    break;
  case mxUINT16_CLASS:
    GETVALUE(uint16_T);
    break;
  case mxINT32_CLASS:
    GETVALUE(int32_T);
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    GETVALUE(int64_T);
    break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt(("Parameter " + paramName + " has unknown type.").c_str());
    break;
  default:
    mexErrMsgTxt(("Parameter " + paramName + " has invalid type.").c_str());
    break;
  }

#undef GETVALUE

  return param;
}

template <class ParamType, class ParamValueType>
ParamType MatlabImportFilter::GetRowVectorArgument(unsigned int idx, 
						   std::string paramName,
						   ParamType def) {

  // if user didn't provide a value, or provided an empty array,
  // return default
  if (idx >= this->args.size() || mxIsEmpty(this->args[idx])) {
    return def;
  }

  // check for null pointer
  if (this->args[idx] == NULL) {
    mexErrMsgTxt(("Parameter " + paramName + " provided, but NULL pointer.").c_str());
  }

  // check that we have a row vector, numeric or boolean
  if (mxGetM(this->args[idx]) != 1) {
    mexErrMsgTxt(("Parameter " + paramName + " must be a row vector.").c_str());
  }
  if (!mxIsNumeric(this->args[idx]) && !mxIsLogical(this->args[idx])) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must be a numeric or logical vector.").c_str());
  }

  // the syntax of this function without specifying a row is the same as specifying row 0
  return MatlabImportFilter::GetRowVectorArgument<ParamType, ParamValueType>(idx, 0, paramName, def);

}

// function to read a static 3-vector (a vector with 3 elements that
// can only be created using the constructor) from one row of a 2D
// input matrix with 3 columns
//
// ParamType: a class of objects that are constructed like
//            x(-2, 0.3, 5), e.g. CGAL::Point_3< CGAL::Simple_cartesian<double> >
// idx: parameter index
// row: row index (C++ convention row = 0, 1, 2, ..., nrows-1)
template <class ParamType>
ParamType 
MatlabImportFilter::GetStaticVector3Argument(unsigned int idx, 
					     mwIndex row, 
					     std::string paramName,
					     ParamType def) {
  
  // if user didn't provide a value, or provided an empty array,
  // return default
  if (idx >= this->args.size() || mxIsEmpty(this->args[idx])) {
    return def;
  }
  
  // check for null pointer
  if (this->args[idx] == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " provided, but NULL pointer.").c_str());
  }

  // check that we have a 2D matrix, numeric or boolean
  if (mxGetNumberOfDimensions(this->args[idx]) > 2) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must be a 2D matrix.").c_str());
  }
  if (!mxIsNumeric(this->args[idx]) && !mxIsLogical(this->args[idx])) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must be a numeric or logical matrix.").c_str());
  }

  // matrix dimensions
  mwSize nrows = mxGetM(this->args[idx]);
  mwSize ncols = mxGetN(this->args[idx]);

  // check that we have 3 columns
  if (ncols != 3) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must have 3 columns.").c_str());
  }

  // check that row index is within range
  if (row < 0 || row >= nrows) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + ": row index out of bounds.").c_str());
  }

  // input matrix type
  mxClassID inputVoxelClassId = mxGetClassID(this->args[idx]);
  
  // parameter to be returned
  ParamType param(mxGetNaN(), mxGetNaN(), mxGetNaN());

  // macro to make the code in the switch statement cleaner
  // col is the column index
  // i is the linearised matrix index
#define GETVALUE(Tx)							\
  {									\
    Tx *valuep = (Tx *)mxGetData(this->args[idx]);			\
    param = ParamType((Tx)valuep[row],					\
		      (Tx)valuep[row + nrows],				\
		      (Tx)valuep[row + nrows + nrows]);			\
      }
  
  // cast the class type provided by Matlab to the type requested by
  // the user
  switch(inputVoxelClassId)  { 
  case mxLOGICAL_CLASS:
    GETVALUE(mxLogical);
    break;
  case mxDOUBLE_CLASS:
    GETVALUE(double);
    break;
  case mxSINGLE_CLASS:
    GETVALUE(float);
    break;
  case mxINT8_CLASS:
    GETVALUE(int8_T);
    break;
  case mxUINT8_CLASS:
    GETVALUE(uint8_T);
    break;
  case mxINT16_CLASS:
    GETVALUE(int16_T);
    break;
  case mxUINT16_CLASS:
    GETVALUE(uint16_T);
    break;
  case mxINT32_CLASS:
    GETVALUE(int32_T);
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    GETVALUE(int64_T);
    break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt(("Parameter " + paramName + " has unknown type.").c_str());
    break;
  default:
    mexErrMsgTxt(("Parameter " + paramName + " has invalid type.").c_str());
    break;
  }

#undef GETVALUE

  // return row from input matrix
  return param;
  
}

// function to get an input argument that is an image
template <class TPixel, unsigned int VImageDimension>
typename itk::Image<TPixel, VImageDimension>::Pointer
MatlabImportFilter::GetImageArgument(unsigned int idx, std::string paramName) {
  
  // note that:
  //
  // 1) in ITK we have X,Y,Z indices, while in Matlab we have R,C,S
  // (row, column, slice)
  //
  // 2) matrices in ITK are read by columns, while in Matlab
  // they are read by rows 
  //
  // So imagine we have this (2, 3) matrix in Matlab
  //
  //   a b   |
  //   c d   | y-axis (resolution 1.0)
  //   e f   |
  //   ---
  //   x-axis (resolution 0.5)
  //
  //   size = [3 2 1]
  //
  // The C-style array is going to be (reading by rows)
  //
  //   im = [a c e b d f]
  //
  // ITK is going to read by colums, thinking that the size is
  //
  //   size = [sx sy sz] = [3 2 1]
  //
  //   a c e   |
  //   b d f   | y-axis (resolution 0.5)
  //   -----
  //   x-axis (resolution 1.0)
  //
  // Note that the matrix has been transposed, but this is not a
  // problem, because the resolution values have been "transposed"
  // too
  //
  // Having the matrix transposed may make us feel a bit uneasy, but
  // it has the advantage that Matlab and ITK can use the same C-style
  // array, without having to rearrange its elements

  // instantiate the filter that will act as an interface between the
  // Matlab image array and the ITK filter
  typedef itk::ImportImageFilter<TPixel, VImageDimension> ImportImageFilterType;
  typename ImportImageFilterType::RegionType region;
  typename ImportImageFilterType::SizeType size;
  typename ImportImageFilterType::IndexType start;
  typename ImportImageFilterType::SpacingType spacing;
  typename ImportImageFilterType::OriginType origin;
  typename ImportImageFilterType::Pointer importFilter = ImportImageFilterType::New();

  // create a header for the image
  MatlabImageHeader imageHeader(this->args[idx], paramName);

  // convert image header parameters to a format that can be passed to
  // the import filter
  for (size_t i = 0; i < imageHeader.GetNumberOfDimensions(); ++i) {
    start[i] = 0;
    size[i] = imageHeader.size[i];
    spacing[i] = imageHeader.spacing[i];
    origin[i] = imageHeader.origin[i];
  }
  
  // set image metainformation
  region.SetIndex(start);
  region.SetSize(size);
  importFilter->SetRegion(region);
  importFilter->SetSpacing(spacing);
  importFilter->SetOrigin(origin);

  // get pointer to input segmentation mask
  const TPixel *im = (TPixel *)mxGetData(imageHeader.data);

  // pass pointer to Matlab image to the import filter, and tell it to
  // NOT attempt to delete the memory when it's destructor is
  // called. This is important, because the input image still has to
  // live in Matlab's memory after running the filter
  const bool importFilterWillOwnTheBuffer = false;
  importFilter->SetImportPointer(const_cast<TPixel *>(im),
				 mxGetNumberOfElements(this->args[idx]),
				 importFilterWillOwnTheBuffer);

  // actually import the image
  importFilter->Update();

  // succesful exit
  return importFilter->GetOutput();

}

// function to read a matrix where each row is a static 3-vector
template <class ParamType>
std::vector<ParamType>
MatlabImportFilter::GetVectorOfStaticVector3Argument(unsigned int idx, 
						     std::string paramName,
						     std::vector<ParamType> def) {

  // if user didn't provide a value, or provided an empty array,
  // return default
  if (idx >= this->args.size() || mxIsEmpty(this->args[idx])) {
    return def;
  }
  
  // check for null pointer
  if (this->args[idx] == NULL) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " provided, but NULL pointer.").c_str());
  }

  // check that we have a 2D matrix, numeric or boolean
  if (mxGetNumberOfDimensions(this->args[idx]) > 2) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must be a 2D matrix.").c_str());
  }
  if (!mxIsNumeric(this->args[idx]) && !mxIsLogical(this->args[idx])) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must be a numeric or logical matrix.").c_str());
  }

  // matrix dimensions
  mwSize nrows = mxGetM(this->args[idx]);
  mwSize ncols = mxGetN(this->args[idx]);

  // check that we have 3 columns
  if (ncols != 3) {
    mexErrMsgTxt(("Parameter " + paramName 
		  + " must have 3 columns.").c_str());
  }

  // input matrix type
  mxClassID inputVoxelClassId = mxGetClassID(this->args[idx]);
  
  // parameter to be returned
  std::vector<ParamType> param;

  // macro to make the code in the switch statement cleaner
  // col is the column index
  // i is the linearised matrix index
#define GETVALUE(Tx)							\
  {									\
    Tx *valuep = (Tx *)mxGetData(this->args[idx]);			\
    for (mwIndex row = 0; row < nrows; ++row) {				\
      param.push_back(ParamType((Tx)valuep[row],			\
				(Tx)valuep[row + nrows],		\
				(Tx)valuep[row + nrows + nrows]));	\
    }									\
  }
  
  // cast the class type provided by Matlab to the type requested by
  // the user
  switch(inputVoxelClassId)  { 
  case mxLOGICAL_CLASS:
    GETVALUE(mxLogical);
    break;
  case mxDOUBLE_CLASS:
    GETVALUE(double);
    break;
  case mxSINGLE_CLASS:
    GETVALUE(float);
    break;
  case mxINT8_CLASS:
    GETVALUE(int8_T);
    break;
  case mxUINT8_CLASS:
    GETVALUE(uint8_T);
    break;
  case mxINT16_CLASS:
    GETVALUE(int16_T);
    break;
  case mxUINT16_CLASS:
    GETVALUE(uint16_T);
    break;
  case mxINT32_CLASS:
    GETVALUE(int32_T);
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    GETVALUE(int64_T);
    break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt(("Parameter " + paramName + " has unknown type.").c_str());
    break;
  default:
    mexErrMsgTxt(("Parameter " + paramName + " has invalid type.").c_str());
    break;
  }

#undef GETVALUE

  // return row from input matrix
  return param;

}


#endif /* MATLABIMPORTFILTER_HXX */
