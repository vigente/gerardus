/*
 * MatlabImportFilter.hxx
 *
 * Class to provide an interface to import data from Matlab mxArrays
 * into ITK
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012-2013 University of Oxford
  * Version: 0.8.1
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

/* C++ headers */

/* ITK headers */
#include "itkImportImageFilter.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>

/* Gerardus headers */
#include "GerardusCommon.h"
#include "VectorWrapper.h"
#include "MatlabImageHeader.h"
#include "MatlabImportFilter.h"

// constructor
MatlabImportFilter::MatlabImportFilter() {
  // left intentionally empty
}

// destructor
MatlabImportFilter::~MatlabImportFilter() {
  // left intentionally empty
}

// function to import into this class the array with the arguments
// provided by Matlab
void MatlabImportFilter::ConnectToMatlabFunctionInput(int _nrhs, const mxArray *_prhs[]) {
  this->nrhs = _nrhs;
  this->prhs = _prhs;
}

// get number of elements in the prhs list of input arguments
unsigned int MatlabImportFilter::GetNumberOfArguments() {
  return this->nrhs;
}

// function to get direct pointers to the Matlab input arguments
const mxArray* MatlabImportFilter::GetPrhsArgument(int idx) {
  if ((idx >= 0) && (idx < this->nrhs)) {
    return this->prhs[idx];
  } else {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:IndexOutOfRange", 
		      "Index of prhs argument is out of range");
  }

  // function will never reach this point, but we need to have a
  // return to avoid a compiler warning
  mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:AssertionFail", 
		    "GetPrhsArgument() should have returned before getting here");
  return NULL;
}

// function to check that number of prhs arguments is within
// certain limits
void MatlabImportFilter::CheckNumberOfArguments(int min, int max) {
  if (this->nrhs < min) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      "Not enough input arguments");
  }
  if (this->nrhs > max) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      "Too many input arguments");
  }
}

// Functions to register an input at the import filter. 
MatlabImportFilter::MatlabInputPointer
MatlabImportFilter::RegisterInput(int pos, std::string name) {
  
  // if the user has put something in this input, even if it's an
  // empty array
  if (pos < this->nrhs) {
    // then we can register the input with its pointer
    return this->RegisterInput(this->prhs[pos], name);
  } else {
    // there's no pointer for this input, so we pass a NULL
    // pointer. The input will be registered, but will apear as
    // input.isProvided=false
    return this->RegisterInput((mxArray *)NULL, name);
  }

}

MatlabImportFilter::MatlabInputPointer
MatlabImportFilter::RegisterInput(const mxArray *pm, std::string name) {

  MatlabImportFilter::MatlabInput input;

  if (pm == NULL) {
    // the pointer to the input can be NULL. For example, this comes
    // from a mxGetField() of a non-existent struct field. Throwing an
    // error here would force the user to always provide all possible
    // fields. Instead, we just flag the input as unavailable and exit
    input.pm = NULL;
    input.name = name;
    input.isProvided = false;
  } else {
    // assign main variables of the input: its name and the memory
    // address it is stored at
    input.pm = pm;
    input.name = name;
    input.isProvided = (pm != NULL) && !mxIsEmpty(pm);
  }

  // insert the new input at the beginning of the list of registered inputs
  MatlabImportFilter::MatlabInputPointer it;
  try {
    it = this->inputsList.insert(this->inputsList.begin(), input);
  } catch (std::exception& e) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:RegisterInput", 
		      ("Input " + name + ": Cannot insert into list of registered inputs\n" + e.what()).c_str());
  }
  return it;

}

// Function to register a field from a struct input at the import
// filter. Basically, register input.field.
//
// structInput:
//   pointer to an already registered input.
//
// field:
//   name of the field we want to register.
//
// returns:
//   a class of type MatlabInput, defined above
MatlabImportFilter::MatlabInputPointer 
MatlabImportFilter::RegisterStructFieldInput(MatlabImportFilter::MatlabInputPointer structInput, 
					     std::string field) {

  // make sure that input is a struct
  if (!mxIsStruct(structInput->pm)) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:InputType", 
		      ("Input " + structInput->name + " must be a struct").c_str());
  }

  // get pointer to the field argument
  const mxArray *fieldarg = mxGetField(structInput->pm, 0, field.c_str());

  // register the field argument as an input
  return RegisterInput(fieldarg, (structInput->name + "." + field).c_str());
  
}

// function to get a pointer to a registered Matlab input by
// providing its name
MatlabImportFilter::MatlabInputPointer 
MatlabImportFilter::GetRegisteredInput(std::string name) {
  
  // iterator to iterate through the list of registered inputs
  std::list<MatlabInput>::iterator it;

  // look for the requested input
  for (it = this->inputsList.begin(); it != this->inputsList.end(); ++it) {
    if (it->name == name) {
      return it;
    }
  }

  // if the input has not been found, we throw an error (We make a
  // mandatory requirement that users only try to GetRegisteredInput()
  // after RegisterInput())
  mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:UnregisteredInputRequested", 
		    ("GetRegisteredInput(): Input " + name 
		     + " requested, but it has not been registered").c_str());

  // the program will never reach this point, but we need a return
  // statement to avoid a warning with the compiler
  return it;

}


// function to get the size of a Matlab array. It simplifies having
// to run mxGetNumberOfDimensions() and mxGetDimensions(), and then
// casting the result into e.g. itk::Size to pass it to ITK
template <class VectorValueType, class VectorType>
VectorType MatlabImportFilter::ReadMatlabArraySize(MatlabImportFilter::MatlabInputPointer input,
						   VectorType def){

  // template VectorSize will be ignored in this syntax, so we assign
  // to it an arbitrary value=0
  return MatlabImportFilter::ReadMatlabArraySize<VectorValueType, VectorType, 0>(input, def);

}

template <class VectorValueType, class VectorType, unsigned int VectorSize>
VectorType MatlabImportFilter::ReadMatlabArraySize(MatlabImportFilter::MatlabInputPointer input,
						   VectorType def){

  // if user didn't provide a value, or provided an empty array, return the default
  if (!input->isProvided) {
    return def;
  }

  // wrap VectorType output into a VectorWrapper, so that we don't
  // need to write different code here for each different output
  // vector
  VectorWrapper<VectorValueType, VectorType, void> sizeWrap;
  return sizeWrap.ReadSize(input->pm, input->name);

}

// function to get the half-size of a Matlab array. Some ITK filters
// request the "half-size" (called radius) of a Matlab array,
// instead of its size. By "half-size" we mean the length of the side to
// the left or right of the central pixel. For example, an array
// with size=[3, 7] has a half-size or radius=[1, 3]. I.e. 
// size = 2 * halfsize + 1
template <class VectorValueType, class VectorType>
VectorType MatlabImportFilter::ReadMatlabArrayHalfSize(MatlabImportFilter::MatlabInputPointer input,
						       VectorType def){

  return MatlabImportFilter::ReadMatlabArrayHalfSize<VectorValueType, VectorType, 0>(input, def);

}

template <class VectorValueType, class VectorType, unsigned int VectorSize>
VectorType MatlabImportFilter::ReadMatlabArrayHalfSize(MatlabImportFilter::MatlabInputPointer input,
						       VectorType def){

  // if user didn't provide a value, or provided an empty array, return the default
  if (!input->isProvided) {
    return def;
  }

  // wrap VectorType output into a VectorWrapper, so that we don't
  // need to write different code here for each different output
  // vector
  VectorWrapper<VectorValueType, VectorType, void> sizeWrap;
  return sizeWrap.ReadHalfSize(input->pm, input->name);

}

// function to get the value of input arguments that are strings
std::string MatlabImportFilter::ReadStringFromMatlab(MatlabImportFilter::MatlabInputPointer input,
						     std::string def) {
  
  // if user didn't provide a value, or provided an empty array, return the default
  if (!input->isProvided) {
    return def;
  }

  // if user provided a value, check that it's a string
  if (!mxIsChar(input->pm)) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " must be a string").c_str());
  }
  
  // import string
  return mxArrayToString(input->pm);

}

// function to get the value of input parameters that are numeric
// scalars from the array of input arguments
template <class ParamType>
ParamType MatlabImportFilter::ReadScalarFromMatlab(MatlabImportFilter::MatlabInputPointer input,
						   ParamType def) {
  return MatlabImportFilter::ReadScalarFromMatlab<ParamType>(input, 0, 0, def);
}

// function to get one scalar value from an input argument that is a matrix
template <class ParamType>
ParamType MatlabImportFilter::ReadScalarFromMatlab(MatlabImportFilter::MatlabInputPointer input,
						   mwIndex row, mwIndex col, ParamType def) {

  // if user didn't provide a value, or provided an empty array, return the default
  if (!input->isProvided) {
    return def;
  }
  
  // check for null pointer
  if (input->pm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:AssertionFail", 
		      ("Input " + input->name + " flagged as provided, but pointer to input is NULL").c_str());
  }

  // if user provided a parameter, check that it's a scalar, whether
  // in numeric or logical form
  if (!mxIsNumeric(input->pm) && !mxIsLogical(input->pm)) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " must be of scalar or logical type").c_str());
  }
  
  // get size of input matrix
  mwSize nrows = mxGetM(input->pm);
  mwSize ncols = mxGetN(input->pm);

  // check that requested row and column are within the matrix range
  if (row < 0 || row >= nrows) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + ": row index out of bounds").c_str());
  }
  if (col < 0 || col >= ncols) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + ": column index out of bounds").c_str());
  }

  // output
  ParamType value = 0;
  
  // input image type
  mxClassID inputVoxelClassId = mxGetClassID(input->pm);
  
  // macro to make the code in the switch statement cleaner
#define GETVALUE(Tx)						\
  {								\
    Tx *valuep = (Tx *)mxGetData(input->pm);		\
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
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " has unknown type").c_str());
    break;
  default:
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " has invalid type").c_str());
    break;
  }
  
#undef GETVALUE
  
  return value;
}

// function to read a row from a Matlab 2D matrix into a C++
// "vector". By "vector" we mean a C++ class that is vector-like,
// e.g. std::vector, CGAL::Point_3 or ITK::Size.
//
// Read the help of the VectorWrapper class defined in VectorWrapper.h for
// a list of supported vector-like types.
//
// Note that you don't need to worry about the type of the scalars in
// Matlab. The type will be automatically detected and cast to the
// vector element type.
template <class VectorValueType, class VectorType>
VectorType MatlabImportFilter::ReadRowVectorFromMatlab(MatlabImportFilter::MatlabInputPointer input, 
						       mwIndex row, VectorType def) {
  
  // if user didn't provide a value, or provided an empty array,
  // return default
  if (!input->isProvided) {
    return def;
  }
  
  // check for null pointer
  if (input->pm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:AssertionFail", 
		      ("Input " + input->name + " flagged as provided, but pointer to input is NULL").c_str());
  }

  // check that we have a 2D matrix, numeric or boolean
  if (mxGetNumberOfDimensions(input->pm) > 2) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " cannot have more than 2 dimensions").c_str());
  }
  if (!mxIsNumeric(input->pm) && !mxIsLogical(input->pm)) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " must be numeric or logical").c_str());
  }

  // input matrix type
  mxClassID inputVoxelClassId = mxGetClassID(input->pm);
  
  // cast the class type provided by Matlab to the type requested by
  // the user
  switch(inputVoxelClassId)  { 
  case mxLOGICAL_CLASS:
    {VectorWrapper<VectorValueType, VectorType, mxLogical> paramWrap;
      return paramWrap.ReadRowVector(input->pm, row, input->name);}
    break;
  case mxDOUBLE_CLASS:
    {VectorWrapper<VectorValueType, VectorType, double> paramWrap;
      return paramWrap.ReadRowVector(input->pm, row, input->name);}
    break;
  case mxSINGLE_CLASS:
    {VectorWrapper<VectorValueType, VectorType, float> paramWrap;
      return paramWrap.ReadRowVector(input->pm, row, input->name);}
    break;
  case mxINT8_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int8_T> paramWrap;
      return paramWrap.ReadRowVector(input->pm, row, input->name);}
    break;
  case mxUINT8_CLASS:
    {VectorWrapper<VectorValueType, VectorType, uint8_T> paramWrap;
      return paramWrap.ReadRowVector(input->pm, row, input->name);}
    break;
  case mxINT16_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int16_T> paramWrap;
      return paramWrap.ReadRowVector(input->pm, row, input->name);}
    break;
  case mxUINT16_CLASS:
    {VectorWrapper<VectorValueType, VectorType, uint16_T> paramWrap;
      return paramWrap.ReadRowVector(input->pm, row, input->name);}
    break;
  case mxINT32_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int32_T> paramWrap;
      return paramWrap.ReadRowVector(input->pm, row, input->name);}
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    // Note: mxINT64_CLASS causes compilation errors in Windows 64 bit, even though it's fine in linux
    //{VectorWrapper<VectorValueType, VectorType, int64_T> paramWrap;
    //  return paramWrap.ReadRowVector(input->pm, row, input->name);}
    //break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " has unknown type").c_str());
    break;
  default:
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " has invalid type").c_str());
    break;
  }
  
  // we should never get here, but we need to provide a return, else
  // the compiler will give a warning
  mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:AssertionFail", 
		    "ReadRowVectorFromMatlab() should have returned before getting here");
  return def;

}

// particular case in which the input matrix must be a row vector
template <class VectorValueType, class VectorType>
VectorType MatlabImportFilter::ReadRowVectorFromMatlab(MatlabImportFilter::MatlabInputPointer input,
						       VectorType def) {

  // if user didn't provide a value, or provided an empty array,
  // return default
  if (!input->isProvided) {
    return def;
  }

  // check that we have a numeric or boolean row vector
  if (mxGetM(input->pm) != 1) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " must be a row vector").c_str());
  }
  if (!mxIsNumeric(input->pm) && !mxIsLogical(input->pm)) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " must be a numeric or logical vector").c_str());
  }

  // the syntax of this function without specifying a row is the same as specifying row 0
  return MatlabImportFilter::ReadRowVectorFromMatlab<VectorValueType, VectorType>(input, 0, def);

}

// function to read a Matlab 2D matrix row by row. It returns the
// matrix as an std::vector of rows. Each row is read as a C++
// "vector". By "vector" we mean a C++ class that is vector-like,
// e.g. std::vector, CGAL::Point_3 or ITK::Size.
//
// Read the help of the VectorWrapper class defined in VectorWrapper.h for
// a list of supported vector-like types.
//
// Note that you don't need to worry about the type of the scalars in
// Matlab. The type will be automatically detected and cast to the
// vector element type.
//
// VectorValueType is the type of each element in the "vector".
// VectorType      is the type of the "vector" itself
template <class VectorValueType, class VectorType>
std::vector<VectorType> 
MatlabImportFilter::ReadVectorOfVectorsFromMatlab(MatlabImportFilter::MatlabInputPointer input,
						  std::vector<VectorType> def) {

  // if user didn't provide a value, or provided an empty array,
  // return default
  if (!input->isProvided) {
    return def;
  }
  
  // check for null pointer
  if (input->pm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:AssertionFail", 
		      ("Input " + input->name + " flagged as provided, but pointer to input is NULL").c_str());
  }

  // check that we have a 2D matrix, numeric or boolean
  if (mxGetNumberOfDimensions(input->pm) > 2) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " cannot have more than 2 dimensions").c_str());
  }
  if (!mxIsNumeric(input->pm) && !mxIsLogical(input->pm)) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " must be numeric or logical").c_str());
  }

  // number of rows in the matrix
  mwSize nrows = mxGetM(input->pm);

  // output vector where we are going to put all the rows from the matrix
  std::vector<VectorType> v(nrows);

  // input matrix type
  mxClassID inputVoxelClassId = mxGetClassID(input->pm);
  
  // cast the class type provided by Matlab to the type requested by
  // the user
  switch(inputVoxelClassId)  { 
  case mxLOGICAL_CLASS:
    {VectorWrapper<VectorValueType, VectorType, mxLogical> paramWrap;
      for (mwIndex row = 0; row < nrows; ++row) {
	v[row] = paramWrap.ReadRowVector(input->pm, row, input->name);
      }
    }
    break;
  case mxDOUBLE_CLASS:
    {VectorWrapper<VectorValueType, VectorType, double> paramWrap;
      for (mwIndex row = 0; row < nrows; ++row) {
	v[row] = paramWrap.ReadRowVector(input->pm, row, input->name);
      }
    }
    break;
  case mxSINGLE_CLASS:
    {VectorWrapper<VectorValueType, VectorType, float> paramWrap;
      for (mwIndex row = 0; row < nrows; ++row) {
	v[row] = paramWrap.ReadRowVector(input->pm, row, input->name);
      }
    }
    break;
  case mxINT8_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int8_T> paramWrap;
      for (mwIndex row = 0; row < nrows; ++row) {
	v[row] = paramWrap.ReadRowVector(input->pm, row, input->name);
      }
    }
    break;
  case mxUINT8_CLASS:
    {VectorWrapper<VectorValueType, VectorType, uint8_T> paramWrap;
      for (mwIndex row = 0; row < nrows; ++row) {
	v[row] = paramWrap.ReadRowVector(input->pm, row, input->name);
      }
    }
    break;
  case mxINT16_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int16_T> paramWrap;
      for (mwIndex row = 0; row < nrows; ++row) {
	v[row] = paramWrap.ReadRowVector(input->pm, row, input->name);
      }
    }
    break;
  case mxUINT16_CLASS:
    {VectorWrapper<VectorValueType, VectorType, uint16_T> paramWrap;
      for (mwIndex row = 0; row < nrows; ++row) {
	v[row] = paramWrap.ReadRowVector(input->pm, row, input->name);
      }
    }
    break;
  case mxINT32_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int32_T> paramWrap;
      for (mwIndex row = 0; row < nrows; ++row) {
	v[row] = paramWrap.ReadRowVector(input->pm, row, input->name);
      }
    }
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int64_T> paramWrap;
      for (mwIndex row = 0; row < nrows; ++row) {
	v[row] = paramWrap.ReadRowVector(input->pm, row, input->name);
      }
    }
    break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " has unknown type").c_str());
    break;
  default:
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:BadInputFormat", 
		      ("Input " + input->name + " has invalid type").c_str());
    break;
  }

  // return the vector of rows, now that it has been populated
  return v;

}

// function to read a Matlab array into a vector. This is the
// equivalent to the linearising operator A(:) in Matlab
template <class VectorValueType, class VectorType>
VectorType
MatlabImportFilter::ReadArrayAsVectorFromMatlab(MatlabImportFilter::MatlabInputPointer input,
						VectorType def) {
  
  // if user didn't provide a value, or provided an empty array,
  // return default
  if (!input->isProvided) {
    return def;
  }
  
  // check for null pointer
  if (input->pm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabImportFilter:AssertionFail", 
		      ("Input " + input->name + " flagged as provided, but pointer to input is NULL").c_str());
  }

  // check that we have a numerical or logical array
  if (!mxIsNumeric(input->pm) && !mxIsLogical(input->pm)) {
    mexErrMsgTxt(("Input " + input->name 
		  + " must be a numeric or logical array.").c_str());
  }

  // input matrix type
  mxClassID inputVoxelClassId = mxGetClassID(input->pm);
  
  // cast the class type provided by Matlab to the type requested by
  // the user
  switch(inputVoxelClassId)  { 
  case mxLOGICAL_CLASS:
    {VectorWrapper<VectorValueType, VectorType, mxLogical> paramWrap;
      return paramWrap.ReadArrayAsVector(input->pm, input->name);}
    break;
  case mxDOUBLE_CLASS:
    {VectorWrapper<VectorValueType, VectorType, double> paramWrap;
      return paramWrap.ReadArrayAsVector(input->pm, input->name);}
    break;
  case mxSINGLE_CLASS:
    {VectorWrapper<VectorValueType, VectorType, float> paramWrap;
      return paramWrap.ReadArrayAsVector(input->pm, input->name);}
    break;
  case mxINT8_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int8_T> paramWrap;
      return paramWrap.ReadArrayAsVector(input->pm, input->name);}
    break;
  case mxUINT8_CLASS:
    {VectorWrapper<VectorValueType, VectorType, uint8_T> paramWrap;
      return paramWrap.ReadArrayAsVector(input->pm, input->name);}
    break;
  case mxINT16_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int16_T> paramWrap;
      return paramWrap.ReadArrayAsVector(input->pm, input->name);}
    break;
  case mxUINT16_CLASS:
    {VectorWrapper<VectorValueType, VectorType, uint16_T> paramWrap;
      return paramWrap.ReadArrayAsVector(input->pm, input->name);}
    break;
  case mxINT32_CLASS:
    {VectorWrapper<VectorValueType, VectorType, int32_T> paramWrap;
      return paramWrap.ReadArrayAsVector(input->pm, input->name);}
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    // Note: mxINT64_CLASS causes compilation errors in Windows 64 bit, even though it's fine in linux
    //{VectorWrapper<VectorValueType, VectorType, int64_T> paramWrap;
    //  return paramWrap.ReadArrayAsVector(input->pm, input->name);}
    //break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt(("Input " + input->name + " has unknown type.").c_str());
    break;
  default:
    mexErrMsgTxt(("Input " + input->name + " has invalid type.").c_str());
    break;
  }

  // we should never get here, but we need to provide a return, else
  // the compiler will give a warning
  mexErrMsgTxt("MatlabImportFilter::ReadArrayAsVectorFromMatlab: This function should have returned before getting here");
  return def;

}

// function to get an input argument that is an image
template <class TPixel, unsigned int VImageDimension>
typename itk::Image<TPixel, VImageDimension>::Pointer
MatlabImportFilter::GetImagePointerFromMatlab(MatlabImportFilter::MatlabInputPointer input) {
  
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
  MatlabImageHeader imageHeader(input->pm, input->name);

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
				 mxGetNumberOfElements(input->pm),
				 importFilterWillOwnTheBuffer);

  // actually import the image
  importFilter->Update();

  // succesful exit
  return importFilter->GetOutput();

}

// function to get a pointer to a Matlab image, and the
// metainformation, in a CGAL::_image format
_image*
MatlabImportFilter::ReadCgalImageFromMatlab(MatlabInputPointer input) {

  // get metainformation of the input image
  MatlabImageHeader imHeader(input->pm, input->name);
  
  if (imHeader.size.size() != 3) {
    mexErrMsgTxt(("Input " + input->name + " must be a 3D image.").c_str());
  }

  // to create an image, CGAL requires as input arguments "word size",
  // "kind of image word" and "image word sign"
  int wordSize = 0;
  WORD_KIND wordKind = WK_UNKNOWN;
  SIGN wordSign = SGN_UNKNOWN;
  switch (imHeader.type) {
  case mxLOGICAL_CLASS:
    wordSize = sizeof(mxLogical);
    wordKind = WK_FLOAT;
    wordSign = SGN_SIGNED;
    break;
  case mxDOUBLE_CLASS:
    wordSize = sizeof(double);
    wordKind = WK_FLOAT;
    wordSign = SGN_UNKNOWN;
    break;
  case mxSINGLE_CLASS:
    wordSize = sizeof(float);
    wordKind = WK_FLOAT;
    wordSign = SGN_UNKNOWN;
    break;
  case mxINT8_CLASS:
    wordSize = sizeof(int8_T);
    wordKind = WK_FIXED;
    wordSign = SGN_SIGNED;
    break;
  case mxUINT8_CLASS:
    wordSize = sizeof(uint8_T);
    wordKind = WK_FIXED;
    wordSign = SGN_UNSIGNED;
    break;
  case mxINT16_CLASS:
    wordSize = sizeof(int16_T);
    wordKind = WK_FIXED;
    wordSign = SGN_SIGNED;
    break;
  case mxUINT16_CLASS:
    wordSize = sizeof(uint16_T);
    wordKind = WK_FIXED;
    wordSign = SGN_UNSIGNED;
    break;
  case mxINT32_CLASS:
    wordSize = sizeof(int32_T);
    wordKind = WK_FIXED;
    wordSign = SGN_SIGNED;
    break;
  // case mxUINT32_CLASS:
  //   wordSize = sizeof();
  //   wordKind = WK_FIXED;
  //   wordSign = SGN_UNSIGNED;
  //   break;
  case mxINT64_CLASS:
    wordSize = sizeof(int64_T);
    wordKind = WK_FIXED;
    wordSign = SGN_SIGNED;
    break;
  // case mxUINT64_CLASS:
  //   wordSize = sizeof();
  //   wordKind = WK_FIXED;
  //   wordSign = SGN_UNSIGNED;
  //   break;
  default:
    mexErrMsgTxt(("Input " + input->name + " has invalid type.").c_str());
  }

  // convert image header from Gerardus to CGAL format
  // _createImage in include/CGAL/ImageIO.h
  //
  // Important: While in Matlab x->cols, y->rows, in CGAL y->cols,
  // x->rows, so we need to swap the row, col dimensions stored in
  // imHeader. That is why below we have e.g. x <-> imHeader.size[0],
  // y <-> imHeader.size[1]
  _image *im = _createImage(
			    imHeader.size[0],     // number of rows, image x dimension
			    imHeader.size[1],     // number of cols, image y dimension
			    imHeader.size[2],     // number of slices, image z dimension
			    1,                    // image vectorial dimension (1 = scalar voxels)
			    imHeader.spacing[0],  // image voxel size in x dimension
			    imHeader.spacing[1],  // image voxel size in y dimension
			    imHeader.spacing[2],  // image voxel size in z dimension
			    wordSize,             // image word size in bytes
			    wordKind,             // image word kind
			    wordSign              // image word sign
			    );

  // set image offset
  im->tx = imHeader.origin[0];
  im->ty = imHeader.origin[1];
  im->tz = imHeader.origin[2];

  // duplicate the input Matlab buffer with the image. The reason for
  // this is that when a variable "_image im" is wrapped by
  // Image_3(im), the program will attempt to free im when all
  // references to im have disappeared. This can be seen in file
  // src/CGAL_ImageIO/Image_3.cpp, function
  // Image_3::private_read(_image* im). im becomes a shared pointer
  //
  // image_ptr = Image_shared_ptr(im, Image_deleter());
  //
  // so when all references to im disappear, then Image_deleter() will
  // try to free up the memory, including im->data. This causes Matlab
  // to crash, because im->data now points to a const array managed by
  // Matlab
  switch (imHeader.type) {
  case mxLOGICAL_CLASS:
    {
      bool *p = (bool *)mxGetData(imHeader.data);
      im->data = new bool [mxGetNumberOfElements(imHeader.data)];
      std::copy(p, p + mxGetNumberOfElements(imHeader.data), (bool *)im->data);
      break;
    }
  case mxDOUBLE_CLASS:
    {
      double *p = (double *)mxGetData(imHeader.data);
      im->data = new double [mxGetNumberOfElements(imHeader.data)];
      std::copy(p, p + mxGetNumberOfElements(imHeader.data), (double *)im->data);
      break;
    }
  case mxSINGLE_CLASS:
    {
      float *p = (float *)mxGetData(imHeader.data);
      im->data = new float [mxGetNumberOfElements(imHeader.data)];
      std::copy(p, p + mxGetNumberOfElements(imHeader.data), (float *)im->data);
      break;
    }
  case mxINT8_CLASS:
    {
      int8_T *p = (int8_T *)mxGetData(imHeader.data);
      im->data = new int8_T [mxGetNumberOfElements(imHeader.data)];
      std::copy(p, p + mxGetNumberOfElements(imHeader.data), (int8_T *)im->data);
      break;
    }
  case mxUINT8_CLASS:
    {
      uint8_T *p = (uint8_T *)mxGetData(imHeader.data);
      im->data = new uint8_T [mxGetNumberOfElements(imHeader.data)];
      std::copy(p, p + mxGetNumberOfElements(imHeader.data), (uint8_T *)im->data);
      break;
    }
  case mxINT16_CLASS:
    {
      int16_T *p = (int16_T *)mxGetData(imHeader.data);
      im->data = new int16_T [mxGetNumberOfElements(imHeader.data)];
      std::copy(p, p + mxGetNumberOfElements(imHeader.data), (int16_T *)im->data);
      break;
    }
  case mxUINT16_CLASS:
    {
      uint16_T *p = (uint16_T *)mxGetData(imHeader.data);
      im->data = new uint16_T [mxGetNumberOfElements(imHeader.data)];
      std::copy(p, p + mxGetNumberOfElements(imHeader.data), (uint16_T *)im->data);
      break;
    }
  case mxINT32_CLASS:
    {
      int32_T *p = (int32_T *)mxGetData(imHeader.data);
      im->data = new int32_T [mxGetNumberOfElements(imHeader.data)];
      std::copy(p, p + mxGetNumberOfElements(imHeader.data), (int32_T *)im->data);
      break;
    }
  // case mxUINT32_CLASS:
  //   break;
  case mxINT64_CLASS:
    {
      int64_T *p = (int64_T *)mxGetData(imHeader.data);
      im->data = new int64_T [mxGetNumberOfElements(imHeader.data)];
      std::copy(p, p + mxGetNumberOfElements(imHeader.data), (int64_T *)im->data);
      break;
    }
  // case mxUINT64_CLASS:
  //   break;
  default:
    mexErrMsgTxt(("Input " + input->name + " has invalid type.").c_str());
  }

  // not sure what these variables do, but if we don't initialize
  // them, they'll have random values. They don't seem to do anything,
  // but it'd better to get a consistent error anyway that can be
  // debugged rather than one that will happen now and then
  im->spm_offset = 0;
  im->spm_scale = 0;

  // // DEBUG:
  // std::cout << "---------------------------" << std::endl;
  // std::cout << "Metainformation for _image:" << std::endl;
  // std::cout << "---------------------------" << std::endl;
  // std::cout << "Within this block, x <-> rows, y <-> cols, z <-> slices" << std::endl;
  // std::cout << "while in Matlab,   x <-> cols, y <-> rows, z <-> slices" << std::endl;
  // std::cout << "*data = " << im->data << std::endl;
  // std::cout << "xdim = " << im->xdim << std::endl;
  // std::cout << "ydim = " << im->ydim << std::endl;
  // std::cout << "zdim = " << im->zdim << std::endl;
  // std::cout << "vectorial dimension = " << im->vdim << std::endl;
  // std::cout << "voxel size in x = " << im->vx << std::endl;
  // std::cout << "voxel size in y = " << im->vy << std::endl;
  // std::cout << "voxel size in z = " << im->vz << std::endl;
  // std::cout << "offset in x = " << im->tx << std::endl;
  // std::cout << "offset in y = " << im->ty << std::endl;
  // std::cout << "offset in z = " << im->tz << std::endl;
  // std::cout << "rotation in x = " << im->rx << std::endl;
  // std::cout << "rotation in y = " << im->ry << std::endl;
  // std::cout << "rotation in z = " << im->rz << std::endl;
  // std::cout << "centre in x = " << im->cx << std::endl;
  // std::cout << "centre in y = " << im->cy << std::endl;
  // std::cout << "centre in z = " << im->cz << std::endl;
  // std::cout << "spm offset = " << im->spm_offset << std::endl;
  // std::cout << "spm scale = " << im->spm_scale << std::endl;
  // std::cout << "word size = " << im->wdim << " bytes" << std::endl;
  // std::cout << "image format to use for I/0 = " << im->imageFormat << std::endl;
  // std::cout << "data buffer vectors are interlaced or non interlaced = " << im->vectMode << std::endl;
  // std::cout << "image words kind = " << im->wordKind << std::endl;
  // std::cout << "image words sign = " << im->sign << std::endl;
  // std::cout << "kind of image file descriptor = " << im->openMode << std::endl;
  // std::cout << "written words endianness = " << im->endianness << std::endl;
  // std::cout << "kind of image data encoding = " << im->dataMode << std::endl;
  // std::cout << "---------------------------" << std::endl;

  return im;

}

// function to swap the values of X and Y in a vector of vectors.
template <class ValueType, class OutsideVectorType>
void 
MatlabImportFilter::SwapXYInVectorOfVectors(OutsideVectorType &vv, mwSize len) {
  
  // loop the outside vector
  for (mwIndex i = 0; i < len; ++i) {
    ValueType temp = vv[i][0];
    vv[i][0] = vv[i][1];
    vv[i][1] = temp;
  }
}

// function to swap the values of X and Y in a vector.
template <class ValueType, class VectorType>
void 
MatlabImportFilter::SwapXYInVector(VectorType &v) {
  
  ValueType temp = v[0];
  v[0] = v[1];
  v[1] = temp;
}

#endif /* MATLABIMPORTFILTER_HXX */
