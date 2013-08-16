/*
 * MatlabExportFilter.hxx
 *
 * Class to provide an interface to graft ITK images onto Matlab
 * outputs
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012-2013 University of Oxford
  * Version: 0.4.0
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

#ifndef MATLABEXPORTFILTER_HXX
#define MATLABEXPORTFILTER_HXX

/* Boost headers */
#include <boost/lexical_cast.hpp> // doesn't need linking

/* ITK headers */
#include "itkImageRegionConstIterator.h"

/* Gerardus headers */
#include "MatlabExportFilter.h"

// constructor
MatlabExportFilter::MatlabExportFilter() {
  this->plhs = NULL;
  this->nlhs = 0;
}

// destructor
// this class is only an interface to handle pointer, so we must not
// attempt to delete the argument list provided by Matlab
MatlabExportFilter::~MatlabExportFilter() {}

// get number of elements in the list of arguments
int MatlabExportFilter::GetNumberOfOutputArguments() {
  return this->nlhs;
}

// function to check that number of input arguments is within
// certain limits
void MatlabExportFilter::CheckNumberOfArguments(int min, int max) {
  if (this->nlhs < min) {
    mexErrMsgTxt("Not enough output arguments");
  }
  if (this->nlhs > max) {
    mexErrMsgTxt("Too many output arguments");
  }
}

// Functions to register an output at the export filter. 
MatlabExportFilter::MatlabOutputPointer
MatlabExportFilter::RegisterOutput(int pos, std::string name) {
  
  return this->RegisterOutput(this->plhs, pos, name);

}
MatlabExportFilter::MatlabOutputPointer
MatlabExportFilter::RegisterOutput(mxArray **base, int pos, std::string name) {

  if (base == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:NullPointer", 
		      ("Output " + name + " cannot be registered at a NULL address").c_str());
  }

  // assign main variables of the output: its name and the memory
  // address it should be stored at
  MatlabExportFilter::MatlabOutput out;
  out.pm = &(base[pos]);
  out.name = name;

  // assign the flags to the output

  // is this output directly in the Matlab function's plhs output array?
  out.isTopLevel = (base == this->plhs);

  // if its at the top level, is it also the first output? This output
  // always needs to be allocated, even if the user didn't request it
  // (if not requested, we can return an empty matrix, to save time
  // and memory)
  out.isTopLevelFirst = out.isTopLevel && (pos == 0);


  // outputs that are regular output variables, e.g. y = myfun(x),
  // are requested if the user put them in the output vector. For
  // example:
  // [y,z] = myfun(x); // y and z have been requested
  // y = myfun(x);     // y has been requested, z has not been requested
  //
  // outputs that are not directly in the function output array are
  // always considered to have been requested by the user. For
  // example, an output variable is y (x = myfun(x)), and this is a
  // cell array. Output y is directly in the output array. But
  // output y{2} is not directly in the function array. So we assume
  // that as long as y is requested, then y{2} is requested too
  if (out.isTopLevel) {
    out.isRequested = (pos < this->GetNumberOfOutputArguments());
  } else {
    out.isRequested = true;
  }

  // insert the new output at the beginning of the list
  MatlabExportFilter::MatlabOutputPointer it;
  try {
    it = this->outputsList.insert(this->outputsList.begin(), out);
  } catch (std::exception& e) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:InsertionIntoList", 
		      ("Output " + name + ": Cannot insert into outputs list\n" + e.what()).c_str());
  }
  return it;

}

// Function to allocate memory for a column vector in Matlab, and get the
// data pointer back. 
template<class TData>
TData *
MatlabExportFilter::AllocateColumnVectorInMatlab(MatlabExportFilter::MatlabOutputPointer output, 
						 mwSize len) {
  
  // vector with matrix dimensions
  std::vector<mwSize> size;
  size.push_back(len);
  size.push_back(1);

  return this->AllocateNDArrayInMatlab<TData>(output, size);

}

// Function to allocate memory for a row vector in Matlab, and get the
// data pointer back. 
template<class TData>
TData *
MatlabExportFilter::AllocateRowVectorInMatlab(MatlabExportFilter::MatlabOutputPointer output, 
					      mwSize len) {

  // vector with matrix dimensions
  std::vector<mwSize> size;
  size.push_back(1);
  size.push_back(len);

  return this->AllocateNDArrayInMatlab<TData>(output, size);

}

// Function to allocate memory for a matrix in Matlab, and get the
// data pointer back. 
template<class TData>
TData *
MatlabExportFilter::AllocateMatrixInMatlab(MatlabExportFilter::MatlabOutputPointer output, 
					   mwSize nrows, mwSize ncols) {

  // vector with matrix dimensions
  std::vector<mwSize> size;
  size.push_back(nrows);
  size.push_back(ncols);

  return this->AllocateNDArrayInMatlab<TData>(output, size);

}

// Function to allocate memory for an N-dimensional array in Matlab, and get the
// data pointer back. 
template<class TData>
TData *
MatlabExportFilter::AllocateNDArrayInMatlab(MatlabExportFilter::MatlabOutputPointer output, 
					    std::vector<mwSize> size) {

  // get the Matlab class ID for the element type we need
  mxClassID outputClassId = convertCppDataTypeToMatlabCassId<TData>();

  // check whether there's already memory allocated in the output, and
  // in that case, de-allocate it
  if (*output->pm != NULL) {
    mxDestroyArray(*output->pm);
  }
  // create output matrix for Matlab's result
  mwSize ndim = size.size();
  mwSize dims[ndim];
  for (size_t i = 0; i < ndim; ++i) {
    dims[i] = size[i];
  }

  *output->pm = (mxArray *)mxCreateNumericArray(ndim, dims, outputClassId, mxREAL);
  if (*output->pm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:MemoryAllocation", 
		      ("Cannot allocate memory for output " + output->name).c_str());
  }
  
  // pointer to the Matlab output buffer
  TData *buffer =  (TData *)mxGetData(*output->pm);
  if(buffer == NULL && !mxIsEmpty(*output->pm)) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:MemoryAccess", 
		      ("Cannot get pointer to allocated memory for output " + output->name).c_str());
  }

  return buffer;

}

// Function to allocate memory for a column vector in Matlab, and get the
// data pointer back. 
template<class TData>
TData *
MatlabExportFilter::AllocateColumnVectorInCellInMatlab(MatlabExportFilter::MatlabOutputPointer output, 
						       int pos, mwSize len) {
  
  // vector with matrix dimensions
  std::vector<mwSize> size;
  size.push_back(len);
  size.push_back(1);

  return this->AllocateNDArrayInCellInMatlab<TData>(output, pos, size);

}

// Function to allocate memory for a row vector in Matlab, and get the
// data pointer back. 
template<class TData>
TData *
MatlabExportFilter::AllocateRowVectorInCellInMatlab(MatlabExportFilter::MatlabOutputPointer output, 
						    int pos, mwSize len) {

  // vector with matrix dimensions
  std::vector<mwSize> size;
  size.push_back(1);
  size.push_back(len);

  return this->AllocateNDArrayInCellInMatlab<TData>(output, pos, size);

}

// Function to allocate memory for a matrix in a cell in Matlab, and
// get the data pointer back.
template<class TData>
TData *
MatlabExportFilter::AllocateMatrixInCellInMatlab(MatlabExportFilter::MatlabOutputPointer output, 
						 int pos, mwSize nrows, mwSize ncols) {

  // vector with matrix dimensions
  std::vector<mwSize> size;
  size.push_back(nrows);
  size.push_back(ncols);

  return this->AllocateNDArrayInCellInMatlab<TData>(output, pos, size);

}

// Function to allocate memory for an N-dimensional array in a cell in
// Matlab, and get the data pointer back.
template<class TData>
TData *
MatlabExportFilter::AllocateNDArrayInCellInMatlab(MatlabExportFilter::MatlabOutputPointer output, 
						  int pos, std::vector<mwSize> size) {

  // get the Matlab class ID for the element type we need
  mxClassID outputClassId = convertCppDataTypeToMatlabCassId<TData>();

  // check whether there's already memory allocated in the cell, and
  // in that case, de-allocate it
  mxArray *cell = mxGetCell(*output->pm, pos);
  if (cell != NULL) {
    mxDestroyArray(cell);
  }

  // allocate memory for the new array
  mwSize ndim = size.size();
  mwSize dims[ndim];
  for (size_t i = 0; i < ndim; ++i) {
    dims[i] = size[i];
  }

  cell = mxCreateNumericArray(ndim, dims, outputClassId, mxREAL);
  if (cell == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:MemoryAllocation", 
		      ("Cannot allocate memory for output " + output->name + "{" 
		       + boost::lexical_cast<std::string>(pos) + "}").c_str());
  }

  // place the new array into the cell array
  mxSetCell(*output->pm, pos, cell);

  // pointer to the Matlab output buffer. If the array created in the
  // cell is empty, mxGetData will return a NULL pointer. Do not treat
  // this case as an error
  TData *buffer =  (TData *)mxGetData(cell);
  if(buffer == NULL && !mxIsEmpty(cell)) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:MemoryAccess", 
		      ("Cannot get pointer to allocated memory for output " + output->name + "{" 
		       + boost::lexical_cast<std::string>(pos) + "}").c_str());
  }

  return buffer;

}

// Function to create an empty output in Matlab.
void 
MatlabExportFilter::CopyEmptyArrayToMatlab(MatlabOutputPointer output) {
  
  // if we are asked to copy the data to an output argument that the
  // user has not requested, we avoid wasting time and memory, and
  // simply exit. The only exception is for the first output argument,
  // plhs[0]. Even if the user has not requested any output arguments,
  // e.g. my_mex_function([1 4.4 2]); we need to allocate memory for
  // an empty matrix at plhs[0], otherwise Matlab will give an "One or
  // more output arguments not assigned during call to..." error
  if (!output->isRequested) {
    if (output->isTopLevelFirst) {
      *output->pm = mxCreateDoubleMatrix(0, 0, mxREAL);
      return;
    } else {
      return;
    }
  }

  // the output has been requested, so we always have to create the
  // empty matrix
  *output->pm = mxCreateDoubleMatrix(0, 0, mxREAL);

}

// function to allocate memory in Matlab and copy an ITK filter
// output to this buffer. In principle, it's better to use
// GraftItkImageOntoMatlab() than CopyItkImageToMatlab(), because
// then ITK and Matlab share the same memory buffer, but the former
// approach does not work with some filters
//
// size is a vector with the dimensions of the output image in
// Matlab. For example, for a 256x200x512 image, size = {256, 200, 512}
//
// size is the same for vector or scalar images. 
template <class TPixel, unsigned int VectorDimension, class TVector>
void
MatlabExportFilter::GraftItkImageOntoMatlab(typename itk::DataObject::Pointer image,
					    std::vector<unsigned int> size,
					    int idx, std::string paramName) {

  if (size.size() != VectorDimension) {
    mexErrMsgTxt(("MatlabExportFilter: Output image " + paramName 
		  + ": wrong length for size vector").c_str());
  }

  // get the Matlab class ID for the element type we need
  mxClassID outputVoxelClassId = convertCppDataTypeToMatlabCassId<TPixel>();

  // dimensions for the output array
  typedef typename itk::Image<TVector, VectorDimension> OutputImageType;
  OutputImageType *pOutput = 
    dynamic_cast<OutputImageType *>(image.GetPointer());
  if (pOutput == NULL) {
    mexErrMsgTxt("MatlabExportFilter: Cannot get pointer to filter output");
  }

  // if the output image is a vector image, extend the size vector
  if (!TypesAreEqual<TPixel, TVector>::value) {
    size.insert(size.begin(), VectorDimension);
  }

  mwSize ndim = size.size();
  mwSize *dims;
  dims = (mwSize *)malloc(ndim * sizeof(mwSize));
  if (dims == NULL) {
	  mexErrMsgTxt("MatlabExportFilter: Cannot allocate memory for dims vector");
  }
  bool isEmptyMatrix = false;
  for (mwSize i = 0; i < ndim; ++i) {
    dims[i] = size[i];
    isEmptyMatrix |= (size[i] == 0);
  }

  // create output matrix for Matlab's result
  if (isEmptyMatrix) {
    this->plhs[idx] = mxCreateDoubleMatrix(0, 0, mxREAL);
  } else {
    this->plhs[idx] = (mxArray *)mxCreateNumericArray(ndim, dims,
  						      outputVoxelClassId,
 						      mxREAL);
  }
  if (this->plhs[idx] == NULL) {
    mexErrMsgTxt("MatlabExportFilter: Cannot allocate memory for output matrix");
  }
  
  // pointer to the Matlab output buffer
  TVector *buffer =  (TVector *)mxGetData(this->plhs[idx]);
  if(buffer == NULL) {
    mexErrMsgTxt("MatlabExportFilter: Cannot get pointer to allocated memory for output matrix");
  }

  // impersonate the data buffer in the filter with the Matlab output
  // buffer
  //
  // note that SetImportPointer() does not create a memory leak,
  // because at this point the output has size 0 (it has not been
  // allocated yet). After running SetImportPointer(), the output has
  // size>0, which means that the filter will see that this output's
  // memory has been allocated already, and won't need to do it itself
  if (pOutput->GetPixelContainer()->Size() != 0) {
    mexWarnMsgTxt("Memory leak, ITK output has size>0 before grafting it onto Matlab");
  }
  const bool filterWillDeleteTheBuffer = false;
  pOutput->GetPixelContainer()->SetImportPointer(buffer,
  						 mxGetNumberOfElements(this->plhs[idx]),
  						 filterWillDeleteTheBuffer);
  
  return;

}


// function to allocate memory in Matlab and copy an ITK filter
// output to this buffer. In principle, it's better to use
// GraftItkImageOntoMatlab() than CopyItkImageToMatlab(), because
// then ITK and Matlab share the same memory buffer, but the former
// approach does not work with some filters
//
// size is a vector with the dimensions of the output image in
// Matlab. For example, for a 256x200x512 image, size = {256, 200, 512}
//
// size is the same for vector or scalar images. 
template <class TPixel, unsigned int VectorDimension, class TVector>
void
MatlabExportFilter::CopyItkImageToMatlab(typename itk::DataObject::Pointer image,
					 std::vector<unsigned int> size,
					 int idx, std::string paramName) {

  if (size.size() != VectorDimension) {
    mexErrMsgTxt(("MatlabExportFilter: Output image " + paramName 
		  + ": wrong length for size vector").c_str());
  }

  // pointer to the filter output
  typedef typename itk::Image<TVector, VectorDimension> OutputImageType;
  OutputImageType *pOutput = 
    dynamic_cast<OutputImageType *>(image.GetPointer());
  if (pOutput == NULL) {
    mexErrMsgTxt("MatlabExportFilter: Cannot get pointer to filter output");
  }

  // if the output image is a vector image, extend the size vector
  if (!TypesAreEqual<TPixel, TVector>::value) {
    size.insert(size.begin(), VectorDimension);
  }

  // get the Matlab class ID for the element type we need
  mxClassID outputVoxelClassId = convertCppDataTypeToMatlabCassId<TPixel>();

  mwSize ndim = size.size();
  mwSize *dims;
  dims = (mwSize *)malloc(ndim * sizeof(mwSize));
  if (dims == NULL) {
	  mexErrMsgTxt("MatlabExportFilter: Cannot allocate memory for dims vector");
  }
  bool isEmptyMatrix = false;
  for (mwSize i = 0; i < ndim; ++i) {
    dims[i] = size[i];
    isEmptyMatrix |= (size[i] == 0);
  }

  // create output matrix for Matlab's result
  if (isEmptyMatrix) {
    this->plhs[idx] = mxCreateDoubleMatrix(0, 0, mxREAL);
  } else {
    this->plhs[idx] = (mxArray *)mxCreateNumericArray(ndim, dims,
  						      outputVoxelClassId,
 						      mxREAL);
  }
  if (this->plhs[idx] == NULL) {
    mexErrMsgTxt("MatlabExportFilter: Cannot allocate memory for output matrix");
  }
  
  // pointer to the Matlab output buffer
  TVector *buffer =  (TVector *)mxGetData(this->plhs[idx]);
  if(buffer == NULL) {
    mexErrMsgTxt("MatlabExportFilter: Cannot get pointer to allocated memory for output matrix");
  }

  // copy ITK filter output to Matlab buffer
  // pOutput->GetBufferPointer()
  typedef typename itk::ImageRegionConstIterator<OutputImageType> IteratorType;
  IteratorType it(pOutput, pOutput->GetLargestPossibleRegion());

  it.GoToBegin();
  while(!it.IsAtEnd()) {
    *buffer = it.Get();
    ++it;
    ++buffer;
  }

  return;
}

// function to allocate memory on the Matlab side and copy a vector
// of scalars from the C++ side.
template <class TPixel, class TVector>
void 
MatlabExportFilter::CopyVectorOfScalarsToMatlab(MatlabExportFilter::MatlabOutputPointer output,
						TVector v, mwSize size) {

  // if we are asked to copy the data to an output argument that the
  // user has not requested, we avoid wasting time and memory, and
  // simply exit. The only exception is for the first output argument,
  // plhs[0]. Even if the user has not requested any output arguments,
  // e.g. my_mex_function([1 4.4 2]); we need to allocate memory for
  // an empty matrix at plhs[0], otherwise Matlab will give an "One or
  // more output arguments not assigned during call to..." error
  if (!output->isRequested) {
    if (output->isTopLevelFirst) {
      this->CopyEmptyArrayToMatlab(output);
      return;
    } else {
      return;
    }
  }

  // the output is requested, so we are going to create it

  // get the Matlab class ID for the element type we need
  mxClassID outputVoxelClassId = convertCppDataTypeToMatlabCassId<TPixel>();

  // create output matrix for Matlab's result
  mwSize ndim = 2;
  mwSize dims[2] = {size, 1};
  *output->pm = (mxArray *)mxCreateNumericArray(ndim, dims,
						outputVoxelClassId,
						mxREAL);
  if (*output->pm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:CopyVectorOfScalarsToMatlab:MemoryAllocation", 
		      ("Cannot allocate memory for output " + output->name).c_str());
  }
  
  // pointer to the Matlab output buffer
  TPixel *buffer =  (TPixel *)mxGetData(*output->pm);
  if(buffer == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:CopyVectorOfScalarsToMatlab:MemoryAccess", 
		      ("Cannot get pointer to allocated memory for output matrix" + output->name).c_str());
  }

  // copy vector to the output
  for (mwIndex row = 0; row < size; ++row) {
    buffer[row] = v[row];
  }

}

// function to allocate memory on the Matlab side and copy a vector
// of vectors from the C++ side.
template <class TPixel, class TInsideVector, class TOutsideVector>
void 
MatlabExportFilter::CopyVectorOfVectorsToMatlab(MatlabExportFilter::MatlabOutputPointer output,
						TOutsideVector v, 
						mwSize outsideSize, mwSize insideSize) {

  // if we are asked to copy the data to an output argument that the
  // user has not requested, we avoid wasting time and memory, and
  // simply exit. The only exception is for the first output argument,
  // plhs[0]. Even if the user has not requested any output arguments,
  // e.g. my_mex_function([1 4.4 2]); we need to allocate memory for
  // an empty matrix at plhs[0], otherwise Matlab will give an "One or
  // more output arguments not assigned during call to..." error
  if (!output->isRequested) {
    if (output->isTopLevelFirst) {
      this->CopyEmptyArrayToMatlab(output);
      return;
    } else {
      return;
    }
  }

  // the output is requested, so we are going to create it

  // get the Matlab class ID for the element type we need
  mxClassID outputVoxelClassId = convertCppDataTypeToMatlabCassId<TPixel>();

  // allocate memory for the output
  mwSize ndim = 2;
  mwSize dims[2] = {outsideSize, insideSize};
  *output->pm = (mxArray *)mxCreateNumericArray(ndim, dims,
						outputVoxelClassId,
						mxREAL);

  if (*output->pm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:MemoryAllocation", 
		      "Cannot allocate memory for output matrix");
  }
  
  // pointer to the Matlab output buffer
  TPixel *buffer =  (TPixel *)mxGetData(*output->pm);
  if(buffer == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:MemoryAccess", 
		      "Cannot get pointer to allocated memory for output matrix");
  }

  // copy vector to the output
  for (mwIndex row = 0; row < outsideSize; ++row) {
    for (mwIndex col = 0; col < insideSize; ++col) {
      buffer[col * outsideSize + row] = v[row][col];
    }
  }

}

#endif /* MATLABEXPORTFILTER_HXX */
