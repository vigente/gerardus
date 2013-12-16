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
  * Version: 0.6.1
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

// get number of elements in the list of plhs arguments
int MatlabExportFilter::GetNumberOfArguments() {
  return this->nlhs;
}

// function to check that number of plhs arguments is within
// certain limits
void MatlabExportFilter::CheckNumberOfArguments(int min, int max) {
  if (this->nlhs < min) {
    mexErrMsgTxt("Not enough output arguments");
  }
  if (this->nlhs > max) {
    mexErrMsgTxt("Too many output arguments");
  }
}

// function to import into this class the array with the arguments
// provided by Matlab
void MatlabExportFilter::ConnectToMatlabFunctionOutput(int _nlhs, mxArray *_plhs[]) {
  this->nlhs = _nlhs;
  this->plhs = _plhs;
}

// Functions to register an output at the export filter. 
MatlabExportFilter::MatlabOutputPointer
MatlabExportFilter::RegisterOutput(int pos, std::string name) {
  
  // assign main variables of the output: its name and the memory
  // address it should be stored at
  MatlabExportFilter::MatlabOutput out;
  out.ppm = &(this->plhs[pos]);
  out.name = name;

  // outputs that are regular output variables, e.g. y = myfun(x),
  // are requested if the user put them in the output vector. For
  // example:
  // [y,z] = myfun(x); // y and z have been requested
  // y = myfun(x);     // y has been requested, z has not been requested
  out.isRequested = (pos < this->nlhs);

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
  if (*output->ppm != NULL) {
    mxDestroyArray(*output->ppm);
  }
  // create output matrix for Matlab's result
  mwSize ndim = size.size();
  if (ndim > 0) { // non-empty matrix
    mwSize *dims = new mwSize[ndim];
    for (size_t i = 0; i < ndim; ++i) {
      dims[i] = size[i];
    }

    *output->ppm = (mxArray *)mxCreateNumericArray(ndim, dims, outputClassId, mxREAL);
  } else { // empty matrix
	*output->ppm = (mxArray *)mxCreateDoubleMatrix(0, 0, mxREAL);
  }
  if (*output->ppm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:MemoryAllocation", 
		      ("Cannot allocate memory for output " + output->name).c_str());
  }
  
  // pointer to the Matlab output buffer
  TData *buffer =  (TData *)mxGetData(*output->ppm);
  if(buffer == NULL && !mxIsEmpty(*output->ppm)) {
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
  mxArray *cell = mxGetCell(*output->ppm, pos);
  if (cell != NULL) {
    mxDestroyArray(cell);
  }

  // allocate memory for the new array
  mwSize ndim = size.size();
  if (ndim > 0) {
	mwSize *dims = new mwSize[ndim];
    for (size_t i = 0; i < ndim; ++i) {
      dims[i] = size[i];
    }

    cell = mxCreateNumericArray(ndim, dims, outputClassId, mxREAL);
  } else {
	cell = mxCreateDoubleMatrix(0, 0, mxREAL);
  }
  if (cell == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:MemoryAllocation", 
		      ("Cannot allocate memory for output " + output->name + "{" 
		       + boost::lexical_cast<std::string>(pos) + "}").c_str());
  }

  // place the new array into the cell array
  mxSetCell(*output->ppm, pos, cell);

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
  
  // only copy the output argument if the user has requested it at the
  // output, to avoid wasting time otherwise
  if (output->isRequested) {
    *output->ppm = mxCreateDoubleMatrix(0, 0, mxREAL);
  }

}

// Function to allocate memory in Matlab and hijack it to be used as
// an ITK filter output. Syntax for images with vector voxels
template <class TPixel, unsigned int VectorDimension, class TVector>
void
MatlabExportFilter::GraftItkImageOntoMatlab(MatlabOutputPointer output, 
					    typename itk::DataObject::Pointer image,
					    std::vector<mwSize> size) {

  // if the user hasn't requested the output from Matlab, then we
  // don't need to graft the ITK image onto Matlab
  if (!output->isRequested) {
    return;
  }

  // check that the number of dimensions is the same in the ITK image
  // and in the Matlab image
  if (size.size() != VectorDimension) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:GraftItkImageOntoMatlab:WrongDataFormat", 
		      ("Number of dimensions in output image" + output->name
		       + " and image to graft disagree").c_str());
  }

  // dimensions for the output array
  typedef typename itk::Image<TVector, VectorDimension> OutputImageType;
  OutputImageType *pOutput = 
    dynamic_cast<OutputImageType *>(image.GetPointer());
  if (pOutput == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:GraftItkImageOntoMatlab:MemoryAccess", 
  		      ("Cannot get pointer to allocated memory for output image" + output->name).c_str());
  }

  // size tells us the size of the image. But if we have vector
  // voxels, the Matlab output will need to be 4-D
  if (!TypesAreEqual<TPixel, TVector>::value) {
    size.insert(size.begin(), VectorDimension);
  }

  // allocate memory for the 3-D or 4-D image in Matlab, and get a pointer to the buffer
  //
  // this is a rather ugly and dangerous hack. But
  // AllocateNDArrayInMatlab() will not accept a generalised TVector,
  // because it only accepts basic types in the template (int, bool,
  // etc). At the same time, SetImportPointer() below requires a
  // TVector pointer, otherwise it gives a compilation problem. We
  // solve this assuming that the TVector is simply a concatenation of
  // TPixel, so we can do a reinterpret_cast(), but this can cause
  // leaks and segfaults if this is not the case.
  TPixel *buffer =  this->AllocateNDArrayInMatlab<TPixel>(output, size);
  TVector *buffer2 = reinterpret_cast<TVector *>(buffer);

  // impersonate the data buffer in the filter with the Matlab output
  // buffer
  //
  // note that SetImportPointer() does not create a memory leak,
  // because at this point the output has size 0 (it has not been
  // allocated yet). After running SetImportPointer(), the output has
  // size>0, which means that the filter will see that this output's
  // memory has been allocated already, and won't need to do it itself
  if (pOutput->GetPixelContainer()->Size() != 0) {
    mexWarnMsgIdAndTxt("Gerardus:MatlabExportFilter:GraftItkImageOntoMatlab:MemoryLeak",
		       ("Memory leak, ITK output has size>0 before grafting it onto Matlab output " 
			+ output->name).c_str());
  }
  const bool filterWillDeleteTheBuffer = false;
  pOutput->GetPixelContainer()->SetImportPointer(buffer2,
  						 mxGetNumberOfElements(*output->ppm),
  						 filterWillDeleteTheBuffer);
  
}

// Function to allocate memory in Matlab and hijack it to be used as
// an ITK filter output. Syntax for images with plain scalar voxels
// (as opposed to vector voxels)
template <class TPixel, unsigned int VectorDimension>
void 
MatlabExportFilter::GraftItkImageOntoMatlab(MatlabOutputPointer output, 
					    typename itk::DataObject::Pointer image, 
					    std::vector<mwSize> size) {
  GraftItkImageOntoMatlab<TPixel, VectorDimension, TPixel>(output, image, size);
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
template <class TPixel, unsigned int VectorDimension>
void 
MatlabExportFilter::CopyItkImageToMatlab(MatlabOutputPointer output, 
					 typename itk::DataObject::Pointer image, 
					 std::vector<mwSize> size) {
  CopyItkImageToMatlab<TPixel, VectorDimension, TPixel>(output, image, size);
}

template <class TPixel, unsigned int VectorDimension, class TVector>
void
MatlabExportFilter::CopyItkImageToMatlab(MatlabOutputPointer output, 
					 typename itk::DataObject::Pointer image, 
					 std::vector<mwSize> size) {

  // if we are asked to copy the data to an output argument that the
  // user has not requested, we avoid wasting time and memory, and
  // simply exit
  if (!output->isRequested) {
    return;
  }

  // check that the number of dimensions is the same in the ITK image
  // and in the Matlab image
  if (size.size() != VectorDimension) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:GraftItkImageOntoMatlab:WrongDataFormat", 
		      ("Number of dimensions in output image" + output->name
		       + " and image to graft disagree").c_str());
  }

  // pointer to the filter output
  typedef typename itk::Image<TVector, VectorDimension> OutputImageType;
  OutputImageType *pOutput = 
    dynamic_cast<OutputImageType *>(image.GetPointer());
  if (pOutput == NULL) {
    mexErrMsgTxt("MatlabExportFilter: Cannot get pointer to filter output");
  }

  // size tells us the size of the image. But if we have vector
  // voxels, the Matlab output will need to be 4-D
  if (!TypesAreEqual<TPixel, TVector>::value) {
    size.insert(size.begin(), VectorDimension);
  }

  // allocate memory for the 3-D or 4-D image in Matlab, and get a pointer to the buffer
  TVector *buffer =  this->AllocateNDArrayInMatlab<TVector>(output, size);

  // copy ITK filter output to Matlab buffer
  typedef typename itk::ImageRegionConstIterator<OutputImageType> IteratorType;
  IteratorType it(pOutput, pOutput->GetLargestPossibleRegion());

  it.GoToBegin();
  while(!it.IsAtEnd()) {
    *buffer = it.Get();
    ++it;
    ++buffer;
  }

}

// function to allocate memory on the Matlab side and copy a vector
// of scalars from the C++ side.
template <class TPixel, class TVector>
void 
MatlabExportFilter::CopyVectorOfScalarsToMatlab(MatlabExportFilter::MatlabOutputPointer output,
						TVector v, mwSize size) {

  // if we are asked to copy the data to an output argument that the
  // user has not requested, we avoid wasting time and memory, and
  // simply exit
  if (!output->isRequested) {
    return;
  }

  // get the Matlab class ID for the element type we need
  mxClassID outputVoxelClassId = convertCppDataTypeToMatlabCassId<TPixel>();

  // create output matrix for Matlab's result
  mwSize ndim = 2;
  mwSize dims[2] = {size, 1};
  *output->ppm = (mxArray *)mxCreateNumericArray(ndim, dims,
						 outputVoxelClassId,
						 mxREAL);
  if (*output->ppm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:CopyVectorOfScalarsToMatlab:MemoryAllocation", 
		      ("Cannot allocate memory for output " + output->name).c_str());
  }
  
  // pointer to the Matlab output buffer
  TPixel *buffer =  (TPixel *)mxGetData(*output->ppm);
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
template <class TPixel, class TOutsideVector>
void 
MatlabExportFilter::CopyVectorOfVectorsToMatlab(MatlabExportFilter::MatlabOutputPointer output,
						TOutsideVector v, 
						mwSize outsideSize, mwSize insideSize) {

  // if we are asked to copy the data to an output argument that the
  // user has not requested, we avoid wasting time and memory, and
  // simply exit
  if (!output->isRequested) {
    return;
  }

  // get the Matlab class ID for the element type we need
  mxClassID outputVoxelClassId = convertCppDataTypeToMatlabCassId<TPixel>();

  // allocate memory for the output
  mwSize ndim = 2;
  mwSize dims[2] = {outsideSize, insideSize};
  *output->ppm = (mxArray *)mxCreateNumericArray(ndim, dims,
						 outputVoxelClassId,
						 mxREAL);

  if (*output->ppm == NULL) {
    mexErrMsgIdAndTxt("Gerardus:MatlabExportFilter:MemoryAllocation", 
		      "Cannot allocate memory for output matrix");
  }
  
  // pointer to the Matlab output buffer
  TPixel *buffer =  (TPixel *)mxGetData(*output->ppm);
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
