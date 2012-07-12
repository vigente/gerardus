/*
 * MatlabExportFilter.hxx
 *
 * Class to provide an interface to graft ITK images onto Matlab
 * outputs
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012 University of Oxford
  * Version: 0.1.0
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

/* ITK headers */

/* Gerardus headers */
#include "MatlabExportFilter.h"

// constructor
MatlabExportFilter::MatlabExportFilter() {
  this->args = NULL;
  this->numArgs = 0;
}

// destructor
// this class is only an interface to handle pointer, so we must not
// attempt to delete the argument list provided by Matlab
MatlabExportFilter::~MatlabExportFilter() {}

// function to check that number of input arguments is within
// certain limits
void MatlabExportFilter::CheckNumberOfArguments(unsigned int min, unsigned int max) {
  if (this->numArgs < min) {
    mexErrMsgTxt("Not enough output arguments");
  }
  if (this->numArgs > max) {
    mexErrMsgTxt("Too many output arguments");
  }
}

// function to allocate memory in Matlab and hijack it to be used as
// an ITK filter output
//
// size is a vector with the dimensions of the output image in
// Matlab. For example, for a 256x200x512 image, size = {256, 200, 512}
//
// size is the same for vector or scalar images. 
template <class TPixel, unsigned int VectorDimension, class TVector>
void
MatlabExportFilter::GraftItkImageOntoMatlab(typename itk::DataObject::Pointer image,
					    std::vector<unsigned int> size,
					    unsigned int idx, std::string paramName) {

  if (size.size() != VectorDimension) {
    mexErrMsgTxt(("Output image " + paramName + ": wrong length for size vector").c_str());
  }

  // convert output data type to output class ID
  mxClassID outputVoxelClassId = mxUNKNOWN_CLASS;
  if (TypeIsBool<TPixel>::value) {
    outputVoxelClassId = mxLOGICAL_CLASS;
  } else if (TypeIsUint8<TPixel>::value) {
    outputVoxelClassId = mxUINT8_CLASS;
  } else if (TypeIsInt8<TPixel>::value) {
    outputVoxelClassId = mxINT8_CLASS;
  } else if (TypeIsUint16<TPixel>::value) {
    outputVoxelClassId = mxUINT16_CLASS;
  } else if (TypeIsInt16<TPixel>::value) {
    outputVoxelClassId = mxINT16_CLASS;
  } else if (TypeIsInt32<TPixel>::value) {
    outputVoxelClassId = mxINT32_CLASS;
  } else if (TypeIsInt64<TPixel>::value) {
    outputVoxelClassId = mxINT64_CLASS;
    // a signed long can be 32- or 64-bit depending on the
    // architecture. 64-bit is the safest option
  } else if (TypeIsSignedLong<TPixel>::value) {
    outputVoxelClassId = mxINT64_CLASS;
  } else if (TypeIsFloat<TPixel>::value) {
    outputVoxelClassId = mxSINGLE_CLASS;
  } else if (TypeIsDouble<TPixel>::value) {
    outputVoxelClassId = mxDOUBLE_CLASS;
  } else {
    mexErrMsgTxt("Assertion fail: Unrecognised output data type");
  }

  // dimensions for the output array
  typedef typename itk::Image<TVector, VectorDimension> OutputImageType;
  OutputImageType *pOutput = 
    dynamic_cast<OutputImageType *>(image.GetPointer());
  if (pOutput == NULL) {
    mexErrMsgTxt("Cannot get pointer to filter output");
  }

  // if the output image is a vector image, extend the size vector
  if (!TypesAreEqual<TPixel, TVector>::value) {
    size.insert(size.begin(), VectorDimension);
  }

  mwSize ndim = size.size();
  mwSize dims[ndim];
  bool isEmptyMatrix = false;
  for (mwSize i = 0; i < ndim; ++i) {
    dims[i] = size[i];
    isEmptyMatrix |= (size[i] == 0);
  }

  // create output matrix for Matlab's result
  if (isEmptyMatrix) {
    this->args[idx] = mxCreateDoubleMatrix(0, 0, mxREAL);
  } else {
    this->args[idx] = (mxArray *)mxCreateNumericArray(ndim, dims,
  						      outputVoxelClassId,
 						      mxREAL);
  }
  if (this->args[idx] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output matrix");
  }
  
  // pointer to the Matlab output buffer
  TVector *imOutp =  (TVector *)mxGetData(this->args[idx]);
  if(imOutp == NULL) {
    mexErrMsgTxt("Cannot get pointer to allocated memory for output matrix");
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
  pOutput->GetPixelContainer()->SetImportPointer(imOutp,
  						 mxGetNumberOfElements(this->args[idx]),
  						 filterWillDeleteTheBuffer);
  
  return;

}

#endif /* MATLABEXPORTFILTER_HXX */
