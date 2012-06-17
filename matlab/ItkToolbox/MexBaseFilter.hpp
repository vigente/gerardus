/* 
 * MexBaseFilter.hpp
 *
 * MexBaseFilter<InVoxelType, OutVoxelType>: This is where
 * the code to actually run the filter on the image lives.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.7.4
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

#ifndef MEXBASEFILTER_HPP
#define MEXBASEFILTER_HPP

/* ITK headers */
#include "itkImportImageFilter.h"
#include "itkImageToImageFilter.h"
#include "itkVectorImage.h"

/* Gerardus headers */
#include "GerardusCommon.hpp"
#include "NrrdImage.hpp"

/* 
 * MexBaseFilter
 */
template <class InVoxelType, class OutVoxelType>
class MexBaseFilter {
private:

  typedef double TScalarType; // data type for scalars
  
protected:

  typedef typename itk::Image<InVoxelType, Dimension> InImageType;
  typedef typename itk::Image<OutVoxelType, Dimension> OutImageType;
  typedef typename itk::ImageToImageFilter<InImageType, OutImageType> ImageToImageFilterType;
  typedef itk::ImportImageFilter<InVoxelType, Dimension> ImportFilterType;
  typename ImportFilterType::Pointer importFilter;
  NrrdImage nrrd;
  int nargout;
  mxArray** argOut;
  mwSize nparam;
  const mxArray** argParam;
  typename ImageToImageFilterType::Pointer filter;

  // get the value of input parameters that are numeric scalars from
  // the array of input arguments
  template <class ParamType>
  ParamType GetScalarParamValue(std::string paramName,
				mwIndex idx, ParamType def);

  // allocate memory for a Matlab output buffer
  template <class OutputType>
  void MallocMatlabOutputBuffer(unsigned int idx, int vectorSize);

  // make the filter use a Matlab buffer for one of its outputs. This
  // is the fastest and least memory-consuming way of working with the
  // filter outputs, because results don't need to be duplicated by
  // copying them to a separate Matlab array.
  template <class OutputType>
  void GraftMatlabOutputBufferIntoItkFilterOutput(unsigned int idx);

public:

  // constructor. We assume that the first two input arguments are
  // filter type and image (common to all filters), and the rest are
  // user-provided parameters (filter-dependent), e.g. dilation radius
  MexBaseFilter(const NrrdImage &_nrrd, int _nargout, mxArray** _argOut,
		const int _nargin, const mxArray** _argIn);
  
  // destructor to deallocate any memory that the Matlab garbage
  // collector won't free automatically. The destructor has to be
  // virtual so that derived classes can add their own local memory
  // deallocation steps
  virtual ~MexBaseFilter();

  // check numer of outputs requested by the user. By default, the
  // function provides 0 or 1 output (the filtered image), but this
  // method can be overriden in child filters with more outputs
  virtual void CheckNumberOfOutputs();

  // put the pointer to the Matlab image buffer into an
  // itk::ImportImageFilter. We can feed the ImportImageFilter directly
  // to the filter. This way, we don't need to copy the buffer
  // to an itk::Image, which saves time, and avoids duplicating the
  // image in memory
  void GraftMatlabInputBufferIntoItkImportFilter();

  // filter setup code common to all filters: pass image to the
  // filter, allocate memory for the Matlab output, and graft the
  // Matlab output into the filter output
  void FilterBasicSetup();

  // by default, this method doesn't do anything, but can be overriden
  // when a child filter needs some extra setput steps (e.g. passing
  // parameters)
  virtual void FilterAdvancedSetup();

  // filter the image
  void RunFilter();

  // by default, this method doesn't do anything, but can be overriden
  // when a child filter provides other outputs apart from the
  // filtered image
  virtual void ExportOtherFilterOutputsToMatlab();

};

// BaseFilter cannot be invoked by the user, but declaring these
// static strings is necessary when we use EXCLUDEFILTER with
// derived filters
template <>
class MexBaseFilter< std::string, std::string > {
public:
  
  static const std::string longname;
  static const std::string shortname;
  
};

/*
 * Templated member functions. They need to go here in the header
 * file, because if they go in the .cpp, then they need instance
 * specialization for every <InVoxelType, OutVoxelType, ParamType>
 *
 */

// function to get the value of input parameters that are numeric
// scalars from the array of input arguments
template <class InVoxelType, class OutVoxelType>
template <class ParamType>
ParamType MexBaseFilter<InVoxelType, 
			OutVoxelType>::GetScalarParamValue(std::string paramName,
							mwIndex idx, 
							ParamType def) {
  
  // if user didn't provide a value, or provided an empty array, return the default
  if (idx >= this->nparam || mxIsEmpty(this->argParam[idx])) {
    return def;
  }

  // if user provided a value, check that it's a scalar, whether in
  // numeric or logical form
  if ((!mxIsNumeric(this->argParam[idx])
       && !mxIsLogical(this->argParam[idx]))
      || mxGetNumberOfElements(this->argParam[idx]) != 1) {
    mexErrMsgTxt((std::string("Parameter ") + paramName 
		  + std::string(" must be a scalar")).c_str());
  }
  
  // output
  ParamType value = 0;
  
  // input image type
  mxClassID inputVoxelClassId = mxGetClassID(this->argParam[idx]);

  // macro to make the code in the switch statement cleaner
#define GETVALUE(Tx, x, Ty, y)						\
  if (x) {								\
    Tx *valuep = (Tx *)mxGetData(x);					\
    y = (Ty)valuep[0];							\
  } else {								\
    mexErrMsgTxt((std::string("Parameter ") + paramName			\
		  + std::string(" couldn't be read")).c_str());		\
}
  
    // cast the class type provided by Matlab to the type requested by
    // the user
  switch(inputVoxelClassId)  { 
  case mxLOGICAL_CLASS:
    GETVALUE(mxLogical, this->argParam[idx], ParamType, value);
    break;
  case mxDOUBLE_CLASS:
    GETVALUE(double, this->argParam[idx], ParamType, value);
    break;
  case mxSINGLE_CLASS:
    GETVALUE(float, this->argParam[idx], ParamType, value);
    break;
  case mxINT8_CLASS:
    GETVALUE(int8_T, this->argParam[idx], ParamType, value);
    break;
  case mxUINT8_CLASS:
    GETVALUE(uint8_T, this->argParam[idx], ParamType, value);
    break;
  case mxINT16_CLASS:
    GETVALUE(int16_T, this->argParam[idx], ParamType, value);
    break;
  case mxUINT16_CLASS:
    GETVALUE(uint16_T, this->argParam[idx], ParamType, value);
    break;
  case mxINT32_CLASS:
    GETVALUE(int32_T, this->argParam[idx], ParamType, value);
    break;
    // case mxUINT32_CLASS:
    //   break;
  case mxINT64_CLASS:
    GETVALUE(int64_T, this->argParam[idx], ParamType, value);
    break;
    // case mxUINT64_CLASS:
    //   break;
  case mxUNKNOWN_CLASS:
    mexErrMsgTxt("Input parameter has unknown type.");
    break;
  default:
    mexErrMsgTxt("Input parameter has invalid type.");
      break;
  }

#undef GETVALUE
  
  return value;
}


// allocate memory for a Matlab output buffer
template <class InVoxelType, class OutVoxelType>
template <class OutputType>
void MexBaseFilter<InVoxelType, 
		   OutVoxelType>::MallocMatlabOutputBuffer(unsigned int idx, int vectorSize=1) {

  // convert output data type to output class ID
  mxClassID outputVoxelClassId = mxUNKNOWN_CLASS;
  if (TypeIsBool<OutputType>::value) {
    outputVoxelClassId = mxLOGICAL_CLASS;
  } else if (TypeIsUint8<OutputType>::value) {
    outputVoxelClassId = mxUINT8_CLASS;
  } else if (TypeIsInt8<OutputType>::value) {
    outputVoxelClassId = mxINT8_CLASS;
  } else if (TypeIsUint16<OutputType>::value) {
    outputVoxelClassId = mxUINT16_CLASS;
  } else if (TypeIsInt16<OutputType>::value) {
    outputVoxelClassId = mxINT16_CLASS;
  } else if (TypeIsInt32<OutputType>::value) {
    outputVoxelClassId = mxINT32_CLASS;
  } else if (TypeIsInt64<OutputType>::value) {
    outputVoxelClassId = mxINT64_CLASS;
    // a signed long can be 32- or 64-bit depending on the
    // architecture. 64-bit is the safest option
  } else if (TypeIsSignedLong<OutputType>::value) {
    outputVoxelClassId = mxINT64_CLASS;
  } else if (TypeIsFloat<OutputType>::value) {
    outputVoxelClassId = mxSINGLE_CLASS;
  } else if (TypeIsDouble<OutputType>::value) {
    outputVoxelClassId = mxDOUBLE_CLASS;
  } else {
    mexErrMsgTxt("Assertion fail: Unrecognised output data type");
  }

  // dimensions for the output array. For vector images, we have an
  // extra dimension apart from the image size itself
  mwSize dims[this->nrrd.getNdim() + 1]; // allocate enough space even
                                         // if vectorSize==1 and we don't use it
  mwSize ndim;
  if (vectorSize > 1) {
    dims[0] = vectorSize;
    for (mwSize i = 0; i < this->nrrd.getNdim(); ++i) {
      dims[i+1] = this->nrrd.getDims()[i];
    }
    ndim = this->nrrd.getNdim() + 1;
  } else{
    for (mwSize i = 0; i < this->nrrd.getNdim(); ++i) {
      dims[i] = this->nrrd.getDims()[i];
    }
    ndim = this->nrrd.getNdim();
  }

  // create output matrix for Matlab's result
  if (this->nrrd.getR() == 0 || this->nrrd.getC() == 0) {
    this->argOut[idx] = mxCreateDoubleMatrix(0, 0, mxREAL);
  } else {
    this->argOut[idx] = (mxArray *)mxCreateNumericArray(ndim, 
							dims,
							outputVoxelClassId,
							mxREAL);
  }
  if (this->argOut[0] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output matrix");
  }

  return;

}

// filter setup code common to all filters: pass image to the
// filter, allocate memory for the Matlab output, and graft the
// Matlab output into the filter output
template <class InVoxelType, class OutVoxelType>
void MexBaseFilter<InVoxelType, OutVoxelType>::FilterBasicSetup() {

  // pass input image to filter
  this->filter->SetInput(this->importFilter->GetOutput());

  // link the filter main output (the filtered image) to the first
  // Matlab output buffer
  this->template MallocMatlabOutputBuffer<OutVoxelType>(0);
  this->template GraftMatlabOutputBufferIntoItkFilterOutput<OutVoxelType>(0);


}

// function to make the filter use a Matlab buffer for one of its
// outputs. The OutpuType is the type of the elements in that output,
// that may be different from the OutVoxelType
template <class InVoxelType, class OutVoxelType>
template <class OutputType>
void MexBaseFilter<InVoxelType, 
		   OutVoxelType>::GraftMatlabOutputBufferIntoItkFilterOutput(unsigned int idx) {

  // pointer to the Matlab output buffer
  OutputType *imOutp =  (OutputType *)mxGetData(this->argOut[idx]);
  if (imOutp == NULL) {
    mexErrMsgTxt("Output buffer memory has not been allocated before trying to graft it");
  }

  // impersonate the data buffer in the filter with the Matlab output
  // buffer
  const bool filterWillDeleteTheBuffer = false;
  typedef typename itk::Image<OutputType, Dimension> ScalarImageType;
  ScalarImageType *ps;

  ps = dynamic_cast<ScalarImageType *>(this->filter->GetOutputs()[idx].GetPointer());
  if (ps == NULL) {
    mexErrMsgTxt("Cannot get pointer to filter scalar output");
  }
  ps->GetPixelContainer()->SetImportPointer(imOutp,
					    mxGetNumberOfElements(this->nrrd.getData()),
					    filterWillDeleteTheBuffer);
  
  return;

}



#endif /* MEXBASEFILTER_HPP */
