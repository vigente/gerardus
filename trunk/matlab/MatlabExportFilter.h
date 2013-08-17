/*
 * MatlabExportFilter.h
 *
 * Class to provide an interface to pass C++ data to Matlab
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

#ifndef MATLABEXPORTFILTER_H
#define MATLABEXPORTFILTER_H

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <vector>

/* ITK headers */
#include "itkSmartPointer.h"

/* Gerardus headers */
#include "GerardusCommon.h"

// class to interface with the Matlab MEX function outputs
class MatlabExportFilter: public itk::Object {

 public:

  // struct to encapsulate each of the outputs to Matlab
  struct MatlabOutput {
    mxArray **pm;     // pointer to the Matlab output (itself a pointer)
    std::string name; // name of the output for error/debug messages
    bool isRequested; // flag: has the user requested this output?
    bool isTopLevel;  // flag: is this output in the plhs array?
    bool isTopLevelFirst; // flag: is this the first top level output?

    // note: if the output is top level, we need to know whether it's
    // the first one, because in that case we need to allocate an
    // empty matrix in Matlab for it, even if the user has not
    // requested it
  };
  
  typedef std::list<MatlabOutput>::iterator MatlabOutputPointer;

 private:
  
  // these are the variables provided by the MEX API for the output of
  // the function
  int             nlhs;
  mxArray         **plhs;

  // list of the outputs registered at this exporter
  std::list<MatlabOutput> outputsList;
  
 protected:
  
  MatlabExportFilter();
  ~MatlabExportFilter();
  
 public:
  
  // standard class typedefs
  typedef MatlabExportFilter                Self;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;

  // method for creation through the object factory
  itkNewMacro(Self);

  // run-time type information (and related methods)
  itkTypeMacro(MatlabExportFilter, Object);

  // function to import into this class the array with the arguments
  // provided by Matlab
  void ConnectToMatlabFunctionOutput(int _nlhs, mxArray *_plhs[]) {
    this->nlhs = _nlhs;
    this->plhs = _plhs;
  }

  // get number of elements in the list of arguments
  int GetNumberOfOutputArguments();

  // function to check that number of input arguments is within
  // certain limits
  void CheckNumberOfArguments(int min, int max);

  // Function to register an output at the export filter. 
  //
  // Registration basically means "this output in Matlab is going to
  // correspond to X". Once an output has been registered, it can be
  // passed to the Copy or Graft methods to actually allocate the
  // memory, connect to a filter, copy the data over, etc.
  //
  // We provide two syntaxes. The first one registers the output on
  // the Matlab output array. The second one registers the output on
  // any array (e.g. an output argument that is a cell array)
  //
  // pos: position index within the base array
  // base: base array
  //
  // returns: a class of type MatlabOutput, defined above
  MatlabOutputPointer RegisterOutput(int pos, std::string name);
  MatlabOutputPointer RegisterOutput(mxArray **base, int pos, std::string name);

  // Functions to allocate memory for vectors, matrices and
  // N-dimensional arrays in Matlab, and get the data pointer back. 
  //
  // This is basically mxCreateNumericArray() + mxGetData() + some checks
  //
  // These functions are necessary because not always are we going to
  // have C++ vectors that we can copy to Matlab. Sometimes, we have
  // to generate the values of the output vector one by one.
  //
  // output: pointer to a registered output
  //
  // len:    number of elements in the vector
  //
  // nrows:  number of rows of allocated matrix
  //
  // ncols:  number of columns of allocated matrix
  //
  // size:   vector with number of rows, cols, slices, etc of allocated N-dimensional array
  //
  // returns: pointer to the data buffer. If the allocated array has
  //          size zero, the returned pointer is NULL
  template<class TData>
    TData *AllocateColumnVectorInMatlab(MatlabOutputPointer output, mwSize len);
  template<class TData>
    TData *AllocateRowVectorInMatlab(MatlabOutputPointer output, mwSize len);
  template<class TData>
    TData *AllocateMatrixInMatlab(MatlabOutputPointer output, mwSize nrows, mwSize ncols);
  template<class TData>
    TData *AllocateNDArrayInMatlab(MatlabOutputPointer output, std::vector<mwSize> size);

  // Functions to allocate memory for vectors, matrices and
  // N-dimensional arrays within a cell of a cell array in Matlab, and
  // get the data pointer back.
  //
  // output:  pointer a registered output with the cell array
  //
  // pos:     index of the cell within the cell array
  //
  // len:     number of elements in the vector
  //
  // nrows:   number of rows of allocated matrix
  //
  // ncols:   number of columns of allocated matrix
  //
  // size:    vector with number of rows, cols, slices, etc of allocated N-dimensional array
  //
  // returns: pointer to the data buffer. If the allocated array has
  //          size zero, the returned pointer is NULL
  template<class TData>
    TData *AllocateColumnVectorInCellInMatlab(MatlabOutputPointer output, int pos, mwSize len);
  template<class TData>
    TData *AllocateRowVectorInCellInMatlab(MatlabOutputPointer output, int pos, mwSize len);
  template<class TData>
    TData *AllocateMatrixInCellInMatlab(MatlabOutputPointer output, int pos, mwSize nrows, mwSize ncols);
  template<class TData>
    TData *AllocateNDArrayInCellInMatlab(MatlabOutputPointer output, int pos, std::vector<mwSize> size);

  // Function to create an empty output in Matlab.
  void CopyEmptyArrayToMatlab(MatlabOutputPointer output);

  // Function to allocate memory in Matlab and hijack it to be used as
  // an ITK filter output. This only works with those filters that
  // create their own output. For other filters, you cannot graft the
  // output; instead, use CopyItkImageOntoMatlab() after running the
  // filter
  //
  // size is a vector with the dimensions of the output image in
  // Matlab. For example, for a 256x200x512 image, size = {256, 200, 512}
  //
  // size is the same for vector or scalar images.
  //
  // There are two syntaxes:
  //
  // * Syntax for images with plain scalar voxels (as opposed to vector voxels)
  //
  // * Syntax for images with vector voxels
  template <class TPixel, unsigned int VectorDimension>
    void GraftItkImageOntoMatlab(MatlabOutputPointer output, 
				 typename itk::DataObject::Pointer image, std::vector<mwSize> size);

  template <class TPixel, unsigned int VectorDimension, class TVector>
    void GraftItkImageOntoMatlab(MatlabOutputPointer output, 
				 typename itk::DataObject::Pointer image, std::vector<mwSize> size);

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
    void CopyItkImageToMatlab(MatlabOutputPointer output, 
			      typename itk::DataObject::Pointer image, std::vector<mwSize> size);

  template <class TPixel, unsigned int VectorDimension>
    void CopyItkImageToMatlab(MatlabOutputPointer output, 
			      typename itk::DataObject::Pointer image, std::vector<mwSize> size);

  // function to allocate memory on the Matlab side and copy a vector
  // of scalars from the C++ side. It can export any C++ class that
  // accepts operator[], e.g. v[i].
  //
  // output: pointer to a MatlabOutput that has been registered with the exporter.
  //
  // v:    vector we want to export to Matlab
  //
  // size: length of the v vector.
  //
  // Having size as a parameter here saves us from having to
  // do a specialization of this function for every class, as
  // different classes provide their size with different methods,
  // e.g. v.size(), v.Size(), v.length(), etc.
  template <class TPixel, class TVector>
    void CopyVectorOfScalarsToMatlab(MatlabExportFilter::MatlabOutputPointer output,
				     TVector v, mwSize size);

  // function to allocate memory on the Matlab side and copy a vector
  // of vectors from the C++ side. Both vectors must accept
  // operator[], e.g. w = v[i], w[j]. The inside vectors will be
  // stacked by rows. E.g., if you have a vector of Points, each point
  // will be copied to a row.
  //
  // output: pointer to a MatlabOutput that has been registered with the exporter.
  //
  // v:     vector of vectors we want to export to Matlab. This should
  //        generally be an std::vector or any other class that implements
  //        operator[]. Examples on how to get vectors:
  //        - ITK:
  //          Coordinates of points:
  //          const std::vector<PointType> x = registration->GetMovingPointSet()->GetPoints()->CastToSTLConstContainer();
  //          Data associated to points (e.g. intensity, colour):
  //          registration->GetMovingPointSet()->GetDataPoints()->CastToSTLConstContainer();
  //
  // outsideSize: length of the outside vector = Number of rows in the Matlab output
  //
  // insideSize:  length of the inside vector = Number of cols in the Matlab output
  //
  // Having these lengths as parameters here saves us from having to
  // do a specialization of this function for every class, as
  // different classes provide their size with different methods,
  // e.g. v.size(), v.Size(), v.length(), etc.
  template <class TPixel, class TInsideVector, class TOutsideVector>
    void CopyVectorOfVectorsToMatlab(MatlabExportFilter::MatlabOutputPointer output,
				     TOutsideVector v, 
				     mwSize outsideSize, mwSize insideSize);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "MatlabExportFilter.hxx"
#endif
#endif /* MATLABEXPORTFILTER_H */
