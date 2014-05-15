/*
 * MatlabExportFilter.h
 *
 * Class to provide an interface to pass C++ data to Matlab
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012-2013 University of Oxford
  * Version: 0.5.0
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
    // it ain't pretty having an mxArray** here, but having just an
    // mxArray* won't work. If we have "mxArray *pm" and do
    //   out->pm = plhs[3];
    //   out->pm = mxCreateNumericArray(...);
    // this will allocate memory but not put it in the Matlab output.
    // For that, we need an "mxArray **ppm" and do
    //   out->ppm = &(plhs[3]);
    //   *out->ppm = mxCreateNumericArray(...);
    // For outputs that don't go in the plhs array, e.g. outputs that go in a cell or struct field, 
    mxArray **ppm;    // pointer to Matlab MEX output
    std::string name; // name of the output for error/debug messages
    bool isRequested; // flag: has the user requested this output?
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
  void ConnectToMatlabFunctionOutput(int _nlhs, mxArray *_plhs[]);

  // get number of elements in the list of plhs arguments
  int GetNumberOfArguments();

  // function to check that number of plhs arguments is within
  // certain limits
  void CheckNumberOfArguments(int min, int max);

  // Function to register an output at the export filter. 
  //
  // Registration basically means "this output in Matlab is going to
  // be called X". Once an output has been registered, it can be
  // passed to the Copy or Graft methods to actually allocate the
  // memory, connect to a filter, copy the data over, etc.
  //
  // pos: position index within the base plhs array
  //
  // returns: a class of type MatlabOutput, defined above
  MatlabOutputPointer RegisterOutput(int pos, std::string name);

  // TODO: RegisterCellInOutput(MatlabOutputPointer output, int pos, std::string name)
  // TODO: RegisterStructFieldInOutput(MatlabOutputPointer output, std::string field)
  // So that we can have outputs nested into other outputs. Nested
  // outputs will inherit the isRequested status of their parent.
  //
  // Note that we already have Allocate*InCellInMatlab() methods
  // below, but we do this directly, without previous registration.

  // Functions to allocate memory for vectors, matrices and
  // N-dimensional arrays in Matlab, and get the data pointer back. 
  //
  // This is basically mxCreateNumericArray() + mxGetData() + some checks
  //
  // These functions are necessary because not always are we going to
  // have C++ vectors that we can copy to Matlab. Sometimes, we have
  // to generate the values of the output vector one by one.
  //
  // Note: Memory is always allocated, regardless of whether the user
  // has requested the output or not. It's up to the main program to
  // check before allocating, e.g.
  //   if (outYI->isRequested) {
  //     double *yi = matlabExport->AllocateColumnVectorInMatlab<TScalarType>(outYI, 10);
  //     for (int i=0; i<10; ++i) {
  //       yi[i] = i;
  //     }
  //   }
  // Why? Because the alternative is to return NULL if the user has
  // not requested the output. Then the main program has to check
  // anyway for a NULL pointer instead of isRequested, and this is
  // prone to create segfaults. If we always allocate the memory, and
  // the user forgets to check isRequested, in the worst case we waste
  // memory and time, but won't get segfaults.
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
  template <class TPixel, class TOutsideVector>
    void CopyVectorOfVectorsToMatlab(MatlabExportFilter::MatlabOutputPointer output,
				     TOutsideVector v, 
				     mwSize outsideSize, mwSize insideSize);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "MatlabExportFilter.hxx"
#endif
#endif /* MATLABEXPORTFILTER_H */
