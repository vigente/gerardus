/*
 * MatlabExportFilter.h
 *
 * Class to provide an interface to pass C++ data to Matlab
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

#ifndef MATLABEXPORTFILTER_H
#define MATLABEXPORTFILTER_H

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <vector>

/* ITK headers */
#include "itkSmartPointer.h"

/* Gerardus common functions */
#include "GerardusCommon.h"

class MatlabExportFilter: public itk::Object {

private:

  unsigned int    numArgs;
  mxArray         **args;

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
  void SetMatlabArgumentsPointer(int _nlhs, mxArray *_plhs[]) {
    this->numArgs = _nlhs;
    this->args = _plhs;
  }

  // get number of elements in the list of arguments
  unsigned int GetNumberOfArguments() {
    return this->numArgs;
  }

  // function to check that number of input arguments is within
  // certain limits
  void CheckNumberOfArguments(unsigned int min, unsigned int max);

  // function to allocate memory in Matlab and hijack it to be used as
  // an ITK filter output. This only works with those filters that
  // create their own output. For other filters, you cannot graft the
  // output; instead, use CopyItkImageOntoMatlab() after running the
  // filter
  //
  // size is a vector with the dimensions of the output image in
  // Matlab. For example, for a 256x200x512 image, size = {256, 200, 512}
  //
  // size is the same for vector or scalar images. 
  template <class TPixel, unsigned int VectorDimension, class TVector>
    void GraftItkImageOntoMatlab(typename itk::DataObject::Pointer image, 
				 std::vector<unsigned int> size,
				 unsigned int idx, std::string paramName);
  template <class TPixel, unsigned int VectorDimension>
    void GraftItkImageOntoMatlab(typename itk::DataObject::Pointer image, 
				 std::vector<unsigned int> size,
				 unsigned int idx, std::string paramName) {
    GraftItkImageOntoMatlab<TPixel, VectorDimension, TPixel>(image, size, 
							     idx, paramName);
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
    void CopyItkImageToMatlab(typename itk::DataObject::Pointer image, 
			      std::vector<unsigned int> size,
			      unsigned int idx, std::string paramName);
  template <class TPixel, unsigned int VectorDimension>
    void CopyItkImageToMatlab(typename itk::DataObject::Pointer image, 
			      std::vector<unsigned int> size,
			      unsigned int idx, std::string paramName) {
    CopyItkImageToMatlab<TPixel, VectorDimension, TPixel>(image, size, 
							  idx, paramName);
  }

  // function to allocate memory on the Matlab side and copy a vector
  // of scalars from the C++ side. It can export any C++ class that
  // accepts operator[], e.g. v[i].
  //
  // idx:  index of the Matlab output
  // v:    vector we want to export to Matlab
  // size: length of the v vector. Having size as a parameter here
  //       saves us from having to do a specialization of this function for
  //       every class, as different classes provide their size with
  //       different methods, e.g. v.size(), v.Size(), v.length(), etc.
  template <class TPixel, class TVector>
    void CopyVectorOfScalarsToMatlab(unsigned int idx, std::string paramName, 
				     TVector v, mwSize size);

  // function to allocate memory on the Matlab side and copy a vector
  // of vectors from the C++ side. Both vectors must accept
  // operator[], e.g. w = v[i], w[j]. The inside vectors will be
  // stacked by rows. E.g., if you have a vector of Points, each point
  // will be copied to a row.
  //
  // idx:   index of the Matlab output
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
  //        Having these lengths as parameters here
  //        saves us from having to do a specialization of this function for
  //        every class, as different classes provide their size with
  //        different methods, e.g. v.size(), v.Size(), v.length(), etc.
  template <class TPixel, class TInsideVector, class TOutsideVector>
    void CopyVectorOfVectorsToMatlab(unsigned int idx, std::string paramName,
				     TOutsideVector v, 
				     mwSize outsideSize, mwSize insideSize);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "MatlabExportFilter.hxx"
#endif
#endif /* MATLABEXPORTFILTER_H */
