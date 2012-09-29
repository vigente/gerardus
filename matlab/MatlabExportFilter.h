/*
 * MatlabExportFilter.h
 *
 * Class to provide an interface to pass C++ data to Matlab
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012 University of Oxford
  * Version: 0.1.1
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
  // an ITK filter output
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

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "MatlabExportFilter.hxx"
#endif
#endif /* MATLABEXPORTFILTER_H */
