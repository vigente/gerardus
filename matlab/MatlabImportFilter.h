/*
 * MatlabImportFilter.h
 *
 * Class to provide an interface to import data from Matlab mxArrays
 * into ITK
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

#ifndef MATLABIMPORTFILTER_H
#define MATLABIMPORTFILTER_H

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>

/* ITK headers */
#include "itkSmartPointer.h"
#include "itkImportImageFilter.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>

class MatlabImportFilter: public itk::Object {

private:

  unsigned int    numArgs;
  const mxArray **args;

protected:

  MatlabImportFilter();
  ~MatlabImportFilter();

public:

  // standard class typedefs
  typedef MatlabImportFilter                Self;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;

  // method for creation through the object factory
  itkNewMacro(Self);

  // run-time type information (and related methods)
  itkTypeMacro(MatlabImportFilter, Object);

  // function to capture the array with the arguments provided by
  // Matlab
  void SetMatlabArgumentsPointer(int _nrhs, const mxArray *_prhs[]) {
    this->numArgs = _nrhs;
    this->args = _prhs;
    if (this->args == NULL) {
      mexErrMsgTxt("Pointer to arguments array is NULL");
    }
  }

  // get number of elements in the list of arguments
  unsigned int GetNumberOfArguments() {
    return this->numArgs;
  }

  // function to get direct pointers to the Matlab input arguments
  //
  // idx: parameter index
  const mxArray *GetArg(unsigned int idx) {
    if ((idx >= 0) && (idx < this->numArgs)) {
      return this->args[idx];
    } else {
      mexErrMsgTxt("Index out of range of Matlab input arguments");
    }
    return NULL;
  }

  // function to check that number of input arguments is within
  // certain limits
  void CheckNumberOfArguments(unsigned int min, unsigned int max);

  // function to get the value of input arguments that are strings
  //
  // idx: parameter index
  // def: value returned by default if argument is empty or not provided
  std::string GetStringArgument(unsigned int idx,
				std::string paramName,
				std::string def);

  // function to get the value of input arguments that are numeric
  // scalars from the array of input arguments
  //
  // idx: parameter index
  // def: value returned by default if argument is empty or not provided
  template <class ParamType>
  ParamType GetScalarArgument(unsigned int idx, 
			      std::string paramName,
			      ParamType def);

  // function to get the an input argument that is a vector of scalars
  // from the array of input arguments
  //
  // idx: parameter index
  // def: value returned by default if argument is empty or not provided
  template <class ParamType, class ParamValueType>
  ParamType GetVectorArgument(unsigned int idx, 
			      std::string paramName,
			      ParamType def);

  // function to get an input argument that is an image. This function
  // returns an itk::ImportImageFilter, which can be used wherever an
  // itk:Image is required, without having to duplicate the Matlab
  // buffer
  //
  // idx: parameter index
  template <class TPixel, unsigned int VImageDimension>
  typename itk::Image<TPixel, VImageDimension>::Pointer
  GetImageArgument(unsigned int idx, std::string paramName);

  // function to read a static 3-vector (a vector with 3 elements that
  // can only be created using the constructor) from one row of a 2D
  // input matrix with 3 columns
  //
  // ParamType: a class of objects that are constructed like
  //            x(-2, 0.3, 5), e.g. CGAL::Point_3< CGAL::Simple_cartesian<double> >
  // idx: parameter index
  // row: row index (C++ convention row = 0, 1, 2, ..., n-1)
  template <class ParamType>
    ParamType GetStaticVector3Argument(unsigned int idx, 
				       mwIndex row, 
				       std::string paramName,
				       ParamType def);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "MatlabImportFilter.hxx"
#endif
#endif /* MATLABIMPORTFILTER_H */
