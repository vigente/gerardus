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

  // input arguments registered with this import interface
  std::vector<const mxArray*> args;

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
    if (_prhs == NULL) {
      mexErrMsgTxt("Pointer to arguments array is NULL");
    }
    for (mwIndex i = 0; i < (mwIndex)_nrhs; i++) {
      this->args.push_back(_prhs[i]);
    }
  }

  // function to register an additional array to import Matlab data
  // from, but that is not directly a top level input argument. For
  // example, if we pass the argument {[1 2], [0 -3]}, the
  // corresponding top level element in the input argument list is a
  // cell array, but we have no way to access the vectors [1 2] and 
  // [0 -3]. With this function, we can register these vectors in the
  // MatlabImport interface, and then access them as if they were
  // normal top level input arguments.
  //
  // returns: an index to identify the newly registered input
  //          array. This index can then be used to access the array, e.g. with
  //          GetStringArgument(idx, ...)
  size_t SetAdditionalMatlabArgumentPointer(const mxArray *mp) {
    if (mp == NULL) {
      mexErrMsgTxt("Pointer to argument is NULL");
    }
    this->args.push_back(mp);
    return this->args.size() - 1;
  }

  // get number of elements in the list of arguments
  unsigned int GetNumberOfArguments() {
    return this->args.size();
  }

  // function to get direct pointers to the Matlab input arguments
  //
  // idx: parameter index
  const mxArray *GetArg(unsigned int idx) {
    if ((idx >= 0) && (idx < this->args.size())) {
      return this->args[idx];
    } else {
      mexErrMsgTxt("Index out of range of Matlab input arguments");
    }
    return NULL;
  }

  // function to check that number of input arguments is within
  // certain limits
  void CheckNumberOfArguments(unsigned int min, unsigned int max);

  // function to get the size of a Matlab array. It simplifies having
  // to run mxGetNumberOfDimensions() and mxGetDimensions(), and then
  // casting the result into e.g. itk::Size to pass it to ITK
  template <class VectorValueType, class VectorType>
    VectorType GetArraySize(unsigned int idx, 
			    std::string paramName,
			    VectorType def);

  // function to get the half-size of a Matlab array. Some ITK filters
  // request the "half-size" (called radius) of a Matlab array,
  // instead of its size. By "half-size" we mean the length of the side to
  // the left or right of the central pixel. For example, an array
  // with size=[3, 7] has a half-size or radius=[1, 3]. I.e. 
  // size = 2 * radius + 1
  template <class VectorValueType, class VectorType>
    VectorType GetArrayHalfSize(unsigned int idx, 
				std::string paramName,
				VectorType def);

  // function to get the value of input arguments that are strings
  //
  // idx: parameter index
  // def: value returned by default if argument is empty or not provided
  std::string GetStringArgument(unsigned int idx,
				std::string paramName,
				std::string def);

  // function to get the value of an input argument that is a numeric
  // scalar
  //
  // idx: parameter index
  // def: value returned by default if argument is empty or not provided
  template <class ParamType>
  ParamType GetScalarArgument(unsigned int idx, 
			      std::string paramName,
			      ParamType def);

  // function to get one scalar value from an input argument that is a matrix
  //
  // idx: parameter index
  // row: matrix row index of the scalar
  // col: matrix column index of the scalar
  // def: value returned by default if argument is empty or not provided
  template <class ParamType>
  ParamType GetScalarArgument(unsigned int idx,
			      mwIndex row,
			      mwIndex col,
			      std::string paramName,
			      ParamType def);

  // function to get an input argument as a vector of scalars. The
  // argument itself can be a row vector, or a 2D matrix. In the latter
  // case, the user has to select one of the rows of the matrix
  //
  // idx: parameter index within the list of Matlab input arguments
  // row: row index in the input 2D matrix
  // def: value returned by default if argument is empty or not provided
  template <class VectorValueType, class VectorType>
    VectorType GetRowVectorArgument(unsigned int idx, 
				   mwIndex row,
				   std::string paramName,
				   VectorType def);
  template <class VectorValueType, class VectorType>
    VectorType GetRowVectorArgument(unsigned int idx, 
				   std::string paramName,
				   VectorType def);

  // function to read a Matlab 2D matrix row by row. It returns the
  // matrix as a vector of rows. Each row is read as a C++ "vector". By
  // "vector" we mean a C++ class that is vector-like, e.g. std::vector,
  // CGAL::Point_3 or ITK::Size.
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
    GetMatrixAsVectorOfRowVectorsArgument(unsigned int idx, 
					  std::string paramName,
					  std::vector<VectorType> def);
  
  // function to read a Matlab array into a vector. This is the
  // equivalent to A(:) in Matlab
  template <class VectorValueType, class VectorType>
    VectorType
    GetArrayArgumentAsVector(unsigned int idx, 
			     std::string paramName,
			     VectorType def);
  
  // function to get an input argument that is an image. This function
  // returns an itk::ImportImageFilter, which can be used wherever an
  // itk:Image is required, without having to duplicate the Matlab
  // buffer
  //
  // idx: parameter index
  template <class TPixel, unsigned int VImageDimension>
    typename itk::Image<TPixel, VImageDimension>::Pointer
    GetImageArgument(unsigned int idx, std::string paramName);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "MatlabImportFilter.hxx"
#endif

#endif /* MATLABIMPORTFILTER_H */
