/*
 * MatlabImportFilter.h
 *
 * Class to provide an interface to import data from Matlab mxArrays
 * into ITK
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012-2013 University of Oxford
  * Version: 0.6.4
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

  // these are the variables provided by the MEX API for the input of
  // the function (not currently used, as until we add a MatlabInput
  // struct, we are using the vector args below)
  const mxArray **prhs;
  int     nrhs;

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

  // function to register an array of mxArray input arguments with the
  // import filter ("input" means arguments that were passed to the
  // Matlab function). After registration, the Matlab input arguments
  // are available to the filter. The filter can then be used to read
  // scalars, strings, vectors, matrices, fields, cells, etc. from the
  // registered arguments
  //
  // nrhs: number of input arguments to register
  // prhs: array of input arguments
  //
  // returns: registration index of the first element that is registered
  int ConnectToMatlabFunctionInput(int nrhs, const mxArray *prhs[]);

  // as above with ConnectToMatlabFunctionInput(), but for
  // a single input argument
  int RegisterInputArgumentFromMatlab(const mxArray *arg);

  // as above with RegisterInputArgumentFromMatlab(), but for a single
  // input argument that is a struct field, e.g. x.alpha.
  //
  // if the field doesn't exist, or cannot be registered for whatever
  // reason, return -1. The reason is that in general, we are not
  // going to be sure whether the user provided the field or not in a
  // struct. If we threw an error, the interface with the user becomes
  // too rigid. The way we do it, if the field is not provided, we
  // just get an index=-1, that can be easily ignored or replaced by a
  // default value at a later stage
  //
  // arg:   pointer to the input argument struct
  // index: the input argument can be 
  int RegisterInputFieldArgumentFromMatlab(const mxArray *arg, 
					   const char *fieldname);
  // alternatively, we can use the index for already registered input
  // arguments
  int RegisterInputFieldArgumentFromMatlab(int index, 
					   const char *fieldname);

  // get number of elements in the list of arguments
  unsigned int GetNumberOfRegisteredArguments() {
    return this->args.size();
  }

  // function to get direct pointers to the Matlab input arguments
  //
  // idx: parameter index
  const mxArray *GetRegisteredArgument(int idx);

  // function to check that number of input arguments is within
  // certain limits
  void CheckNumberOfArguments(unsigned int min, unsigned int max);

  // function to get the size of a Matlab array. It simplifies having
  // to run mxGetNumberOfDimensions() and mxGetDimensions(), and then
  // casting the result into e.g. itk::Size to pass it to ITK
  template <class VectorValueType, class VectorType>
    VectorType ReadMatlabArraySize(int idx, 
			    std::string paramName,
			    VectorType def);

  template <class VectorValueType, class VectorType, unsigned int VectorSize>
    VectorType ReadMatlabArraySize(int idx, 
			    std::string paramName,
			    VectorType def);

  // function to get the half-size of a Matlab array. Some ITK filters
  // request the "half-size" (called radius) of a Matlab array,
  // instead of its size. By "half-size" we mean the length of the side to
  // the left or right of the central pixel. For example, an array
  // with size=[3, 7] has a half-size or radius=[1, 3]. I.e. 
  // size = 2 * radius + 1
  template <class VectorValueType, class VectorType>
    VectorType ReadMatlabArrayHalfSize(int idx, 
				std::string paramName,
				VectorType def);

  template <class VectorValueType, class VectorType, unsigned int VectorSize>
    VectorType ReadMatlabArrayHalfSize(int idx, 
				std::string paramName,
				VectorType def);

  // function to get the value of input arguments that are strings
  //
  // idx: parameter index
  // def: value returned by default if argument is empty or not provided
  std::string ReadStringFromMatlab(int idx,
				std::string paramName,
				std::string def);

  // function to get the value of an input argument that is a numeric
  // scalar
  //
  // idx: parameter index
  // def: value returned by default if argument is empty or not provided
  template <class ParamType>
  ParamType ReadScalarFromMatlab(int idx, 
			      std::string paramName,
			      ParamType def);

  // function to get one scalar value from an input argument that is a matrix
  //
  // idx: parameter index
  // row: matrix row index of the scalar
  // col: matrix column index of the scalar
  // def: value returned by default if argument is empty or not provided
  template <class ParamType>
  ParamType ReadScalarFromMatlab(int idx,
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
  template <class VectorValueType, class VectorType, unsigned int VectorSize>
    VectorType ReadRowVectorFromMatlab(int idx, 
				   mwIndex row,
				   std::string paramName,
				   VectorType def);
  template <class VectorValueType, class VectorType, unsigned int VectorSize>
    VectorType ReadRowVectorFromMatlab(int idx, 
				   std::string paramName,
				   VectorType def);

  template <class VectorValueType, class VectorType>
    VectorType ReadRowVectorFromMatlab(int idx, 
				   mwIndex row,
				   std::string paramName,
				   VectorType def);
  template <class VectorValueType, class VectorType>
    VectorType ReadRowVectorFromMatlab(int idx, 
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
    ReadVectorOfVectorsFromMatlab(int idx, 
					  std::string paramName,
					  std::vector<VectorType> def);
  
  // function to read a Matlab array into a vector. This is the
  // equivalent to A(:) in Matlab
  template <class VectorValueType, class VectorType>
    VectorType
    ReadArrayAsVectorFromMatlab(int idx, 
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
    GetImagePointerFromMatlab(int idx, std::string paramName);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "MatlabImportFilter.hxx"
#endif

#endif /* MATLABIMPORTFILTER_H */
