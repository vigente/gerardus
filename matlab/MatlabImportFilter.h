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
  * Version: 0.8.0
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
//#include <CGAL/Image_3.h>
#include <CGAL/ImageIO.h>

class MatlabImportFilter: public itk::Object {

 public:

  // struct to encapsulate each of the inputs to Matlab
  struct MatlabInput {
    const mxArray *pm;  // Matlab MEX input
    std::string name;   // name of the input for error/debug messages
    bool isProvided;    // flag: has the user provided this input?
  };
  
  typedef std::list<MatlabInput>::iterator MatlabInputPointer;

private:

  // these are the variables provided by the MEX API for the input of
  // the function (not currently used, as until we add a MatlabInput
  // struct, we are using the vector args below)
  const mxArray **prhs;
  int     nrhs;

  // list of the inputs registered at this importer
  std::list<MatlabInput> inputsList;
  
protected:

  MatlabImportFilter();
  ~MatlabImportFilter();

public:

  // standard class typedefs
  typedef MatlabImportFilter                Self;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;

  // method for creation through the object factory.
  itkNewMacro(Self);

  // run-time type information (and related methods)
  itkTypeMacro(MatlabImportFilter, Object);

  // function to import into this class the array with the arguments
  // provided by Matlab.
  void ConnectToMatlabFunctionInput(int _nrhs, const mxArray *_prhs[]);

  // get number of elements in the prhs list of input arguments
  unsigned int GetNumberOfArguments();

  // function to get direct pointers to the Matlab input arguments.
  //
  // idx: parameter index
  const mxArray *GetPrhsArgument(int idx);

  // function to check that number of prhs arguments is within
  // certain limits
  void CheckNumberOfArguments(int min, int max);

  // Function to register an input at the import filter.
  //
  // Registration basically means "this input in Matlab is going to
  // correspond to X". Once an input has been registered, it can be
  // passed to the Read methods to copy the data over, etc.
  //
  // pos: 
  //   position index within the base array prhs
  //
  // pm:
  //   direct pointer to an input. pm == base[pos] in the Matlab
  //   array, but this syntax allows to register child inputs, e.g. a
  //   field in a struct
  //
  // name:
  //   name to associate to the input for debugging purposes
  //
  // returns:
  //   a class of type MatlabInput, defined above
  //
  // syntax 1: the input is in the Matlab default array, so we only
  //           need to provide the position of the particular input we want to
  //           register.
  //
  // syntax 2: valid for single inputs. This syntax is useful to
  // register e.g. a struct field or a cell element within a cell array.
  MatlabInputPointer RegisterInput(int pos, std::string name);
  MatlabInputPointer RegisterInput(const mxArray *pm, std::string name);

  // Function to register a field from a struct input at the import filter.
  //
  // structInput:
  //   pointer to an already registered input, that must be of type struct.
  //
  // field:
  //   name of the field we want to register.
  //
  // returns:
  //   a class of type MatlabInput, defined above
  MatlabInputPointer RegisterStructFieldInput(MatlabInputPointer structInput,
					      std::string field);

  // function to get a pointer to a registered Matlab input by
  // providing its name.
  //
  // If the user runs this method on a name that has not been
  // registered, the function will throw an error and exit Matlab.
  MatlabInputPointer GetRegisteredInput(std::string name);

  // function to get the size of a Matlab array. It simplifies having
  // to run mxGetNumberOfDimensions() and mxGetDimensions(), and then
  // casting the result into e.g. itk::Size to pass it to ITK.
  //
  // input:
  //   pointer to a registered input
  //
  // def:
  //   default value to return if the user has not provided an input
  template <class VectorValueType, class VectorType>
    VectorType ReadMatlabArraySize(MatlabInputPointer input,
				   VectorType def);

  template <class VectorValueType, class VectorType, unsigned int VectorSize>
    VectorType ReadMatlabArraySize(MatlabInputPointer input,
			    VectorType def);

  // function to get the half-size of a Matlab array. Some ITK filters
  // request the "half-size" (called radius) of a Matlab array,
  // instead of its size. By "half-size" we mean the length of the side to
  // the left or right of the central pixel. For example, an array
  // with size=[3, 7] has a half-size or radius=[1, 3]. I.e. 
  // size = 2 * radius + 1.
  //
  // input:
  //   pointer to a registered input
  //
  // def:
  //   default value to return if the user has not provided an input
  template <class VectorValueType, class VectorType>
    VectorType ReadMatlabArrayHalfSize(MatlabInputPointer input,
				VectorType def);

  template <class VectorValueType, class VectorType, unsigned int VectorSize>
    VectorType ReadMatlabArrayHalfSize(MatlabInputPointer input,
				       VectorType def);

  // function to get the value of input arguments that are strings.
  //
  // input:
  //   pointer to a registered input
  //
  // def:
  //   default value to return if the user has not provided an input
  std::string ReadStringFromMatlab(MatlabInputPointer input,
				   std::string def);

  // function to get the value of an input argument that is a numeric
  // scalar.
  //
  // input:
  //   pointer to a registered input
  //
  // def:
  //   default value to return if the user has not provided an input
  template <class ParamType>
  ParamType ReadScalarFromMatlab(MatlabInputPointer input,
				 ParamType def);

  // function to get one scalar value from an input argument that is a matrix.
  //
  // input:
  //   pointer to a registered input
  //
  // row: 
  //   matrix row index of the scalar
  //
  // col: 
  //   matrix column index of the scalar
  //
  // def:
  //   default value to return if the user has not provided an input
  template <class ParamType>
  ParamType ReadScalarFromMatlab(MatlabInputPointer input,
				 mwIndex row, mwIndex col, ParamType def);

  // function to get an input argument as a vector of scalars. The
  // argument itself can be a row vector, or a 2D matrix. In the latter
  // case, the user has to select one of the rows of the matrix.
  //
  // input:
  //   pointer to a registered input
  //
  // row: 
  //   matrix row index of the scalar
  //
  // col: 
  //   matrix column index of the scalar
  //
  // def:
  //   default value to return if the user has not provided an input
  template <class VectorValueType, class VectorType>
    VectorType ReadRowVectorFromMatlab(MatlabInputPointer input,
				       mwIndex row, VectorType def);
  template <class VectorValueType, class VectorType>
    VectorType ReadRowVectorFromMatlab(MatlabInputPointer input,
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
  //
  // input:
  //   pointer to a registered input
  //
  // def:
  //   default value to return if the user has not provided an input
  template <class VectorValueType, class VectorType>
    std::vector<VectorType> 
    ReadVectorOfVectorsFromMatlab(MatlabInputPointer input,
				  std::vector<VectorType> def);

  // function to read a Matlab array into a vector. This is the
  // equivalent to A(:) in Matlab.
  //
  // input:
  //   pointer to a registered input
  //
  // def:
  //   default value to return if the user has not provided an input
  template <class VectorValueType, class VectorType>
    VectorType
    ReadArrayAsVectorFromMatlab(MatlabInputPointer input,
				VectorType def);

  // function to get an input argument that is an image. This function
  // returns an itk::ImportImageFilter, which can be used wherever an
  // itk:Image is required, without having to duplicate the Matlab
  // buffer.
  //
  // input:
  //   pointer to a registered input
  template <class TPixel, unsigned int VImageDimension>
    typename itk::Image<TPixel, VImageDimension>::Pointer
    GetImagePointerFromMatlab(MatlabInputPointer input);

  // function to read a Matlab 3D array into a CGAL::_image*.
  //
  // Because CGAL will try to free() the array memory, Matlab would
  // crash if we make the CGAL::_image* image point directly to the
  // Matlab buffer. Thus, unfortunately, this function needs to make a
  // duplicate of the Matlab array with the image, and feed that
  // duplicate to CGAL::_image*.
  //
  // Note: While in Matlab x->cols, y->rows, in CGAL y->cols,
  // x->rows, so we need to swap the row, col dimensions stored in
  // data objects designed for Matlab, e.g. MatlabImageHeader. That's
  // why the implementation of this function has e.g. x <-> imHeader.size[0],
  // y <-> imHeader.size[1]. Normally, when we deal with Matlab code,
  // we would swap that relationship.
  // 
  // input:
  //   pointer to a registered input
  _image*
    ReadCgalImageFromMatlab(MatlabInputPointer input);

  // function to swap the values of X and Y in a vector of vectors. That is,
  //
  // [1 4 7]     [4 1 7]
  // [2 5 8]  -> [5 2 8]
  // [3 6 9]     [6 3 9]
  //
  // Matlab stores matrices as (rows <-> y, cols <-> x), while ITK
  // uses (rows <-> x, cols <-> y). Often, this is not a problem
  // because when we read from Matlab what happens is that ITK
  // receives the image transposed in X and Y, it processes the image
  // and transposes it back to Matlab. However, if we import a list of
  // point coordinates, [x y z], to apply to the image, after
  // importing from Matlab into ITK we'll need to swap the coordinates
  // as [y x z] before we can apply them correctly to the image.
  //
  // Note that we only need to swap if the data has the [x y z]
  // format, not if it has the [r c s] format.
  //
  // @param vv: Vector of vectors, e.g. std::vector<itk::Point>
  // @param len: length of the outside vector
  template <class ValueType, class OutsideVectorType>
    void SwapXYInVectorOfVectors(OutsideVectorType &vv, mwSize len);

  // function to swap the values of X and Y in a vector.
  //
  // See the description above in SwapXYInVectorOfVectors().
  template <class ValueType, class VectorType>
    void SwapXYInVector(VectorType &v);

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "MatlabImportFilter.hxx"
#endif

#endif /* MATLABIMPORTFILTER_H */
