/* 
 * BaseFilter.hpp
 *
 * BaseFilter<InVoxelType, OutVoxelType, FilterType>: This is where
 * the code to actually run the filter on the image lives.
 *
 * Instead of having a function (e.g. runFilter), we have the code in
 * the constructor of class FilterFactory.
 *
 * The reason is that template explicit specialization is only
 * possible in classes, not in functions. We need explicit
 * specialization to prevent the compiler from compiling certain
 * input/output image data types for some filters that don't accept
 * them.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
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

#ifndef BASEFILTER_HPP
#define BASEFILTER_HPP

/*
 * Global variables
 */
static const unsigned int Dimension = 3; // volume data dimension
                                         // (3D volume)

/*
 * Block of functions to allow testing of template types
 */
template< class T >
struct TypeIsBool
{ static const bool value = false; };

template<>
struct TypeIsBool< bool >
{ static const bool value = true; };

template< class T >
struct TypeIsUint8
{ static const bool value = false; };

template<>
struct TypeIsUint8< uint8_T >
{ static const bool value = true; };

template< class T >
struct TypeIsUint16
{ static const bool value = false; };

template<>
struct TypeIsUint16< uint16_T >
{ static const bool value = true; };

template< class T >
struct TypeIsFloat
{ static const bool value = false; };

template<>
struct TypeIsFloat< float >
{ static const bool value = true; };

template< class T >
struct TypeIsDouble
{ static const bool value = false; };

template<>
struct TypeIsDouble< double >
{ static const bool value = true; };

/*
 * Parser functions to map type variables to type templates
 */

void parseFilterTypeToTemplate(char *filter,
			       NrrdImage nrrd,
			       int nargout,
			       mxArray** &argOut);
template <class InVoxelType>
void parseOutputTypeToTemplate(char *filter,
			       NrrdImage nrrd,
			       int nargout,
			       mxArray** &argOut);
void parseInputTypeToTemplate(mxClassID inputVoxelClassId, 
			      char *filter,
			      NrrdImage nrrd,
			      int nargout,
			      mxArray** &argOut);

/*
 * FilterParamFactory: class to pass parameters specific to one filter
 * but not the others
 */

// default: any filter without an explicit specialization
template <class InVoxelType, class OutVoxelType, 
	  class FilterType>
class FilterParamFactory {
public:
  FilterParamFactory(typename FilterType::Pointer filter) {
    // by default, we assume that filters do not need parameters. If a
    // filter needs some specific parameters, or setting any flags, we
    // need to declare a explicit specialization of this class, and
    // put the corresponding code there
    //
    // Hence, this constructor is empty
    ;
  }
};

// SignedMaurerDistanceMapImageFilter: Specific parameters
template <class InVoxelType, class OutVoxelType>
class FilterParamFactory<InVoxelType, OutVoxelType,
			 itk::SignedMaurerDistanceMapImageFilter< 
			   itk::Image<InVoxelType, Dimension>,
			   itk::Image<OutVoxelType, Dimension> > 
			 > {
public:
  FilterParamFactory(typename itk::SignedMaurerDistanceMapImageFilter< 
		     itk::Image<InVoxelType, Dimension>,
		     itk::Image<OutVoxelType, Dimension> >::Pointer filter) {
    // we want the output in real world coordinates by default. If the
    // user wants voxel units, then provide a plain image at the
    // input, or make the spacing in the NRRD struct = [1.0, 1.0, 1.0]
    filter->SetUseImageSpacing(true);
  
    // we want actual Euclidean distances, not squared ones
    filter->SquaredDistanceOff();
  }
};

/* 
 * BaseFilter
 */
template <class InVoxelType, class OutVoxelType, class FilterType>
class BaseFilter {
private:
  int nargout;
public:
  BaseFilter(char *filterType, NrrdImage &nrrd, 
		int _nargout, mxArray** &argOut);
};

#endif /* BASEFILTER_HPP */
