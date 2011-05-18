/*
 * FilterExclusions.cpp
 *
 * This is a file to remove some cluttering from file
 * ItkImFilter.cpp. It contains the template explicit specialization
 * code that prevents the compiler from trying to instantiate certain
 * input/output image types for certain filters.
 *
 * The way ItkImFilter.cpp works is that it figures out the input and
 * output types of the image at run-time, and selects the right
 * templated filter.
 *
 * However, this implies that the compiler has to compile every
 * possible combination of input/output voxel type for every
 * filter. 
 *
 * Some filters do not accept certain types. For example, the
 * BinaryThinningImageFilter3D must have the same input and output
 * types.
 *
 * We solve this problem with template explicit specialization. That
 * is, by default all data type combinations are compiled. If we want
 * to prevent the compiler from compiling a particular case, we
 * have to add a definition of BaseFilter for that particular
 * triplet (input data type, output data type, filter type), and give
 * it an empty constructor. This way, for that particular triplet,
 * the compiler won't try to instantiate the filter.
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
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

#ifndef FILTEREXCLUSIONS_CPP
#define FILTEREXCLUSIONS_CPP

// avoid compiling combinations of input/output data types that are
// not available for some filters.  Read the header help of
// FilterExclusions.cpp for more info
/*
 * BinaryThinningImageFilter3D
 */

typedef itk::BinaryThinningImageFilter3D< itk::Image<bool, Dimension>,
					  itk::Image<uint8_T, Dimension>
					  > ThinningFilter_bool_uint8_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<bool, Dimension>,
					  itk::Image<uint16_T, Dimension>
					  > ThinningFilter_bool_uint16_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<bool, Dimension>,
					  itk::Image<float, Dimension>
					  > ThinningFilter_bool_float;
typedef itk::BinaryThinningImageFilter3D< itk::Image<bool, Dimension>,
					  itk::Image<double, Dimension>
					  > ThinningFilter_bool_double;

typedef itk::BinaryThinningImageFilter3D< itk::Image<uint8_T, Dimension>,
					  itk::Image<bool, Dimension>
					  > ThinningFilter_uint8_T_bool;
typedef itk::BinaryThinningImageFilter3D< itk::Image<uint8_T, Dimension>,
					  itk::Image<uint16_T, Dimension>
					  > ThinningFilter_uint8_T_uint16_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<uint8_T, Dimension>,
					  itk::Image<float, Dimension>
					  > ThinningFilter_uint8_T_float;
typedef itk::BinaryThinningImageFilter3D< itk::Image<uint8_T, Dimension>,
					  itk::Image<double, Dimension>
					  > ThinningFilter_uint8_T_double;

typedef itk::BinaryThinningImageFilter3D< itk::Image<int8_T, Dimension>,
					  itk::Image<bool, Dimension>
					  > ThinningFilter_int8_T_bool;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int8_T, Dimension>,
					  itk::Image<uint8_T, Dimension>
					  > ThinningFilter_int8_T_uint8_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int8_T, Dimension>,
					  itk::Image<uint16_T, Dimension>
					  > ThinningFilter_int8_T_uint16_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int8_T, Dimension>,
					  itk::Image<float, Dimension>
					  > ThinningFilter_int8_T_float;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int8_T, Dimension>,
					  itk::Image<double, Dimension>
					  > ThinningFilter_int8_T_double;

typedef itk::BinaryThinningImageFilter3D< itk::Image<uint16_T, Dimension>,
					  itk::Image<bool, Dimension>
					  > ThinningFilter_uint16_T_bool;
typedef itk::BinaryThinningImageFilter3D< itk::Image<uint16_T, Dimension>,
					  itk::Image<uint8_T, Dimension>
					  > ThinningFilter_uint16_T_uint8_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<uint16_T, Dimension>,
					  itk::Image<float, Dimension>
					  > ThinningFilter_uint16_T_float;
typedef itk::BinaryThinningImageFilter3D< itk::Image<uint16_T, Dimension>,
					  itk::Image<double, Dimension>
					  > ThinningFilter_uint16_T_double;

typedef itk::BinaryThinningImageFilter3D< itk::Image<int16_T, Dimension>,
					  itk::Image<bool, Dimension>
					  > ThinningFilter_int16_T_bool;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int16_T, Dimension>,
					  itk::Image<uint8_T, Dimension>
					  > ThinningFilter_int16_T_uint8_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int16_T, Dimension>,
					  itk::Image<uint16_T, Dimension>
					  > ThinningFilter_int16_T_uint16_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int16_T, Dimension>,
					  itk::Image<float, Dimension>
					  > ThinningFilter_int16_T_float;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int16_T, Dimension>,
					  itk::Image<double, Dimension>
					  > ThinningFilter_int16_T_double;

typedef itk::BinaryThinningImageFilter3D< itk::Image<int32_T, Dimension>,
					  itk::Image<bool, Dimension>
					  > ThinningFilter_int32_T_bool;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int32_T, Dimension>,
					  itk::Image<uint8_T, Dimension>
					  > ThinningFilter_int32_T_uint8_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int32_T, Dimension>,
					  itk::Image<uint16_T, Dimension>
					  > ThinningFilter_int32_T_uint16_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int32_T, Dimension>,
					  itk::Image<float, Dimension>
					  > ThinningFilter_int32_T_float;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int32_T, Dimension>,
					  itk::Image<double, Dimension>
					  > ThinningFilter_int32_T_double;

typedef itk::BinaryThinningImageFilter3D< itk::Image<int64_T, Dimension>,
					  itk::Image<bool, Dimension>
					  > ThinningFilter_int64_T_bool;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int64_T, Dimension>,
					  itk::Image<uint8_T, Dimension>
					  > ThinningFilter_int64_T_uint8_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int64_T, Dimension>,
					  itk::Image<uint16_T, Dimension>
					  > ThinningFilter_int64_T_uint16_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int64_T, Dimension>,
					  itk::Image<float, Dimension>
					  > ThinningFilter_int64_T_float;
typedef itk::BinaryThinningImageFilter3D< itk::Image<int64_T, Dimension>,
					  itk::Image<double, Dimension>
					  > ThinningFilter_int64_T_double;

typedef itk::BinaryThinningImageFilter3D< itk::Image<float, Dimension>,
					  itk::Image<bool, Dimension>
					  > ThinningFilter_float_bool;
typedef itk::BinaryThinningImageFilter3D< itk::Image<float, Dimension>,
					  itk::Image<uint8_T, Dimension>
					  > ThinningFilter_float_uint8_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<float, Dimension>,
					  itk::Image<uint16_T, Dimension>
					  > ThinningFilter_float_uint16_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<float, Dimension>,
					  itk::Image<double, Dimension>
					  > ThinningFilter_float_double;

typedef itk::BinaryThinningImageFilter3D< itk::Image<double, Dimension>,
					  itk::Image<bool, Dimension>
					  > ThinningFilter_double_bool;
typedef itk::BinaryThinningImageFilter3D< itk::Image<double, Dimension>,
					  itk::Image<uint8_T, Dimension>
					  > ThinningFilter_double_uint8_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<double, Dimension>,
					  itk::Image<uint16_T, Dimension>
					  > ThinningFilter_double_uint16_T;
typedef itk::BinaryThinningImageFilter3D< itk::Image<double, Dimension>,
					  itk::Image<float, Dimension>
					  > ThinningFilter_double_float;

template <>
class BaseFilter< bool, uint8_T, ThinningFilter_bool_uint8_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< bool, uint16_T, ThinningFilter_bool_uint16_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< bool, float, ThinningFilter_bool_float >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< bool, double, ThinningFilter_bool_double >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class BaseFilter< uint8_T, bool, ThinningFilter_uint8_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< uint8_T, uint16_T, ThinningFilter_uint8_T_uint16_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< uint8_T, float, ThinningFilter_uint8_T_float >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< uint8_T, double, ThinningFilter_uint8_T_double >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class BaseFilter< int8_T, bool, ThinningFilter_int8_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int8_T, uint8_T, ThinningFilter_int8_T_uint8_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int8_T, uint16_T, ThinningFilter_int8_T_uint16_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int8_T, float, ThinningFilter_int8_T_float >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int8_T, double, ThinningFilter_int8_T_double >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class BaseFilter< uint16_T, bool, ThinningFilter_uint16_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< uint16_T, uint8_T, ThinningFilter_uint16_T_uint8_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< uint16_T, float, ThinningFilter_uint16_T_float >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< uint16_T, double, ThinningFilter_uint16_T_double >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class BaseFilter< int16_T, bool, ThinningFilter_int16_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int16_T, uint8_T, ThinningFilter_int16_T_uint8_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int16_T, uint16_T, ThinningFilter_int16_T_uint16_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int16_T, float, ThinningFilter_int16_T_float >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int16_T, double, ThinningFilter_int16_T_double >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class BaseFilter< int32_T, bool, ThinningFilter_int32_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int32_T, uint8_T, ThinningFilter_int32_T_uint8_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int32_T, uint16_T, ThinningFilter_int32_T_uint16_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int32_T, float, ThinningFilter_int32_T_float >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int32_T, double, ThinningFilter_int32_T_double >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class BaseFilter< int64_T, bool, ThinningFilter_int64_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int64_T, uint8_T, ThinningFilter_int64_T_uint8_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int64_T, uint16_T, ThinningFilter_int64_T_uint16_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int64_T, float, ThinningFilter_int64_T_float >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< int64_T, double, ThinningFilter_int64_T_double >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class BaseFilter< float, bool, ThinningFilter_float_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< float, uint8_T, ThinningFilter_float_uint8_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< float, uint16_T, ThinningFilter_float_uint16_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< float, double, ThinningFilter_float_double >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class BaseFilter< double, bool, ThinningFilter_double_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< double, uint8_T, ThinningFilter_double_uint8_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< double, uint16_T, ThinningFilter_double_uint16_T >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class BaseFilter< double, float, ThinningFilter_double_float >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

/*
 * SignedMaurerDistanceMapImageFilter
 */

typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<bool, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_bool_bool;
typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<uint8_T, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_uint8_T_bool;
typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<int8_T, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_int8_T_bool;
typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<uint16_T, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_uint16_T_bool;
typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<int16_T, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_int16_T_bool;
typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<uint32_T, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_uint32_T_bool;
typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<int32_T, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_int32_T_bool;
typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<int64_T, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_int64_T_bool;
typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<float, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_float_bool;
typedef itk::SignedMaurerDistanceMapImageFilter< 
  itk::Image<double, Dimension>,
  itk::Image<bool, Dimension>
  > MaurerFilter_double_bool;

template <>
class BaseFilter< bool, bool, MaurerFilter_bool_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};
template <>
class BaseFilter< uint8_T, bool, MaurerFilter_uint8_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};
template <>
class BaseFilter< int8_T, bool, MaurerFilter_int8_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};
template <>
class BaseFilter< uint16_T, bool, MaurerFilter_uint16_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};
template <>
class BaseFilter< int16_T, bool, MaurerFilter_int16_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};
template <>
class BaseFilter< uint32_T, bool, MaurerFilter_uint32_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};
template <>
class BaseFilter< int32_T, bool, MaurerFilter_int32_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};
template <>
class BaseFilter< int64_T, bool, MaurerFilter_int64_T_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};
template <>
class BaseFilter< float, bool, MaurerFilter_float_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};
template <>
class BaseFilter< double, bool, MaurerFilter_double_bool >
{
public:
  BaseFilter(char *, NrrdImage, int, mxArray**) {
    mexErrMsgTxt("Invalid input or output image type for itk::SignedMaurerDistanceMapImageFilter");
  }
};

#endif // FILTEREXCLUSIONS_CPP
