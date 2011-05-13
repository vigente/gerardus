/*
 * FilterFactoryExclusions.hpp
 *
 * This is a file to remove some cluttering from file
 * ItkImFilter.cpp. It contains the template partial specialization
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
 * We solve this problem with template partial specialization. That
 * is, by default all data type combinations are compiled. If we want
 * to prevent the compiler from compiling a particular case, we
 * have to add a definition of FilterFactory for that particular
 * triplet (input data type, output data type, filter type), and give
 * it an empty constructor. This way, for that particular triplet,
 * the compiler won't try to instantiate the filter.
 *
 */

#ifndef FILTERFACTORYEXCLUSIONS_HPP
#define FILTERFACTORYEXCLUSIONS_HPP

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
class FilterFactory< bool, uint8_T, ThinningFilter_bool_uint8_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< bool, uint16_T, ThinningFilter_bool_uint16_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< bool, float, ThinningFilter_bool_float >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< bool, double, ThinningFilter_bool_double >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class FilterFactory< uint8_T, bool, ThinningFilter_uint8_T_bool >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< uint8_T, uint16_T, ThinningFilter_uint8_T_uint16_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< uint8_T, float, ThinningFilter_uint8_T_float >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< uint8_T, double, ThinningFilter_uint8_T_double >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class FilterFactory< int8_T, bool, ThinningFilter_int8_T_bool >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int8_T, uint8_T, ThinningFilter_int8_T_uint8_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int8_T, uint16_T, ThinningFilter_int8_T_uint16_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int8_T, float, ThinningFilter_int8_T_float >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int8_T, double, ThinningFilter_int8_T_double >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class FilterFactory< uint16_T, bool, ThinningFilter_uint16_T_bool >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< uint16_T, uint8_T, ThinningFilter_uint16_T_uint8_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< uint16_T, float, ThinningFilter_uint16_T_float >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< uint16_T, double, ThinningFilter_uint16_T_double >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class FilterFactory< int16_T, bool, ThinningFilter_int16_T_bool >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int16_T, uint8_T, ThinningFilter_int16_T_uint8_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int16_T, uint16_T, ThinningFilter_int16_T_uint16_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int16_T, float, ThinningFilter_int16_T_float >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int16_T, double, ThinningFilter_int16_T_double >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class FilterFactory< int32_T, bool, ThinningFilter_int32_T_bool >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int32_T, uint8_T, ThinningFilter_int32_T_uint8_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int32_T, uint16_T, ThinningFilter_int32_T_uint16_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int32_T, float, ThinningFilter_int32_T_float >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int32_T, double, ThinningFilter_int32_T_double >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class FilterFactory< int64_T, bool, ThinningFilter_int64_T_bool >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int64_T, uint8_T, ThinningFilter_int64_T_uint8_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int64_T, uint16_T, ThinningFilter_int64_T_uint16_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int64_T, float, ThinningFilter_int64_T_float >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< int64_T, double, ThinningFilter_int64_T_double >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class FilterFactory< float, bool, ThinningFilter_float_bool >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< float, uint8_T, ThinningFilter_float_uint8_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< float, uint16_T, ThinningFilter_float_uint16_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< float, double, ThinningFilter_float_double >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

template <>
class FilterFactory< double, bool, ThinningFilter_double_bool >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< double, uint8_T, ThinningFilter_double_uint8_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< double, uint16_T, ThinningFilter_double_uint16_T >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};
template <>
class FilterFactory< double, float, ThinningFilter_double_float >
{
public:
  FilterFactory(char *, NrrdImage, mxArray*) {
    mexErrMsgTxt("Invalid input or output image type for itk::BinaryThinningImageFilter3D");
  }
};

#endif // FILTERFACTORYEXCLUSIONS_HPP
