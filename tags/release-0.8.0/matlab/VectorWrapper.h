/*
 * VectorWrapper.h
 *
 * There are classes in ITK, e.g.  that are conceptually vectors, but
 * don't have the same methods as e.g. std::vector. For example,
 * itk::Size<Dimension> has the static method
 * itk::Size<Dimension>::GetDimensionSize, while std::vector has the
 * method v.size().
 *
 * This means that a function in Gerardus would need to be coded
 * twice, in both cases with the same code except for the
 * ::GetDimensionSize (know at compilation time) or .size() (known at
 * run time).
 *
 * Another example. In the CGAL library,
 * CGAL::Point_3<CGAL::Simple_cartesian<double> > is a 3-vector. This
 * vector is different from the other two in that its elements cannot
 * be populated using e.g. v[0] = 1.0; Instead, the vector is
 * populated with the constructor,
 * CGAL::Point_3<CGAL::Simple_cartesian<double> > v 
 *   = CGAL::Point_3<CGAL::Simple_cartesian<double> >(1.0, -2.3, 4.5);
 *
 * To overcome this problem, this wrapper provides a common interface
 * to different vector classes, so that we can interact with all these
 * vector-like classes the same way. By default, we assume that we
 * have a std::vector. For other types, we specialise the
 * VectorWrapper class
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

#ifndef VECTORWRAPPER_H
#define VECTORWRAPPER_H

/* ITK headers */
#include "itkSize.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>

/* Boost headers */
#include <boost/lexical_cast.hpp>

/*
 * By default, VectorWrapper assumes that we want to put Matlab's row data into an
 * std::vector<type>
 */
template<class VectorValueType, class VectorType, class MatlabValueType>
  class VectorWrapper{
 public:
  VectorWrapper() {}

  // read a row from a Matlab matrix
  VectorType ReadRowVector(const mxArray *pm, mwIndex row, std::string paramName);

  // read a whole array into a vector
  VectorType ReadArrayAsVector(const mxArray *pm, std::string paramName);

  // read the argument dimensions into a vector (size = 2 * halfsize + 1)
  VectorType ReadSize(const mxArray *pm, std::string paramName);
  VectorType ReadHalfSize(const mxArray *pm, std::string paramName);

};

/*
 * Partial specialisation if we want to put Matlab's row data into an
 * itk::Size<Dimension>::SizeType vector-like class
 */

// auxiliary functions so that we don't need to rewrite this code in
// every partial specialization

// ItkSizeCommonReadRowVector<Dimension>
template <class MatlabValueType, unsigned int Dimension>
typename itk::Size<Dimension>::SizeType
ItkSizeCommonReadRowVector(const mxArray *pm, mwIndex row, std::string paramName);

// ItkSizeCommonReadSize<Dimension>
template <unsigned int Dimension>
typename itk::Size<Dimension>::SizeType
ItkSizeCommonReadSize(const mxArray *pm, std::string paramName);

// ItkSizeCommonReadHalfSize<Dimension>
template <unsigned int Dimension>
typename itk::Size<Dimension>::SizeType
ItkSizeCommonReadHalfSize(const mxArray *pm, std::string paramName);

// partial specialisations
template<class MatlabValueType>
class VectorWrapper<itk::Size<2>::SizeValueType, itk::Size<2>::SizeType, MatlabValueType>{
 public:
  VectorWrapper() {}
  itk::Size<2>::SizeType ReadRowVector(const mxArray *pm, mwIndex row, std::string paramName) {
    return ItkSizeCommonReadRowVector<MatlabValueType, 2>(pm, row, paramName);
  }
  itk::Size<2>::SizeType ReadSize(const mxArray *pm, std::string paramName) {
    return ItkSizeCommonReadSize<2>(pm, paramName);
  }
  itk::Size<2>::SizeType ReadHalfSize(const mxArray *pm, std::string paramName) {
    return ItkSizeCommonReadHalfSize<2>(pm, paramName);
  }
};

template<class MatlabValueType>
class VectorWrapper<itk::Size<3>::SizeValueType, itk::Size<3>::SizeType, MatlabValueType>{
 public:
  VectorWrapper() {}
  itk::Size<3>::SizeType ReadRowVector(const mxArray *pm, mwIndex row, std::string paramName) {
    return ItkSizeCommonReadRowVector<MatlabValueType, 3>(pm, row, paramName);
  }
  itk::Size<3>::SizeType ReadSize(const mxArray *pm, std::string paramName) {
    return ItkSizeCommonReadSize<3>(pm, paramName);
  }
  itk::Size<3>::SizeType ReadHalfSize(const mxArray *pm, std::string paramName) {
    return ItkSizeCommonReadHalfSize<3>(pm, paramName);
  }
};

template<class MatlabValueType>
class VectorWrapper<itk::Size<4>::SizeValueType, itk::Size<4>::SizeType, MatlabValueType>{
 public:
  VectorWrapper() {}
  itk::Size<4>::SizeType ReadRowVector(const mxArray *pm, mwIndex row, std::string paramName) {
    return ItkSizeCommonReadRowVector<MatlabValueType, 4>(pm, row, paramName);
  }
  itk::Size<4>::SizeType ReadSize(const mxArray *pm, std::string paramName) {
    return ItkSizeCommonReadSize<4>(pm, paramName);
  }
  itk::Size<4>::SizeType ReadHalfSize(const mxArray *pm, std::string paramName) {
    return ItkSizeCommonReadHalfSize<4>(pm, paramName);
  }
};

/*
 * Partial specialisation if we want to put Matlab's row data into an
 * CGAL::Point_3<CGAL::Simple_cartesian<type> > vector-like class
 */

// CgalCommonReadStaticRowVector<VectorType>
//
// auxiliary function so that we don't need to rewrite this code in
// every partial specialization
template <class VectorValueType, class VectorType, class MatlabValueType>
VectorType
CgalCommonReadStaticRowVector(const mxArray *pm, mwIndex row, std::string paramName);

// partial specialisation for CGAL::Point_3<CGAL::Simple_cartesian<double> >
template<class MatlabValueType>
class VectorWrapper<double, CGAL::Point_3<CGAL::Simple_cartesian<double> >, 
  MatlabValueType>{
 public:
  VectorWrapper() {}
  CGAL::Point_3<CGAL::Simple_cartesian<double> >
    ReadRowVector(const mxArray *pm, mwIndex row, std::string paramName) {
    return CgalCommonReadStaticRowVector<double, 
      CGAL::Point_3<CGAL::Simple_cartesian<double> >, MatlabValueType>(pm, row, paramName);
  }
};

// partial specialisation for CGAL::Direction_3<CGAL::Simple_cartesian<double> >
template<class MatlabValueType>
class VectorWrapper<double, CGAL::Direction_3<CGAL::Simple_cartesian<double> >, 
  MatlabValueType>{
 public:
  VectorWrapper() {}
  CGAL::Direction_3<CGAL::Simple_cartesian<double> >
    ReadRowVector(const mxArray *pm, mwIndex row, std::string paramName) {
    return CgalCommonReadStaticRowVector<double, 
      CGAL::Direction_3<CGAL::Simple_cartesian<double> >, MatlabValueType>(pm, row, paramName);
  }
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "VectorWrapper.hxx"
#endif

#endif /* VECTORWRAPPER_H */
