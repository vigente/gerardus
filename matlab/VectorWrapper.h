/* VectorWrapper: there are classes in ITK, e.g.  that are
 * conceptually vectors, but don't have the same methods as
 * e.g. std::vector. For example, itk::Size<Dimension> has the static
 * method itk::Size<Dimension>::GetDimensionSize, while std::vector
 * has the method v.size().
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

#ifndef VECTORWRAPPER_H
#define VECTORWRAPPER_H

/* ITK headers */
#include "itkSize.h"

/* CGAL headers */
#include <CGAL/Simple_cartesian.h>

// std::vector<type>
template<class VectorType, class VectorValueType>
  class VectorWrapper{
 private:
  VectorType *v;
 public:
 VectorWrapper(VectorType &_v): v(&_v) {}
  unsigned int Size() {return v->size();}
  // a std::vector can accomodate any vector length, so this method always returns true
  bool IsCompatibleWithSize(mwSize) {return true;}
  void Resize(mwSize len) {v->resize(len);}
};

// CGAL::Point_3<CGAL::Simple_cartesian<double> >
template<class VectorValueType>
class VectorWrapper<CGAL::Point_3<CGAL::Simple_cartesian<VectorValueType> >,
  VectorValueType> {
 private:
  CGAL::Point_3<CGAL::Simple_cartesian<VectorValueType> > *v;
 public:
  
};

// itk::Size<Dimension>
template<>
class VectorWrapper<itk::Size<1>::SizeType, itk::Size<1>::SizeValueType> {
 private:
  itk::Size<1> *v;
 public:
 VectorWrapper(itk::Size<1> &_v): v(&_v) {}
  unsigned int Size() {return 1;}
  // can only accommodate vectors of length 1
  bool IsCompatibleWithSize(mwSize len) {return (len==1);}
  // itk::Size cannot be resized. If the new length is the same, we
  // just ignore it. If we are trying to resize to a different length, give
  // error message
  void Resize(mwSize len) {
    if (len != 1) 
      mexErrMsgTxt("GerardusCommon: VectorWrapper: Fixed length vector cannot change size");}
};

template<>
class VectorWrapper<itk::Size<2>::SizeType, itk::Size<2>::SizeValueType> {
 private:
  itk::Size<2> *v;
 public:
 VectorWrapper(itk::Size<2> &_v): v(&_v) {}
  unsigned int Size() {return 2;}
  // can only accommodate vectors of length 2
  bool IsCompatibleWithSize(mwSize len) {return (len==2);}
  // itk::Size cannot be resized. If the new length is the same, we
  // just ignore it. If we are trying to resize to a different length, give
  // error message
  void Resize(mwSize len) {
    if (len != 2) 
      mexErrMsgTxt("GerardusCommon: VectorWrapper: Fixed length vector cannot change size");}
};

template<>
class VectorWrapper<itk::Size<3>::SizeType, itk::Size<3>::SizeValueType> {
 private:
  itk::Size<3> *v;
 public:
 VectorWrapper(itk::Size<3> &_v): v(&_v) {}
  unsigned int Size() {return 3;}
  // can only accommodate vectors of length 3
  bool IsCompatibleWithSize(mwSize len) {return (len==3);}
  // itk::Size cannot be resized. If the new length is the same, we
  // just ignore it. If we are trying to resize to a different length, give
  // error message
  void Resize(mwSize len) {
    if (len != 3) 
      mexErrMsgTxt("GerardusCommon: VectorWrapper: Fixed length vector cannot change size");}
};

template<>
class VectorWrapper<itk::Size<4>::SizeType, itk::Size<4>::SizeValueType> {
 private:
  itk::Size<4> *v;
 public:
 VectorWrapper(itk::Size<4> &_v): v(&_v) {}
  unsigned int Size() {return 4;}
  // can only accommodate vectors of length 4
  bool IsCompatibleWithSize(mwSize len) {return (len==4);}
  // itk::Size cannot be resized. If the new length is the same, we
  // just ignore it. If we are trying to resize to a different length, give
  // error message
  void Resize(mwSize len) {
    if (len != 4) 
      mexErrMsgTxt("GerardusCommon: VectorWrapper: Fixed length vector cannot change size");}
};

template<>
class VectorWrapper<itk::Size<5>::SizeType, itk::Size<5>::SizeValueType> {
 private:
  itk::Size<5> *v;
 public:
 VectorWrapper(itk::Size<5> &_v): v(&_v) {}
  unsigned int Size() {return 5;}
  // can only accommodate vectors of length 5
  bool IsCompatibleWithSize(mwSize len) {return (len==5);}
  // itk::Size cannot be resized. If the new length is the same, we
  // just ignore it. If we are trying to resize to a different length, give
  // error message
  void Resize(mwSize len) {
    if (len != 5) 
      mexErrMsgTxt("GerardusCommon: VectorWrapper: Fixed length vector cannot change size");}
};

#endif /* VECTORWRAPPER_H */
