/*
 * GerardusCommon.h
 *
 * Miscellaneous functions of general use.
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2011 University of Oxford
  * Version: 0.9.2
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

#ifndef GERARDUSCOMMON_H
#define GERARDUSCOMMON_H

/* mex headers */
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

/* ITK headers */
#include "itkOffset.h"

/*
 * utIsInterruptPending(): "undocumented MATLAB API implemented in
 * libut.so, libut.dll, and included in the import library
 * libut.lib. To use utIsInterruptPending in a mex-file, one must
 * manually declare bool utIsInterruptPending() because this function
 * is not included in any header files shipped with MATLAB. Since
 * libut.lib, by default, is not linked by mex, one must explicitly
 * tell mex to use libut.lib." -- Wotao Yin, 
 * http://www.caam.rice.edu/~wy1/links/mex_ctrl_c_trick/
 *
 */
#ifdef __cplusplus 
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif

/*
 * ctrlcCheckPoint(): function to check whether the user has pressed
 * Ctrl+C, and if so, terminate execution returning an error message
 * with a hyperlink to the offending function's help, and a hyperlink
 * to the line in the source code file this function was called from
 *
 * It is implemented as a C++ macro to check for the CTRL+C flag, and
 * a call to function ctrlcErrMsgTxt() inside, to throw the error. The
 * reason is that if ctrlcCheckPoint() were a function instead of a
 * macro, this would introduce a function call at every iteration of
 * the loop, which is very expensive. But then we don't want to put
 * the whole error message part inside a macro, it's bug-prone and bad
 * programming practice. And once the CTRL+C has been detected,
 * whether the error message is generated a bit faster or not is not
 * important.
 *
 * In practice, to use this function put a call like this e.g. inside
 * loops that may take for a very long time:
 *
 *    // exit if user pressed Ctrl+C
 *    ctrlcCheckPoint(__FILE__, __LINE__);
 *
 * sourceFile: full path and name of the C++ file that calls this
 *             function. This should usually be the preprocessor
 *             directive __FILE__
 *
 * lineNumber: line number where this function is called from. This
 *             should usually be the preprocessor directive __LINE__
 *
 */
inline
void ctrlcErrMsgTxt(std::string sourceFile, int lineNumber) {
  
  // run from here the following code in the Matlab side:
  //
  // >> path = mfilename('fullpath')
  //
  // this provides the full path and function name of the function
  // that called ctrlcCheckPoint()
  int nlhs = 1; // number of output arguments we expect
  mxArray *plhs[1]; // to store the output argument
  int nrhs = 1; // number of input arguments we are going to pass
  mxArray *prhs[1]; // to store the input argument we are going to pass
  prhs[0] = mxCreateString("fullpath"); // input argument to pass
  if (mexCallMATLAB(nlhs, plhs, nrhs, prhs, "mfilename")) { // run mfilename('fullpath')
    mexErrMsgTxt("ctrlcCheckPoint(): mfilename('fullpath') returned error");
  }
  if (plhs == NULL) {
    mexErrMsgTxt("ctrlcCheckPoint(): mfilename('fullpath') returned NULL array of outputs");
  }
  if (plhs[0] == NULL) {
    mexErrMsgTxt("ctrlcCheckPoint(): mfilename('fullpath') returned NULL output instead of valid path");
  }
  
  // get full path to current function, including function's name
  // (without the file extension)
  char *pathAndName = mxArrayToString(plhs[0]);
  if (pathAndName == NULL) {
    mexErrMsgTxt("ctrlcCheckPoint(): mfilename('fullpath') output cannot be converted to string");
  }
  
  // for some reason, using mexErrMsgTxt() to give this output
  // doesn't work. Instead, we have to give the output to the
  // standar error, and then call mexErrMsgTxt() to terminate
  // execution of the program
  std::cerr << "Operation terminated by user during "
	    << "<a href=\"matlab:helpUtils.errorDocCallback('"
	    << mexFunctionName()
	    << "', '" << pathAndName << ".m', " << lineNumber << ")\">"
	    << mexFunctionName()
	    << "</a> (<a href=\"matlab:opentoline('"
	    << sourceFile
	    << "'," << lineNumber << ",0)\">line " << lineNumber
	    << "</a>)"
	    << std::endl;
  mexErrMsgTxt("");
}

#define ctrlcCheckPoint(sourceFile, lineNumber)				\
  if (utIsInterruptPending()) {						\
    ctrlcErrMsgTxt(sourceFile, lineNumber);				\
  }

/* 
 * CAST2MWSIZE(): macro to cast to mwSize type. This definition is
 *                necessary for ITK v3.20.0 to avoid an error when
 *                trying to compile
 *                itk::FixedArray::operator[](unsigned __int64) for
 *                Windows 64 bit, but * maybe we can remove it when
 *                ITK v4.0.0 is released
 */
#ifdef _WIN64
#define CAST2MWSIZE(x) static_cast<unsigned long>(x)
#else
#define CAST2MWSIZE(x) static_cast<mwSize>(x)
#endif

/*
 * sub2ind(): function that converts r, c, s indices to linear indices
 *            in a 3D array (same as Matlab's function sub2ind(),
 *            although in Matlab indices start at 1, and in C++, they
 *            start at 0)
 *
 * R, C, S: size of the array in rows, columns and slices, respectively
 * rcs: subindices to be converted
 * r, c, s: subindices to be converted
 *
 */
inline
mwIndex sub2ind(mwSize R, mwSize C, mwSize S, std::vector<mwIndex> rcs);
inline
mwIndex sub2ind(mwSize R, mwSize C, mwSize S, itk::Offset<3> rcs);
inline
mwIndex sub2ind(mwSize R, mwSize C, mwSize S, mwIndex r, mwIndex c, mwIndex s);

/*
 * ind2sub(): function that converts linear indices in a 3D array to
 *            r, c, s indices (same as Matlab's function ind2sub(),
 *            although in Matlab indices start at 1, and in C++, they
 *            start at 0)
 *
 * R, C, S: size of the array in rows, columns and slices, respectively
 * rcs: subindices to be converted
 *
 */
inline
std::vector<mwIndex> ind2sub(mwSize R, mwSize C, mwSize S, mwIndex idx);
inline
itk::Offset<3> ind2sub_itkOffset(mwSize R, mwSize C, mwSize S, mwIndex idx);

/*
 * Block of functions to allow testing of template types
 */
template< class T1, class T2 >
struct TypesAreEqual
{ static const bool value = false; };

template< class T >
struct TypesAreEqual< T, T >
{ static const bool value = true; };

template< class T >
struct TypeIsBool
{ static const bool value = false; };

template<>
struct TypeIsBool< mxLogical >
{ static const bool value = true; };

template< class T >
struct TypeIsUint8
{ static const bool value = false; };

template<>
struct TypeIsUint8< uint8_T >
{ static const bool value = true; };

template< class T >
struct TypeIsInt8
{ static const bool value = false; };

template<>
struct TypeIsInt8< int8_T >
{ static const bool value = true; };

template< class T >
struct TypeIsUint16
{ static const bool value = false; };

template<>
struct TypeIsUint16< uint16_T >
{ static const bool value = true; };

template< class T >
struct TypeIsInt16
{ static const bool value = false; };

template<>
struct TypeIsInt16< int16_T >
{ static const bool value = true; };

template< class T >
struct TypeIsInt32
{ static const bool value = false; };

template<>
struct TypeIsInt32< int32_T >
{ static const bool value = true; };

template< class T >
struct TypeIsInt64
{ static const bool value = false; };

template<>
struct TypeIsInt64< int64_T >
{ static const bool value = true; };

template< class T >
struct TypeIsSignedLong
{ static const bool value = false; };

template<>
struct TypeIsSignedLong< signed long >
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

/**
 * function to print the content "of any sequence that supports input
 * iterators, e.g. std::vector, std::deque, std::list, std::string and
 * even plain arrays". Adapted from "Functors with state - 3 (print
 * contents of vector using std::copy)". By Abhishek Padmanabh,
 * http://learningcppisfun.blogspot.co.uk/2007/02/functors-with-state-3-print-contents-of.html
 *
 * SeqType: type of the sequence
 *
 * ostr:           stream to print the sequence to, e.g. std::cout, std::cerr
 * v:              sequence to print
 * delimiter:      delimiter character to print between sequence elements, default " "
 * endchar:        character to print right after last sequence element, default "" (nothing)
 *
 * Example: 
 *
 *   std::vector<double> v(5, 1.0);
 *   printSeq<std::vector<double> >(std::cout, v, " ", "*");
 *
 *   Output:
 *   1 1 1 1 1*
 */
template<typename SeqType>
void printSeq(std::ostream& ostr, 
	      SeqType v, 
	      const std::string& delimiter = " ",
	      const std::string& endchar = ""){
  //print each element but the last followed by the delimiter
  std::copy(v.begin(), v.end()-1, 
	    std::ostream_iterator<typename SeqType::iterator::value_type>(ostr, 
									  delimiter.c_str()));
  //print last element followed by the end character
  std::copy(v.end()-1, v.end(),
	    std::ostream_iterator<typename SeqType::iterator::value_type>(ostr, 
									  endchar.c_str()));
}

/**
 * Debugging class to get a string with the type of a template, e.g.
 *   std::cout << print_T<OffsetType>() << std::endl;
 * returns
 *   OffsetType = std::string print_T() 
 *   [with T = itk::Offset<3u>, std::string = std::basic_string<char>]
 */
template<typename T>
std::string print_T() {
  return __PRETTY_FUNCTION__;
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "GerardusCommon.hxx"
#endif
#endif /* GERARDUSCOMMON_H */
