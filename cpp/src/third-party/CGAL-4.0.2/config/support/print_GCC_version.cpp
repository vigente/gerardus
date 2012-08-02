// Copyright (c) 2005  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Installation/config/support/print_GCC_version.cpp $
// $Id: print_GCC_version.cpp 67093 2012-01-13 11:22:39Z lrineau $
// 
//
// Author(s)     : Fernando Cacciola

// Print out the gcc version

#include <iostream>


//
// Just in case this is called with a non-gcc compiler such as pgCC
//
 
#ifdef __clang_version__
#  define _CGAL_GCC_VERSION "Not GNU/CC (CLANG)"
#else
#  ifndef __GNUC__
#    define _CGAL_GCC_VERSION "Not GNU/CC"
#  endif
#endif
#ifndef _CGAL_GCC_VERSION
#  ifdef __VERSION__
#    define _CGAL_GCC_VERSION __VERSION__
#  else
#    define _CGAL_GCC_VERSION "Unknown version (_CGAL_GCC_VERSION is not defined)"
#  endif
#endif

int main()
{
  std::cout << "version=" << _CGAL_GCC_VERSION << std::endl;
  return 0;
}
