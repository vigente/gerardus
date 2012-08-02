// Copyright (c) 1999,2007  
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Number_types/include/CGAL/sse2.h $
// $Id: sse2.h 67093 2012-01-13 11:22:39Z lrineau $
//
//
// Author(s)     :  Andreas Fabri

#ifndef CGAL_SSE2_H
#define CGAL_SSE2_H

#include <emmintrin.h>

#if defined ( _MSC_VER )
#define CGAL_ALIGN_16  __declspec(align(16))
#elif defined( __GNU__ )
#define  CGAL_ALIGN_16 __attribute__((aligned(16))) 
#endif

#endif // CGAL_SSE2_H
