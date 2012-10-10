// Copyright (c) 2010 GeometryFactory (France). All rights reserved.
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Number_types/include/CGAL/gmp.h $
// $Id: gmp.h 67093 2012-01-13 11:22:39Z lrineau $
// 
// Author: Andreas Fabri

#ifndef CGAL_GMP_H
#define CGAL_GMP_H 1

#include <CGAL/config.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4244 4146) // conversion with loss of data
                                     // warning on - applied on unsigned number
#endif

#include <gmp.h>


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_GMP_H
