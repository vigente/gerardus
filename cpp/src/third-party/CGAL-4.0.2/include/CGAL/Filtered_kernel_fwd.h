// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// Copyright (c) 2009  GeometryFactory Sarl (France)
// All rights reserved.
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Filtered_kernel/include/CGAL/Filtered_kernel_fwd.h $
// $Id: Filtered_kernel_fwd.h 67093 2012-01-13 11:22:39Z lrineau $
// 
//
// Author(s)     : Sylvain Pion, Laurent Rineau

#ifndef CGAL_FILTERED_KERNEL_FWD_H
#define CGAL_FILTERED_KERNEL_FWD_H

#include <CGAL/config.h>

namespace CGAL {

#ifdef CGAL_NO_STATIC_FILTERS
template < typename CK, bool UseStaticFilters = false >
#else
template < typename CK, bool UseStaticFilters = true >
#endif
struct Filtered_kernel;

} // end namespace CGAL

#endif // CGAL_FILTERED_KERNEL_FWD_H
