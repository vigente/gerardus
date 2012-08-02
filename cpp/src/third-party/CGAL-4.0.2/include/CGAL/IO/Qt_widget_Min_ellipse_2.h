// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Qt_widget/include/CGAL/IO/Qt_widget_Min_ellipse_2.h $
// $Id: Qt_widget_Min_ellipse_2.h 67093 2012-01-13 11:22:39Z lrineau $
// 
//
// Author(s)     : Radu Ursu

#ifndef CGAL_QT_WIDGET_MIN_ELLIPSE_2_H
#define CGAL_QT_WIDGET_MIN_ELLIPSE_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Min_ellipse_2.h>

namespace CGAL{

template< class Traits_ >
Qt_widget&
operator<<(Qt_widget &ws,
              const CGAL::Min_ellipse_2<Traits_>& min_ellipse)
{
    typedef CGAL::Min_ellipse_2<Traits_>::Point_iterator  Point_iterator;

    Point_iterator  first( min_ellipse.points_begin());
    Point_iterator  last ( min_ellipse.points_end());
    for ( ; first != last; ++first)
        ws << *first;
    return ws;
}

}//end namespace CGAL

#endif
