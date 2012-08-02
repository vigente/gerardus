// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Boolean_set_operations_2/include/CGAL/Boolean_set_operations_2/Gps_intersection_functor.h $
// $Id: Gps_intersection_functor.h 67117 2012-01-13 18:14:48Z lrineau $
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_INTERSECTION_FUNCTOR_H
#define CGAL_GPS_INTERSECTION_FUNCTOR_H

#include <CGAL/Boolean_set_operations_2/Gps_base_functor.h>

namespace CGAL {

template <class Arrangement_>
class Gps_intersection_functor : public Gps_base_functor<Arrangement_>
{
public:

  typedef Arrangement_                                    Arrangement_2;

  typedef typename Arrangement_2::Face_const_handle       Face_const_handle;
  typedef typename Arrangement_2::Face_handle             Face_handle;


  void create_face (Face_const_handle f1,
                    Face_const_handle f2,
                    Face_handle res_f)
  {
    if(f1->contained() && f2->contained())
    {
      res_f->set_contained(true);
    }
  }
};


} //namespace CGAL

#endif
