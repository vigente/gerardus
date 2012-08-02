// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Apollonius_graph_2/include/CGAL/Apollonius_graph_data_structure_2.h $
// $Id: Apollonius_graph_data_structure_2.h 67117 2012-01-13 18:14:48Z lrineau $
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H
#define CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H

#include <CGAL/Triangulation_data_structure_2.h>


namespace CGAL {

// For backward compatibility

template <class Vb, class Fb>
class Apollonius_graph_data_structure_2
  : public Triangulation_data_structure_2<Vb, Fb>
{};


} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H
