// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Mesh_2/include/CGAL/Mesh_2/Output_stream.h $
// $Id: Output_stream.h 67117 2012-01-13 18:14:48Z lrineau $
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESHES_OUTPUT_STREAM_H

#ifdef CGAL_MESHES_NO_OUTPUT
#  include <CGAL/IO/Verbose_ostream.h>
#  define CGAL_MESHES_OUTPUT_STREAM CGAL::Verbose_ostream()
#else
#  ifndef CGAL_MESHES_OUTPUT_ON_CERR
#    define CGAL_MESHES_OUTPUT_STREAM std::cout
#  else
#    define CGAL_MESHES_OUTPUT_STREAM std::cerr
#  endif
#endif

#endif
