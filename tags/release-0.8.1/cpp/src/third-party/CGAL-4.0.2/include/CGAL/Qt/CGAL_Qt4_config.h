// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/GraphicsView/include/CGAL/Qt/CGAL_Qt4_config.h $
// $Id: CGAL_Qt4_config.h 67141 2012-01-16 10:54:58Z lrineau $
// 
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_QT4_CONFIG_H
#define CGAL_QT4_CONFIG_H

#include <QtCore/qglobal.h>

#ifdef CGAL_Qt4_DLL
#  ifdef CGAL_Qt4_EXPORTS
#    define CGAL_QT4_EXPORT Q_DECL_EXPORT
#  else
#    define CGAL_QT4_EXPORT Q_DECL_IMPORT
#  endif
#else
// empty definition
#  define CGAL_QT4_EXPORT
#endif

#endif // CGAL_QT4_CONFIG_H
