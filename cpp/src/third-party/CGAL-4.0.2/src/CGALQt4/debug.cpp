// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/GraphicsView/src/CGALQt4/debug.cpp $
// $Id: debug.cpp 67117 2012-01-13 18:14:48Z lrineau $
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#include <CGAL/Qt/debug.h>
#include <QDir>

#include <iostream>

namespace CGAL {
namespace Qt {


void traverse_resources(const QString& name, const QString& dirname, int indent)
{
  std::cerr << qPrintable(QString(indent, ' '))
            << qPrintable(name);
  QString fullname = 
    dirname.isEmpty() ?
    name :
    dirname + "/" + name;
  QDir dir(fullname);
  if(dir.exists()) {
    std::cerr << "/\n";
    Q_FOREACH(QString path, dir.entryList())
    {
      traverse_resources(path, fullname, indent + 2);
    }
  }
  else {
    std::cerr << "\n";
  }
}

} // namesapce Qt
} // namespace CGAL
