// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Linear_cell_complex/demo/Linear_cell_complex/Linear_cell_complex_3_demo.cpp $
// $Id: Linear_cell_complex_3_demo.cpp 67427 2012-01-24 18:39:05Z lrineau $
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#include "MainWindow.h"
#include "typedefs.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>

int main(int argc, char** argv)
{
  std::cout<<"Size of dart: "<<sizeof(LCC::Dart)<<std::endl;

  QApplication application(argc,argv);
  
  application.setOrganizationDomain("cgal.org");
  application.setOrganizationName("CNRS and LIRIS' Establishments");
  application.setApplicationName("3D Linear Cell Complex");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  CGAL_Qt4_init_resources(); // that function is in a DLL
  Q_INIT_RESOURCE(Linear_cell_complex_3);
  MainWindow mw;
  mw.show();

  return application.exec();
}
