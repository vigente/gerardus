// Copyright (c) 2010 Inria Lorraine (France). All rights reserved.
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Installation/config/support/print_MPFI_version.cpp $
// $Id: print_MPFI_version.cpp 67093 2012-01-13 11:22:39Z lrineau $
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#include <iostream>
#include <string>
#include <mpfi.h>

int main(){
    std::cout<<"version="<<std::string(mpfi_get_version())<<std::endl;
    return 0;
}
