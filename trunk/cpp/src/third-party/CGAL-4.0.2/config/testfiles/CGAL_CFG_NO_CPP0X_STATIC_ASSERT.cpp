// Copyright (c) 2011  INRIA Saclay (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Installation/config/testfiles/CGAL_CFG_NO_CPP0X_STATIC_ASSERT.cpp $
// $Id: CGAL_CFG_NO_CPP0X_STATIC_ASSERT.cpp 67093 2012-01-13 11:22:39Z lrineau $
//
// Author(s)     : Marc Glisse

//| If a compiler does not support static_assert (from C++0x)
//| CGAL_CFG_NO_CPP0X_STATIC_ASSERT is set. 

int main(){
	static_assert(true,"Everything is fine");
}
