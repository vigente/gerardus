// Copyright (c) 2005  Stanford University (USA).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Kinetic_data_structures/include/CGAL/Polynomial/Root_stack_default_traits.h $
// $Id: Root_stack_default_traits.h 67093 2012-01-13 11:22:39Z lrineau $
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_ROOT_ENUMERATOR_DEFAULT_TRAITS_H
#define CGAL_POLYNOMIAL_ROOT_ENUMERATOR_DEFAULT_TRAITS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Root_stack_traits_base.h>

namespace CGAL { namespace POLYNOMIAL {

template <class Polynomial>
class Root_stack_default_traits: public internal::Root_stack_traits_base<Polynomial>
{
    private:
        typedef internal::Root_stack_traits_base<Polynomial> Base;

    public:
};

} } //namespace CGAL::POLYNOMIAL
#endif
