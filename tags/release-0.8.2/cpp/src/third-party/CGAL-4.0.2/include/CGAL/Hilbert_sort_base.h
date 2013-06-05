// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Spatial_sorting/include/CGAL/Hilbert_sort_base.h $
// $Id: Hilbert_sort_base.h 67093 2012-01-13 11:22:39Z lrineau $
//
// Author(s)     : Christophe Delage

#ifndef CGAL_HILBERT_SORT_BASE_H
#define CGAL_HILBERT_SORT_BASE_H

#include <CGAL/basic.h>
#include <algorithm>

namespace CGAL {

namespace internal {

    template <class RandomAccessIterator, class Cmp>
    RandomAccessIterator
    hilbert_split (RandomAccessIterator begin, RandomAccessIterator end,
                   Cmp cmp = Cmp ())
    {
        if (begin >= end) return begin;

        RandomAccessIterator middle = begin + (end - begin) / 2;
        std::nth_element (begin, middle, end, cmp);
        return middle;
    }
}

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_BASE_H
