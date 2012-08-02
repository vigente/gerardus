// Copyright (c) 2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.0-branch/Installation/config/testfiles/CGAL_CFG_USING_BASE_MEMBER_BUG_2.cpp $
// $Id: CGAL_CFG_USING_BASE_MEMBER_BUG_2.cpp 67093 2012-01-13 11:22:39Z lrineau $
// 
//
// Author(s)     : Sylvain Pion

//| This flag is set if the compiler bugs with some "using Base::Member;" in
//| a derived class.  The workaround is to write a forwarder or not use using.
//| At least SunPro CC 5.3 has this bug where the typical error message is :
//| "Error: The function B<int>::g() has not had a body defined."
//| Note that the subtlely is that the error message does not mention
//| "Member"...
//| This test is updated (hijacked) to detect an issue with sunpro 5.9
//| that instantiates a function twice.

template < class T >
struct A
{
  A() {}

  void f();
};

template < class T >
void A<T>::f(){}

template < class GT >
struct B : A<GT>
{
  using A<GT>::f;

  B() {}

  void g();
};

template < class GT >
void
B<GT>::g()
{
  f();
}

template struct B<int>;

int main()
{
  B<int> b;
  b.g();
  return 0;
}
