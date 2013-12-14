/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996      Roland Wunderling                              */
/*                  1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include "spxdefines.h"
#include "vector_exact.h"
#include "vector.h"
#include "svector.h"

namespace soplex
{
Vector_exact& Vector_exact::operator=(const Vector_exact& vec)
{
   if (this != &vec)
   {
      assert(dim() == vec.dim());

      for( int i = 0; i < dimen; i++ )
      {
         val[i] = vec[i];
      }

      assert(isConsistent());
   }

   return *this;
}

bool operator==(const Vector_exact& vec1, const Vector_exact& vec2)
{
   if( &vec1 == &vec2 )
      return true;
   else if( vec1.dim() != vec2.dim() )
      return false;
   else
   {
      for( int i = 0; i < vec1.dim(); i++ )
      {
         if( vec1[i] != vec2[i] )
            return false;
      }
   }

   return true;
}

Vector_exact& Vector_exact::operator=(const Vector& vec)
{
   assert(dim() == vec.dim());

   for( int i = 0; i < dimen; i++ )
   {
      val[i] = vec[i];
   }

   assert(isConsistent());

   return *this;
}

Vector_exact& Vector_exact::operator=(const SVector& psv)
{
   clear();

   for( int i = 0; i < psv.size(); i++ )
   {
      assert(psv.index(i) < dim());
      val[psv.index(i)] = psv.value(i);
   }

   assert(isConsistent());

   return *this;
}

Vector_exact& Vector_exact::operator+=(const Vector_exact& vec)
{
   assert(dim() == vec.dim());

   for( int i = 0; i < dimen; i++ )
   {
      val[i] += vec[i];
   }

   return *this;
}

Vector_exact& Vector_exact::operator+=(const Vector& vec)
{
   assert(dim() == vec.dim());

   for( int i = 0; i < dim(); i++ )
   {
      val[i] += vec[i];
   }

   return *this;
}

Vector_exact& Vector_exact::operator+=(const SVector& vec)
{
   for( int i = 0; i < vec.size(); i++ )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] += vec.value(i);
   }

   return *this;
}

/* stopped here */

Vector_exact& Vector_exact::operator-=(const Vector_exact& vec)
{
   assert(dim() == vec.dim());

   for( int i = 0; i < dim(); i++ )
   {
      val[i] -= vec[i];
   }

   return *this;
}

Vector_exact& Vector_exact::operator-=(const Vector& vec)
{
   assert(dim() == vec.dim());

   for( int i = 0; i < dim(); i++ )
   {
      val[i] -= vec[i];
   }

   return *this;
}

Vector_exact& Vector_exact::operator-=(const SVector& vec)
{
   for( int i = 0; i < vec.size() ; i++ )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] -= vec.value(i);
   }

   return *this;
}

Vector_exact& Vector_exact::operator*=(MpqReal x)
{
   for( int i = 0; i < dim(); i++ )
   {
      val[i] *= x;
   }

   return *this;
}

MpqReal Vector_exact::maxAbs() const
{
   MpqReal maxi = 0.0;

   for( int i = 0; i < dim(); i++ )
   {
      MpqReal x = abs(val[i]);

      if( x > maxi )
         maxi = x;
   }

   assert(maxi >= 0.0);

   return maxi;
}

MpqReal Vector_exact::minAbs() const
{
   assert(dim() > 0);

   MpqReal mini = abs(val[0]);

   for( int i = 1; i < dim(); i++ )
   {
      MpqReal x = abs(val[i]);

      if( x < mini )
         mini = x;
   }

   assert(mini >= 0.0);

   return mini;
}

Vector_exact& Vector_exact::multAdd(MpqReal x, const Vector_exact& vec)
{
   assert(vec.dim() == dimen);

   for( int i = 0; i < dimen; i++ )
      val[i] += x * vec.val[i];

   return *this;
}

Vector_exact& Vector_exact::multAdd(MpqReal x, const SVector& vec)
{
   for( int i = 0; i < vec.size(); i++ )
   {
      assert(vec.index(i) < dim());
      val[vec.index(i)] += x * vec.value(i);
   }

   return *this;
}

std::ostream& operator<<(std::ostream& s, const Vector_exact& vec)
{
   int i;
   s << '(';
   for (i = 0; i < vec.dim() - 1; ++i)
      s << vec[i] << ", ";
   s << vec[i] << ')';
   return s;
}

bool Vector_exact::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (dim() > 0 && val == 0)
      return MSGinconsistent("Vector_exact");
#endif

   return true;
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
