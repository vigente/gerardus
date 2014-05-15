/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
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
#include "vector.h"
#include "vector_exact.h"
#include "ssvector.h"
#include "svector.h"

namespace soplex
{
Vector& Vector::operator=(const Vector& vec)
{
   if (this != &vec)
   {
      assert(dim() == vec.dim());
      memcpy(val, vec.val, dimen*sizeof(Real));

      assert(isConsistent());
   }
   return *this;
}

Vector& Vector::operator=(const SVector& vec)
{
   clear();
   assign(vec);

   assert(isConsistent());

   return *this;
}

Vector& Vector::operator=(const Vector_exact& vec)
{
   assert(dim() == vec.dim());

   for( int i = 0; i < dim(); i++ )
   {
      val[i] = get_d(vec[i]);
   }

   assert(isConsistent());

   return *this;
}

Vector& Vector::assign(const SVector& psv)
{
   for (int i = psv.size(); i-- > 0;)
      val[psv.index(i)] = psv.value(i);

   assert(isConsistent());
      
   return *this;
}

Vector& Vector::operator+=(const Vector& vec)
{
   assert(dim() == vec.dim());
   for (int i = 0; i < dim(); ++i)
      val[i] += vec[i];
   return *this;
}

Vector& Vector::operator+=(const SVector& vec)
{
   for( int i = vec.size()-1; i >= 0 ; i-- )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] += vec.value(i);
   }
   return *this;
}

Vector& Vector::operator+=(const SSVector& vec)
{
   for( int i = vec.size()-1; i >= 0 ; i-- )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] -= vec.value(i);
   }
   return *this;
}

Vector& Vector::operator-=(const Vector& vec)
{
   assert(dim() == vec.dim());
   for (int i = 0; i < dim(); ++i)
      val[i] -= vec[i];
   return *this;
}

Vector& Vector::operator-=(const SVector& vec)
{
   for( int i = vec.size()-1; i >= 0 ; i-- )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] -= vec.value(i);
   }
   return *this;
}

Vector& Vector::operator-=(const SSVector& vec)
{
   for( int i = vec.size()-1; i >= 0 ; i-- )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] -= vec.value(i);
   }
   return *this;
}

Vector& Vector::operator*=(Real x)
{
   for (int i = 0; i < dim(); ++i)
      val[i] *= x;

   return *this;
}

Real Vector::length() const
{
   return sqrt(length2());
}

Real Vector::length2() const
{
   return *this * *this;
}

Real Vector::maxAbs() const
{
   Real maxi = 0.0;

   for(int i = 0; i < dim(); ++i)
      if (fabs(val[i]) > maxi)
         maxi = fabs(val[i]);

   assert(maxi >= 0.0);

   return maxi;
}

Real Vector::minAbs() const
{
   Real mini = infinity;

   for(int i = 0; i < dim(); ++i)
      if (fabs(val[i]) < mini)
         mini = fabs(val[i]);

   assert(mini >= 0.0);

   return mini;
}

std::ostream& operator<<(std::ostream& s, const Vector& vec)
{
   int i;
   s << '(';
   for (i = 0; i < vec.dim() - 1; ++i)
      s << vec[i] << ", ";
   s << vec[i] << ')';
   return s;
}

Real Vector::operator*(const SVector& v) const
{
   assert(dim() >= v.dim());
   int i;
   Real x = 0;
   for (i = v.size(); i-- > 0;)
      x += val[v.index(i)] * v.value(i);
   return x;
}

bool Vector::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (dim() > 0 && val == 0)
      return MSGinconsistent("Vector");
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
