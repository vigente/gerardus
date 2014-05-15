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
#include <iomanip>
#include <assert.h>

#include "spxdefines.h"
#include "ssvector.h"
#include "svset.h"
#include "exceptions.h"

/**@file ssvector.cpp
 * @todo There is a lot pointer arithmetic done here. It is not clear if
 *       this is an advantage at all. See all the function int() casts.
 * @todo Several operations like maxAbs could setup the vector while
 *       computing the result.
 */
namespace soplex
{
#define MARKER   1e-100

static const Real shortProductFactor = 0.5;

void SSVector::setMax(int newmax)
{
   assert(idx    != 0);
   assert(newmax != 0);
   assert(newmax >= IdxSet::size());

   // len = (newmax < IdxSet::max()) ? IdxSet::max() : newmax;
   len = newmax;

   spx_realloc(idx, len);
}

void SSVector::reDim (int newdim)
{
   for (int i = IdxSet::size() - 1; i >= 0; --i)
      if (index(i) >= newdim)
         remove(i);
   DVector::reDim(newdim);
   setMax(DVector::memSize() + 1);
   assert(isConsistent());
}

void SSVector::reMem(int newsize)
{
   DVector::reSize(newsize);
   assert(isConsistent());
   setMax(DVector::memSize() + 1);
}

void SSVector::clear()
{
   if (isSetup())
   {
      for(int i = 0; i < num; ++i)
         val[idx[i]] = 0.0;
   }
   else
      Vector::clear();

   IdxSet::clear();
   setupStatus = true;
   assert(isConsistent());
}

void SSVector::setValue(int i, Real x)
{
   assert(i >= 0 && i < DVector::dim());

   if (isSetup())
   {
      int n = number(i);

      if (n < 0)
      {
         if (isNotZero(x, epsilon))
            IdxSet::add(1, &i);
      }
      else if (x == 0)
         clearNum(n);
   }
   val[i] = x;

   assert(isConsistent());
}

// #ifdef USE_OLD // old version
// void SSVector::setup()
// {
//    if (!isSetup())
//    {
//       IdxSet::clear();
// 
// // #define      TWO_LOOPS
// #ifdef  TWO_LOOPS
//       int i = 0;
//       int n = 0;
//       int* id = idx;
//       Real* v = val;
//       const Real* end = val + dim();
// 
//       while (v < end)
//       {
//          id[n] = i++;
//          n += (*v++ != 0);
//       }
// 
//       Real x;
//       int* ii = idx;
//       int* last = idx + n;
//       v = val;
// 
//       for (; id < last; ++id)
//       {
//          x = v[*id];
//          if (isNotZero(x, epsilon))
//             *ii++ = *id;
//          else
//             v[*id] = 0;
//       }
//       num = ii - idx;
// 
// #else
// 
//       if (dim() <= 1)
//       {
//          if (dim())
//          {
//             if (isNotZero(*val, epsilon))
//                IdxSet::add(0);
//             else
//                *val = 0;
//          }
//       }
//       else
//       {
//          int* ii = idx;
//          Real* v = val;
//          Real* end = v + dim() - 1;
// 
//          /* setze weissen Elefanten */
//          Real last = *end;
//          *end = MARKER;
// 
//          /* erstes element extra */
//          if (isNotZero(*v, epsilon))
//             *ii++ = 0;
//          else
//             *v = 0.0;
// 
//          for(;;)
//          {
//             while (*++v == 0.0)
//                ;
//             if (isNotZero(*v, epsilon))
//             {
//                *ii++ = int(v - val);
//             }
//             else
//             {
//                *v = 0.0;
//                if (v == end)
//                   break;
//             }
//          }
//          /* fange weissen Elefanten wieder ein */
//          if (isNotZero(last, epsilon))
//          {
//             *v = last;
//             *ii++ = dim() - 1;
//          }
//          else
//             *v = 0;
// 
//          num = int(ii - idx);
//       }
// 
// #endif
// 
//       setupStatus = true;
//       assert(isConsistent());
//    }
// }
// #else // new version, not yet fully tested
void SSVector::setup()
{
   if (!isSetup())
   {
      IdxSet::clear();

      num = 0;

      for(int i = 0; i < dim(); ++i)
      {
         if (val[i] != 0.0)
         {
            if (isZero(val[i], epsilon))
               val[i] = 0.0;
            else
            {
               idx[num] = i;
               num++;
            }
         }
      }
      setupStatus = true;
      assert(isConsistent());
   }
}
// #endif

SSVector& SSVector::operator+=(const Vector& vec)
{
   Vector::operator+=(vec);

   if (isSetup())
   {
      setupStatus = false;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator+=(const SVector& vec)
{
   Vector::operator+=(vec);

   if (isSetup())
   {
      setupStatus = false;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator+=(const SSVector& vec)
{
   for (int i = 0; i < vec.size(); ++i)
      val[vec.index(i)] += vec.value(i);

   if (isSetup())
   {
      setupStatus = false;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator-=(const Vector& vec)
{
   Vector::operator-=(vec);

   if (isSetup())
   {
      setupStatus = false;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator-=(const SVector& vec)
{
   Vector::operator-=(vec);

   if (isSetup())
   {
      setupStatus = false;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator-=(const SSVector& vec)
{
   if (vec.isSetup())
   {
      for (int i = 0; i < vec.size(); ++i)
         val[vec.index(i)] -= vec.value(i);
   }
   else
   {
      Vector::operator-=(Vector(vec));
   }

   if (isSetup())
   {
      setupStatus = false;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator*=(Real x)
{
   assert(isSetup());

   for (int i = 0; i < size(); ++i)
      val[index(i)] *= x;

   assert(isConsistent());

   return *this;
}

Real SSVector::maxAbs() const
{
   if (isSetup())
   {
      Real maxabs = REAL(0.0);

      for(int i = 0; i < num; ++i)
      {
         Real x = fabs(val[idx[i]]);

         if (x > maxabs)
            maxabs = x;
      }
      return maxabs;
   }
   else
      return Vector::maxAbs();
}

Real SSVector::length2() const
{
   Real x = REAL(0.0);

   if (isSetup())
   {
      for(int i = 0; i < num; ++i)
         x += val[idx[i]] * val[idx[i]];
   }
   else
      x = Vector::length2();

   return x;
}

Real SSVector::length() const
{
   return sqrt(length2());
}

#if 0 // buggy and not used
/* @todo check if really not used or if the Vector version is used instead.
 */
SSVector& SSVector::multAdd(Real xx, const SSVector& svec)
{
   if (svec.isSetup())
   {
      if (isSetup())
      {
         int i, j;
         Real x;

         for (i = svec.size() - 1; i >= 0; --i)
         {
            j = svec.index(i);
            if (val[j])
            {
               x = val[j] + xx * svec.value(i);
               if (isNotZero(x, epsilon))
                  val[j] = x;
               else
               {
                  val[j] = 0;
                  for (--i; i >= 0; --i)
                     val[svec.index(i)] += xx * svec.value(i);
                  unSetup();
                  break;
               }
            }
            else
            {
               x = xx * svec.value(i);
               if (isNotZero(x, epsilon))
               {
                  val[j] = x;
                  addIdx(j);
               }
            }
         }
      }
      else
         Vector::multAdd(xx, svec);
   }
   else
   {
      /**@todo this code does not work, because in is never something
       *       added to v. Also the idx will not be setup correctly
       *       Fortunately the whole function seems not to be called
       *       at all. 
       */
      throw SPxInternalCodeException("XSSVEC01 This should never happen.");

      Real y;
      int* ii = idx;
      Real* v = val;
      Real* rv = static_cast<Real*>(svec.val);
      Real* last = rv + svec.dim() - 1;
      Real x = *last;

      *last = MARKER;
      for(;;)
      {
         while (!*rv)
         {
            ++rv;
            ++v;
         }
         y = *rv++ * xx;
         if (isNotZero(y, epsilon))
         {
            *ii++ = int(v - val);
            *v++ = y;
         }
         else if (rv == last)
            break;
         else
            v++;
      }
      *rv = x;

      x *= xx;
      if (isNotZero(x, epsilon))
      {
         *ii++ = int(v - val);
         *v = x;
      }
      num = int(ii - idx);

      setupStatus = true;
   }

   assert(isConsistent());
   return *this;
}
#endif // 0

///@todo SSVector::multAdd() should be rewritten without pointer arithmetic.
SSVector& SSVector::multAdd(Real xx, const SVector& svec)
{
   if (isSetup())
   {
      int i, j;
      Real x;
      Real* v = val;
      int adjust = 0;

      for (i = svec.size() - 1; i >= 0; --i)
      {
         j = svec.index(i);
         if (v[j])
         {
            x = v[j] + xx * svec.value(i);
            if (isNotZero(x, epsilon))
               v[j] = x;
            else
            {
               adjust = 1;
               v[j] = MARKER;
            }
         }
         else
         {
            x = xx * svec.value(i);
            if (isNotZero(x, epsilon))
            {
               v[j] = x;
               addIdx(j);
            }
         }
      }

      if (adjust)
      {
         int* iptr = idx;
         int* iiptr = idx;
         int* endptr = idx + num;
         for (; iptr < endptr; ++iptr)
         {
            x = v[*iptr];
            if (isNotZero(x, epsilon))
               *iiptr++ = *iptr;
            else
               v[*iptr] = 0;
         }
         num = int(iiptr - idx);
      }
   }
   else
      Vector::multAdd(xx, svec);

   assert(isConsistent());
   return *this;
}

SSVector& SSVector::multAdd(Real x, const Vector& vec)
{
   Vector::multAdd(x, vec);

   if (isSetup())
   {
      setupStatus = false;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator=(const SSVector& rhs)
{
   assert(rhs.isConsistent());

   if (this != &rhs)
   {
      clear();
      epsilon = rhs.epsilon;
      setMax(rhs.max());
      DVector::reDim(rhs.dim());

      if (rhs.isSetup())
      {
         IdxSet::operator=(rhs);

         for(int i = 0; i < size(); ++i)
         {
            int j  = index(i);
            val[j] = rhs.val[j];
         }
      }
      else
      {
         num = 0;

         for(int i = 0; i < rhs.dim(); ++i)
         {
            if (isNotZero(rhs.val[i], epsilon))
            {
               val[i]       = rhs.val[i];
               idx[num]     = i;
               num++;
            }
         }
      }
      setupStatus = true;
   }
   assert(isConsistent());

   return *this;
}

// setup rhs and assign to this
void SSVector::setup_and_assign(SSVector& rhs)
{
   clear();
   epsilon = rhs.epsilon;
   setMax(rhs.max());
   DVector::reDim(rhs.dim());

   if (rhs.isSetup())
   {
      IdxSet::operator=(rhs);
      
      for(int i = 0; i < size(); ++i)
      {
         int j  = index(i);
         val[j] = rhs.val[j];
      }
   }
   else
   {
      num = 0;

      for(int i = 0; i < rhs.dim(); ++i)
      {
         if (rhs.val[i] != 0.0)
         {
            if (isNotZero(rhs.val[i], epsilon))
            {
               rhs.idx[num] = i;
               idx[num]     = i;
               val[i]       = rhs.val[i];
               num++;
            }
            else
            {
               rhs.val[i] = 0.0;
            }
         }
      }
      rhs.num         = num;
      rhs.setupStatus = true;
   }
   setupStatus = true;

   assert(rhs.isConsistent());
   assert(isConsistent());
}


SSVector& SSVector::operator=(const SVector& rhs)
{
   clear();
   return assign(rhs);
}

// @todo implementation of SSVector::assign() could be put into operator=()
SSVector& SSVector::assign(const SVector& rhs)
{
   assert(rhs.dim() <= Vector::dim());

   num = 0;

   for(int i = 0; i < rhs.size(); ++i)
   {
      int  k = rhs.index(i);
      Real v = rhs.value(i);

      if (isZero(v, epsilon))
         val[k] = REAL(0.0);
      else
      {
         val[k]     = v;
         idx[num++] = k;
      }
   }
   setupStatus = true;

   assert(isConsistent());

   return *this;
}

SSVector& SSVector::assign2product1(const SVSet& A, const SSVector& x)
{
   assert(x.isSetup());
   assert(x.size() == 1);

   // get the nonzero value of x and the corresponding vector in A:
   const int      nzidx = x.idx[0];
   const Real     nzval = x.val[nzidx];
   const SVector& Ai    = A    [nzidx];

   // compute A[nzidx] * nzval:
   if ( isZero(nzval, epsilon) || Ai.size() == 0 )
      clear();    // this := zero vector
   else 
   {
      num = Ai.size();
      for ( register int j = 0; j < num; j++ )
      {
         const SVector::Element& Aij = Ai.element(j);
         idx[j]       = Aij.idx;
         val[Aij.idx] = nzval * Aij.val;
      }
   }

   assert(isConsistent());
   return *this;
}

SSVector& SSVector::assign2productShort(const SVSet& A, const SSVector& x)
{
   assert(x.isSetup());

   if (x.size() == 0) // x can be setup but have size 0 => this := zero vector
   {
      clear();
      return *this;
   }

   // compute x[0] * A[0]
   int            curidx = x.idx[0];
   const Real     x0     = x.val[curidx];
   const SVector& A0     = A    [curidx];
   int            nonzero_idx = 0;

   num = A0.size();
   if ( isZero(x0, epsilon) || num == 0 )
   {
      // A[0] == 0 or x[0] == 0 => this := zero vector
      clear();
   }
   else 
   {
      for ( register int j = 0; j < num; ++j )
      {
         const SVector::Element& elt     = A0.element(j);
         const Real              product = x0 * elt.val;

         idx[ nonzero_idx ] = elt.idx;
         val[ elt.idx ]     = product; // store the value in any case
         if ( product != 0 )           // not 'isNotZero(product, epsilon)'
            ++nonzero_idx;             // count only non-zero values
      }
   }

   // Compute the other x[i] * A[i] and add them to the existing vector.
   for ( register int i = 1; i < x.size(); ++i )
   {
      curidx                = x.idx[i];
      const Real     xi     = x.val[curidx];
      const SVector& Ai     = A    [curidx];

      // If A[i] == 0 or x[i] == 0, do nothing.
      if ( isNotZero(xi, epsilon) || Ai.size() == 0 )
      {
         // Compute x[i] * A[i] and add it to the existing vector.
         for ( register int j = 0; j < Ai.size(); ++j )
         {
            const SVector::Element& elt  = Ai.element(j);
            idx[ nonzero_idx ]           = elt.idx;
            Real                 oldval  = val[elt.idx];

            // An old value of exactly 0 means the position is still unused.
            // It will be used now (either by a new nonzero or by a MARKER), 
            // so increase the counter. If oldval != 0, we just 
            // change an existing NZ-element, so don't increase the counter.
            if ( oldval == 0 )
               ++nonzero_idx;

            // Add the current product x[i] * A[i][j]; if oldval was
            // MARKER before, it does not hurt because MARKER is really small.
            oldval += xi * elt.val;

            // If the new value is exactly 0, mark the index as used
            // by setting a value which is nearly 0; otherwise, store 
            // the value. Values below epsilon will be removed later.
            if ( oldval == 0 )
               val[ elt.idx ] = MARKER;
            else
               val[ elt.idx ] = oldval;
         }
      }
   }

   // Clean up by shifting all nonzeros (w.r.t. epsilon) to the front of idx, 
   // zeroing all values which are nearly 0, and setting #num# appropriately.
   int nz_counter = 0;
   for ( register int i = 0; i < nonzero_idx; ++i )
   {
      curidx = idx[i];
      if ( isZero( val[curidx], epsilon ) )
         val[curidx] = 0;
      else
      {  
         idx[nz_counter] = curidx;
         ++nz_counter;
      }
      num = nz_counter;
   }

   assert(isConsistent());
   return *this;
}

SSVector& SSVector::assign2productFull(const SVSet& A, const SSVector& x)
{
   assert(x.isSetup());

   if (x.size() == 0) // x can be setup but have size 0 => this := zero vector
   {
      clear();
      return *this;
   }

   bool A_is_zero = true;
   for ( int i = 0; i < x.size(); ++i )
   {
      const int      curidx = x.idx[i];
      const Real     xi     = x.val[curidx];
      const SVector& Ai     = A    [curidx];

      if ( A_is_zero && Ai.size() > 0 )
         A_is_zero = false;

      for ( register int j = 0; j < Ai.size(); ++j )
      {
         const SVector::Element& elt  = Ai.element(j);
         val[ elt.idx ] += xi * elt.val;
      }
   }

   if ( A_is_zero )
      clear();       // case x != 0 but A == 0

   return *this;
}

SSVector& SSVector::assign2product4setup(const SVSet& A, const SSVector& x)
{
   assert(A.num() == x.dim());
   assert(x.isSetup());
   clear();

   if (x.size() == 1)
   {
      assign2product1(A, x);
      setupStatus = true;
   }

   else if (Real(x.size()) * A.memSize() <= shortProductFactor * dim() * A.num()
             && isSetup())
   {
      assign2productShort(A, x);
      setupStatus = true;
   }

   else
   {
      assign2productFull(A, x);
      setupStatus = false;
   }
   assert(isConsistent());

   return *this;
}

SSVector& SSVector::assign2product(const SSVector& x, const SVSet& A)
{
   assert(A.num() == dim());

   Real y;

   clear();

   for (int i = dim(); i-- > 0;)
   {
      y = A[i] * x;

      if (isNotZero(y, epsilon))
      {
         val[i] = y;
         IdxSet::addIdx(i);
      }
   }
   assert(isConsistent());

   return *this;
}

SSVector& SSVector::assign2productAndSetup(const SVSet& A, SSVector& x)
{
   if (x.isSetup())
      return assign2product4setup(A, x);

#if 0   // buggy (missing test for 'svec.size() == 0' and 'x.dim() == 0' )
   int*  xi  = x.idx;
   Real* xv  = x.val;
   Real* end = xv + x.dim() - 1;

   /* setze weissen Elefanten */
   Real lastval = *end;
   *end = MARKER;

   for(;;)
   {
      while (!*xv)
         ++xv;
      if (isNotZero(*xv, epsilon))
      {
         const SVector& svec = A[ *xi++ = int(xv - x.val) ];
         const SVector::Element* elem = &(svec.element(0));
         const SVector::Element* last = elem + svec.size();
         for (; elem < last; ++elem)
            val[elem->idx] += *xv * elem->val;
      }
      else
      {
         *xv = 0;
         if (xv == end)
            break;
      }
      xv++;
   }

   /* fange weissen Elefanten wieder ein */
   if (isNotZero(lastval, epsilon))
   {
      *xv = lastval;
      const SVector& svec = A[ *xi++ = int(xv - x.val) ];
      const SVector::Element* elem = &(svec.element(0));
      const SVector::Element* last = elem + svec.size();
      for (; elem < last; ++elem)
         val[elem->idx] += *xv * elem->val;
   }
   else
      *xv = 0;

   x.num = int(xi - x.idx);

#else

   if (x.dim() == 0) // x == 0 => this := zero vector
      clear();

   // x is not setup, so walk through its value vector
   int nzcount = 0;
   for ( register int i = 0; i < x.dim(); ++i )
   {
      // advance to the next element != 0
      Real& xi = x.val[i];
      if (xi == 0)
         continue;

      // If x[i] is really nonzero, compute A[i] * x[i] and adapt x.idx,
      // otherwise set x[i] to 0.
      if (isNotZero(xi, epsilon))
      {
         x.idx[ nzcount++ ]    = i;
         const SVector& Ai     = A[i];
         const int      Aisize = Ai.size();
         if ( Aisize > 0 ) // otherwise: Ai == 0 => do nothing
         {
            for ( register int j = 0; j < Aisize; ++j )
            {
               const SVector::Element& elt = Ai.element(j);
               val[elt.idx] += xi * elt.val;
            }
         }
      }
      else
      {
         xi = 0;
      }
   }
   x.num = nzcount;

#endif

   x.setupStatus = true;
   setupStatus   = false;
   assert(isConsistent());

   return *this;
}


bool SSVector::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (Vector::dim() > IdxSet::max())
      return MSGinconsistent("SSVector");

   if (Vector::dim() < IdxSet::dim())
      return MSGinconsistent("SSVector");

   if (isSetup())
   {
      for (int i = 0; i < Vector::dim(); ++i)
      {
         int j = number(i);

         if (j < 0 && fabs(val[i]) > 0.0) 
         {
            MSG_ERROR( spxout << "ESSVEC01 i = " << i 
                              << "\tidx = " << j 
                              << "\tval = " << std::setprecision(16) << val[i] 
                              << std::endl; )
            return MSGinconsistent("SSVector");
         }
      }
   }
   return DVector::isConsistent() && IdxSet::isConsistent();
#else
   return true;
#endif
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

