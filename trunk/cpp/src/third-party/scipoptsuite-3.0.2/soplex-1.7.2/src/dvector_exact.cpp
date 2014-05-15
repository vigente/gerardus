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


#include "spxdefines.h"
#include "dvector_exact.h"
#include "spxalloc.h"

namespace soplex
{

void DVector_exact::reSize(int newsize)
{
   reSize(newsize, dim());
}

void DVector_exact::reSize(int newsize, int newdim)
{
   assert(newsize >= newdim);

   if( newsize > memsize )
   {
      MpqReal* olddata = mem;

      mem = new MpqReal[newsize]();
      assert(mem != 0);

      if( dimen > 0 )
      {
         for( int i = 0; i < dimen; i++ )
            mem[i] = olddata[i];

         delete[] olddata;
      }
   }

   val = mem;
   memsize = newsize;
   dimen = newdim;
}

void DVector_exact::reDim(int newdim)
{
   assert(memsize >= 0);

   if( newdim > memsize )
   {
      reSize(int(newdim + 0.2 * memsize));
   }

   for( int i = dimen; i < newdim; i++ )
      mem[i] = 0;

   dimen = newdim;
}

DVector_exact::DVector_exact(const DVector_exact& old)
   : Vector_exact(0, 0)
   , mem(0)
{
   dimen = old.dim();
   memsize = old.memsize;

   mem = new MpqReal[memsize]();
   assert(mem != 0);

   val = mem;
   *this = old;

   assert(DVector_exact::isConsistent());
}

DVector_exact::DVector_exact(const Vector_exact& old)
   : Vector_exact(0, 0)
   , mem(0)
{
   dimen = old.dim();
   memsize = dimen;

   mem = new MpqReal[memsize]();
   assert(mem != 0);

   val = mem;
   *this = old;

   assert(DVector_exact::isConsistent());
}

DVector_exact::DVector_exact(const Vector& old)
   : Vector_exact(0, 0)
   , mem(0)
{
   dimen = old.dim();
   memsize = dimen;

   mem = new MpqReal[memsize]();
   assert(mem != 0);

   val = mem;
   Vector_exact::operator=(old);

   assert(DVector_exact::isConsistent());
}

DVector_exact::DVector_exact(int p_dim)
   : Vector_exact(0, 0)
   , mem(0)
{
   memsize = (p_dim > 0) ? p_dim : 4;

   mem = new MpqReal[memsize]();
   assert(mem != 0);

   val = mem;
   dimen = p_dim;

   assert(DVector_exact::isConsistent());
}

DVector_exact::~DVector_exact()
{
   if( mem != 0 )
      delete [] mem;
}

bool DVector_exact::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if( val != mem || dimen > memsize || dimen < 0 )
      return MSGinconsistent("DVector_exact");

   return Vector_exact::isConsistent();
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
