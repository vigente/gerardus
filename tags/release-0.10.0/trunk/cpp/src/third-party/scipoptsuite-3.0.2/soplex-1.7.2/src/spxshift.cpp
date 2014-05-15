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

//#define DEBUGGING 1

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "spxsolver.h"
#include "spxout.h"

namespace soplex
{
void SPxSolver::shiftFvec()
{
   METHOD( "SPxSolver::shiftFvec()" );

   /* the allowed tolerance is (rep() == COLUMN) ? feastol() : opttol() because theFvec is the primal vector in COLUMN
    * and the dual vector in ROW representation; this is equivalent to entertol()
    */
   Random mult(10.0 * entertol(), 100.0 * entertol());
   Real allow = entertol() - epsilon();

   assert(type() == ENTER);
   assert(allow > 0);

   for (int i = dim() - 1; i >= 0; --i)
   {
      if (theUBbound[i] + allow < (*theFvec)[i])
      {
         MSG_DEBUG( spxout << "DSHIFT08 theUBbound[" << i << "] violated by " << (*theFvec)[i] - theUBbound[i] - allow << std::endl );

         if (theUBbound[i] != theLBbound[i])
            shiftUBbound(i, (*theFvec)[i] + Real(mult));
         else
         {
            shiftUBbound(i, (*theFvec)[i]);
            theLBbound[i] = theUBbound[i];
         }
      }
      else if ((*theFvec)[i] < theLBbound[i] - allow)
      {
         MSG_DEBUG( spxout << "DSHIFT08 theLBbound[" << i << "] violated by " << theLBbound[i] - (*theFvec)[i] - allow << std::endl );

         if (theUBbound[i] != theLBbound[i])
            shiftLBbound(i, (*theFvec)[i] - Real(mult));
         else
         {
            shiftLBbound(i, (*theFvec)[i]);
            theUBbound[i] = theLBbound[i];
         }
      }
   }

#ifndef NDEBUG
   testBounds();
   MSG_DEBUG( spxout << "DSHIFT01 shiftFvec: OK" << std::endl; )
#endif
}

// -----------------------------------------------------------------

/*
    This methods assumes correctly setup vectors |pVec| and |coPvec| and bound
    vectors for leaving simplex. Then it checks all values of |pVec| and
    |coPvec| to obey these bounds and enlarges them if neccessary.
 */
void SPxSolver::shiftPvec()
{
   METHOD( "SPxSolver::shiftPvec()" );

   /* the allowed tolerance is (rep() == ROW) ? feastol() : opttol() because thePvec is the primal vector in ROW and the
    * dual vector in COLUMN representation; this is equivalent to leavetol()
    */
   Random mult(10.0 * leavetol(), 100.0 * leavetol());
   Real allow = leavetol() - epsilon();
   int i, tmp;

   assert(type() == LEAVE);
   assert(allow > 0.0);

   for (i = dim() - 1; i >= 0; --i)
   {
      tmp = !isBasic(coId(i));
      if ((*theCoUbound)[i] + allow <= (*theCoPvec)[i] && tmp)
      {
         if ((*theCoUbound)[i] != (*theCoLbound)[i])
            shiftUCbound(i, (*theCoPvec)[i] + Real(mult));
         else
         {
            shiftUCbound(i, (*theCoPvec)[i]);
            (*theCoLbound)[i] = (*theCoUbound)[i];
         }
      }
      else if ((*theCoLbound)[i] - allow >= (*theCoPvec)[i] && tmp)
      {
         if ((*theCoUbound)[i] != (*theCoLbound)[i])
            shiftLCbound(i, (*theCoPvec)[i] - Real(mult));
         else
         {
            shiftLCbound(i, (*theCoPvec)[i]);
            (*theCoUbound)[i] = (*theCoLbound)[i];
         }
      }
   }

   for (i = coDim() - 1; i >= 0; --i)
   {
      tmp = !isBasic(id(i));
      if ((*theUbound)[i] + allow <= (*thePvec)[i] && tmp)
      {
         if ((*theUbound)[i] != (*theLbound)[i])
            shiftUPbound(i, (*thePvec)[i] + Real(mult));
         else
         {
            shiftUPbound(i, (*thePvec)[i]);
            (*theLbound)[i] = (*theUbound)[i];
         }
      }
      else if ((*theLbound)[i] - allow >= (*thePvec)[i] && tmp)
      {
         if ((*theUbound)[i] != (*theLbound)[i])
            shiftLPbound(i, (*thePvec)[i] - Real(mult));
         else
         {
            shiftLPbound(i, (*thePvec)[i]);
            (*theUbound)[i] = (*theLbound)[i];
         }
      }
   }

#ifndef NDEBUG
   testBounds();
   MSG_DEBUG( spxout << "DSHIFT02 shiftPvec: OK" << std::endl; )
#endif
}
// -----------------------------------------------------------------
void SPxSolver::perturbMin(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   Real eps,
   Real delta,
   int start,
   int incr)
{
   METHOD( "SPxSolver::perturbMin()" );
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const Real* vec = uvec.get_const_ptr();
   const Real* upd = uvec.delta().values();
   const IdxSet& idx = uvec.delta().indices();
   Random mult(10.0 * delta, 100.0 * delta);
   Real x, l, u;
   int i, j;

#ifdef  FULL_SHIFT
   eps = delta;

   for (i = uvec.dim() - start - 1; i >= 0; i -= incr)
   {
      u = p_up[i];
      l = p_low[i];

      if (p_up[i] <= vec[i] + eps)
      {
         p_up[i] = vec[i] + Real(mult);
         theShift += p_up[i] - u;
      }
      if (p_low[i] >= vec[i] - eps)
      {
         p_low[i] = vec[i] - Real(mult);
         theShift -= p_low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = p_up[i];
      l = p_low[i];

      if (x < epsilon())
      {
         if (u != l && vec[i] >= u - eps)
         {
            p_up[i] = vec[i] + Real(mult);
            theShift += p_up[i] - u;
         }
      }
      else if (x > epsilon())
      {
         if (u != l && vec[i] <= l + eps)
         {
            p_low[i] = vec[i] - Real(mult);
            theShift -= p_low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT
}
// -----------------------------------------------------------------
void SPxSolver::perturbMax(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   Real eps,
   Real delta,
   int start,
   int incr) 
{
   METHOD( "SPxSolver::perturbMax()" );
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const Real* vec = uvec.get_const_ptr();
   const Real* upd = uvec.delta().values();
   const IdxSet& idx = uvec.delta().indices();
   Random mult(10.0 * delta, 100.0 * delta);
   Real x, l, u;
   int i, j;

#ifdef  FULL_SHIFT
   eps = delta;
   for (i = uvec.dim() - start - 1; i >= 0; i -= incr)
   {
      u = p_up[i];
      l = p_low[i];
      if (p_up[i] <= vec[i] + eps)
      {
         p_up[i] = vec[i] + Real(mult);
         theShift += p_up[i] - u;
      }
      if (p_low[i] >= vec[i] - eps)
      {
         p_low[i] = vec[i] - Real(mult);
         theShift -= p_low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = p_up[i];
      l = p_low[i];
      if (x > epsilon())
      {
         if (u != l && vec[i] >= u - eps)
         {
            p_up[i] = vec[i] + Real(mult);
            theShift += p_up[i] - u;
         }
      }
      else if (x < epsilon())
      {
         if (u != l && vec[i] <= l + eps)
         {
            p_low[i] = vec[i] - Real(mult);
            theShift -= p_low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT
}

void SPxSolver::perturbMinEnter(void)
{
   METHOD( "SPxSolver::perturbMinEnter()" );
   MSG_DEBUG( spxout << "DSHIFT03 iteration= " << iteration() << ": perturbing " << shift(); )
   fVec().delta().setup();
   perturbMin(fVec(), lbBound(), ubBound(), epsilon(), entertol());
   MSG_DEBUG( spxout << "\t->" << shift() << std::endl; )
}


void SPxSolver::perturbMaxEnter(void)
{
   METHOD( "SPxSolver::perturbMaxEnter()" );
   MSG_DEBUG( spxout << "DSHIFT04 iteration= " << iteration() << ": perturbing " << shift(); )
   fVec().delta().setup();
   perturbMax(fVec(), lbBound(), ubBound(), epsilon(), entertol());
   MSG_DEBUG( spxout << "\t->" << shift() << std::endl; )
}


Real SPxSolver::perturbMin(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   Real eps,
   Real p_delta,
   const SPxBasis::Desc::Status* stat,
   int start,
   int incr) const
{
   METHOD( "SPxSolver::perturbMin()" );
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const Real* vec = uvec.get_const_ptr();
   const Real* upd = uvec.delta().values();
   const IdxSet& idx = uvec.delta().indices();
   Random mult(10*p_delta, 100*p_delta);
   Real x, l, u;
   int i, j;
   Real l_theShift = 0;

#ifdef  FULL_SHIFT
   eps = p_delta;
   for (i = uvec.dim() - start - 1; i >= 0; i -= incr)
   {
      u = p_up[i];
      l = p_low[i];
      if (p_up[i] <= vec[i] + eps && rep()*stat[i] < 0)
      {
         p_up[i] = vec[i] + Real(mult);
         l_theShift += p_up[i] - u;
      }
      if (p_low[i] >= vec[i] - eps && rep()*stat[i] < 0)
      {
         p_low[i] = vec[i] - Real(mult);
         l_theShift -= p_low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = p_up[i];
      l = p_low[i];
      if (x < eps)
      {
         if (u != l && vec[i] >= u - eps && rep()*stat[i] < 0)
         {
            p_up[i] = vec[i] + Real(mult);
            l_theShift += p_up[i] - u;
         }
      }
      else if (x > eps)
      {
         if (u != l && vec[i] <= l + eps && rep()*stat[i] < 0)
         {
            p_low[i] = vec[i] - Real(mult);
            l_theShift -= p_low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT 
   return l_theShift;
}

Real SPxSolver::perturbMax(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   Real eps,
   Real p_delta,
   const SPxBasis::Desc::Status* stat,
   int start,
   int incr) const
{
   METHOD( "SPxSolver::perturbMax()" );
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const Real* vec = uvec.get_const_ptr();
   const Real* upd = uvec.delta().values();
   const IdxSet& idx = uvec.delta().indices();
   Random mult(10*p_delta, 100*p_delta);
   Real x, l, u;
   int i, j;
   Real l_theShift = 0;

#ifdef  FULL_SHIFT
   eps = p_delta;
   for (i = uvec.dim() - start - 1; i >= 0; i -= incr)
   {
      u = p_up[i];
      l = p_low[i];
      if (p_up[i] <= vec[i] + eps && rep()*stat[i] < 0)
      {
         p_up[i] = vec[i] + Real(mult);
         l_theShift += p_up[i] - u;
      }
      if (p_low[i] >= vec[i] - eps && rep()*stat[i] < 0)
      {
         p_low[i] = vec[i] - Real(mult);
         l_theShift -= p_low[i] - l;
      }
   }

#else   // !FULL_SHIFT
   for (j = uvec.delta().size() - start - 1; j >= 0; j -= incr)
   {
      i = idx.index(j);
      x = upd[i];
      u = p_up[i];
      l = p_low[i];
      if (x > eps)
      {
         if (u != l && vec[i] >= u - eps && rep()*stat[i] < 0)
         {
            p_up[i] = vec[i] + Real(mult);
            l_theShift += p_up[i] - u;
         }
      }
      else if (x < eps)
      {
         if (u != l && vec[i] <= l + eps && rep()*stat[i] < 0)
         {
            p_low[i] = vec[i] - Real(mult);
            l_theShift -= p_low[i] - l;
         }
      }
   }
#endif  // !FULL_SHIFT 
   return l_theShift;
}


void SPxSolver::perturbMinLeave(void)
{
   METHOD( "SPxSolver::perturbMinLeave()" );
   MSG_DEBUG( spxout << "DSHIFT05 iteration= " << iteration() << ": perturbing " << shift(); )
   pVec().delta().setup();
   coPvec().delta().setup();
   theShift += perturbMin(pVec(), lpBound(), upBound(), epsilon(), leavetol(),
      desc().status(), 0, 1);
   theShift += perturbMin(coPvec(), lcBound(), ucBound(), epsilon(), leavetol(),
      desc().coStatus(), 0, 1);
   MSG_DEBUG( spxout << "\t->" << shift() << std::endl; )
}


void SPxSolver::perturbMaxLeave(void)
{
   METHOD( "SPxSolver::perturbMaxLeave()" );
   MSG_DEBUG( spxout << "DSHIFT06 iteration= " << iteration() << ": perturbing " << shift(); )
   pVec().delta().setup();
   coPvec().delta().setup();
   theShift += perturbMax(pVec(), lpBound(), upBound(), epsilon(), leavetol(),
      desc().status(), 0, 1);
   theShift += perturbMax(coPvec(), lcBound(), ucBound(), epsilon(), leavetol(),
      desc().coStatus(), 0, 1);
   MSG_DEBUG( spxout << "\t->" << shift() << std::endl; )
}


void SPxSolver::unShift(void)
{
   METHOD( "SPxSolver::unShift()" );
   MSG_INFO3( spxout << "DSHIFT07 = " << "unshifting ..." << std::endl; );

   if (isInitialized())
   {
      int i;
      Real t_up, t_low;
      const SPxBasis::Desc& ds = desc();

      theShift = 0;
      if (type() == ENTER)
      {
         Real eps = entertol();

         if (rep() == COLUMN)
         {
            for (i = dim(); i-- > 0;)
            {
               SPxId l_id = baseId(i);
               int l_num = number(l_id);
               if (l_id.type() == SPxId::ROW_ID)
               {
                  t_up = -lhs(l_num);
                  t_low = -rhs(l_num);
               }
               else
               {
                  assert(l_id.type() == SPxId::COL_ID);
                  t_up = upper(l_num);
                  t_low = lower(l_num);
               }
               if (t_up != t_low)
               {
                  if ((*theFvec)[i] < t_up - eps)
                     theUBbound[i] = t_up;
                  else if ((*theFvec)[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  if ((*theFvec)[i] > t_low + eps)
                     theLBbound[i] = t_low;
                  else if ((*theFvec)[i] < t_low)
                     theShift -= theLBbound[i] - t_low;
               }
               else
               {
                  if (theUBbound[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  else if (theLBbound[i] < t_low)
                     theShift += t_low - theLBbound[i];
               }
            }
            for (i = nRows(); i-- > 0;)
            {
               if (!isBasic(ds.rowStatus(i)))
               {
                  t_up = -lhs(i);
                  t_low = -rhs(i);
                  if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  if (t_low > theLRbound[i])
                     theShift += t_low - theLRbound[i];
               }
            }
            for (i = nCols(); i-- > 0;)
            {
               if (!isBasic(ds.colStatus(i)))
               {
                  t_up = upper(i);
                  t_low = lower(i);
                  if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  if (t_low > theLCbound[i])
                     theShift += t_low - theLCbound[i];
               }
            }
         }
         else
         {
            assert(rep() == ROW);
            for (i = dim(); i-- > 0;)
            {
               SPxId l_id = baseId(i);
               int l_num = number(l_id);
               t_up = t_low = 0;
               if (l_id.type() == SPxId::ROW_ID)
                  clearDualBounds(ds.rowStatus(l_num), t_up, t_low);
               else
                  clearDualBounds(ds.colStatus(l_num), t_up, t_low);
               if (theUBbound[i] != theLBbound[i])
               {
                  if (theUBbound[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  else
                     theShift -= theUBbound[i] - t_up;
               }
               else
               {
                  /* if the basic (primal or dual) variable is fixed (e.g., basis status P_FREE in row representation)
                   * then shiftFvec() and shiftPvec() do not relax the bounds, but shift both, hence they may be outside
                   * of [t_low,t_up] */
                  assert(theLBbound[i] == theUBbound[i] || theUBbound[i] >= t_up);
                  assert(theLBbound[i] == theUBbound[i] || theLBbound[i] <= t_low);

                  if ((*theFvec)[i] < t_up - eps)
                     theUBbound[i] = t_up;
                  else if ((*theFvec)[i] > t_up)
                     theShift += theUBbound[i] - t_up;

                  if ((*theFvec)[i] > t_low + eps)
                     theLBbound[i] = t_low;
                  else if ((*theFvec)[i] < t_low)
                     theShift -= theLBbound[i] - t_low;
               }
            }
            for (i = nRows(); i-- > 0;)
            {
               if (!isBasic(ds.rowStatus(i)))
               {
                  t_up = t_low = 0;
                  clearDualBounds(ds.rowStatus(i), t_up, t_low);
                  if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  if (t_low > theLRbound[i])
                     theShift += t_low - theLRbound[i];
               }
            }
            for (i = nCols(); i-- > 0;)
            {
               if (!isBasic(ds.colStatus(i)))
               {
                  t_up = t_low = 0;
                  clearDualBounds(ds.colStatus(i), t_up, t_low);
                  if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  if (t_low > theLCbound[i])
                     theShift += t_low - theLCbound[i];
               }
            }
         }
      }
      else
      {
         assert(type() == LEAVE);

         Real eps = leavetol();

         if (rep() == COLUMN)
         {
            for (i = nRows(); i-- > 0;)
            {
               t_up = t_low = 0;
               clearDualBounds(ds.rowStatus(i), t_up, t_low);
               if (!isBasic(ds.rowStatus(i)))
               {
                  if ((*theCoPvec)[i] < t_up - eps)
                  {
                     theURbound[i] = t_up;
                     if( t_up == t_low )
                        theLRbound[i] = t_low; // for fixed rows we change both bounds
                  }
                  else
                     theShift += theURbound[i] - t_up;
                  if ((*theCoPvec)[i] > t_low + eps)
                  {
                     theLRbound[i] = t_low;
                     if( t_up == t_low )
                        theURbound[i] = t_up; // for fixed rows we change both bounds
                  }
                  else
                     theShift += t_low - theLRbound[i];
               }
               else if (theURbound[i] > t_up)
                  theShift += theURbound[i] - t_up;
               else if (theLRbound[i] < t_low)
                  theShift += t_low - theLRbound[i];
            }
            for (i = nCols(); i-- > 0;)
            {
               t_up = t_low = -maxObj(i);
               clearDualBounds(ds.colStatus(i), t_low, t_up);
               if (!isBasic(ds.colStatus(i)))
               {
                  if ((*thePvec)[i] < -t_up - eps)
                  {
                     theUCbound[i] = -t_up;
                     if( t_up == t_low )
                        theLCbound[i] = -t_low; // for fixed variables we change both bounds
                  }
                  else
                     theShift += theUCbound[i] - (-t_up);
                  if ((*thePvec)[i] > -t_low + eps)
                  {
                     theLCbound[i] = -t_low;
                     if( t_up == t_low )
                        theUCbound[i] = -t_up; // for fixed variables we change both bounds
                  }
                  else
                     theShift += (-t_low) - theLCbound[i];
               }
               else if (theUCbound[i] > -t_up)
                  theShift += theUCbound[i] - (-t_up);
               else if (theLCbound[i] < -t_low)
                  theShift += (-t_low) - theLCbound[i];
            }
         }
         else
         {
            assert(rep() == ROW);
            for (i = nRows(); i-- > 0;)
            {
               t_up = rhs(i);
               t_low = lhs(i);
               if (t_up == t_low)
               {
                  if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  else
                     theShift += t_low - theLRbound[i];
               }
               else
                  if (!isBasic(ds.rowStatus(i)))
                  {
                     if ((*thePvec)[i] < t_up - eps)
                        theURbound[i] = t_up;
                     else
                        theShift += theURbound[i] - t_up;
                     if ((*thePvec)[i] > t_low + eps)
                        theLRbound[i] = t_low;
                     else
                        theShift += t_low - theLRbound[i];
                  }
                  else if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  else if (theLRbound[i] < t_low)
                     theShift += t_low - theLRbound[i];
            }
            for (i = nCols(); i-- > 0;)
            {
               t_up = upper(i);
               t_low = lower(i);
               if (t_up == t_low)
               {
                  if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  else
                     theShift += t_low - theLCbound[i];
               }
               else
                  if (!isBasic(ds.colStatus(i)))
                  {
                     if ((*theCoPvec)[i] < t_up - eps)
                        theUCbound[i] = t_up;
                     else
                        theShift += theUCbound[i] - t_up;
                     if ((*theCoPvec)[i] > t_low + eps)
                        theLCbound[i] = t_low;
                     else
                        theShift += t_low - theLCbound[i];
                  }
                  else if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  else if (theLCbound[i] < t_low)
                     theShift += t_low - theLCbound[i];
            }
         }
      }
   }
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
