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
#include <stdio.h>

#include "spxdefines.h"
#include "spxfastrt.h"


/*
    Here comes our implementation of the Harris procedure improved by shifting
    bounds. The basic idea is to allow a slight infeasibility within |delta| to
    allow for more freedom when selecting the leaving variable. This freedom
    may then be used for selecting numerically stable variables yielding great
    improvements.

    The algorithms operates in two phases. In a first phase, the maximum value
    |val| is determined, when infeasibility within |delta| is allowed. In the second
    phase, between all variables with value < |val| the one is selected which
    has the largest update vector component. However, this may not
    always yield an improvement. In that case, we shift the variable towards
    infeasibility and retry. This avoids cycling in the shifted LP.
 */

namespace soplex
{

#define MINSTAB         1e-5
#define LOWSTAB         1e-10
#define TRIES           2
#define SHORT           1e-5
#define DELTA_SHIFT     1e-5
#define EPSILON         1e-10

void SPxFastRT::resetTols()
{
   // epsilon = thesolver->epsilon();
   epsilon = EPSILON;
   /*
       if(thesolver->basis().stability() < 1e-4)
           epsilon *= 1e-4 / thesolver->basis().stability();
   */
}

void SPxFastRT::tighten()
{
   /*
   if((delta > 1.99 * DELTA_SHIFT  &&  thesolver->theShift < 1e-4) ||
       (delta > 1e-4   &&  thesolver->theShift > 1e-4))
    */
   // if(delta > 1.99 * DELTA_SHIFT)
   if (fastDelta >= delta + DELTA_SHIFT)
   {
      fastDelta -= DELTA_SHIFT;
      if (fastDelta > 1e-4)
         fastDelta -= 2 * DELTA_SHIFT;
   }

   if (minStab < MINSTAB)
   {
      minStab /= 0.90;
      if (minStab < 1e-6)
         minStab /= 0.90;
   }
}

void SPxFastRT::relax()
{
   minStab *= 0.95;
   fastDelta += 3 * DELTA_SHIFT;
   // delta   += 2 * (thesolver->theShift > delta) * DELTA_SHIFT;
}

Real SPxFastRT::minStability(Real maxabs)
{
   if (maxabs < 1000.0)
      return minStab;
   return maxabs*minStab / 1000.0;
}

/* The code below implements phase 1 of the ratio test. It differs from the description in the
 * Ph.D. thesis page 57 as follows: It uses \f$\delta_i = d_i - s_i - \delta\f$ if \f$d_i > s_i\f$.
 *
 * This change leads to the following behavior of the source code. Consider the first case (x >
 * epsilon, u < infinity): If u - vec[i] <= 0, vec[i] violates the upper bound. In the Harris ratio
 * test, we would compute (u - vec[i] + delta)/upd[i]. The code computes instead delta/upd[i].
 */
int SPxFastRT::maxDelta(
   Real& val,                                /* on return: maximum step length */
   Real& maxabs,                             /* on return: maximum absolute value in upd vector */
   UpdateVector& update,
   const Vector& lowBound,
   const Vector& upBound,
   int start,
   int incr) const
{
   int i, sel;
   Real x, y, max;
   Real u, l;
   bool leaving = m_type == SPxSolver::LEAVE;

   Real mabs = maxabs;

   const Real* up = upBound.get_const_ptr();
   const Real* low = lowBound.get_const_ptr();
   const Real* vec = update.get_const_ptr();
   const Real* upd = update.delta().values();
   const int* idx = update.delta().indexMem();

   sel = -1;
   max = val;

   if (update.delta().isSetup())
   {
      const int* last = idx + update.delta().size();
      for (idx += start; idx < last; idx += incr)
      {
         i = *idx;
         x = upd[i];

         /* in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables */
         if( leaving && ((iscoid && thesolver->isCoBasic(i)) || (!iscoid && thesolver->isBasic(i))) )
            continue;

         if (x > epsilon)
         {
            // @todo check wether mabs should be computed only over bounded vars, i.e., in the if block below
            mabs = (x > mabs) ? x : mabs;
            u = up[i];
            if (u < infinity)
            {
               y = u - vec[i];
               if (y <= 0)
                  x = fastDelta / x;
               else
                  x = (y + fastDelta) / x;

               if (x < max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
         else if (x < -epsilon)
         {
            // @todo check wether mabs should be computed only over bounded vars, i.e., in the if block below
            mabs = (-x > mabs) ? -x : mabs;
            l = low[i];
            if (l > -infinity)
            {
               y = l - vec[i];
               if ( y >= 0 )
                  x = - fastDelta / x;
               else
                  x = ( y - fastDelta ) / x;

               if (x < max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
      }
   }
   else
   {
      /* In this case, the indices of the semi-sparse vector update.delta() are not set up and are filled below. */
      int* l_idx = update.delta().altIndexMem();
      Real* uval = update.delta().altValues();
      const Real* uend = uval + update.dim();
      assert( uval == upd );

      for (i = 0; uval < uend; ++uval, ++i)
      {
         x = *uval;
         if (x)
         {
            if (x >= -epsilon && x <= epsilon)
            {
               *uval = 0.0;
               continue;
            }
            else
               *l_idx++ = i;

            /* in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables */
            if( leaving && ((iscoid && thesolver->isCoBasic(i)) || (!iscoid && thesolver->isBasic(i))) )
               continue;

            if (x > epsilon)
            {
               mabs = (x > mabs) ? x : mabs;
               u = up[i];
               if (u < infinity)
               {
                  y = u - vec[i];
                  if (y <= 0)
                     x = fastDelta / x;
                  else
                     x = (y + fastDelta) / x;

                  if (x < max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
            else if (x < -epsilon)
            {
               mabs = (-x > mabs) ? -x : mabs;
               l = low[i];
               if (l > -infinity)
               {
                  y = l - vec[i];
                  if ( y >= 0 )
                     x = - fastDelta / x;
                  else
                     x = ( y - fastDelta ) / x;

                  if (x < max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
         }
      }
      update.delta().setSize(int(l_idx - update.delta().indexMem()));
      update.delta().forceSetup();
   }

   val = max;
   maxabs = mabs;

   return sel;
}

/* See maxDelta() */
int SPxFastRT::minDelta(
   Real& val,
   Real& maxabs,
   UpdateVector& update,
   const Vector& lowBound,
   const Vector& upBound,
   int start,
   int incr) const
{
   int i, sel;
   Real x, y, max;
   Real u, l;
   bool leaving = m_type == SPxSolver::LEAVE;

   Real mabs = maxabs;

   const Real* up = upBound.get_const_ptr();
   const Real* low = lowBound.get_const_ptr();
   const Real* vec = update.get_const_ptr();
   const Real* upd = update.delta().values();
   const int* idx = update.delta().indexMem();

   sel = -1;
   max = val;

   if (update.delta().isSetup())
   {
      const int* last = idx + update.delta().size();
      for (idx += start; idx < last; idx += incr)
      {
         i = *idx;
         x = upd[i];

         /* in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables */
         if( leaving && ((iscoid && thesolver->isCoBasic(i)) || (!iscoid && thesolver->isBasic(i))) )
            continue;

         if (x > epsilon)
         {
            // @todo check wether mabs should be computed only over bounded vars, i.e., in the if block below
            mabs = (x > mabs) ? x : mabs;
            l = low[i];
            if (l > -infinity)
            {
               y = l - vec[i];
               if ( y >= 0 )
                  x = - fastDelta / x;
               else
                  x = ( y - fastDelta ) / x;

               if (x > max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
         else if (x < -epsilon)
         {
            // @todo check wether mabs should be computed only over bounded vars, i.e., in the if block below
            mabs = (-x > mabs) ? -x : mabs;
            u = up[i];
            if (u < infinity)
            {
               y = u - vec[i];
               if (y <= 0)
                  x = fastDelta / x;
               else
                  x = (y + fastDelta) / x;

               if (x > max)
               {
                  max = x;
                  sel = i;
               }
            }
         }
      }
   }
   else
   {
      /* In this case, the indices of the semi-sparse vector update.delta() are not set up and are filled below. */
      int* l_idx = update.delta().altIndexMem();
      Real* uval = update.delta().altValues();
      const Real* uend = uval + update.dim();
      assert( uval == upd );

      for (i = 0; uval < uend; ++uval, ++i)
      {
         x = *uval;

         if (x)
         {
            if (x >= -epsilon && x <= epsilon)
            {
               *uval = 0.0;
               continue;
            }
            else
               *l_idx++ = i;

            /* in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables */
            if( leaving && ((iscoid && thesolver->isCoBasic(i)) || (!iscoid && thesolver->isBasic(i))) )
               continue;

            if (x > epsilon)
            {
               mabs = (x > mabs) ? x : mabs;
               l = low[i];
               if (l > -infinity)
               {
                  y = l - vec[i];
                  if ( y >= 0 )
                     x = - fastDelta / x;
                  else
                     x = ( y - fastDelta ) / x;

                  if (x > max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
            else if (x < -epsilon)
            {
               mabs = (-x > mabs) ? -x : mabs;
               u = up[i];
               if (u < infinity)
               {
                  y = u - vec[i];
                  if (y <= 0)
                     x = fastDelta / x;
                  else
                     x = (y + fastDelta) / x;

                  if (x > max)
                  {
                     max = x;
                     sel = i;
                  }
               }
            }
         }
      }
      update.delta().setSize(int(l_idx - update.delta().indexMem()));
      update.delta().forceSetup();
   }

   val = max;
   maxabs = mabs;

   return sel;
}

int SPxFastRT::maxDelta(
   Real& val,
   Real& maxabs)
{
   return maxDelta(val, maxabs,
      thesolver->fVec(), thesolver->lbBound(), thesolver->ubBound(), 0, 1);
}

int SPxFastRT::minDelta(
   Real& val,
   Real& maxabs)
{
   return minDelta(val, maxabs,
      thesolver->fVec(), thesolver->lbBound(), thesolver->ubBound(), 0, 1);
}

SPxId SPxFastRT::maxDelta(
   int& nr,
   Real& max,                                /* on return: maximum step length */
   Real& maxabs)                             /* on return: maximum absolute value in delta vector */
{
   /* The following cause side effects on coPvec and pVec - both changes may be needed later in
      maxSelect(). We can therefore not move the first function after the (indp >= 0) check. */
   iscoid = true;
   int indc = maxDelta(max, maxabs,
      thesolver->coPvec(), thesolver->lcBound(), thesolver->ucBound(), 0, 1);
   iscoid = false;
   int indp = maxDelta(max, maxabs,
      thesolver->pVec(), thesolver->lpBound(), thesolver->upBound(), 0, 1);

   if (indp >= 0)
   {
      nr = indp;
      return thesolver->id(indp);
   }
   if (indc >= 0)
   {
      nr = indc;
      return thesolver->coId(indc);
   }
   nr = -1;
   return SPxId();
}

SPxId SPxFastRT::minDelta(
   int& nr,
   Real& max,
   Real& maxabs)
{
   /* The following cause side effects on coPvec and pVec - both changes may be needed later in
      minSelect(). We can therefore not move the first function after the (indp >= 0) check. */
   iscoid = true;
   const int indc = minDelta(max, maxabs,
      thesolver->coPvec(), thesolver->lcBound(), thesolver->ucBound(), 0, 1);
   iscoid = false;
   const int indp = minDelta(max, maxabs,
      thesolver->pVec(), thesolver->lpBound(), thesolver->upBound(), 0, 1);

   if (indp >= 0)
   {
      nr = indp;
      return thesolver->id(indp);
   }
   if (indc >= 0)
   {
      nr = indc;
      return thesolver->coId(indc);
   }
   nr = -1;
   return SPxId();
}

/* \p best returns the minimum update value such that the corresponding value of \p upd.delta() is
 * at least \p stab and the update value is smaller than \p max. If no valid update value has been
 * found \p bestDelta returns the slack to the bound corresponding to the index used for \p best. */
int SPxFastRT::minSelect(
   Real& val,
   Real& stab,
   Real& best,
   Real& bestDelta,
   Real max,
   const UpdateVector& update,
   const Vector& lowBound,
   const Vector& upBound,
   int start,
   int incr) const
{
   int i;
   Real x, y;
   bool leaving = m_type == SPxSolver::LEAVE;

   const Real* up = upBound.get_const_ptr();
   const Real* low = lowBound.get_const_ptr();
   const Real* vec = update.get_const_ptr();
   const Real* upd = update.delta().values();
   const int* idx = update.delta().indexMem();
   const int* last = idx + update.delta().size();

   int nr = -1;
   int bestNr = -1;

   for (idx += start; idx < last; idx += incr)
   {
      i = *idx;
      x = upd[i];

      // in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables
      if( leaving && ((iscoid && thesolver->isCoBasic(i)) || (!iscoid && thesolver->isBasic(i))) )
         continue;

      if (x > stab)
      {
         y = (low[i] - vec[i]) / x;

         if (y >= max)
         {
            val = y;
            nr = i;
            stab = x;
         }
         else if (y < best)
         {
            best = y;
            bestNr = i;
         }
      }
      else if (x < -stab)
      {
         y = (up[i] - vec[i]) / x;
         if (y >= max)
         {
            val = y;
            nr = i;
            stab = -x;
         }
         else if (y < best)
         {
            best = y;
            bestNr = i;
         }
      }
   }

   if (nr < 0 && bestNr > 0)
   {
      if (upd[bestNr] < 0)
         bestDelta = up[bestNr] - vec[bestNr];
      else
         bestDelta = vec[bestNr] - low[bestNr];
   }
   return nr;
}

/* See minSelect() */
int SPxFastRT::maxSelect(
   Real& val,
   Real& stab,
   Real& best,
   Real& bestDelta,
   Real max,
   const UpdateVector& update,
   const Vector& lowBound,
   const Vector& upBound,
   int start,
   int incr) const
{
   int i;
   Real x, y;
   bool leaving = m_type == SPxSolver::LEAVE;

   const Real* up = upBound.get_const_ptr();
   const Real* low = lowBound.get_const_ptr();
   const Real* vec = update.get_const_ptr();
   const Real* upd = update.delta().values();
   const int* idx = update.delta().indexMem();
   const int* last = idx + update.delta().size();

   int nr = -1;
   int bestNr = -1;

   for (idx += start; idx < last; idx += incr)
   {
      i = *idx;
      x = upd[i];

      // in the dual algorithm, bound flips cannot happen, hence we only consider nonbasic variables
      if( leaving && ((iscoid && thesolver->isCoBasic(i)) || (!iscoid && thesolver->isBasic(i))) )
         continue;

      if (x > stab)
      {
         y = (up[i] - vec[i]) / x;
         if (y <= max)
         {
            val = y;
            nr = i;
            stab = x;
         }
         else if (y > best)
         {
            best = y;
            bestNr = i;
         }
      }
      else if (x < -stab)
      {
         y = (low[i] - vec[i]) / x;
         if (y <= max)
         {
            val = y;
            nr = i;
            stab = -x;
         }
         else if (y > best)
         {
            best = y;
            bestNr = i;
         }
      }
   }

   if (nr < 0 && bestNr > 0)
   {
      if (upd[bestNr] > 0)
         bestDelta = up[bestNr] - vec[bestNr];
      else
         bestDelta = vec[bestNr] - low[bestNr];
   }

   return nr;
}


int SPxFastRT::maxSelect(
   Real& val,
   Real& stab,
   Real& bestDelta,
   Real max)
{
   Real best = -infinity;
   bestDelta = 0.0;
   return maxSelect(val, stab, best, bestDelta, max,
      thesolver->fVec(), thesolver->lbBound(), thesolver->ubBound(),  0, 1);
}

SPxId SPxFastRT::maxSelect(
   int& nr,
   Real& val,
   Real& stab,
   Real& bestDelta,
   Real max
)
{
   int indp, indc;
   Real best = -infinity;
   bestDelta = 0.0;
   iscoid = true;
   indc = maxSelect(val, stab, best, bestDelta, max,
      thesolver->coPvec(), thesolver->lcBound(), thesolver->ucBound(), 0, 1);
   iscoid = false;
   indp = maxSelect(val, stab, best, bestDelta, max,
      thesolver->pVec(), thesolver->lpBound(), thesolver->upBound(), 0, 1);

   if (indp >= 0)
   {
      nr = indp;
      return thesolver->id(indp);
   }
   if (indc >= 0)
   {
      nr = indc;
      return thesolver->coId(indc);
   }
   nr = -1;
   return SPxId();
}

int SPxFastRT::minSelect(
   Real& val,
   Real& stab,
   Real& bestDelta,
   Real max)
{
   Real best = infinity;
   bestDelta = 0.0;
   return minSelect(val, stab, best, bestDelta, max,
      thesolver->fVec(), thesolver->lbBound(), thesolver->ubBound(), 0, 1);
}

SPxId SPxFastRT::minSelect(
   int& nr,
   Real& val,
   Real& stab,
   Real& bestDelta,
   Real max)
{
   Real best = infinity;
   bestDelta = 0.0;
   iscoid = true;
   int indc = minSelect(val, stab, best, bestDelta, max,
      thesolver->coPvec(), thesolver->lcBound(), thesolver->ucBound(), 0, 1);
   iscoid = false;
   int indp = minSelect(val, stab, best, bestDelta, max,
      thesolver->pVec(), thesolver->lpBound(), thesolver->upBound(), 0, 1);

   if (indp >= 0)
   {
      nr = indp;
      return thesolver->id(indp);
   }
   if (indc >= 0)
   {
      nr = indc;
      return thesolver->coId(indc);
   }
   nr = -1;
   return SPxId();
}


bool SPxFastRT::maxShortLeave(Real& sel, int leave, Real maxabs)
{
   assert(leave >= 0);
   assert(maxabs >= 0);

   sel = thesolver->fVec().delta()[leave];

   if (sel > maxabs*SHORT)
   {
      sel = (thesolver->ubBound()[leave] - thesolver->fVec()[leave]) / sel;
      return true;
   }

   if (sel < -maxabs*SHORT)
   {
      sel = (thesolver->lbBound()[leave] - thesolver->fVec()[leave]) / sel;
      return true;
   }

   return false;
}

bool SPxFastRT::minShortLeave(Real& sel, int leave, Real maxabs)
{
   assert(leave >= 0);
   assert(maxabs >= 0);

   sel = thesolver->fVec().delta()[leave];

   if (sel > maxabs*SHORT)
   {
      sel = (thesolver->lbBound()[leave] - thesolver->fVec()[leave]) / sel;
      return true;
   }

   if ( sel < -maxabs*SHORT)
   {
      sel = (thesolver->ubBound()[leave] - thesolver->fVec()[leave]) / sel;
      return true;
   }

   return false;
}

bool SPxFastRT::maxReLeave(Real& sel, int leave, Real maxabs)
{
   UpdateVector& vec = thesolver->fVec();
   Vector& low = thesolver->lbBound();
   Vector& up = thesolver->ubBound();

   if (leave < 0)
      return true;

   if (up[leave] > low[leave])
   {
      Real x = vec.delta()[leave];

      if (sel < -fastDelta / maxabs)
      {
         sel = 0.0;
         if (x < 0.0)
            thesolver->shiftLBbound(leave, vec[leave]);
         else
            thesolver->shiftUBbound(leave, vec[leave]);
      }
   }
   else
   {
      sel = 0.0;
      thesolver->shiftLBbound(leave, vec[leave]);
      thesolver->shiftUBbound(leave, vec[leave]);
   }

   return false;
}

bool SPxFastRT::minReLeave(Real& sel, int leave, Real maxabs)
{
   UpdateVector& vec = thesolver->fVec();
   Vector& low = thesolver->lbBound();
   Vector& up = thesolver->ubBound();

   if (leave < 0)
      return true;

   if (up[leave] > low[leave])
   {
      Real x = vec.delta()[leave];

      if (sel > fastDelta / maxabs)
      {
         if (x > 0.0)
         {
            thesolver->theShift += low[leave];
            sel = 0.0;
            low[leave] = vec[leave] + sel * x;
            thesolver->theShift -= low[leave];
         }
         else
         {
            thesolver->theShift -= up[leave];
            sel = 0.0;
            up[leave] = vec[leave] + sel * x;
            thesolver->theShift += up[leave];
         }
      }
   }
   else
   {
      sel = 0.0;
      if (vec[leave] < low[leave])
         thesolver->theShift += low[leave] - vec[leave];
      else
         thesolver->theShift += vec[leave] - up[leave];
      low[leave] = up[leave] = vec[leave];
   }

   return false;
}

int SPxFastRT::selectLeave(Real& val, SPxId)
{
   Real maxabs, max, sel;
   int leave = -1;
   int cnt = 0;

   assert( m_type == SPxSolver::ENTER );

   resetTols();

   if (val > epsilon)
   {
      do
      {
         // phase 1:
         max = val;
         maxabs = 0.0;
         leave = maxDelta(max, maxabs);
         if (max == val)
            return -1;

         if (!maxShortLeave(sel, leave, maxabs))
         {
            // phase 2:
            Real stab, bestDelta;
            stab = 100.0 * minStability(maxabs);
            leave = maxSelect(sel, stab, bestDelta, max);
            if (bestDelta < DELTA_SHIFT*TRIES)
               cnt++;
            else
               cnt += TRIES;
         }
         if (!maxReLeave(sel, leave, maxabs))
            break;
         relax();
      }
      while (cnt < TRIES);
   }
   else if (val < -epsilon)
   {
      do
      {
         max = val;
         maxabs = 0;
         leave = minDelta(max, maxabs);

         assert(leave < 0 || !(thesolver->baseId(leave).isSPxColId()) || thesolver->desc().colStatus(thesolver->number(SPxColId(thesolver->baseId(leave)))) != SPxBasis::Desc::P_FIXED);

         if (max == val)
            return -1;

         if (!minShortLeave(sel, leave, maxabs))
         {
            // phase 2:
            Real stab, bestDelta;
            stab = 100.0 * minStability(maxabs);
            leave = minSelect(sel, stab, bestDelta, max);

            assert(leave < 0 || !(thesolver->baseId(leave).isSPxColId()) || thesolver->desc().colStatus(thesolver->number(SPxColId(thesolver->baseId(leave)))) != SPxBasis::Desc::P_FIXED);

            if (bestDelta < DELTA_SHIFT*TRIES)
               cnt++;
            else
               cnt += TRIES;
         }
         if (!minReLeave(sel, leave, maxabs))
            break;
         relax();
      }
      while (cnt < TRIES);
   }
   else
      return -1;

   MSG_DEBUG(
      if (leave >= 0)
         spxout
           << "DFSTRT01 "
           << thesolver->basis().iteration() << "("
           << std::setprecision(6) << thesolver->value() << ","
           << std::setprecision(2) << thesolver->basis().stability() << "):"
           << leave << "\t"
           << std::setprecision(4) << sel << " "
           << std::setprecision(4) << thesolver->fVec().delta()[leave] << " "
           << std::setprecision(6) << maxabs
           << std::endl;
      else
         spxout << "DFSTRT02 " << thesolver->basis().iteration()
                << ": skipping instable pivot" << std::endl;
   )

   if (leave >= 0 || minStab > 2*solver()->epsilon())
   {
      val = sel;
      if (leave >= 0)
         tighten();
   }

   assert(leave < 0 || !(thesolver->baseId(leave).isSPxColId()) || thesolver->desc().colStatus(thesolver->number(SPxColId(thesolver->baseId(leave)))) != SPxBasis::Desc::P_FIXED);

   return leave;
}


bool SPxFastRT::maxReEnter(
   Real& sel,
   Real maxabs,
   const SPxId& id,
   int nr)
{
   Real x, d;
   Vector* up;
   Vector* low;

   UpdateVector& pvec = thesolver->pVec();
   SSVector& pupd = thesolver->pVec().delta();
   Vector& upb = thesolver->upBound();
   Vector& lpb = thesolver->lpBound();
   UpdateVector& cvec = thesolver->coPvec();
   SSVector& cupd = thesolver->coPvec().delta();
   Vector& ucb = thesolver->ucBound();
   Vector& lcb = thesolver->lcBound();

   if (thesolver->isCoId(id))
   {
      if (thesolver->isCoBasic(nr))
      {
         cupd.clearIdx(nr);
         return true;
      }

      x = cvec[nr];
      d = cupd[nr];
      up = &ucb;
      low = &lcb;

      if (d < 0.0)
         sel = (lcb[nr] - cvec[nr]) / d;
      else
         sel = (ucb[nr] - cvec[nr]) / d;
   }
   else if (thesolver->isId(id))
   {
      pvec[nr] = thesolver->vector(nr) * cvec;
      if (thesolver->isBasic(nr))
      {
         pupd.clearIdx(nr);
         return true;
      }

      x = pvec[nr];
      d = pupd[nr];
      up = &upb;
      low = &lpb;

      if (d < 0.0)
         sel = (lpb[nr] - pvec[nr]) / d;
      else
         sel = (upb[nr] - pvec[nr]) / d;
   }
   else
      return true;

   if ((*up)[nr] != (*low)[nr])
   {
      if (sel < -fastDelta / maxabs)
      {
         if (d > 0.0)
         {
            thesolver->theShift -= (*up)[nr];
            sel = 0.0;
            (*up)[nr] = x + sel * d;
            thesolver->theShift += (*up)[nr];
         }
         else
         {
            thesolver->theShift += (*low)[nr];
            sel = 0.0;
            (*low)[nr] = x + sel * d;
            thesolver->theShift -= (*low)[nr];
         }
      }
   }
   else
   {
      sel = 0.0;
      if (x > (*up)[nr])
         thesolver->theShift += x - (*up)[nr];
      else
         thesolver->theShift += (*low)[nr] - x;
      (*up)[nr] = (*low)[nr] = x;
   }

   return false;
}

bool SPxFastRT::minReEnter(
   Real& sel,
   Real maxabs,
   const SPxId& id,
   int nr)
{
   Real x, d;
   Vector* up;
   Vector* low;

   UpdateVector& pvec = thesolver->pVec();
   SSVector& pupd = thesolver->pVec().delta();
   Vector& upb = thesolver->upBound();
   Vector& lpb = thesolver->lpBound();
   UpdateVector& cvec = thesolver->coPvec();
   SSVector& cupd = thesolver->coPvec().delta();
   Vector& ucb = thesolver->ucBound();
   Vector& lcb = thesolver->lcBound();

   if (thesolver->isCoId(id))
   {
      if (thesolver->isCoBasic(nr))
      {
         cupd.clearIdx(nr);
         return true;
      }
      x = cvec[nr];
      d = cupd[nr];
      up = &ucb;
      low = &lcb;
      if (d > 0.0)
         sel = (thesolver->lcBound()[nr] - cvec[nr]) / d;
      else
         sel = (thesolver->ucBound()[nr] - cvec[nr]) / d;
   }

   else if (thesolver->isId(id))
   {
      pvec[nr] = thesolver->vector(nr) * cvec;
      if (thesolver->isBasic(nr))
      {
         pupd.clearIdx(nr);
         return true;
      }
      x = pvec[nr];
      d = pupd[nr];
      up = &upb;
      low = &lpb;
      if (d > 0.0)
         sel = (thesolver->lpBound()[nr] - pvec[nr]) / d;
      else
         sel = (thesolver->upBound()[nr] - pvec[nr]) / d;
   }

   else
      return true;

   if ((*up)[nr] != (*low)[nr])
   {
      if (sel > fastDelta / maxabs)
      {
         if (d < 0.0)
         {
            thesolver->theShift -= (*up)[nr];
            sel = 0.0;
            (*up)[nr] = x + sel * d;
            thesolver->theShift += (*up)[nr];
         }
         else
         {
            thesolver->theShift += (*low)[nr];
            sel = 0.0;
            (*low)[nr] = x + sel * d;
            thesolver->theShift -= (*low)[nr];
         }
      }
   }
   else
   {
      sel = 0.0;
      if (x > (*up)[nr])
         thesolver->theShift += x - (*up)[nr];
      else
         thesolver->theShift += (*low)[nr] - x;
      (*up)[nr] = (*low)[nr] = x;
   }

   return false;
}

bool SPxFastRT::shortEnter(
   const SPxId& enterId,
   int nr,
   Real max,
   Real maxabs) const
{
   if (thesolver->isCoId(enterId))
   {
      if (max != 0.0)
      {
         Real x = thesolver->coPvec().delta()[nr];
         if (x < maxabs * SHORT && -x < maxabs * SHORT)
            return false;
      }
      return true;
   }
   else if (thesolver->isId(enterId))
   {
      if (max != 0.0)
      {
         Real x = thesolver->pVec().delta()[nr];
         if (x < maxabs * SHORT && -x < maxabs * SHORT)
            return false;
      }
      return true;
   }

   return false;
}

SPxId SPxFastRT::selectEnter(Real& val, int)
{
   SPxId enterId;
   Real max, sel;
   Real maxabs = 0.0;
   int nr;
   int cnt = 0;

   assert( m_type == SPxSolver::LEAVE );

   // force instable pivot iff true (see explanation in leave.cpp and spxsolve.cpp)
   bool instable = solver()->instableLeave;
   Real lowstab = LOWSTAB;
   assert(!instable || solver()->instableLeaveNum >= 0);

   resetTols();
   sel = 0.0;

   if (val > epsilon)
   {
      do
      {
         maxabs = 0.0;
         max = val;

         enterId = maxDelta(nr, max, maxabs);
         if (!enterId.isValid())
            return enterId;
         assert(max >= 0.0);

         if (!shortEnter(enterId, nr, max, maxabs))
         {
            Real bestDelta, stab;

            stab = minStability(maxabs);

            // force instable pivot iff instable is true (see explanation in leave.cpp and spxsolve.cpp)
            if (instable)
            {
               enterId = maxSelect(nr, sel, lowstab, bestDelta, max);
            }
            else
            {
               enterId = maxSelect(nr, sel, stab, bestDelta, max);
            }
            if (bestDelta < DELTA_SHIFT*TRIES)
               cnt++;
            else
               cnt += TRIES;
         }
         if (!maxReEnter(sel, maxabs, enterId, nr))
            break;
         relax();
      }
      while (cnt < TRIES);
   }
   else if (val < -epsilon)
   {
      do
      {
         maxabs = 0.0;
         max = val;
         enterId = minDelta(nr, max, maxabs);
         if (!enterId.isValid())
            return enterId;
         assert(max <= 0.0);

         if (!shortEnter(enterId, nr, max, maxabs))
         {
            Real bestDelta, stab;

            stab = minStability(maxabs);

            // force instable pivot iff instable is true (see explanation in leave.cpp and spxsolve.cpp)
            if (instable)
            {
               enterId = minSelect(nr, sel, lowstab, bestDelta, max);
            }
            else
            {
               enterId = minSelect(nr, sel, stab, bestDelta, max);
            }
            if (bestDelta < DELTA_SHIFT*TRIES)
               cnt++;
            else
               cnt += TRIES;
         }
         if (!minReEnter(sel, maxabs, enterId, nr))
            break;
         relax();
      }
      while (cnt < TRIES);
   }

   MSG_DEBUG(
      if (enterId.isValid())
      {
         Real x;
         if (thesolver->isCoId(enterId))
            x = thesolver->coPvec().delta()[ thesolver->number(enterId) ];
         else
            x = thesolver->pVec().delta()[ thesolver->number(enterId) ];
         spxout << "DFSTRT03 " << thesolver->basis().iteration() << ": "
                << sel << '\t' << x << " (" << maxabs << ")" << std::endl;
      }
      else
         spxout << "DFSTRT04 " << thesolver->basis().iteration()
                << ": skipping instable pivot" << std::endl;
   )


   if (enterId.isValid() || minStab > 2*epsilon)
   {
      val = sel;
      if (enterId.isValid())
         tighten();
   }

   return enterId;
}

void SPxFastRT::load(SPxSolver* spx)
{
   thesolver = spx;
   setType(spx->type());
}

void SPxFastRT::setType(SPxSolver::Type type)
{
   m_type = type;

   minStab = MINSTAB;
   fastDelta = delta;
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
