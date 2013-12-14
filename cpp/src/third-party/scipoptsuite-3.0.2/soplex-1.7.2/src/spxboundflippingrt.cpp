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

// #define DEBUGGING

#include <assert.h>
#include "spxdefines.h"
#include "spxboundflippingrt.h"
#include "sorter.h"
#include "spxsolver.h"
#include "spxout.h"
#include "spxid.h"

namespace soplex
{

#define MINSTAB          1e-5
#define LOWSTAB          1e-10
#define MAX_RELAX_COUNT  2
#define LONGSTEP_FREQ    100
#define MIN_LONGSTEP     1e-5


/** perform necessary bound flips to restore dual feasibility */
void SPxBoundFlippingRT::flipAndUpdate(
   int&                  usedBp              /**< number of bounds that should be flipped */
   )
{
   assert(thesolver->rep() == SPxSolver::COLUMN);

   int skipped;

   updPrimRhs.setup();
   updPrimRhs.reDim(thesolver->dim());
   updPrimVec.reDim(thesolver->dim());
   updPrimRhs.clear();
   updPrimVec.clear();

   skipped = 0;
   for( int i = 0; i < usedBp; ++i )
   {
      int idx;
      idx = breakpoints[i].idx;
      if( idx < 0 )
      {
         ++skipped;
         continue;
      }
      Real range;
      Real upper;
      Real lower;
      SPxBasis::Desc::Status stat;
      SPxBasis::Desc& ds = thesolver->basis().desc();
      range = 0;
      if( breakpoints[i].src == PVEC )
      {
         stat = ds.status(idx);
         upper = thesolver->upper(idx);
         lower = thesolver->lower(idx);
         switch( stat )
         {
            case SPxBasis::Desc::P_ON_UPPER :
               ds.status(idx) = SPxBasis::Desc::P_ON_LOWER;
               range = lower - upper;
               assert((*thesolver->theLbound)[idx] == -infinity);
               (*thesolver->theLbound)[idx] = (*thesolver->theUbound)[idx];
               (*thesolver->theUbound)[idx] = infinity;
               break;
            case SPxBasis::Desc::P_ON_LOWER :
               ds.status(idx) = SPxBasis::Desc::P_ON_UPPER;
               range = upper - lower;
               assert((*thesolver->theUbound)[idx] == infinity);
               (*thesolver->theUbound)[idx] = (*thesolver->theLbound)[idx];
               (*thesolver->theLbound)[idx] = -infinity;
               break;
            default :
               ++skipped;
               MSG_WARNING( spxout << "PVEC unexpected status: " << stat
                                   << " index: " << idx
                                   << " val: " << thesolver->pVec()[idx]
                                   << " upd: " << thesolver->pVec().delta()[idx]
                                   << " lower: " << lower
                                   << " upper: " << upper
                                   << " bp.val: " << breakpoints[i].val
                                   << std::endl; )
         }
         MSG_DEBUG( spxout << "PVEC flipped from: " << stat
                           << " index: " << idx
                           << " val: " << thesolver->pVec()[idx]
                           << " upd: " << thesolver->pVec().delta()[idx]
                           << " lower: " << lower
                           << " upper: " << upper
                           << " bp.val: " << breakpoints[i].val
                           << std::endl; )
         assert(fabs(range) < 1e20);
         updPrimRhs.multAdd(range, thesolver->vector(idx));
      }
      else if( breakpoints[i].src == COPVEC )
      {
         stat = ds.coStatus(idx);
         upper = thesolver->rhs(idx);
         lower = thesolver->lhs(idx);
         switch( stat )
         {
            case SPxBasis::Desc::P_ON_UPPER :
               ds.coStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
               range = lower - upper;
               assert((*thesolver->theCoUbound)[idx] == infinity);
               (*thesolver->theCoUbound)[idx] = -(*thesolver->theCoLbound)[idx];
               (*thesolver->theCoLbound)[idx] = -infinity;
               break;
            case SPxBasis::Desc::P_ON_LOWER :
               ds.coStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
               range = upper - lower;
               assert((*thesolver->theCoLbound)[idx] == -infinity);
               (*thesolver->theCoLbound)[idx] = -(*thesolver->theCoUbound)[idx];
               (*thesolver->theCoUbound)[idx] = infinity;
               break;
            default :
               ++skipped;
               MSG_WARNING( spxout << "COPVEC unexpected status: " << stat
                                   << " index: " << idx
                                   << " val: " << thesolver->coPvec()[idx]
                                   << " upd: " << thesolver->coPvec().delta()[idx]
                                   << " lower: " << lower
                                   << " upper: " << upper
                                   << " bp.val: " << breakpoints[i].val
                                   << std::endl; )
         }
         MSG_DEBUG( spxout << "COPVEC flipped from: " << stat
                           << " index: " << idx
                           << " val: " << thesolver->coPvec()[idx]
                           << " upd: " << thesolver->coPvec().delta()[idx]
                           << " lower: " << lower
                           << " upper: " << upper
                           << " bp.val: " << breakpoints[i].val
                           << std::endl; )
         assert(fabs(range) < 1e20);
         updPrimRhs.setValue(idx, updPrimRhs[idx] - range);
      }
   }
   usedBp -= skipped;
   if( usedBp > 0 )
   {
      thesolver->primRhs -= updPrimRhs;
      thesolver->setup4solve2(&updPrimVec, &updPrimRhs);
   }

   return;
}

/** store all available pivots/breakpoints in an array (positive pivot search direction) */
void SPxBoundFlippingRT::collectBreakpointsMax(
   int&                  nBp,                /**< number of found breakpoints so far */
   int&                  minIdx,             /**< index to current minimal breakpoint */
   const int*            idx,                /**< pointer to indices of current vector */
   int                   nnz,                /**< number of nonzeros in current vector */
   const Real*           upd,                /**< pointer to update values of current vector */
   const Real*           vec,                /**< pointer to values of current vector */
   const Real*           upp,                /**< pointer to upper bound/rhs of current vector */
   const Real*           low,                /**< pointer to lower bound/lhs of current vector */
   BreakpointSource      src                 /**< type of vector (pVec or coPvec)*/
   )
{
   Real minVal;
   Real curVal;
   const int* last;

   minVal = ( nBp == 0 ) ? infinity : breakpoints[minIdx].val;

   last = idx + nnz;

   for( ; idx < last; ++idx )
   {
      int i = *idx;
      Real x = upd[i];
      if( x > epsilon )
      {
         if( upp[i] < infinity )
         {
            Real y = upp[i] - vec[i];
            curVal = (y <= 0) ? fastDelta / x : (y + fastDelta) / x;
            assert(curVal > 0);

            breakpoints[nBp].idx = i;
            breakpoints[nBp].src = src;
            breakpoints[nBp].val = curVal;

            if( curVal < minVal )
            {
               minVal = curVal;
               minIdx = nBp;
            }

            nBp++;
         }
      }
      else if( x < -epsilon )
      {
         if (low[i] > -infinity)
         {
            Real y = low[i] - vec[i];
            curVal = (y >= 0) ? -fastDelta / x : (y - fastDelta) / x;
            assert(curVal > 0);

            breakpoints[nBp].idx = i;
            breakpoints[nBp].src = src;
            breakpoints[nBp].val = curVal;

            if( curVal < minVal )
            {
               minVal = curVal;
               minIdx = nBp;
            }

            nBp++;
         }
      }
      if( nBp >= breakpoints.size() )
         breakpoints.reSize(nBp * 2);
   }

   return;
}

/** store all available pivots/breakpoints in an array (negative pivot search direction) */
void SPxBoundFlippingRT::collectBreakpointsMin(
   int&                  nBp,                /**< number of found breakpoints so far */
   int&                  minIdx,             /**< index to current minimal breakpoint */
   const int*            idx,                /**< pointer to indices of current vector */
   int                   nnz,                /**< number of nonzeros in current vector */
   const Real*           upd,                /**< pointer to update values of current vector */
   const Real*           vec,                /**< pointer to values of current vector */
   const Real*           upp,                /**< pointer to upper bound/rhs of current vector */
   const Real*           low,                /**< pointer to lower bound/lhs of current vector */
   BreakpointSource      src                 /**< type of vector (pVec or coPvec)*/
   )
{
   Real minVal;
   Real curVal;
   const int* last;

   minVal = ( nBp == 0 ) ? infinity : breakpoints[minIdx].val;

   last = idx + nnz;

   for( ; idx < last; ++idx )
   {
      int i = *idx;
      Real x = upd[i];
      if( x > epsilon )
      {
         if( low[i] > -infinity )
         {
            Real y = low[i] - vec[i];

            curVal = (y >= 0) ? fastDelta / x : (fastDelta - y) / x;
            assert(curVal > 0);

            breakpoints[nBp].idx = i;
            breakpoints[nBp].src = src;
            breakpoints[nBp].val = curVal;

            if( curVal < minVal )
            {
               minVal = curVal;
               minIdx = nBp;
            }

            nBp++;
         }
      }
      else if( x < -epsilon )
      {
         if (upp[i] < infinity)
         {
            Real y = upp[i] - vec[i];
            curVal = (y <= 0) ? -fastDelta / x : -(y + fastDelta) / x;
            assert(curVal > 0);

            breakpoints[nBp].idx = i;
            breakpoints[nBp].src = src;
            breakpoints[nBp].val = curVal;

            if( curVal < minVal )
            {
               minVal = curVal;
               minIdx = nBp;
            }

            nBp++;
         }
      }
      if( nBp >= breakpoints.size() )
         breakpoints.reSize(nBp * 2);
   }
   return;
}

/** get values for entering index and perform shifts if necessary */
bool SPxBoundFlippingRT::getData(
   Real&                 val,
   SPxId&                enterId,
   int                   idx,
   Real                  stab,
   Real                  degeneps,
   const Real*           upd,
   const Real*           vec,
   const Real*           low,
   const Real*           upp,
   BreakpointSource      src,
   Real                  max
   )
{
   if( src == PVEC )
   {
      thesolver->pVec()[idx] = thesolver->vector(idx) * thesolver->coPvec();
      Real x = upd[idx];
      // skip breakpoint if it is too small
      if( fabs(x) < stab )
      {
         return false;
      }
      enterId = thesolver->id(idx);
      val = (max * x > 0) ? upp[idx] : low[idx];
      val = (val - vec[idx]) / x;
      if( upp[idx] == low[idx] )
      {
         val = 0.0;
         if( vec[idx] > upp[idx] )
            thesolver->theShift += vec[idx] - upp[idx];
         else
            thesolver->theShift += low[idx] - vec[idx];
         thesolver->upBound()[idx] = thesolver->lpBound()[idx] = vec[idx];
      }
      else if( (max > 0 && val < -degeneps) || (max < 0 && val > degeneps) )
      {
         val = 0.0;
         if( max * x > 0 )
            thesolver->shiftUPbound(idx, vec[idx]);
         else
            thesolver->shiftLPbound(idx, vec[idx]);
      }
   }
   else // breakpoints[usedBp].src == COPVEC
   {
      Real x = upd[idx];
      if( fabs(x) < stab )
      {
         return false;
      }
      enterId = thesolver->coId(idx);
      val = (max * x > 0.0) ? upp[idx] : low[idx];
      val = (val - vec[idx]) / x;
      if( upp[idx] == low[idx] )
      {
         val = 0.0;
         if( vec[idx] > upp[idx] )
            thesolver->theShift += vec[idx] - upp[idx];
         else
            thesolver->theShift += low[idx] - vec[idx];
         thesolver->ucBound()[idx] = thesolver->lcBound()[idx] = vec[idx];
      }
      else if( (max > 0 && val < -degeneps) || (max < 0 && val > degeneps) )
      {
         val = 0.0;
         if( max * x > 0 )
            thesolver->shiftUCbound(idx, vec[idx]);
         else
            thesolver->shiftLCbound(idx, vec[idx]);
      }
   }
   return true;
}

/** determine entering variable */
SPxId SPxBoundFlippingRT::selectEnter(
   Real&                 val,
   int                   leaveIdx
   )
{
   assert( m_type == SPxSolver::LEAVE );
   assert(thesolver->boundflips == 0);

   // reset the history and try again to do some long steps
   if( thesolver->leaveCount % LONGSTEP_FREQ == 0 )
   {
      MSG_DEBUG( spxout << "ILSTEP06 resetting long step history" << std::endl; )
      flipPotential = 1;
   }
   if( !enableLongsteps || thesolver->rep() == SPxSolver::ROW || flipPotential < 0.01 )
   {
      MSG_DEBUG( spxout << "ILSTEP07 switching to fast ratio test" << std::endl; )
      return SPxFastRT::selectEnter(val, leaveIdx);
   }
   const Real*  pvec = thesolver->pVec().get_const_ptr();
   const Real*  pupd = thesolver->pVec().delta().values();
   const int*   pidx = thesolver->pVec().delta().indexMem();
   int          pupdnnz = thesolver->pVec().delta().size();
   const Real*  lpb  = thesolver->lpBound().get_const_ptr();
   const Real*  upb  = thesolver->upBound().get_const_ptr();

   const Real*  cvec = thesolver->coPvec().get_const_ptr();
   const Real*  cupd = thesolver->coPvec().delta().values();
   const int*   cidx = thesolver->coPvec().delta().indexMem();
   int          cupdnnz = thesolver->coPvec().delta().size();
   const Real*  lcb  = thesolver->lcBound().get_const_ptr();
   const Real*  ucb  = thesolver->ucBound().get_const_ptr();

   resetTols();

   Real max;

   // index in breakpoint array of minimal value (i.e. choice of normal RT)
   int minIdx;

   // temporary breakpoint data structure to make swaps possible
   Breakpoint tmp;

   // most stable pivot value in candidate set
   Real moststable;

   // initialize invalid enterId
   SPxId enterId;

   // slope of objective function improvement
   Real slope;

   // number of found breakpoints
   int nBp;

   // number of latest skipable breakpoint
   int usedBp;

   Real degeneps;
   Real stab;
   bool instable;

   max = val;
   val = 0.0;
   moststable = 0.0;
   nBp = 0;
   minIdx = -1;

   // get breakpoints and and determine the index of the minimal value
   if( max > 0 )
   {
      collectBreakpointsMax(nBp, minIdx, pidx, pupdnnz, pupd, pvec, upb, lpb, PVEC);
      collectBreakpointsMax(nBp, minIdx, cidx, cupdnnz, cupd, cvec, ucb, lcb, COPVEC);
   }
   else
   {
      collectBreakpointsMin(nBp, minIdx, pidx, pupdnnz, pupd, pvec, upb, lpb, PVEC);
      collectBreakpointsMin(nBp, minIdx, cidx, cupdnnz, cupd, cvec, ucb, lcb, COPVEC);
   }

   if( nBp == 0 )
      return enterId;

   assert(minIdx >= 0);

   // swap smallest breakpoint to the front to skip the sorting phase if no bound flip is possible
   tmp = breakpoints[minIdx];
   breakpoints[minIdx] = breakpoints[0];
   breakpoints[0] = tmp;

   // compute initial slope
   slope = fabs(thesolver->fTest()[leaveIdx]);
   if( slope == 0 )
   {
      // this may only happen if SoPlex decides to make an instable pivot
      assert(thesolver->instableLeaveNum >= 0);
      // restore original slope
      slope = fabs(thesolver->instableLeaveVal);
   }

   // set up structures for the quicksort implementation
   BreakpointCompare compare;
   compare.entry = breakpoints.get_const_ptr();

   // pointer to end of sorted part of breakpoints
   int sorted = 0;
   // minimum number of entries that are supposed to be sorted by partial sort
   int sortsize = 4;

   // get all skipable breakpoints
   for( usedBp = 0; usedBp < nBp && slope > 0; ++usedBp)
   {
      // sort breakpoints only partially to save time
      if( usedBp > sorted )
      {
         sorted = sorter_qsortPart(breakpoints.get_ptr(), compare, sorted + 1, nBp, sortsize);
      }
      int i = breakpoints[usedBp].idx;
      // compute new slope
      if( breakpoints[usedBp].src == PVEC )
      {
         if( thesolver->isBasic(i) )
         {
            // mark basic indices
            breakpoints[usedBp].idx = -1;
            thesolver->pVec().delta().clearIdx(i);
         }
         else
         {
            Real absupd = fabs(pupd[i]);
            slope -= (thesolver->upper(i) * absupd) - (thesolver->lower(i) * absupd);
            // get most stable pivot
            if( absupd > moststable )
               moststable = absupd;
         }
      }
      else
      {
         assert(breakpoints[usedBp].src == COPVEC);
         if( thesolver->isCoBasic(i) )
         {
            // mark basic indices
            breakpoints[usedBp].idx = -1;
            thesolver->coPvec().delta().clearIdx(i);
         }
         else
         {
            Real absupd = fabs(cupd[i]);
            slope -= (thesolver->rhs(i) * absupd) - (thesolver->lhs(i) * absupd);
            if( absupd > moststable )
               moststable = absupd;
         }
      }
   }
   --usedBp;
   assert(usedBp >= 0);

   // check for unboundedness/infeasibility
   if( slope > delta && usedBp >= nBp - 1 )
   {
      MSG_DEBUG( spxout << "DLSTEP02 " << thesolver->basis().iteration()
                        << ": unboundedness in ratio test" << std::endl; )
      flipPotential *= 0.95;
      val = max;
      return SPxFastRT::selectEnter(val, leaveIdx);
   }

   MSG_DEBUG( spxout << "DLSTEP01 "
                     << thesolver->basis().iteration()
                     << ": number of flip candidates: "
                     << usedBp
                     << std::endl; )

   // try to get a more stable pivot by looking at those with similar step length
   int stableBp;              // index to walk over additional breakpoints (after slope change)
   int bestBp = -1;           // breakpoints index with best possible stability
   Real bestDelta = breakpoints[usedBp].val;  // best step length (after bound flips)
   for( stableBp = usedBp + 1; stableBp < nBp; ++stableBp )
   {
      Real stableDelta;
      // get next breakpoints in increasing order
      if( stableBp > sorted )
      {
         sorted = sorter_qsortPart(breakpoints.get_ptr(), compare, sorted + 1, nBp, sortsize);
      }
      int idx = breakpoints[stableBp].idx;
      if( breakpoints[stableBp].src == PVEC )
      {
         if( thesolver->isBasic(idx) )
         {
            // mark basic indices
            // TODO this is probably unnecessary since we do not look at these breakpoints anymore
            breakpoints[stableBp].idx = -1;
            thesolver->pVec().delta().clearIdx(idx);
            continue;
         }
         thesolver->pVec()[idx] = thesolver->vector(idx) * thesolver->coPvec();
         Real x = pupd[idx];
         stableDelta = (x > 0.0) ? upb[idx] : lpb[idx];
         stableDelta = (stableDelta - pvec[idx]) / x;

         if( stableDelta <= bestDelta)
         {
            if( fabs(x) > moststable )
            {
               moststable = fabs(x);
               bestBp = stableBp;
            }
         }
         // stop searching if the step length is too big
         else if( stableDelta > 2 * bestDelta )
            break;
      }
      else
      {
         if( thesolver->isCoBasic(idx) )
         {
            // mark basic indices
            breakpoints[usedBp].idx = -1;
            thesolver->coPvec().delta().clearIdx(idx);
            continue;
         }
         Real x = cupd[idx];
         stableDelta = (x > 0.0) ? ucb[idx] : lcb[idx];
         stableDelta = (stableDelta - cvec[idx]) / x;

         if( stableDelta <= bestDelta)
         {
            if( fabs(x) > moststable )
            {
               moststable = fabs(x);
               bestBp = stableBp;
            }
         }
         // TODO check whether 2 is a good factor
         // stop searching if the step length is too big
         else if( stableDelta > 2 * bestDelta )
            break;
      }
   }

   degeneps = fastDelta / moststable;  /* as in SPxFastRT */
   // get stability requirements
   instable = thesolver->instableLeave;
   assert(!instable || thesolver->instableLeaveNum >= 0);
   stab = instable ? LOWSTAB : SPxFastRT::minStability(moststable);

   bool foundStable= false;

   if( bestBp > -1 )
   {
      // we found a more stable pivot
      if( moststable > stab )
      {
         // stability requirements are satisfied
         int idx = breakpoints[bestBp].idx;
         if( breakpoints[bestBp].src == PVEC )
            foundStable = getData(val, enterId, idx, stab, degeneps, pupd, pvec, lpb, upb, PVEC, max);
         else
            foundStable = getData(val, enterId, idx, stab, degeneps, cupd, cvec, lcb, ucb, COPVEC, max);
      }
   }

#if 0
   // do not make long steps if the gain in the dual objective is too small, except to avoid degenerate steps
   if( usedBp > 0 && breakpoints[0].val > epsilon && breakpoints[usedBp].val - breakpoints[0].val < MIN_LONGSTEP * breakpoints[0].val )
   {
      MSG_DEBUG( spxout << "DLSTEP03 "
                        << thesolver->basis().iteration()
                        << ": bound flip gain is too small"
                        << std::endl; )

      // ensure that the first breakpoint is nonbasic
      while( breakpoints[usedBp].idx < 0 && usedBp < nBp )
         ++usedBp;
      // @todo make sure that the selected pivot element is stable
   }
#endif

   // scan pivot candidates from back to front and stop as soon as a good one is found
   // @todo select the moststable pivot or one that is comparably stable
//    stab = (moststable > stab) ? moststable : stab;

   while( !foundStable && usedBp >= 0 )
   {
      int idx = breakpoints[usedBp].idx;

      // only look for non-basic variables
      if( idx >= 0 )
      {
         if( breakpoints[usedBp].src == PVEC )
            foundStable = getData(val, enterId, idx, stab, degeneps, pupd, pvec, lpb, upb, PVEC, max);
         else
            foundStable = getData(val, enterId, idx, stab, degeneps, cupd, cvec, lcb, ucb, COPVEC, max);
      }
      --usedBp;
   }

   if( !foundStable )
   {
//       assert(usedBp < 0);
      assert(!enterId.isValid());
      if( relax_count < MAX_RELAX_COUNT )
      {
         MSG_DEBUG( spxout << "DLSTEP04 "
                           << thesolver->basis().iteration()
                           << ": no valid enterId found - relaxing..."
                           << std::endl; )
         relax();
         ++relax_count;
         // restore original value
         val = max;
         // try again with relaxed delta
         return SPxBoundFlippingRT::selectEnter(val, leaveIdx);
      }
      else
      {
         MSG_DEBUG( spxout << "DLSTEP05 "
                           << thesolver->basis().iteration()
                           << " no valid enterId found - breaking..."
                           << std::endl; )
         assert(!enterId.isValid());
         return enterId;
      }
   }
   else
   {
//       assert(usedBp >= 0);
      relax_count = 0;
      tighten();
   }

   // flip bounds of skipped breakpoints only if a nondegenerate step is to be performed
   if( usedBp > 0 && fabs(breakpoints[usedBp].val) > fastDelta )
   {
      flipAndUpdate(usedBp);
      thesolver->boundflips = usedBp;
   }
   else
      thesolver->boundflips = 0;

   //TODO incorporate the ratio between usedBp, nBp and dim/coDim
   //     to get an idea of effort and speed

   // estimate wether long steps may be possible in future iterations
   flipPotential *= (usedBp + 0.95);

   MSG_DEBUG( spxout << "DLSTEP06 "
                     << thesolver->basis().iteration()
                     << ": selected Id: "
                     << enterId
                     << " number of candidates: "
                     << nBp
                     << std::endl; )
   return enterId;
}

int SPxBoundFlippingRT::selectLeave(
   Real&                 max,
   SPxId                 enterId
   )
{
   return SPxFastRT::selectLeave(max, enterId);
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
