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

#include "spxdefines.h"
#include "spxdevexpr.h"

#define DEVEX_REFINETOL 2.0

namespace soplex
{

void SPxDevexPR::load(SPxSolver* base)
{
   thesolver = base;
   setRep(base->rep());
   assert(isConsistent());
}

bool SPxDevexPR::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (thesolver != 0)
      if (penalty.dim() != thesolver->coDim()
           || coPenalty.dim() != thesolver->dim())
         return MSGinconsistent("SPxDevexPR");
#endif

   return true;
}

void SPxDevexPR::init(SPxSolver::Type tp)
{
   int i;
   if (tp == SPxSolver::ENTER)
   {
      for (i = penalty.dim(); --i >= 0;)
         penalty[i] = 2;
      for (i = coPenalty.dim(); --i >= 0;)
         coPenalty[i] = 2;
   }
   else
   {
      for (i = coPenalty.dim(); --i >= 0;)
         coPenalty[i] = 1;
   }
   assert(isConsistent());
}

void SPxDevexPR::setType(SPxSolver::Type tp)
{
   init(tp);
   refined = false;
}

/**@todo suspicious: Shouldn't the relation between dim, coDim, Vecs, 
 *       and CoVecs be influenced by the representation ?
 */
void SPxDevexPR::setRep(SPxSolver::Representation)
{
   if (thesolver != 0)
   {
      addedVecs(thesolver->coDim());
      addedCoVecs(thesolver->dim());
      assert(isConsistent());
   }
}

int SPxDevexPR::selectLeave()
{
   int retid;
   Real val;

#ifdef PARTIAL_PRICING
   retid = selectLeavePart(val, theeps);
#else
   if (thesolver->sparsePricingLeave)
      retid = selectLeaveSparse(val, theeps);
   else
      retid = selectLeaveX(val, theeps);
#endif
   if ( retid < 0 && !refined )
   {
      refined = true;
      MSG_INFO3( spxout << "WDEVEX02 trying refinement step..\n"; )
      retid = selectLeaveX(val, theeps/DEVEX_REFINETOL);
   }

   return retid;
}

int SPxDevexPR::selectLeaveX(Real& best, Real feastol, int start, int incr)
{
   Real x;

   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* cpen = coPenalty.get_const_ptr();
   Real bstX = 0;
   int bstI = -1;
   int end = coPenalty.dim();

   for (; start < end; start += incr)
   {
      if (fTest[start] < -feastol)
      {
         x = fTest[start] * fTest[start] / cpen[start];
         if (x > bstX)
         {
            bstX = x;
            bstI = start;
            last = cpen[start];
         }
      }
   }
   best = bstX;
   return bstI;
}

int SPxDevexPR::selectLeavePart(Real& best, Real feastol)
{
   Real x;

   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* cpen = coPenalty.get_const_ptr();
   Real bstX = 0;
   int bstI = -1;
   int dim = coPenalty.dim();
   int count = 0;
   int oldstartpricing = startpricing;
   int end = oldstartpricing + IMPROVEMENT_STEPLENGTH;

   for (int i = oldstartpricing; i < dim; ++i)
   {
      Real fTesti = fTest[i];
      if (fTesti < -feastol)
      {
         Real cpeni = cpen[i];
         x = fTesti * fTesti / cpeni;
         if (x > bstX * IMPROVEMENT_THRESHOLD)
         {
            bstX = x;
            bstI = i;
            last = cpeni;
            end = i + IMPROVEMENT_STEPLENGTH;
//             if (count == 0)
               startpricing = (i + 1) % dim;
            ++count;
         }
      }
      if (i >= end && count >= MAX_PRICING_CANDIDATES)
         break;
   }

   if (end < dim && count >= MAX_PRICING_CANDIDATES)
   {
      best = bstX;
      return bstI;
   }
   else
   {
      end -= dim;
   }

   for (int i = 0; i < oldstartpricing; ++i)
   {
      Real fTesti = fTest[i];
      if (fTesti < -feastol)
      {
         Real cpeni = cpen[i];
         x = fTesti * fTesti / cpeni;
         if (x > bstX * IMPROVEMENT_THRESHOLD)
         {
            bstX = x;
            bstI = i;
            last = cpeni;
            end = i + IMPROVEMENT_STEPLENGTH;
//             if (count == 0)
               startpricing = (i + 1) % dim;
            ++count;
         }
      }
      if (i >= end && count >= MAX_PRICING_CANDIDATES)
         break;
   }

   best = bstX;
   return bstI;
}

int SPxDevexPR::selectLeaveSparse(Real& best, Real feastol)
{
   Real x;

   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* cpen = coPenalty.get_const_ptr();
   Real bstX = 0;
   int bstI = -1;
   int idx = -1;
   Real fTesti;
   Real coPeni;

   for (int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      fTesti = fTest[idx];
      if (fTesti < -feastol)
      {
         coPeni = cpen[idx];
         x = fTesti * fTesti / coPeni;
         if (x > bstX)
         {
            bstX = x;
            bstI = idx;
            last = coPeni;
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);

         assert(thesolver->isInfeasible[idx]);
         thesolver->isInfeasible[idx] = false;
      }
   }
   best = bstX;
   return bstI;
}


void SPxDevexPR::left4(int n, SPxId id)
{
   left4X(n, id, 0, 1);
}

void SPxDevexPR::left4X(int n, const SPxId& id, int start, int incr)
{
   if (id.isValid())
   {
      int i, j;
      Real x;
      const Real* rhoVec = thesolver->fVec().delta().values();
      Real rhov_1 = 1 / rhoVec[n];
      Real beta_q = thesolver->coPvec().delta().length2() * rhov_1 * rhov_1;

#ifndef NDEBUG
      if (fabs(rhoVec[n]) < theeps)
      {
         MSG_ERROR( spxout << "WDEVEX01: rhoVec = "
                           << rhoVec[n] << " with smaller absolute value than theeps = " << theeps << std::endl; )
      }
#endif  // NDEBUG

      //  Update #coPenalty# vector
      const IdxSet& rhoIdx = thesolver->fVec().idx();
      int len = thesolver->fVec().idx().size();
      for (i = len - 1 - start; i >= 0; i -= incr)
      {
         j = rhoIdx.index(i);
         x = rhoVec[j] * rhoVec[j] * beta_q;
         // if(x > coPenalty[j])
         coPenalty[j] += x;
      }

      coPenalty[n] = beta_q;
   }
}



SPxId SPxDevexPR::selectEnter()
{
   assert(thesolver != 0);

   SPxId enterId;

   enterId = selectEnterX(theeps);

   if( !enterId.isValid() && !refined )
   {
      refined = true;
      MSG_INFO3( spxout << "WDEVEX02 trying refinement step..\n"; )
      enterId = selectEnterX(theeps/DEVEX_REFINETOL);
   }

   return enterId;
}

// choose the best entering index among columns and rows but prefer sparsity
SPxId SPxDevexPR::selectEnterX(Real tol)
{
   SPxId enterId;
   SPxId enterIdCo;
   Real best;
   Real bestCo;

   best = 0;
   bestCo = 0;
   enterId = (thesolver->sparsePricingEnter) ? selectEnterSparseDim(best, tol) : selectEnterDenseDim(best, tol);
   enterIdCo = (thesolver->sparsePricingEnterCo) ? selectEnterSparseCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);

   // prefer slack indices to reduce nonzeros in basis matrix
   if( enterId.isValid() && (best > SPARSITY_TRADEOFF * bestCo || !enterIdCo.isValid()) )
      return enterId;
   else
      return enterIdCo;
}

SPxId SPxDevexPR::selectEnterSparseDim(Real& best, Real feastol)
{
   const Real* cTest = thesolver->coTest().get_const_ptr();
   const Real* cpen = coPenalty.get_const_ptr();
   int end = coPenalty.dim();
   int enterIdx = -1;
   int idx;
   Real coTesti;
   Real coPeni;
   Real x;

   assert(end == thesolver->coTest().dim());
   for(int i = thesolver->infeasibilities.size() -1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      coTesti = cTest[idx];
      if (coTesti < -feastol)
      {
         coPeni = cpen[idx];
         x = coTesti * coTesti / coPeni;
         if (x > best)
         {
            best = x;
            enterIdx = idx;
            last = cpen[idx];
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         assert(thesolver->isInfeasible[idx]);
         thesolver->isInfeasible[idx] = false;
      }
   }
   if (enterIdx >= 0)
      return thesolver->coId(enterIdx);

   SPxId none;
   return none;
}


SPxId SPxDevexPR::selectEnterSparseCoDim(Real& best, Real feastol)
{
   const Real* test = thesolver->test().get_const_ptr();
   const Real* pen = penalty.get_const_ptr();
   int end = penalty.dim();
   int enterIdx = -1;
   int idx;
   Real testi;
   Real peni;
   Real x;

   assert(end == thesolver->test().dim());
   for (int i = thesolver->infeasibilitiesCo.size() -1; i >= 0; --i)
   {
      idx = thesolver->infeasibilitiesCo.index(i);
      testi = test[idx];
      if (testi < -feastol)
      {
         peni = pen[idx];
         x = testi * testi / peni;
         if (x > best)
         {
            best = x;
            enterIdx = idx;
            last = pen[idx];
         }
      }
      else
      {
         thesolver->infeasibilitiesCo.remove(i);
         assert(thesolver->isInfeasibleCo[idx]);
         thesolver->isInfeasibleCo[idx] = false;
      }
   }

   if (enterIdx >= 0)
      return thesolver->id(enterIdx);

   SPxId none;
   return none;
}


SPxId SPxDevexPR::selectEnterDenseDim(Real& best, Real feastol, int start, int incr)
{
   const Real* cTest = thesolver->coTest().get_const_ptr();
   const Real* cpen = coPenalty.get_const_ptr();
   int end = coPenalty.dim();
   int enterIdx = -1;
   Real x;

   assert(end == thesolver->coTest().dim());
   for (; start < end; start += incr)
   {
      if (cTest[start] < -feastol)
      {
         x = cTest[start] * cTest[start] / cpen[start];
         if (x > best)
         {
            best = x;
            enterIdx = start;
            last = cpen[start];
         }
      }
   }

   if (enterIdx >= 0)
      return thesolver->coId(enterIdx);

   SPxId none;
   return none;
}


SPxId SPxDevexPR::selectEnterDenseCoDim(Real& best, Real feastol, int start, int incr)
{
   const Real* test = thesolver->test().get_const_ptr();
   const Real* pen = penalty.get_const_ptr();
   int end = penalty.dim();
   int enterIdx = -1;
   Real x;

   assert(end == thesolver->test().dim());
   for (; start < end; start += incr)
   {
      if (test[start] < -feastol)
      {
         x = test[start] * test[start] / pen[start];
         if (x > best)
         {
            best = x;
            enterIdx = start;
            last = pen[start];
         }
      }
   }

   if (enterIdx >= 0)
      return thesolver->id(enterIdx);

   SPxId none;
   return none;
}


void SPxDevexPR::entered4(SPxId id, int n)
{
   entered4X(id, n, 0, 1, 0, 1);
}

/**@todo suspicious: the pricer should be informed, that variable id 
    has entered the basis at position n, but the id is not used here 
    (this is true for all pricers)
*/
void SPxDevexPR::entered4X(SPxId /*id*/, int n,
   int start1, int incr1, int start2, int incr2)
{
   if (n >= 0 && n < thesolver->dim())
   {
      const Real* pVec = thesolver->pVec().delta().values();
      const IdxSet& pIdx = thesolver->pVec().idx();
      const Real* coPvec = thesolver->coPvec().delta().values();
      const IdxSet& coPidx = thesolver->coPvec().idx();
      Real xi_p = 1 / thesolver->fVec().delta()[n];
      int i, j;

      assert(thesolver->fVec().delta()[n] > thesolver->epsilon()
              || thesolver->fVec().delta()[n] < -thesolver->epsilon());

      xi_p = xi_p * xi_p * last;

      for (j = coPidx.size() - 1 - start1; j >= 0; j -= incr1)
      {
         i = coPidx.index(j);
         coPenalty[i] += xi_p * coPvec[i] * coPvec[i];
         if (coPenalty[i] <= 1 || coPenalty[i] > 1e+6)
         {
            init(SPxSolver::ENTER);
            return;
         }
      }

      for (j = pIdx.size() - 1 - start2; j >= 0; j -= incr2)
      {
         i = pIdx.index(j);
         penalty[i] += xi_p * pVec[i] * pVec[i];
         if (penalty[i] <= 1 || penalty[i] > 1e+6)
         {
            init(SPxSolver::ENTER);
            return;
         }
      }
   }
}

void SPxDevexPR::addedVecs (int n)
{
   int initval = (thesolver->type() == SPxSolver::ENTER) ? 2 : 1;
   n = penalty.dim();
   penalty.reDim (thesolver->coDim());
   for (int i = penalty.dim()-1; i >= n; --i )
      penalty[i] = initval;
}

void SPxDevexPR::addedCoVecs(int n)
{
   int initval = (thesolver->type() == SPxSolver::ENTER) ? 2 : 1;
   n = coPenalty.dim();
   coPenalty.reDim(thesolver->dim());
   for (int i = coPenalty.dim()-1; i >= n; --i)
      coPenalty[i] = initval;
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
