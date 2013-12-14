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

// #define DEBUGGING 1

/*      \SubSection{Updating the Basis for Entering Variables}
 */
#include <assert.h>

#include "spxdefines.h"
#include "soplex.h"
#include "spxratiotester.h"
#include "spxout.h"
#include "exceptions.h"

namespace soplex
{

/*
In the entering simplex algorithms (i.e. iteratively a vector is selected to
\em enter the simplex basis as in the dual rowwise and primal columnwise case)
let $A$ denote the current basis, $x$ and entering vector and $f$ the
feasibility vector. For a feasible basis $l \le f \le u$ holds.  
For the rowwise case $f$ is obtained by solving $f^T = c^T A^{-1}$, 
wherease in columnwisecase $f = A^{-1} b$.
 
Let us further consider the rowwise case. Exchanging $x$ with the $i$-th
vector of $A$ yields

\begin{equation}\label{update.eq}
    A^{(i)} = E_i A \hbox{, with } E_i = I + e_i (x^T A^{-1} - e_i^T).
\end{equation}

With $E_i^{-1} = I + e_i \frac{e_i^T - \delta^T}{\delta}$, 
$\delta^T = x^T A^{-1}$ one gets the new feasibility vector

\begin{eqnarray*}
        (f^{(i)})^T
    &=& c^T (A^{(i)})^{-1}      \\
    &=& c^T A^{-1} + c^T A^{-1} e_i \frac{e_i^T - \delta^T}{\delta_i} \\
    &=& f^T + \frac{f_i}{\delta_i} e_i^T - \frac{f_i}{\delta_i} \delta^T. \\
\end{eqnarray*}

The selection of the leaving vector $i^*$ for the basis must ensure, that for
all $j \ne i^*$ $f^{(i^*)}_j$ remains within its bounds $l_j$ and $u_j$.
 */


/*
    Testing all values of |pVec| against its bounds. If $i$, say, is violated
    the violation is saved as negative value in |theTest[i]|.
 */
Real SPxSolver::test(int i, SPxBasis::Desc::Status stat) const
{
   METHOD( "SPxSolver::test()" );
   assert(type() == ENTER);
   assert(!isBasic(stat));

   Real x;

   switch (stat)
   {
   case SPxBasis::Desc::D_FREE:
   case SPxBasis::Desc::D_ON_BOTH:
      assert(rep() == ROW);
      x = (*thePvec)[i] - lhs(i);
      if (x < 0)
         return x;
      // no break: next is else case
      //lint -fallthrough
   case SPxBasis::Desc::D_ON_LOWER:
      assert(rep() == ROW);
      return rhs(i) - (*thePvec)[i];
   case SPxBasis::Desc::D_ON_UPPER:
      assert(rep() == ROW);
      return (*thePvec)[i] - lhs(i);

   case SPxBasis::Desc::P_ON_UPPER:
      assert(rep() == COLUMN);
      return maxObj(i) - (*thePvec)[i];
   case SPxBasis::Desc::P_ON_LOWER:
      assert(rep() == COLUMN);
      return (*thePvec)[i] - maxObj(i);
   case SPxBasis::Desc::P_FREE :
      x = maxObj(i) - (*thePvec)[i];
      return (x < 0) ? x : -x;

   default:
      return 0;
   }
}

void SPxSolver::computeTest()
{
   METHOD( "SPxSolver::computeTest()" );

   const SPxBasis::Desc& ds = desc();
   Real pricingTol = leavetol();
   infeasibilitiesCo.clear();
   int ninfeasibilities = 0;

   for(int i = 0; i < coDim(); ++i)
   {
      SPxBasis::Desc::Status stat = ds.status(i);

      if(isBasic(stat))
         theTest[i] = 0.0;
      else
      {
         theTest[i] = test(i, stat);

         if( remainingRoundsEnterCo == 0 )
         {
            if( theTest[i] < -pricingTol )
            {
               assert(infeasibilitiesCo.size() < infeasibilitiesCo.max());
               infeasibilitiesCo.addIdx(i);
               isInfeasibleCo[i] = true;
               ++ninfeasibilities;
            }
            else
               isInfeasibleCo[i] = false;
            if( ninfeasibilities > sparsityThresholdEnterCo )
            {
               MSG_INFO2( spxout << "IENTER04 too many infeasibilities for sparse pricing"
                                 << std::endl; )
               remainingRoundsEnterCo = DENSEROUNDS;
               sparsePricingEnterCo = false;
               ninfeasibilities = 0;
            }
         }
      }
   }
   if( ninfeasibilities == 0 && !sparsePricingEnterCo )
      --remainingRoundsEnterCo;
   else if( ninfeasibilities <= sparsityThresholdEnterCo && !sparsePricingEnterCo )
   {
      std::streamsize prec = spxout.precision();
      MSG_INFO2( spxout << "IENTER03 sparse pricing active, "
                        << "sparsity: "
                        << std::setw(6) << std::fixed << std::setprecision(4)
                        << (Real) ninfeasibilities/coDim()
                        << std::scientific << std::setprecision(int(prec))
                        << std::endl; )
      sparsePricingEnterCo = true;
   }
}

Real SPxSolver::computePvec(int i)
{
   METHOD( "SPxSolver::computePvec()" );

   return (*thePvec)[i] = vector(i) * (*theCoPvec);
}

Real SPxSolver::computeTest(int i)
{
   METHOD( "SPxSolver::computeTest()" );
   SPxBasis::Desc::Status stat = desc().status(i);
   if (isBasic(stat))
      return theTest[i] = 0;
   else
      return theTest[i] = test(i, stat);
}

/*
    Testing all values of #coPvec# against its bounds. If $i$, say, is violated
    the violation is saved as negative value in |theCoTest[i]|.
 */
Real SPxSolver::coTest(int i, SPxBasis::Desc::Status stat) const
{
   METHOD( "SPxSolver::coTest()" );
   assert(type() == ENTER);
   assert(!isBasic(stat));

   Real x;

   switch (stat)
   {
   case SPxBasis::Desc::D_FREE:
   case SPxBasis::Desc::D_ON_BOTH :
      assert(rep() == ROW);
      x = (*theCoPvec)[i] - SPxLP::lower(i);
      if (x < 0)
         return x;
      // no break: next is else case
      //lint -fallthrough 
   case SPxBasis::Desc::D_ON_LOWER:
      assert(rep() == ROW);
      return SPxLP::upper(i) - (*theCoPvec)[i];
   case SPxBasis::Desc::D_ON_UPPER:
      assert(rep() == ROW);
      return (*theCoPvec)[i] - SPxLP::lower(i);

   case SPxBasis::Desc::P_ON_UPPER:
      assert(rep() == COLUMN);
      return (*theCoPvec)[i] - 0;             // slacks !
   case SPxBasis::Desc::P_ON_LOWER:
      assert(rep() == COLUMN);
      return 0 - (*theCoPvec)[i];             // slacks !

   default:
      return 0;
   }
}

void SPxSolver::computeCoTest()
{
   METHOD( "SPxSolver::computeCoTest()" );
   int i;
   Real pricingTol = leavetol();
   infeasibilities.clear();
   int ninfeasibilities = 0;
   const SPxBasis::Desc& ds = desc();

   for (i = dim() - 1; i >= 0; --i)
   {
      SPxBasis::Desc::Status stat = ds.coStatus(i);
      if (isBasic(stat))
         theCoTest[i] = 0;
      else
      {
         theCoTest[i] = coTest(i, stat);
         if( remainingRoundsEnter == 0 )
         {
            if( theCoTest[i] < -pricingTol )
            {
               assert(infeasibilities.size() < infeasibilities.max());
               infeasibilities.addIdx(i);
               isInfeasible[i] = true;
               ++ninfeasibilities;
            }
            else
               isInfeasible[i] = false;
            if( ninfeasibilities > sparsityThresholdEnter )
            {
               MSG_INFO2( spxout << "IENTER06 too many infeasibilities for sparse pricing"
                                 << std::endl; )
               remainingRoundsEnter = DENSEROUNDS;
               sparsePricingEnter = false;
               ninfeasibilities = 0;
            }
         }
      }
   }
   if( ninfeasibilities == 0 && !sparsePricingEnter )
      --remainingRoundsEnter;
   else if( ninfeasibilities <= sparsityThresholdEnter && !sparsePricingEnter )
   {
      MSG_INFO2( spxout << "IENTER05 sparse pricing active, "
                        << "sparsity: "
                        << std::setw(6) << std::fixed << std::setprecision(4)
                        << (Real) ninfeasibilities/dim()
                        << std::endl; )
      sparsePricingEnter = true;
   }
}


/*
    The following methods require propersy initialized vectors |fVec| and
    #coPvec#.
 */
void SPxSolver::updateTest()
{
   METHOD( "SPxSolver::updateTest()" );
   thePvec->delta().setup();

   const IdxSet& idx = thePvec->idx();
   const SPxBasis::Desc& ds = desc();
   Real pricingTol = epsilon();

   int i;
   for (i = idx.size() - 1; i >= 0; --i)
   {
      int j = idx.index(i);
      SPxBasis::Desc::Status stat = ds.status(j);
      if (!isBasic(stat))
      {
         theTest[j] = test(j, stat);

         if( sparsePricingEnterCo && theTest[j] < -pricingTol )
         {
            assert(remainingRoundsEnterCo == 0);
            if( !isInfeasibleCo[j] )
            {
               infeasibilitiesCo.addIdx(j);
               isInfeasibleCo[j] = true;
            }
         }
      }
      else
         theTest[j] = 0;
   }
}

void SPxSolver::updateCoTest()
{
   METHOD( "SPxSolver::updateCoTest()" );
   theCoPvec->delta().setup();

   const IdxSet& idx = theCoPvec->idx();
   const SPxBasis::Desc& ds = desc();
   Real pricingTol = epsilon();

   int i;
   for (i = idx.size() - 1; i >= 0; --i)
   {
      int j = idx.index(i);
      SPxBasis::Desc::Status stat = ds.coStatus(j);
      if (!isBasic(stat))
      {
         theCoTest[j] = coTest(j, stat);

         if( sparsePricingEnter && theCoTest[j] < -pricingTol )
         {
            assert(remainingRoundsEnter == 0);
            if( !isInfeasible[j] )
            {
               infeasibilities.addIdx(j);
               isInfeasible[j] = true;
            }
         }
      }
      else
         theCoTest[j] = 0;
   }
}



/*  \Section{Compute statistics on entering variable}
    Here is a list of variables relevant when including |Id| to the basis.
    They are computed by |computeEnterStats()|.
 */
void SPxSolver::getEnterVals
(
   SPxId enterId,
   Real& enterTest,
   Real& enterUB,
   Real& enterLB,
   Real& enterVal,
   Real& enterMax,
   Real& enterPric,
   SPxBasis::Desc::Status& enterStat,
   Real& enterRO
)
{
   METHOD( "SPxSolver::getEnterVals()" );
   int enterIdx;
   SPxBasis::Desc& ds = desc();

   if (enterId.isSPxColId())
   {
      enterIdx = number(SPxColId(enterId));
      enterStat = ds.colStatus(enterIdx);
      assert(!isBasic(enterStat));

      /*      For an #Id# to enter the basis we better recompute the Test value.
       */
      if (rep() == COLUMN)
      {
         computePvec(enterIdx);
         enterTest = computeTest(enterIdx);
         theTest[enterIdx] = 0;
      }
      else
      {
         enterTest = coTest()[enterIdx];
         theCoTest[enterIdx] = 0;
      }

      switch (enterStat)
      {
         // primal/columnwise cases:
      case SPxBasis::Desc::P_ON_UPPER :
         assert( rep() == COLUMN );
         enterUB = theUCbound[enterIdx];
         enterLB = theLCbound[enterIdx];
         enterVal = enterUB;
         enterMax = enterLB - enterUB;
         enterPric = (*thePvec)[enterIdx];
         enterRO = maxObj(enterIdx);
         if( enterLB <= -infinity )
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_ON_LOWER;
         else if( EQ( enterLB, enterUB ) )
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_FREE;
         else
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_ON_BOTH;
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert( rep() == COLUMN );
         enterUB = theUCbound[enterIdx];
         enterLB = theLCbound[enterIdx];
         enterVal = enterLB;
         enterMax = enterUB - enterLB;
         enterPric = (*thePvec)[enterIdx];
         enterRO = maxObj(enterIdx);
         if( enterUB >= infinity )
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_ON_UPPER;
         else if( EQ( enterLB, enterUB ) )
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_FREE;
         else
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_ON_BOTH;
         break;
      case SPxBasis::Desc::P_FREE :
         assert( rep() == COLUMN );
         enterUB = theUCbound[enterIdx];
         enterLB = theLCbound[enterIdx];
         enterVal = 0;
         enterPric = (*thePvec)[enterIdx];
         enterRO = maxObj(enterIdx);
         ds.colStatus(enterIdx) = SPxBasis::Desc::D_UNDEFINED;
         enterMax = (enterRO - enterPric > 0) ? infinity : -infinity;
         break;

         // dual/rowwise cases:
      case SPxBasis::Desc::D_ON_UPPER :
         assert( rep() == ROW );
         assert(theUCbound[enterIdx] < infinity);
         enterUB = theUCbound[enterIdx];
         enterLB = -infinity;
         enterMax = -infinity;
         enterVal = enterUB;
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = SPxLP::lower(enterIdx);
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert( rep() == ROW );
         assert(theLCbound[enterIdx] > -infinity);
         enterLB = theLCbound[enterIdx];
         enterUB = infinity;
         enterMax = infinity;
         enterVal = enterLB;
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = SPxLP::upper(enterIdx);
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case SPxBasis::Desc::D_FREE:
         assert( rep() == ROW );
         assert(SPxLP::lower(enterIdx) == SPxLP::upper(enterIdx));
         enterUB = infinity;
         enterLB = -infinity;
         enterVal = 0;
         enterRO = SPxLP::upper(enterIdx);
         enterPric = (*theCoPvec)[enterIdx];
         if (enterPric > enterRO)
            enterMax = infinity;
         else
            enterMax = -infinity;
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert( rep() == ROW );
         enterPric = (*theCoPvec)[enterIdx];
         if (enterPric > SPxLP::upper(enterIdx))
         {
            enterLB = theLCbound[enterIdx];
            enterUB = infinity;
            enterMax = infinity;
            enterVal = enterLB;
            enterRO = SPxLP::upper(enterIdx);
            ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         }
         else
         {
            enterUB = theUCbound[enterIdx];
            enterVal = enterUB;
            enterRO = SPxLP::lower(enterIdx);
            enterLB = -infinity;
            enterMax = -infinity;
            ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         }
         break;
      default:
         throw SPxInternalCodeException("XENTER01 This should never happen.");
      }
      MSG_DEBUG( spxout << "DENTER03 SPxSolver::getEnterVals() : col " << enterIdx
                        << ": " << enterStat
                        << " -> " << ds.colStatus(enterIdx)
                        << std::endl; )
   }

   else
   {
      assert(enterId.isSPxRowId());
      enterIdx = number(SPxRowId(enterId));
      enterStat = ds.rowStatus(enterIdx);
      assert(!isBasic(enterStat));

      /*      For an #Id# to enter the basis we better recompute the Test value.
       */
      if (rep() == ROW)
      {
         computePvec(enterIdx);
         enterTest = computeTest(enterIdx);
         theTest[enterIdx] = 0;
      }
      else
      {
         enterTest = coTest()[enterIdx];
         theCoTest[enterIdx] = 0;
      }

      switch (enterStat)
      {
         // primal/columnwise cases:
      case SPxBasis::Desc::P_ON_UPPER :
         assert( rep() == COLUMN );
         enterUB = theURbound[enterIdx];
         enterLB = theLRbound[enterIdx];
         enterVal = enterLB;
         enterMax = enterUB - enterLB;
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = 0;
         if( enterUB >= infinity )
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_ON_LOWER;
         else if( EQ( enterLB, enterUB ) )
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_FREE;
         else
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_ON_BOTH;
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert( rep() == COLUMN );
         enterUB = theURbound[enterIdx];
         enterLB = theLRbound[enterIdx];
         enterVal = enterUB;
         enterMax = enterLB - enterUB;
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = 0;
         if( enterLB <= -infinity )
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_ON_UPPER;
         else if( EQ( enterLB, enterUB ) )
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_FREE;
         else
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_ON_BOTH;
         break;
      case SPxBasis::Desc::P_FREE :
         assert( rep() == COLUMN );
#if 1
         throw SPxInternalCodeException("XENTER02 This should never happen.");
#else
         MSG_ERROR( spxout << "EENTER99 ERROR: not yet debugged!" << std::endl; )
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = 0;
         ds.rowStatus(enterIdx) = SPxBasis::Desc::D_UNDEFINED;
#endif
         break;

         // dual/rowwise cases:
      case SPxBasis::Desc::D_ON_UPPER :
         assert( rep() == ROW );
         assert(theURbound[enterIdx] < infinity);
         enterUB = theURbound[enterIdx];
         enterLB = -infinity;
         enterVal = enterUB;
         enterMax = -infinity;
         enterPric = (*thePvec)[enterIdx];
         enterRO = lhs(enterIdx);
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert( rep() == ROW );
         assert(theLRbound[enterIdx] > -infinity);
         enterLB = theLRbound[enterIdx];
         enterUB = infinity;
         enterVal = enterLB;
         enterMax = infinity;
         enterPric = (*thePvec)[enterIdx];
         enterRO = rhs(enterIdx);
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case SPxBasis::Desc::D_FREE:
         assert( rep() == ROW );
         assert(rhs(enterIdx) == lhs(enterIdx));
         enterUB = infinity;
         enterLB = -infinity;
         enterVal = 0;
         enterPric = (*thePvec)[enterIdx];
         enterRO = rhs(enterIdx);
         enterMax = (enterPric > enterRO) ? infinity : -infinity;
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert( rep() == ROW );
         enterPric = (*thePvec)[enterIdx];
         if (enterPric > rhs(enterIdx))
         {
            enterLB = theLRbound[enterIdx];
            enterVal = enterLB;
            enterUB = infinity;
            enterMax = infinity;
            enterRO = rhs(enterIdx);
            ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         }
         else
         {
            enterUB = theURbound[enterIdx];
            enterVal = enterUB;
            enterLB = -infinity;
            enterMax = -infinity;
            enterRO = lhs(enterIdx);
            ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         }
         break;

      default:
         throw SPxInternalCodeException("XENTER03 This should never happen.");
      }
      MSG_DEBUG( spxout << "DENTER05 SPxSolver::getEnterVals() : row " 
                        << enterIdx << ": " << enterStat
                        << " -> " << ds.rowStatus(enterIdx)
                        << std::endl; )
   }
}

/*      process leaving variable
 */
void SPxSolver::getEnterVals2
(
   int leaveIdx,
   Real enterMax,
   Real& leavebound
)
{
   METHOD( "SPxSolver::getEnterVals2()" );
   int idx;
   SPxBasis::Desc& ds = desc();
   SPxId leftId = baseId(leaveIdx);

   if (leftId.isSPxRowId())
   {
      idx = number(SPxRowId(leftId));
      SPxBasis::Desc::Status leaveStat = ds.rowStatus(idx);

      switch (leaveStat)
      {
      case SPxBasis::Desc::P_FIXED :
         assert(rep() == ROW);
         throw SPxInternalCodeException("XENTER04 This should never happen.");
         break;
      case SPxBasis::Desc::P_ON_UPPER :
         assert(rep() == ROW);
         leavebound = theLBbound[leaveIdx];
         theLRbound[idx] = leavebound;
         ds.rowStatus(idx) = dualRowStatus(idx);
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert(rep() == ROW);
         leavebound = theUBbound[leaveIdx];
         theURbound[idx] = leavebound;
         ds.rowStatus(idx) = dualRowStatus(idx);
         break;
      case SPxBasis::Desc::P_FREE :
         assert(rep() == ROW);
#if 1
         throw SPxInternalCodeException("XENTER05 This should never happen.");
#else
         MSG_ERROR( spxout << "EENTER98 ERROR: not yet debugged!" << std::endl; )
         if ((*theCoPvec)[leaveIdx] - theLBbound[leaveIdx] <
              theUBbound[leaveIdx] - (*theCoPvec)[leaveIdx])
         {
            leavebound = theLBbound[leaveIdx];
            theLRbound[idx] = leavebound;
         }
         else
         {
            leavebound = theUBbound[leaveIdx];
            theURbound[idx] = leavebound;
         }
         ds.rowStatus(idx) = SPxBasis::Desc::D_UNDEFINED;
#endif
         break;
         // primal/columnwise cases:
      case SPxBasis::Desc::D_UNDEFINED :
         assert(rep() == COLUMN);
         throw SPxInternalCodeException("XENTER06 This should never happen.");
         break;
      case SPxBasis::Desc::D_FREE :
         assert(rep() == COLUMN);
         if (theFvec->delta()[leaveIdx] * enterMax < 0)
            leavebound = theUBbound[leaveIdx];
         else
            leavebound = theLBbound[leaveIdx];
         theLRbound[idx] = leavebound;
         theURbound[idx] = leavebound;
         ds.rowStatus(idx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         assert(rep() == COLUMN);
         leavebound = theUBbound[leaveIdx];
         theURbound[idx] = leavebound;
         ds.rowStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert(rep() == COLUMN);
         leavebound = theLBbound[leaveIdx];
         theLRbound[idx] = leavebound;
         ds.rowStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert(rep() == COLUMN);
         if (enterMax * theFvec->delta()[leaveIdx] < 0)
         {
            leavebound = theUBbound[leaveIdx];
            theURbound[idx] = leavebound;
            ds.rowStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         }
         else
         {
            leavebound = theLBbound[leaveIdx];
            theLRbound[idx] = leavebound;
            ds.rowStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         }
         break;

      default:
         throw SPxInternalCodeException("XENTER07 This should never happen.");
      }
      MSG_DEBUG( spxout << "DENTER06 SPxSolver::getEnterVals2(): row " 
                        << idx << ": " << leaveStat
                        << " -> " << ds.rowStatus(idx)
                        << std::endl; )
   }

   else
   {
      assert(leftId.isSPxColId());
      idx = number(SPxColId(leftId));
      SPxBasis::Desc::Status leaveStat = ds.colStatus(idx);

      switch (leaveStat)
      {
      case SPxBasis::Desc::P_ON_UPPER :
         assert(rep() == ROW);
         leavebound = theLBbound[leaveIdx];
         theLCbound[idx] = leavebound;
         ds.colStatus(idx) = dualColStatus(idx);
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert(rep() == ROW);
         leavebound = theUBbound[leaveIdx];
         theUCbound[idx] = leavebound;
         ds.colStatus(idx) = dualColStatus(idx);
         break;
      case SPxBasis::Desc::P_FREE :
         assert(rep() == ROW);
         if (theFvec->delta()[leaveIdx] * enterMax > 0)
         {
            leavebound = theLBbound[leaveIdx];
            theLCbound[idx] = leavebound;
         }
         else
         {
            leavebound = theUBbound[leaveIdx];
            theUCbound[idx] = leavebound;
         }
         ds.colStatus(idx) = SPxBasis::Desc::D_UNDEFINED;
         break;
      case SPxBasis::Desc::P_FIXED:
         assert(rep() == ROW);
         throw SPxInternalCodeException("XENTER08 This should never happen.");
         break;
         // primal/columnwise cases:
      case SPxBasis::Desc::D_FREE :
         assert(rep() == COLUMN);
         if (theFvec->delta()[leaveIdx] * enterMax > 0)
            leavebound = theLBbound[leaveIdx];
         else
            leavebound = theUBbound[leaveIdx];
         theUCbound[idx] =
            theLCbound[idx] = leavebound;
         ds.colStatus(idx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         assert(rep() == COLUMN);
         leavebound = theLBbound[leaveIdx];
         theLCbound[idx] = leavebound;
         ds.colStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert(rep() == COLUMN);
         leavebound = theUBbound[leaveIdx];
         theUCbound[idx] = leavebound;
         ds.colStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
      case SPxBasis::Desc::D_UNDEFINED :
         assert(rep() == COLUMN);
         if (enterMax * theFvec->delta()[leaveIdx] < 0)
         {
            leavebound = theUBbound[leaveIdx];
            theUCbound[idx] = leavebound;
            ds.colStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         }
         else
         {
            leavebound = theLBbound[leaveIdx];
            theLCbound[idx] = leavebound;
            ds.colStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         }
         break;

      default:
         throw SPxInternalCodeException("XENTER09 This should never happen.");
      }
      MSG_DEBUG( spxout << "DENTER07 SPxSolver::getEnterVals2(): col " 
                        << idx << ": " << leaveStat
                        << " -> " << ds.colStatus(idx)
                        << std::endl; )
   }
}


void
SPxSolver::ungetEnterVal(
   SPxId enterId,
   SPxBasis::Desc::Status enterStat,
   Real leaveVal,
   const SVector& vec
)
{
   METHOD( "SPxSolver::ungetEnterVal()" );
   int enterIdx;
   SPxBasis::Desc& ds = desc();

   if (enterId.isSPxColId())
   {
      enterIdx = number(SPxColId(enterId));
      if (enterStat == SPxBasis::Desc::P_ON_UPPER)
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
      else
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
      theFrhs->multAdd(leaveVal, vec);
   }
   else
   {
      enterIdx = number(SPxRowId(enterId));
      assert(enterId.isSPxRowId());
      if (enterStat == SPxBasis::Desc::P_ON_UPPER)
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
      else
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
      (*theFrhs)[enterIdx] += leaveVal;
   }
   if (isId(enterId))
      theTest[enterIdx] = 0;
   else
      theCoTest[enterIdx] = 0;
}

void SPxSolver::rejectEnter(
   SPxId enterId,
   Real enterTest,
   SPxBasis::Desc::Status enterStat
)
{
   METHOD( "SPxSolver::rejectEnter()" );
   int enterIdx = number(enterId);
   if (isId(enterId))
   {
      theTest[enterIdx] = enterTest;
      desc().status(enterIdx) = enterStat;
   }
   else
   {
      theCoTest[enterIdx] = enterTest;
      desc().coStatus(enterIdx) = enterStat;
   }
}

bool SPxSolver::enter(SPxId& enterId)
{
   METHOD( "SPxSolver::enter()" );
   assert(enterId.isValid());
   assert(type() == ENTER);
   assert(initialized);

   SPxId none;          // invalid id used when enter fails
   Real enterTest;      // correct test value of entering var
   Real enterUB;        // upper bound of entering variable
   Real enterLB;        // lower bound of entering variable
   Real enterVal;       // current value of entering variable
   Real enterMax;       // maximum value for entering shift
   Real enterPric;      // priced value of entering variable
   SPxBasis::Desc::Status enterStat;      // status of entering variable
   Real enterRO;        // rhs/obj of entering variable
   const SVector* enterVec = enterVector(enterId);

   getEnterVals(enterId, enterTest, enterUB, enterLB,
      enterVal, enterMax, enterPric, enterStat, enterRO);

   if (enterTest > -epsilon())
   {
      rejectEnter(enterId, enterTest, enterStat);
      change(-1, none, 0);

      MSG_DEBUG( spxout << "DENTER08 rejecting false enter pivot" << std::endl; )

      return true;
   }

   /*  Before performing the actual basis update, we must determine, how this
       is to be accomplished.
    */
   // BH 2005-11-15: Obviously solve4update() is only called if theFvec.delta()
   // is setup (i.e. the indices of the NZEs are stored within it) and there are
   // 0 NZEs (???).
   // In that case theFvec->delta() is set such that
   //   Base * theFvec->delta() = enterVec 
   if (theFvec->delta().isSetup() && theFvec->delta().size() == 0)
      SPxBasis::solve4update(theFvec->delta(), *enterVec);
#ifdef ENABLE_ADDITIONAL_CHECKS
   else
   {
      // BH 2005-11-29: This code block seems to check the assertion
      //   || Base * theFvec->delta() - enterVec ||_2 <= entertol()
      DVector tmp(dim());
      // BH 2005-11-15: This cast is necessary since SSVector inherits protected from DVector.
      tmp = reinterpret_cast<DVector&>(theFvec->delta());
      multBaseWith(tmp);
      tmp -= *enterVec;
      if (tmp.length() > entertol()) {
         // This happens frequently and does usually not hurt, so print these
         // warnings only with verbose level INFO2 and higher.
         MSG_INFO2( spxout << "WENTER09 fVec updated error = " 
                              << tmp.length() << std::endl; )
      }
   }
#endif  // ENABLE_ADDITIONAL_CHECKS

   if (m_numCycle > m_maxCycle)
   {
      if (-enterMax > 0)
         perturbMaxEnter();
      else
         perturbMinEnter();
   }

   Real leaveVal = -enterMax;

   int leaveIdx = theratiotester->selectLeave(leaveVal, enterId);

   /* in row representation, fixed columns and rows should not leave the basis */
   assert(leaveIdx < 0 || !baseId(leaveIdx).isSPxColId() || desc().colStatus(number(SPxColId(baseId(leaveIdx)))) != SPxBasis::Desc::P_FIXED);
   assert(leaveIdx < 0 || !baseId(leaveIdx).isSPxRowId() || desc().rowStatus(number(SPxRowId(baseId(leaveIdx)))) != SPxBasis::Desc::P_FIXED);

   /*
       We now tried to find a variable to leave the basis. If one has been
       found, a regular basis update is to be performed.
    */
   if (leaveIdx >= 0)
   {
      if (fabs(leaveVal) < entertol())
      {
         if (theUBbound[leaveIdx] != theLBbound[leaveIdx] 
            && enterStat != Desc::P_FREE && enterStat != Desc::D_FREE) 
            m_numCycle++;
      }
      else
         m_numCycle /= 2;

      // setup for updating the copricing vector
      if (coSolveVector2)
         SPxBasis::coSolve(theCoPvec->delta(), *coSolveVector2, unitVecs[leaveIdx], *coSolveVector2rhs);
      else
         SPxBasis::coSolve(theCoPvec->delta(), unitVecs[leaveIdx]);

      (*theCoPrhs)[leaveIdx] = enterRO;
      theCoPvec->value() = (enterRO - enterPric) / theFvec->delta()[leaveIdx];

      if (theCoPvec->value() > epsilon() || theCoPvec->value() < -epsilon())
      {
         if (pricing() == FULL)
         {
            thePvec->value() = theCoPvec->value();
            setupPupdate();
         }
         doPupdate();
      }
      assert(thePvec->isConsistent());
      assert(theCoPvec->isConsistent());

      assert(!baseId(leaveIdx).isSPxRowId() || desc().rowStatus(number(SPxRowId(baseId(leaveIdx)))) != SPxBasis::Desc::P_FIXED);
      assert(!baseId(leaveIdx).isSPxColId() || desc().colStatus(number(SPxColId(baseId(leaveIdx)))) != SPxBasis::Desc::P_FIXED);

      Real leavebound;             // bound on which leaving variable moves
      try
      {
         getEnterVals2(leaveIdx, enterMax, leavebound);
      }
      catch( SPxException F )
      {
         rejectEnter(enterId, enterTest, enterStat);
         change(-1, none, 0);
         throw F;
      }

      //  process entering variable
      theUBbound[leaveIdx] = enterUB;
      theLBbound[leaveIdx] = enterLB;

      //  compute tests:
      updateCoTest();
      if (pricing() == FULL)
         updateTest();


      // update feasibility vectors
      theFvec->value() = leaveVal;
      theFvec->update();
      (*theFvec)[leaveIdx] = enterVal - leaveVal;

      if (leavebound > epsilon() || leavebound < -epsilon())
         theFrhs->multAdd(-leavebound, baseVec(leaveIdx));

      if (enterVal > epsilon() || enterVal < -epsilon())
         theFrhs->multAdd(enterVal, *enterVec);


      //  change basis matrix
      change(leaveIdx, enterId, enterVec, &(theFvec->delta()));
   }
   /*  No leaving vector could be found that would yield a stable pivot step.
    */
   else if (leaveVal != -enterMax)
   {
      rejectEnter(enterId, REAL(0.01) * enterTest - REAL(2.0) * leavetol(), enterStat);
      change(-1, none, 0);
   }
   /*  No leaving vector has been selected from the basis. However, if the
       shift amount for |fVec| is bounded, we are in the case, that the
       entering variable is moved from one bound to its other, before any of
       the basis feasibility variables reaches their bound. This may only
       happen in primal/columnwise case with upper and lower bounds on
       variables.
    */
   else if (leaveVal < infinity && leaveVal > -infinity)
   {
      assert(rep() == COLUMN);
      assert(leaveVal == -enterMax);

      change(-1, enterId, enterVec);

      theFvec->value() = leaveVal;
      theFvec->update();

      ungetEnterVal(enterId, enterStat, leaveVal, *enterVec);
   }
   /*  No variable could be selected to leave the basis and even the entering
       variable is unbounded --- this is a failure.  
    */
   else
   {
      /* The following line originally was in the "lastUpdate() > 1" case;
         we need it in the INFEASIBLE/UNBOUNDED case, too, to have the
         basis descriptor at the correct size. 
       */
      rejectEnter(enterId, enterTest, enterStat);
      change(-1, none, 0);

      if (lastUpdate() > 1)
      {
         MSG_INFO3( spxout << "IENTER01 factorization triggered in "
                              << "enter() for feasibility test" << std::endl; )
         factorize();

         /* after a factorization, the entering column/row might not be infeasible or suboptimal anymore, hence we do
          * not try to call leave(leaveIdx), but rather return to the main solving loop and call the pricer again
          */
         return true;
      }

      MSG_INFO3( spxout << "IENTER02 unboundness/infeasiblity found in "
                           << "enter()" << std::endl; )

      if (rep() == ROW)
      {
         Real sign;

         dualFarkas.clear();
         dualFarkas.setMax(fVec().delta().size() + 1);
         sign = (leaveVal > 0 ? -1.0 : 1.0);

         for( int j = 0; j < fVec().delta().size(); ++j )
         {
            SPxId id = baseId(fVec().idx().index(j));

            if( id.isSPxRowId() )
               dualFarkas.add(number(SPxRowId(id)), sign * fVec().delta().value(j));
         }

         if( enterId.isSPxRowId() )
            dualFarkas.add(number(SPxRowId(enterId)), -sign);

         setBasisStatus(SPxBasis::INFEASIBLE);
      }
      /**@todo if shift() is not zero, we must not conclude primal unboundedness */
      else
      {
         Real sign;

         primalRay.clear();
         primalRay.setMax(fVec().delta().size() + 1);
         sign = leaveVal > 0 ? 1.0 : -1.0;

         for( int j = 0; j < fVec().delta().size(); ++j )
         {
            SPxId i = baseId(fVec().idx().index(j));

            if( i.isSPxColId() )
               primalRay.add(number(SPxColId(i)), sign*fVec().delta().value(j));
         }

         if( enterId.isSPxColId() )
            primalRay.add(number(SPxColId(enterId)), -sign);

         setBasisStatus(SPxBasis::UNBOUNDED);
      }
      return false;
   }
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
