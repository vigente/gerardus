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

namespace soplex
{
/** Setting up the feasiblity bound for normal primal variables is
    straightforward. However, slack variables need some more details
    on how we treate them. This is slightly different from usual
    textbook versions. Let \f$l_i \le A_i^T x \le u_i\f$. This can be
    transformed to \f$A_i^Tx + s_i = 0\f$, with \f$-u_i \le s_i \le
    -l_i\f$. Hence, with this definition of slack variables \f$s_i\f$, we
    can directly use vectors \f$l\f$ and \f$u\f$ as feasibility bounds.  
 */
void SPxSolver::setPrimalBounds()
{
   METHOD( "SPxSolver::setPrimalBounds()" );

   theUCbound = SPxLP::upper();
   theLCbound = SPxLP::lower();

   if (rep() == ROW)
   {
      theURbound = rhs();
      theLRbound = lhs();
   }
   else
   {
      theURbound = lhs();
      theLRbound = rhs();
      theURbound *= -1.0;
      theLRbound *= -1.0;
   }
}

/** Seting up the basis for dual simplex requires to install upper and lower
    feasibility bounds for dual variables (|Lbound| and |Ubound|). Here is a
    list of how these must be set for inequalities of type \f$l \le a^Tx \le u\f$:
 
    \f[
        \begin{tabular}{cccc}
            $l$         &       $u$     & |Lbound|      & |Ubound|      \\
        \hline
        $-\infty=l$     & $u=\infty$    &       0       &       0       \\
        $-\infty<l$     & $u=\infty$    &       0       & $\infty$      \\
        $-\infty=l$     & $u<\infty$    & $-\infty$     &       0       \\
        \multicolumn{2}{c}{
        $-\infty<l \ne u<\infty$}       &       0       &       0       \\
        \multicolumn{2}{c}{
        $-\infty<l  =  u<\infty$}       & $-\infty$     & $\infty$      \\
        \end{tabular}
   \f]
   The case \f$l = -\infty\f$, \f$u = \infty\f$ occurs for unbounded primal variables.
   Such must be treated differently from the general case.
 
   Given possible upper and lower bounds to a dual variable with |Status stat|,
   this function clears the bounds according to |stat| by setting them to
   \f$\infty\f$ or \f$-\infty\f$, respectively.
 */
void SPxSolver::clearDualBounds(
   SPxBasis::Desc::Status stat,
   Real&                  upp,
   Real&                  lw) const
{
   METHOD( "SPxSolver::clearDualBounds()" );

   switch (stat)
   {
   case SPxBasis::Desc::P_ON_UPPER + SPxBasis::Desc::P_ON_LOWER :
   case SPxBasis::Desc::D_FREE :
      upp = infinity;
      lw  = -infinity;
      break;
   case SPxBasis::Desc::P_ON_UPPER :
   case SPxBasis::Desc::D_ON_LOWER :
      upp = infinity;
      break;
   case SPxBasis::Desc::P_ON_LOWER :
   case SPxBasis::Desc::D_ON_UPPER :
      lw  = -infinity;
      break;

   default:
      break;
   }
}

void SPxSolver::setDualColBounds()
{
   METHOD( "SPxSolver::setDualColBounds()" );

   assert(rep() == COLUMN);

   const SPxBasis::Desc& ds = desc();

   int i;

   for(i = 0; i < nRows(); ++i)
   {
      theURbound[i] = 0.0;
      theLRbound[i] = 0.0;

      clearDualBounds(ds.rowStatus(i), theURbound[i], theLRbound[i]);
   }

   for(i = 0; i < nCols(); ++i)
   {
      theUCbound[i] = -maxObj(i);
      theLCbound[i] = -maxObj(i);

      // exchanged ...                 due to definition of slacks!
      clearDualBounds(ds.colStatus(i), theLCbound[i], theUCbound[i]);

      theUCbound[i] *= -1.0;
      theLCbound[i] *= -1.0;
   }
}

void SPxSolver::setDualRowBounds()
{
   METHOD( "SPxSolver::setDualRowBounds()" );

   assert(rep() == ROW);

   int i;

   for(i = 0; i < nRows(); ++i)
   {
      theURbound[i] = 0.0;
      theLRbound[i] = 0.0;

      clearDualBounds(dualRowStatus(i), theURbound[i], theLRbound[i]);
   }

   for(i = 0; i < nCols(); ++i)
   {
      theUCbound[i] = 0.0;
      theLCbound[i] = 0.0;

      clearDualBounds(dualColStatus(i), theUCbound[i], theLCbound[i]);
   }
}


/** This sets up the bounds for basic variables for entering simplex algorithms.
    It requires, that all upper lower feasibility bounds have allready been
    setup. Method |setEnterBound4Row(i, n)| does so for the |i|-th basis
    variable being row index |n|. Equivalently, method
    |setEnterBound4Col(i, n)| does so for the |i|-th basis variable being
    column index |n|.
 */
void SPxSolver::setEnterBound4Row(int i, int n)
{
   METHOD( "SPxSolver::setEnterBound4Row()" );
   assert(baseId(i).isSPxRowId());
   assert(number(SPxRowId(baseId(i))) == n);
   switch (desc().rowStatus(n))
   {
   case SPxBasis::Desc::P_ON_LOWER :
      theLBbound[i] = -infinity;
      theUBbound[i] = theURbound[n];
      break;
   case SPxBasis::Desc::P_ON_UPPER :
      theLBbound[i] = theLRbound[n];
      theUBbound[i] = infinity;
      break;

   default:
      theUBbound[i] = theURbound[n];
      theLBbound[i] = theLRbound[n];
      break;
   }
}

void SPxSolver::setEnterBound4Col(int i, int n)
{
   METHOD( "SPxSolver::setEnterBound4Col()" );
   assert(baseId(i).isSPxColId());
   assert(number(SPxColId(baseId(i))) == n);
   switch (desc().colStatus(n))
   {
   case SPxBasis::Desc::P_ON_LOWER :
      theLBbound[i] = -infinity;
      theUBbound[i] = theUCbound[n];
      break;
   case SPxBasis::Desc::P_ON_UPPER :
      theLBbound[i] = theLCbound[n];
      theUBbound[i] = infinity;
      break;

   default:
      theUBbound[i] = theUCbound[n];
      theLBbound[i] = theLCbound[n];
      break;
   }
}

void SPxSolver::setEnterBounds()
{
   METHOD( "SPxSolver::setEnterBounds()" );

   for (int i = 0; i < dim(); ++i)
   {
      SPxId base_id = baseId(i);

      if (base_id.isSPxRowId())
         setEnterBound4Row(i, number(SPxRowId(base_id)));
      else
         setEnterBound4Col(i, number(SPxColId(base_id)));
   }
}


/** This sets up the bounds for basic variables for leaving simplex algorithms.
    While method |setLeaveBound4Row(i,n)| does so for the |i|-th basic variable
    being the |n|-th row, |setLeaveBound4Col(i,n)| does so for the |i|-th basic
    variable being the |n|-th column.
 */
void SPxSolver::setLeaveBound4Row(int i, int n)
{
   METHOD( "SPxSolver::setLeaveBound4Row()" );
   assert(baseId(i).isSPxRowId());
   assert(number(SPxRowId(baseId(i))) == n);
   switch (desc().rowStatus(n))
   {
   case SPxBasis::Desc::P_ON_LOWER :
      theLBbound[i] = -infinity;
      theUBbound[i] = 0.0;
      break;
   case SPxBasis::Desc::P_ON_UPPER :
      theLBbound[i] = 0.0;
      theUBbound[i] = infinity;
      break;
   case SPxBasis::Desc::P_ON_UPPER + SPxBasis::Desc::P_ON_LOWER :
      theLBbound[i] = -infinity;
      theUBbound[i] = infinity;
      break;
   case SPxBasis::Desc::P_FREE :
      theLBbound[i] = 0.0;
      theUBbound[i] = 0.0;
      break;

   default:
      assert(rep() == COLUMN);
      theLBbound[i] = -rhs(n);                // slacks !
      theUBbound[i] = -lhs(n);                // slacks !
      break;
   }
}

void SPxSolver::setLeaveBound4Col(int i, int n)
{
   METHOD( "SPxSolver::setLeaveBound4Col()" );

   assert(baseId(i).isSPxColId());
   assert(number(SPxColId(baseId(i))) == n);

   switch (desc().colStatus(n))
   {
   case SPxBasis::Desc::P_ON_LOWER :
      theLBbound[i] = -infinity;
      theUBbound[i] = 0;
      break;
   case SPxBasis::Desc::P_ON_UPPER :
      theLBbound[i] = 0;
      theUBbound[i] = infinity;
      break;
   case SPxBasis::Desc::P_FIXED :
      theLBbound[i] = -infinity;
      theUBbound[i] = infinity;
      break;
   case SPxBasis::Desc::P_FREE :
      theLBbound[i] = theUBbound[i] = 0;
      break;

   default:
      theUBbound[i] = SPxLP::upper(n);
      theLBbound[i] = SPxLP::lower(n);
      break;
   }
}

void SPxSolver::setLeaveBounds()
{
   METHOD( "SPxSolver::setLeaveBounds()" );

   for (int i = 0; i < dim(); ++i)
   {
      SPxId base_id = baseId(i);

      if (base_id.isSPxRowId())
         setLeaveBound4Row(i, number(SPxRowId(base_id)));
      else
         setLeaveBound4Col(i, number(SPxColId(base_id)));
   }
}

void SPxSolver::testBounds() const
{
   METHOD( "SPxSolver::testBounds()" );

   if (type() == ENTER)
   {
      Real viol_max = (1 + iterCount) * entertol();
      int nlinesprinted = 0;

      for(int i = 0; i < dim(); ++i )
      {
         // Minor bound violations happen frequently, so print these
         // warnings only with verbose level INFO2 and higher.
         if ((*theFvec)[i] > theUBbound[i] + viol_max)  //@ &&  theUBbound[i] != theLBbound[i])
         {
            MSG_INFO2( spxout << "WBOUND01 Invalid upper enter bound " << i 
                            << " viol_max: " << viol_max
                            << " Fvec: " << (*theFvec)[i] 
                            << " UBbound: "<< theUBbound[i] << std::endl; )
            nlinesprinted++;
         }
         if ((*theFvec)[i] < theLBbound[i] - viol_max)  //@ &&  theUBbound[i] != theLBbound[i])
         {
            MSG_INFO2( spxout << "WBOUND02 Invalid lower enter bound " << i 
                            << " viol_max: " << viol_max
                            << " Fvec: " << (*theFvec)[i] 
                            << " LBbound: "<< theLBbound[i] << std::endl; )
            nlinesprinted++;
         }
         if( nlinesprinted >= 3 )
         {
            MSG_INFO2( spxout << "WBOUND10 suppressing further warnings of type WBOUND{01,02} in this round" << std::endl );
            break;
         }
      }
   }
   else
   {
      Real viol_max = (1 + iterCount) * leavetol();
      int nlinesprinted = 0;

      for(int i = 0; i < dim(); ++i )
      {
         if ((*theCoPvec)[i] > (*theCoUbound)[i] + viol_max) // && (*theCoUbound)[i] != (*theCoLbound)[i])
         {
            MSG_INFO2( spxout << "WBOUND03 Invalid upper cobound " << i
                            << " viol_max: " << viol_max
                            << " CoPvec: " << (*theCoPvec)[i]
                            << " CoUbound: "<< (*theCoUbound)[i] << std::endl; )
            nlinesprinted++;
         }
         if ((*theCoPvec)[i] < (*theCoLbound)[i] - viol_max) // && (*theCoUbound)[i] != (*theCoLbound)[i])
         {
            MSG_INFO2( spxout << "WBOUND04 Invalid lower cobound " << i 
                            << " viol_max: " << viol_max
                            << " CoPvec: " << (*theCoPvec )[i]
                            << " CoLbound: " << (*theCoLbound)[i] 
                            << std::endl; )
            nlinesprinted++;
         }
         if( nlinesprinted >= 3 )
         {
            MSG_INFO2( spxout << "WBOUND11 suppressing further warnings of type WBOUND{03,04} in this round" << std::endl );
            break;
         }
      }

      nlinesprinted = 0;
      for(int i = 0; i < coDim(); ++i )
      {
         if ((*thePvec)[i] > (*theUbound)[i] + viol_max)  // &&  (*theUbound)[i] != (*theLbound)[i])
         {
            MSG_INFO2( spxout << "WBOUND05 Invalid upper bound " << i 
                            << " viol_max: " << viol_max
                            << " Pvec: " << (*thePvec)[i]
                            << " Ubound: " << (*theUbound)[i] << std::endl; )
            nlinesprinted++;
         }
         if ((*thePvec)[i] < (*theLbound)[i] - viol_max)  // &&  (*theUbound)[i] != (*theLbound)[i])
         {
            MSG_INFO2( spxout << "WBOUND06 Invalid lower bound " << i 
                            << " viol_max: " << viol_max
                            << " Pvec: " << (*thePvec)[i]
                            << " Lbound: " << (*theLbound)[i] << std::endl; )
            nlinesprinted++;
         }
         if( nlinesprinted >= 3 )
         {
            MSG_INFO2( spxout << "WBOUND12 suppressing further warnings of type WBOUND{05,06} in this round" << std::endl );
            break;
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
