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
#include "spxpricer.h"
#include "spxratiotester.h"
#include "exceptions.h"

namespace soplex
{

#if 0
void SPxSolver::localAddRows(int start)
{
   METHOD( "SPxSolver::localAddRows()" );
   assert( start <= SPxLP::nRows() );

   /**@todo This method seems to be called, to update
    *       theFvec, theFrhs, ..., but a resolve after
    *       adding a row results in a failure.
    *       To fix this, we call unInit() so that init() is called before solving
    *       in spxsolve.cpp:solve(). In init(), the
    *       vectors are set up, so there is no need
    *       to update them here.
    */
   if( start == SPxLP::nRows() )
      return;

   const SPxBasis::Desc& ds = desc();

   if (type() == ENTER)
   {
      if (rep() == COLUMN)
      {
         int i;
         for (i = start; i < SPxLP::nRows(); ++i)
         {
            theURbound[i] = -lhs(i);
            theLRbound[i] = -rhs(i);
            setEnterBound4Row(i, i);
            computeEnterCoPrhs4Row(i, i);
            // init #theFrhs[i]#:
            Real& v_rhs = (*theFrhs)[i];
            const SVector& row = rowVector(i); 
            v_rhs = 0;
            for (int j = row.size() - 1; j >= 0; --j)
            {
               int idx = row.index(j);
               switch (ds.colStatus(idx))
               {
               case Desc::P_ON_UPPER:
                  v_rhs += row.value(j) * theUCbound[idx];
                  break;
               case Desc::P_ON_LOWER:
               case Desc::P_FIXED:
                  v_rhs += row.value(j) * theLCbound[idx];
                  break;
               default:
                  break;
               }
            }
         }
         SPxBasis::solve (*theFvec, *theFrhs);
         SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
         for (i = start; i < SPxLP::nRows(); ++i)
         {
            /* we need to compare with tolerance (rep() == COLUMN) ? feastol() : opttol() because theFvec is the primal
             * vector in COLUMN and the dual vector in ROW representation; this is equivalent to entertol(); this also
             * fits because we are within the "type() == ENTER" case
             */
            if (theUBbound[i] + entertol() < (*theFvec)[i])
               shiftUBbound(i, (*theFvec)[i]);
            else if ((*theFvec)[i] < theLBbound[i] - entertol())
               shiftLBbound(i, (*theFvec)[i]);
         }
         computePvec();
         computeCoTest();
         computeTest();
      }
      else
      {
         assert(rep() == ROW);
         for (int i = start; i < SPxLP::nRows(); ++i)
         {
            theURbound[i] = theLRbound[i] = 0;
            clearDualBounds(dualRowStatus(i),
                             theURbound[i], theLRbound[i]);
            (*thePvec)[i] = vector(i) * (*theCoPvec);
            theTest[i] = test(i, ds.status(i));
         }
      }
   }
   else
   {
      assert(type() == LEAVE);
      if (rep() == ROW)
      {
         for (int i = start; i < SPxLP::nRows(); ++i)
         {
            theURbound[i] = rhs(i);
            theLRbound[i] = lhs(i);
            (*thePvec)[i] = vector(i) * (*theCoPvec);

            /* we need to compare with tolerance (rep() == ROW) ? feastol() : opttol() because thePvec is the primal
             * vector in ROW and the dual vector in COLUMN representation; this is equivalent to leavetol(); this also
             * fits because we are within the "type() == LEAVE" case
             */
            if (theURbound[i] + leavetol() < (*thePvec)[i])
               shiftUPbound(i, (*thePvec)[i]);
            else if ((*thePvec)[i] < theLRbound[i] - leavetol())
               shiftLPbound(i, (*thePvec)[i]);
         }
      }
      else
      {
         assert(rep() == COLUMN);
         int i;
         for (i = start; i < SPxLP::nRows(); ++i)
         {
            theURbound[i] = theLRbound[i] = 0;
            clearDualBounds(ds.rowStatus(i),
                             theURbound[i], theLRbound[i]);
            setLeaveBound4Row(i, i);
            computeLeaveCoPrhs4Row(i, i);
            // init #theFrhs[i]#
            Real& v_rhs = (*theFrhs)[i];
            const SVector& row = rowVector(i); 
            v_rhs = 0;
            for (int j = row.size() - 1; j >= 0; --j)
            {
               int idx = row.index(j);
               switch (ds.colStatus(idx))
               {
               case Desc::P_ON_UPPER:
                  v_rhs += row.value(j) * SPxLP::upper(idx);
                  break;
               case Desc::P_ON_LOWER:
               case Desc::P_FIXED:
                  v_rhs += row.value(j) * SPxLP::lower(idx);
                  break;
               default:
                  break;
               }
            }
         }
         SPxBasis::solve (*theFvec, *theFrhs);
         SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
         for (i = start; i < SPxLP::nRows(); ++i)
         {
            if ((*theFvec)[i] > theUBbound[i])
               theCoTest[i] = theUBbound[i] - (*theFvec)[i];
            else
               theCoTest[i] = (*theFvec)[i] - theLBbound[i];
         }
      }
   }
}

void SPxSolver::addedRows(int n)
{
   METHOD( "SPxSolver::addedRows()" );

   SPxLP::addedRows(n);

   if( n > 0 )
   {
      reDim();
      
      if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
      {
         SPxBasis::addedRows(n);

         if (isInitialized())
         {
            localAddRows(nRows() - n);

            assert(thepricer != 0);

            if (rep() == ROW)
               thepricer->addedVecs(n);
            else
               thepricer->addedCoVecs(n);            
         }
      }
   }

   /* we must not assert consistency here, since addedCols() might be still necessary to obtain a consistent basis */
}
#endif //0

void SPxSolver::addedRows(int n)
{
   METHOD( "SPxSolver::addedRows()" );

   if( n > 0 )
   {
      SPxLP::addedRows(n);

      unInit();
      reDim();

      if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
         SPxBasis::addedRows(n);
   }

   /* we must not assert consistency here, since addedCols() might be still necessary to obtain a consistent basis */
}

#if 0
void SPxSolver::localAddCols(int start)
{
   METHOD( "SPxSolver::localAddCols()" );
   assert( start <= SPxLP::nCols() );

   /**@todo This method seems to be called, to update
    *       theFvec, theFrhs, ..., but a resolve after
    *       adding a row results in a failure.
    *       To fix this, we call unIinit() so that init() is called before solving
    *       in spxsolve.cpp:solve(). In init(), the
    *       vectors are set up, so there is no need
    *       to update them here.
    */
   if( start == SPxLP::nCols() )
      return;

   const SPxBasis::Desc& ds = desc();

   if (type() == ENTER)
   {
      if (rep() == COLUMN)
      {
         int reSolve = 0;
         int i;
         Real x;
         for (i = start; i < SPxLP::nCols(); ++i)
         {
            (*thePvec)[i] = vector(i) * (*theCoPvec);
            theTest[i] = test(i, ds.colStatus(i));
            theUCbound[i] = SPxLP::upper(i);
            theLCbound[i] = SPxLP::lower(i);
            switch (ds.colStatus(i))
            {
            case SPxBasis::Desc::P_ON_LOWER + SPxBasis::Desc::P_ON_UPPER :
               assert(SPxLP::lower(i) == SPxLP::upper(i));
               /*FALLTHROUGH*/
            case SPxBasis::Desc::P_ON_UPPER :
               x = SPxLP::upper(i);
               break;
            case SPxBasis::Desc::P_ON_LOWER :
               x = SPxLP::lower(i);
               break;
            default:
               x = 0;
               break;
            }
            if (x)
            {
               theFrhs->multAdd(-x, vector(i));
               reSolve = 1;
            }
         }
         if (reSolve)
         {
            SPxBasis::solve(*theFvec, *theFrhs);
            shiftFvec();
         }
      }
      else
      {
         int i;
         for (i = start; i < SPxLP::nCols(); ++i)
         {
            theUCbound[i] = theLCbound[i] = 0;
            (*theFrhs)[i] = SPxLP::spxSense() * SPxLP::obj(i);
            clearDualBounds(ds.colStatus(i),
                             theUCbound[i], theLCbound[i]);
            setEnterBound4Col(i, i);
            computeEnterCoPrhs4Col(i, i);
         }
         SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
         computePvec();
         computeCoTest();
         computeTest();
         SPxBasis::solve(*theFvec, *theFrhs);
         for (i = start; i < SPxLP::nCols(); ++i)
         {
            /* we need to compare with tolerance (rep() == COLUMN) ? feastol() : opttol() because theFvec is the primal
             * vector in COLUMN and the dual vector in ROW representation; this is equivalent to entertol(); this also
             * fits because we are within the "type() == ENTER" case
             */
            if (theUBbound[i] + entertol() < (*theFvec)[i])
               shiftUBbound(i, (*theFvec)[i]);
            if ((*theFvec)[i] < theLBbound[i] - entertol())
               shiftLBbound(i, (*theFvec)[i]);
         }
      }
   }
   else
   {
      if (rep() == ROW)
      {
         int i;
         for (i = start; i < SPxLP::nCols(); ++i)
         {
            theUCbound[i] = SPxLP::upper(i);
            theLCbound[i] = SPxLP::lower(i);
            (*theFrhs)[i] = SPxLP::spxSense() * SPxLP::obj(i);
            setLeaveBound4Col(i, i);
            computeLeaveCoPrhs4Col(i, i);
         }
         SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
         computePvec();
         //          shiftPvec();
         SPxBasis::solve(*theFvec, *theFrhs);
         for (i = start; i < SPxLP::nCols(); ++i)
         {
            if ((*theFvec)[i] > theUBbound[i])
               theCoTest[i] = theUBbound[i] - (*theFvec)[i];
            else
               theCoTest[i] = (*theFvec)[i] - theLBbound[i];
         }
      }
      else
      {
         Real x;
         int i;
         int reSolve = 0;
         for (i = start; i < SPxLP::nCols(); ++i)
         {
            theUCbound[i] = theLCbound[i] = -maxObj(i);
            clearDualBounds(ds.colStatus(i),
                             theLCbound[i], theUCbound[i]);
            theUCbound[i] *= -1;
            theLCbound[i] *= -1;

            (*thePvec)[i] = vector(i) * (*theCoPvec);

            /* we need to compare with tolerance (rep() == ROW) ? feastol() : opttol() because thePvec is the primal
             * vector in ROW and the dual vector in COLUMN representation; this is equivalent to leavetol(); this also
             * fits because we are within the "type() == LEAVE" case
             */
            if (theUCbound[i] + leavetol() < (*thePvec)[i])
               shiftUPbound(i, (*thePvec)[i]);
            if (theLCbound[i] - leavetol() > (*thePvec)[i])
               shiftLPbound(i, (*thePvec)[i]);

            switch (ds.colStatus(i))
            {
            case SPxBasis::Desc::P_ON_LOWER + SPxBasis::Desc::P_ON_UPPER :
               assert(SPxLP::lower(i) == SPxLP::upper(i));
               /*FALLTHROUGH*/
            case SPxBasis::Desc::P_ON_UPPER :
               x = SPxLP::upper(i);
               break;
            case SPxBasis::Desc::P_ON_LOWER :
               x = SPxLP::lower(i);
               break;
            default:
               x = 0;
               break;
            }
            if (x)
            {
               theFrhs->multAdd(-x, vector(i));
               reSolve = 1;
            }
         }
         if (reSolve)
         {
            SPxBasis::solve(*theFvec, *theFrhs);
            computeFtest();
         }
      }
   }
}

void SPxSolver::addedCols(int n)
{
   METHOD( "SPxSolver::addedCols()" );
   SPxLP::addedCols(n);

   if( n > 0 )
   {
      reDim();
      
      if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
      {
         SPxBasis::addedCols(n);
         if (isInitialized())
         {
            localAddCols(nCols() - n);
            assert(thepricer != 0);
            if (rep() == COLUMN)
               thepricer->addedVecs(n);
            else
               thepricer->addedCoVecs(n);
         }
      }
   }

   /* we must not assert consistency here, since addedRows() might be still necessary to obtain a consistent basis */
}
#endif //0

void SPxSolver::addedCols(int n)
{
   METHOD( "SPxSolver::addedCols()" );

   if( n > 0 )
   {
      SPxLP::addedCols(n);

      unInit();
      reDim();

      if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
         SPxBasis::addedCols(n);
   }

   /* we must not assert consistency here, since addedRows() might be still necessary to obtain a consistent basis */
}
   
void SPxSolver::doRemoveRow(int i)
{
   METHOD( "SPxSolver::doRemoveRow()" );

   SPxLP::doRemoveRow(i);

   unInit();

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      removedRow(i);

#if 0
      if (isInitialized())
      {
         int n = SPxLP::nRows();

         theURbound[i] = theURbound[n];
         theLRbound[i] = theLRbound[n];

         if (rep() == ROW)
         {
            (*thePvec)[i] = (*thePvec)[n];
            if (type() == ENTER)
               theTest[i] = theTest[n];
            reDim();
            assert(thepricer != 0);
            thepricer->removedVec(i);
         }
         else
         {
            unInit();
         }
      }
#endif // 0

      switch (SPxBasis::status())
      {
      case SPxBasis::DUAL:
      case SPxBasis::INFEASIBLE:
         setBasisStatus(SPxBasis::REGULAR);
         break;
      case SPxBasis::OPTIMAL:
         setBasisStatus(SPxBasis::PRIMAL);
         break;
      default:
         break;
      }
   }
}

void SPxSolver::doRemoveRows(int perm[])
{
   METHOD( "SPxSolver::doRemoveRows()" );

   SPxLP::doRemoveRows(perm);

   unInit();

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      removedRows(perm);
#if 0
      if (isInitialized())
      {
         int n = SPxLP::nRows();

         if (rep() == ROW)
         {
            if (type() == ENTER)
            {
               for (int i = 0; i < n; ++i)
                  if (perm[i] >= 0)
                  {
                     theURbound[perm[i]] = theURbound[i];
                     theLRbound[perm[i]] = theLRbound[i];
                     (*thePvec)[perm[i]] = (*thePvec)[i];
                     theTest[perm[i]] = theTest[i];
                  }
            }
            else
            {
               for (int i = 0; i < n; ++i)
                  if (perm[i] >= 0)
                  {
                     theURbound[perm[i]] = theURbound[i];
                     theLRbound[perm[i]] = theLRbound[i];
                     (*thePvec)[perm[i]] = (*thePvec)[i];
                  }
            }
            assert(thepricer != 0);
            thepricer->removedVecs(perm);
            reDim();
         }
         else
         {
            unInit();
         }
      }
#endif
      switch (SPxBasis::status())
      {
      case SPxBasis::DUAL:
      case SPxBasis::INFEASIBLE:
         setBasisStatus(SPxBasis::REGULAR);
         break;
      case SPxBasis::OPTIMAL:
         setBasisStatus(SPxBasis::PRIMAL);
         break;
      default:
         break;
      }
   }
}

void SPxSolver::doRemoveCol(int i)
{
   METHOD( "SPxSolver::doRemoveCol()" );

   SPxLP::doRemoveCol(i);

   unInit();

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      removedCol(i);

#if 0
      if (isInitialized())
      {
         int n = SPxLP::nCols();

         theUCbound[i] = theUCbound[n];
         theLCbound[i] = theLCbound[n];
         if (rep() == COLUMN)
         {
            (*thePvec)[i] = (*thePvec)[n];
            if (type() == ENTER)
               theTest[i] = theTest[n];
            assert(thepricer != 0);
            thepricer->removedVec(i);
            reDim();
         }
         else
         {
            unInit();
         }
      }
#endif //0
      switch (SPxBasis::status())
      {
      case SPxBasis::PRIMAL:
      case SPxBasis::UNBOUNDED:
         setBasisStatus(SPxBasis::REGULAR);
         break;
      case SPxBasis::OPTIMAL:
         setBasisStatus(SPxBasis::DUAL);
         break;
      default:
         break;
      }
   }
}

void SPxSolver::doRemoveCols(int perm[])
{
   METHOD( "SPxSolver::doRemoveCols()" );

   SPxLP::doRemoveCols(perm);

   unInit();

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      removedCols(perm);

#if 0
      if (isInitialized())
      {
         int n = SPxLP::nCols();

         if (rep() == COLUMN)
         {
            if (type() == ENTER)
            {
               for (int i = 0; i < n; ++i)
                  if (perm[i] >= 0)
                  {
                     theUCbound[perm[i]] = theUCbound[i];
                     theLCbound[perm[i]] = theLCbound[i];
                     (*thePvec)[perm[i]] = (*thePvec)[i];
                     theTest[perm[i]] = theTest[i];
                  }
            }
            else
            {
               for (int i = 0; i < n; ++i)
                  if (perm[i] >= 0)
                  {
                     theUCbound[perm[i]] = theUCbound[i];
                     theLCbound[perm[i]] = theLCbound[i];
                     (*thePvec)[perm[i]] = (*thePvec)[i];
                  }
            }
            assert(thepricer != 0);
            thepricer->removedVecs(perm);
            reDim();
         }
         else
         {
            unInit();
         }
      }
#endif //0
      switch (SPxBasis::status())
      {
      case SPxBasis::PRIMAL:
      case SPxBasis::UNBOUNDED:
         setBasisStatus(SPxBasis::REGULAR);
         break;
      case SPxBasis::OPTIMAL:
         setBasisStatus(SPxBasis::DUAL);
         break;
      default:
         break;
      }
   }
}

void SPxSolver::changeObj(const Vector& newObj)
{
   METHOD( "SPxSolver::changeObj()" );

   SPxLP::changeObj(newObj);

   /**@todo Factorization remains valid, we do not need a reDim()
    *       pricing vectors should be recomputed.
    */
   unInit();
}

void SPxSolver::changeObj(int i, Real newVal)
{
   METHOD( "SPxSolver::changeObj()" );

   SPxLP::changeObj(i, newVal);

   /**@todo Factorization remains valid, we do not need a reDim()
    *       pricing vectors should be recomputed.
    */
   unInit();
}

static void changeLowerStatus(
   SPxBasis::Desc::Status& stat,
   Real                    newLower,
   Real                    upper,
   const SPxBasis&         basis,
   int                     i)
{
   MSG_DEBUG( spxout << "DCHANG01 changeLowerStatus(): col " << i
                     << "[" << newLower << ":" << upper << "] " << stat; )

   switch (stat)
   {
   case SPxBasis::Desc::P_ON_LOWER:
      if (newLower <= -infinity)
         stat = (upper >= infinity) ? SPxBasis::Desc::P_FREE : SPxBasis::Desc::P_ON_UPPER;
      else if (newLower == upper)
         stat = SPxBasis::Desc::P_FIXED;
      break;
   case SPxBasis::Desc::P_ON_UPPER:
      if (newLower == upper)
         stat = SPxBasis::Desc::P_FIXED;
      break;
   case SPxBasis::Desc::P_FREE:
      if (newLower > -infinity)
         stat = SPxBasis::Desc::P_ON_LOWER;
      break;
   case SPxBasis::Desc::P_FIXED:
      if (newLower != upper)
         stat = SPxBasis::Desc::P_ON_UPPER;
      break;
   case SPxBasis::Desc::D_FREE:
   case SPxBasis::Desc::D_ON_UPPER:
   case SPxBasis::Desc::D_ON_LOWER:
   case SPxBasis::Desc::D_ON_BOTH:
   case SPxBasis::Desc::D_UNDEFINED:
      stat = basis.dualColStatus(i);
      break;
   default:
      throw SPxInternalCodeException("XCHANG01 This should never happen.");
   }
   MSG_DEBUG( spxout << " -> " << stat << std::endl; )
}

void SPxSolver::changeLower(const Vector& newLower)
{
   METHOD( "SPxSolver::changeLower()" );

   SPxLP::changeLower(newLower);

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      for (int i = 0; i < newLower.dim(); ++i)
         changeLowerStatus(desc().colStatus(i), newLower[i], upper(i), *this, i);

      unInit();
   }
}

void SPxSolver::changeLower(int i, Real newLower)
{
   METHOD( "SPxSolver::changeLower()" );

   if (newLower != lower(i))
   {
      // This has to be done before calling changeLowerStatus() because that is calling
      // basis.dualColStatus() which calls lower() and needs the changed value.
      SPxLP::changeLower(i, newLower);

      if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
      {
         changeLowerStatus(desc().colStatus(i), newLower, upper(i), *this, i);
         unInit();
      }
   }
}

static void changeUpperStatus(
   SPxBasis::Desc::Status& stat,
   Real                    newUpper,
   Real                    lower,
   const SPxBasis&         basis,
   int                     i)
{
   MSG_DEBUG( spxout << "DCHANG02 changeUpperStatus(): col " << i
                     << "[" << lower << ":" << newUpper << "] " << stat; )

   switch (stat)
   {
   case SPxBasis::Desc::P_ON_LOWER:
      if (newUpper == lower)
         stat = SPxBasis::Desc::P_FIXED;
      break;
   case SPxBasis::Desc::P_ON_UPPER:
      if (newUpper >= infinity)
         stat = (lower <= -infinity) ? SPxBasis::Desc::P_FREE : SPxBasis::Desc::P_ON_LOWER;
      else if (newUpper == lower)
         stat = SPxBasis::Desc::P_FIXED;
      break;
   case SPxBasis::Desc::P_FREE:
      if (newUpper < infinity)
         stat = SPxBasis::Desc::P_ON_UPPER;
      break;
   case SPxBasis::Desc::P_FIXED:
      if (newUpper != lower)
         stat = SPxBasis::Desc::P_ON_LOWER;
      break;
   case SPxBasis::Desc::D_FREE:
   case SPxBasis::Desc::D_ON_UPPER:
   case SPxBasis::Desc::D_ON_LOWER:
   case SPxBasis::Desc::D_ON_BOTH:
   case SPxBasis::Desc::D_UNDEFINED:
      stat = basis.dualColStatus(i);
      break;
   default:
      throw SPxInternalCodeException("XCHANG02 This should never happen.");
   }
   MSG_DEBUG( spxout << " -> " << stat << std::endl; );
}

void SPxSolver::changeUpper(const Vector& newUpper)
{
   METHOD( "SPxSolver::changeUpper()" );

   SPxLP::changeUpper(newUpper);

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      for (int i = 0; i < newUpper.dim(); ++i)
         changeUpperStatus(desc().colStatus(i), newUpper[i], lower(i), *this, i);

      unInit();
   }
}

void SPxSolver::changeUpper(int i, Real newUpper)
{
   METHOD( "SPxSolver::changeUpper()" );

   if (newUpper != upper(i))
   {
      SPxLP::changeUpper(i, newUpper);

      if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
      {
         changeUpperStatus(desc().colStatus(i), newUpper, lower(i), *this, i);
         unInit();
      }
   }
}

void SPxSolver::changeBounds(const Vector& newLower, const Vector& newUpper)
{
   METHOD( "SPxSolver::changeBounds()" );

   changeLower(newLower);
   changeUpper(newUpper);
}

void SPxSolver::changeBounds(int i, Real newLower, Real newUpper)
{
   METHOD( "SPxSolver::changeBounds()" );

   changeLower(i, newLower);
   changeUpper(i, newUpper);
}

/**@todo Change Lhs/Rhs Status the same way as changeBounds
 */
static void changeLhsStatus(
   SPxBasis::Desc::Status& stat,
   Real newLhs,
   Real rhs,
   const SPxBasis& basis,
   int i)
{
   MSG_DEBUG( spxout << "DCHANG03 changeLhsStatus()  : row " << i
                     << ": " << stat; )
   switch (stat)
   {
   case SPxBasis::Desc::P_ON_LOWER:
      if (newLhs <= -infinity)
         stat = (rhs >= infinity) ? SPxBasis::Desc::P_FREE : SPxBasis::Desc::P_ON_UPPER;
      else if (newLhs == rhs)
         stat = SPxBasis::Desc::P_FIXED;
      break;
   case SPxBasis::Desc::P_ON_UPPER:
      if (newLhs == rhs)
         stat = SPxBasis::Desc::P_FIXED;
      break;
   case SPxBasis::Desc::P_FREE:
      if (newLhs > -infinity)
         stat = SPxBasis::Desc::P_ON_LOWER;
      break;
   case SPxBasis::Desc::P_FIXED:
      if (newLhs != rhs)
         stat = SPxBasis::Desc::P_ON_UPPER;
      break;
   case SPxBasis::Desc::D_FREE:
   case SPxBasis::Desc::D_ON_UPPER:
   case SPxBasis::Desc::D_ON_LOWER:
   case SPxBasis::Desc::D_ON_BOTH:
   case SPxBasis::Desc::D_UNDEFINED:
      stat = basis.dualRowStatus(i);
      break;
   default:
      throw SPxInternalCodeException("XCHANG03 This should never happen.");
   }
   MSG_DEBUG( spxout << " -> " << stat << std::endl; )
}

void SPxSolver::changeLhs(const Vector& newLhs)
{
   METHOD( "SPxSolver::changeLhs()" );

   SPxLP::changeLhs(newLhs);

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      for (int i = 0; i < nRows(); ++i)
         changeLhsStatus(desc().rowStatus(i), newLhs[i], rhs(i), *this, i);

      unInit();
   }
}

void SPxSolver::changeLhs(int i, Real newLhs)
{
   METHOD( "SPxSolver::changeLhs()" );

   SPxLP::changeLhs(i, newLhs);

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      changeLhsStatus(desc().rowStatus(i), newLhs, rhs(i), *this, i);
      unInit();
   }
}

static void changeRhsStatus(
   SPxBasis::Desc::Status& stat,
   Real newRhs,
   Real lhs,
   const SPxBasis& basis,
   int i)
{
   MSG_DEBUG( spxout << "DCHANG04 changeRhsStatus()  : row " << i
                     << ": " << stat; )
   switch (stat)
   {
   case SPxBasis::Desc::P_ON_UPPER:
      if (newRhs >= infinity)
         stat = (lhs <= -infinity) ? SPxBasis::Desc::P_FREE : SPxBasis::Desc::P_ON_LOWER;
      else if (newRhs == lhs)
         stat = SPxBasis::Desc::P_FIXED;
      break;
   case SPxBasis::Desc::P_ON_LOWER:
      if (newRhs == lhs)
         stat = SPxBasis::Desc::P_FIXED;
      break;
   case SPxBasis::Desc::P_FREE:
      if (newRhs < infinity)
         stat = SPxBasis::Desc::P_ON_UPPER;
      break;
   case SPxBasis::Desc::P_FIXED:
      if (newRhs != lhs)
         stat = SPxBasis::Desc::P_ON_LOWER;
      break;
   case SPxBasis::Desc::D_FREE:
   case SPxBasis::Desc::D_ON_UPPER:
   case SPxBasis::Desc::D_ON_LOWER:
   case SPxBasis::Desc::D_ON_BOTH:
   case SPxBasis::Desc::D_UNDEFINED:
      stat = basis.dualRowStatus(i);
      break;
   default:
      throw SPxInternalCodeException("XCHANG04 This should never happen.");
   }
   MSG_DEBUG( spxout << " -> " << stat << std::endl; )
}


void SPxSolver::changeRhs(const Vector& newRhs)
{
   METHOD( "SPxSolver::changeRhs()" );

   SPxLP::changeRhs(newRhs);

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      for (int i = 0; i < nRows(); ++i)
         changeRhsStatus(desc().rowStatus(i), newRhs[i], lhs(i), *this, i);
      unInit();
   }
}

void SPxSolver::changeRhs(int i, Real newRhs)
{
   METHOD( "SPxSolver::changeRhs()" );

   SPxLP::changeRhs(i, newRhs);

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      changeRhsStatus(desc().rowStatus(i), newRhs, lhs(i), *this, i);
      unInit();
   }
}

void SPxSolver::changeRange(const Vector& newLhs, const Vector& newRhs)
{
   METHOD( "SPxSolver::changeRange()" );
   SPxLP::changeLhs(newLhs);
   SPxLP::changeRhs(newRhs);
   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      for (int i = nRows() - 1; i >= 0; --i)
      {
         changeLhsStatus(desc().rowStatus(i), newLhs[i], rhs(i), *this, i);
         changeRhsStatus(desc().rowStatus(i), newRhs[i], lhs(i), *this, i);
      }
      unInit();
   }
}

void SPxSolver::changeRange(int i, Real newLhs, Real newRhs)
{
   METHOD( "SPxSolver::changeRange()" );

   SPxLP::changeLhs(i, newLhs);
   SPxLP::changeRhs(i, newRhs);

   if (SPxBasis::status() > SPxBasis::NO_PROBLEM)
   {
      changeLhsStatus(desc().rowStatus(i), newLhs, rhs(i), *this, i);
      changeRhsStatus(desc().rowStatus(i), newRhs, lhs(i), *this, i);
      unInit();
   }
}

void SPxSolver::changeRow(int i, const LPRow& newRow)
{
   METHOD( "SPxSolver::changeRow()" );

   SPxLP::changeRow(i, newRow);
   if ( SPxBasis::status() > SPxBasis::NO_PROBLEM )
      SPxBasis::changedRow( i );
   unInit();
}

void SPxSolver::changeCol(int i, const LPCol& newCol)
{
   METHOD( "SPxSolver::changeCol()" );

   SPxLP::changeCol(i, newCol);
   if ( SPxBasis::status() > SPxBasis::NO_PROBLEM )
      SPxBasis::changedCol( i );
   unInit();
}

void SPxSolver::changeElement(int i, int j, Real val)
{
   METHOD( "SPxSolver::changeElement()" );

   SPxLP::changeElement(i, j, val);
   if ( SPxBasis::status() > SPxBasis::NO_PROBLEM )
      SPxBasis::changedElement( i, j );
   unInit();
}

void SPxSolver::changeSense(SPxSense sns)
{
   METHOD( "SPxSolver::changeSense()" );

   SPxLP::changeSense(sns);
   unInit();
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
