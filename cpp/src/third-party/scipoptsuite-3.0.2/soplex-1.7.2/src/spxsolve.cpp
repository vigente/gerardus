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

//#define DEBUGGING

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "mpqreal.h"
#include "spxsolver.h"
#include "spxpricer.h"
#include "spxratiotester.h"
#include "spxdefaultrt.h"
#include "spxstarter.h"
#include "spxout.h"
#include "exceptions.h"

#define MAXCYCLES 400
#define MAXSTALLS 10000
#define MAXSTALLRECOVERS 10
#define MAXREFACPIVOTS 10

namespace soplex
{

/// Interval for displaying iteration information.
long iterationInterval = 100;

/**@todo check separately for ENTER and LEAVE algorithm */
bool SPxSolver::precisionReached(Real& newpricertol) const
{
   Real maxViolRedCost;
   Real sumViolRedCost;
   Real maxViolBounds;
   Real sumViolBounds;
   Real maxViolConst;
   Real sumViolConst;

   qualRedCostViolation(maxViolRedCost, sumViolRedCost);
   qualBoundViolation(maxViolBounds, sumViolBounds);
   qualConstraintViolation(maxViolConst, sumViolConst);

   // is the solution good enough ?
   bool reached = maxViolRedCost < opttol() && maxViolBounds < feastol() && maxViolConst < feastol();

   if (!reached)
   {
      newpricertol = thepricer->epsilon() / 10.0;

      MSG_INFO3( spxout << "ISOLVE71 "
                           << "Precision not reached: Pricer tolerance = "
                           << thepricer->epsilon()
                           << " new tolerance = " << newpricertol
                           << std::endl
                           << " maxViolRedCost= " << maxViolRedCost
                           << " maxViolBounds= " << maxViolBounds
                           << " maxViolConst= " << maxViolConst
                           << std::endl
                           << " sumViolRedCost= " << sumViolRedCost
                           << " sumViolBounds= " << sumViolBounds
                           << " sumViolConst= " << sumViolConst
                           << std::endl; );
   }
   return reached;
}

SPxSolver::Status SPxSolver::fpsolve()
{
   METHOD( "SPxSolver::fpsolve()" );

   SPxId enterId;
   int   leaveNum;
   int   loopCount = 0;
   Real  minShift = infinity;
   int   cycleCount = 0;
   bool  priced = false;
   Real  lastDelta = 1;

   /* store the last (primal or dual) feasible objective value to recover/abort in case of stalling */
   Real  stallRefValue;
   Real  stallRefShift;
   int   stallRefIter;
   int   stallNumRecovers;

   if (dim() <= 0 && coDim() <= 0) // no problem loaded
   {
      m_status = NO_PROBLEM;
      throw SPxStatusException("XSOLVE01 No Problem loaded");
   }

   if (slinSolver() == 0) // linear system solver is required.
   {
      m_status = NO_SOLVER;
      throw SPxStatusException("XSOLVE02 No Solver loaded");
   }
   if (thepricer == 0) // pricer is required.
   {
      m_status = NO_PRICER;
      throw SPxStatusException("XSOLVE03 No Pricer loaded");
   }
   if (theratiotester == 0) // ratiotester is required.
   {
      m_status = NO_RATIOTESTER;
      throw SPxStatusException("XSOLVE04 No RatioTester loaded");
   }

   m_numCycle = 0;
   iterCount  = 0;
   if (!isInitialized())
   {
      /*
      if(SPxBasis::status() <= NO_PROBLEM)
          SPxBasis::load(this);
       */
      /**@todo != REGULAR is not enough. Also OPTIMAL/DUAL/PRIMAL should
       * be tested and acted accordingly.
       */
      if (thestarter != 0 && status() != REGULAR)  // no basis and no starter.
         thestarter->generate(*this);              // generate start basis.

      init();

      // Inna/Tobi: init might fail, if the basis is singular
      if( !isInitialized() )
      {
         assert(SPxBasis::status() == SPxBasis::SINGULAR);
         m_status = UNKNOWN;
         return status();
      }
   }

   //setType(type());

   if (!matrixIsSetup)
      SPxBasis::load(this);

   //factorized = false;

   assert(thepricer->solver()      == this);
   assert(theratiotester->solver() == this);

   // maybe this should be done in init() ?
   thepricer->setType(type());
   theratiotester->setType(type());

   MSG_INFO3(
      spxout << "ISOLVE72 starting value = " << value() << std::endl;
      spxout << "ISOLVE73 starting shift = " << shift() << std::endl; 
   )

   if (SPxBasis::status() == SPxBasis::OPTIMAL)
      setBasisStatus(SPxBasis::REGULAR);

   m_status   = RUNNING;
   bool stop  = terminate();
   leaveCount = 0;
   enterCount = 0;
   boundflips = 0;
   totalboundflips = 0;

   stallNumRecovers = 0;

   /* if we run into a singular basis, we will retry from regulardesc with tighter tolerance in the ratio test */
   SPxSolver::Type tightenedtype = type();
   bool tightened = false;

   while (!stop)
   {
      const SPxBasis::Desc regulardesc = desc();

      // we need to reset these pointers to avoid unnecessary/wrong solves in leave() or enter()
      solveVector2 = 0;
      solveVector3 = 0;
      coSolveVector2 = 0;

      try
      {

      if (type() == ENTER)
      {
         int enterCycleCount = 0;
         int enterFacPivotCount = 0;

         stallRefIter = iteration()-1;
         stallRefShift = shift();
         stallRefValue = value();

         /* in the entering algorithm, entertol() should be maintained by the ratio test and leavetol() should be
          * reached by the pricer
          */
         Real maxpricertol = leavetol();
         Real minpricertol = 0.01 * maxpricertol;

         thepricer->setEpsilon(maxpricertol);
         priced = false;

         // to avoid shifts we restrict tolerances in the ratio test
         if( loopCount > 0 )
         {
            lastDelta = (lastDelta < entertol()) ? lastDelta : entertol();
            lastDelta *= 0.01;
            theratiotester->setDelta(lastDelta);
            assert(theratiotester->getDelta() > 0);
            MSG_DEBUG( spxout << "decreased delta for ratiotest to: " << theratiotester->getDelta() << std::endl; )
         }
         else
         {
            lastDelta = 1;
            theratiotester->setDelta(entertol());
         }

         do
         {
            MSG_INFO3(
               if( iteration() % iterationInterval == 0 )
                  spxout << "ISOLVE74 Enter iteration: " << iteration()
                         << ", Value = " << value()
                         << ", Shift = " << shift() << std::endl;
            )
            enterId = thepricer->selectEnter();

            if (!enterId.isValid())
            {
               // we are not infeasible and have no shift
               if (  shift() <= epsilon()
                  && ( SPxBasis::status() == SPxBasis::REGULAR 
                     || SPxBasis::status() == SPxBasis::DUAL 
                     || SPxBasis::status() == SPxBasis::PRIMAL))
               {
                  Real newpricertol = minpricertol;

                  // is the solution good enough ?
                  // max three times reduced
                  if ((thepricer->epsilon() > minpricertol) && !precisionReached(newpricertol))
                  {  // no!
                     // we reduce the pricer tolerance. Note that if the pricer does not find a candiate
                     // with the reduced tolerance, we quit, regardless of the violations.
                     if (newpricertol < minpricertol)
                        newpricertol = minpricertol;

                     thepricer->setEpsilon(newpricertol);

                     MSG_INFO2( spxout << "ISOLVE75 Setting pricer tolerance = "
                                          << thepricer->epsilon()
                                          << std::endl; )
                  }
                  // solution seems good, no check whether we are precise enough
                  else if (lastUpdate() == 0)
                  {
                     priced = true;
                     break;
                  }
                  // We have an iterationlimit and everything looks good? Then stop!
                  // 6 is just a number picked.
                  else if (maxIters > 0 && lastUpdate() < 6)
                  {
                     priced = true;
                     break;
                  }
               }
               MSG_INFO3( spxout << "ISOLVE76 solve(enter) triggers refactorization" << std::endl; )

               // if the factorization is not fresh, we better refactorize and call the pricer again; however, this can
               // create cycling, so it is performed only a limited number of times per ENTER round
               if( lastUpdate() > 0 && enterFacPivotCount < MAXREFACPIVOTS )
               {
                  factorize();

                  // if the factorization was found out to be singular, we have to quit
                  if( SPxBasis::status() < SPxBasis::REGULAR )
                  {
                     MSG_ERROR( spxout << "ESOLVE09 something wrong with factorization, Basis status: " << SPxBasis::status() << std::endl; )
                     stop = true;
                     break;
                  }

                  // call pricer again
                  enterId = thepricer->selectEnter();

                  // count how often the pricer has found something only after refactorizing
                  if( enterId.isValid() )
                     enterFacPivotCount++;
               }

               if( !enterId.isValid() )
               {
                  priced = true;
                  break;
               }
            }

            /* check if we have iterations left */
            if (maxIters >= 0 && iterations() >= maxIters)
            {
               MSG_INFO2( spxout << "ISOLVE53e Maximum number of iterations (" << maxIters
                                 << ") reached" << std::endl; )
               m_status = ABORT_ITER;
               stop = true;
               break;
            }

            enter(enterId);
            assert((testBounds(), 1));
            thepricer->entered4(lastEntered(), lastIndex());
            stop = terminate();
            clearUpdateVecs();
            if( lastEntered().isValid() ) /* either a successful pivot was performed or a nonbasic variable flipped to the other bound */
            {
               enterCount++;
               enterCycleCount = 0;
            }
            else
            {
               enterCycleCount++;
               if( enterCycleCount > MAXCYCLES )
               {
                  MSG_INFO2( spxout << "ISOLVE77 Abort solving due to cycling in "
                                       << "entering algorithm" << std::endl; );
                  m_status = ABORT_CYCLING;
                  stop = true;
               }
            }

            /* check every MAXSTALLS iterations whether shift and objective value have not changed */
            if( (iteration() - stallRefIter) % MAXSTALLS == 0 )
            {
               if( fabs(value() - stallRefValue) <= epsilon() && fabs(shift() - stallRefShift) <= epsilon() )
               {
                  if( stallNumRecovers < MAXSTALLRECOVERS )
                  {
                     /* try to recover by unshifting/switching algorithm up to MAXSTALLRECOVERS times (just a number picked) */
                     MSG_INFO3( spxout << "ISOLVE21 Stalling detected - trying to recover by switching to LEAVING algorithm." << std::endl; )

                     ++stallNumRecovers;
                     break;
                  }
                  else
                  {
                     /* giving up */
                     MSG_INFO2( spxout << "ISOLVE22 Abort solving due to stalling in entering algorithm." << std::endl; );

                     m_status = ABORT_CYCLING;
                     stop = true;
                  }
               }
               else
               {
                  /* merely update reference values */
                  stallRefIter = iteration()-1;
                  stallRefShift = shift();
                  stallRefValue = value();
               }
            }

            //@ assert(isConsistent());
         }
         while (!stop);

         MSG_INFO3(
            spxout << "ISOLVE78 Enter finished. iteration: " << iteration() 
                   << ", value: " << value()
                   << ", shift: " << shift()
                   << ", epsilon: " << epsilon()
                   << ", feastol: " << feastol()
                   << ", opttol: " << opttol()
                   << std::endl
                   << "ISOLVE56 stop: " << stop
                   << ", basis status: " << SPxBasis::status() << " (" << int(SPxBasis::status()) << ")"
                   << ", solver status: " << m_status << " (" << int(m_status) << ")" << std::endl;
         )

         if (!stop)
         {
            /**@todo technically it would be ok to finish already when (priced && maxinfeas + shift() <= entertol()) is
             *  satisfied; maybe at least in the case when SoPlex keeps jumping back between ENTER and LEAVE always
             *  shifting (looping), we may relax this condition here;
             *  note also that unShift may increase shift() slightly due to roundoff errors
             */
            if (shift() <= epsilon())
            {
               // factorize();
               unShift();

               Real maxinfeas = maxInfeas();

               MSG_INFO3(
                  spxout << "ISOLVE79 maxInfeas: " << maxinfeas
                         << ", shift: " << shift()
                         << ", entertol: " << entertol() << std::endl;
               )

               if (priced && maxinfeas + shift() <= entertol())
               {
                  setBasisStatus(SPxBasis::OPTIMAL);
                  m_status = OPTIMAL;
                  break;
               }
            }
            setType(LEAVE);
            init();
            thepricer->setType(type());
            theratiotester->setType(type());
         }
      }
      else
      {
         assert(type() == LEAVE);

         int leaveCycleCount = 0;
         int leaveFacPivotCount = 0;

         instableLeaveNum = -1;
         instableLeave = false;

         stallRefIter = iteration()-1;
         stallRefShift = shift();
         stallRefValue = value();

         /* in the leaving algorithm, leavetol() should be maintained by the ratio test and entertol() should be reached
          * by the pricer
          */
         Real maxpricertol = entertol();
         Real minpricertol = 0.01 * maxpricertol;

         thepricer->setEpsilon(maxpricertol);
         priced = false;

         // to avoid shifts we restrict tolerances in the ratio test
         if( loopCount > 0 )
         {
            lastDelta = (lastDelta < leavetol()) ? lastDelta : leavetol();
            lastDelta *= 0.01;
            theratiotester->setDelta(lastDelta);
            assert(theratiotester->getDelta() > 0);
            MSG_DEBUG( spxout << "decreased delta for ratiotest to: " << theratiotester->getDelta() << std::endl; )
         }
         else
         {
            lastDelta = 1;
            theratiotester->setDelta(leavetol());
         }

         do
         {
            MSG_INFO3(
               if( iteration() % iterationInterval == 0 )
                  spxout << "ISOLVE80 Leave Iteration: " << iteration()
                         << ", Value = " << value()
                         << ", Shift = " << shift() << std::endl;
            )
            
            leaveNum = thepricer->selectLeave();

            if (leaveNum < 0 && instableLeaveNum >= 0 && lastUpdate() == 0)
            {
               /* no leaving variable was found, but because of instableLeaveNum >= 0 we know
                  that this is due to the scaling of theCoTest[...]. Thus, we use 
                  instableLeaveNum and SPxFastRT::selectEnter shall accept even an instable
                  entering variable. */
               MSG_INFO3(
                  spxout << "ISOLVE98 Trying instable leave iteration" << std::endl;
               )
            
               leaveNum = instableLeaveNum;
               instableLeave = true;
            }
            else
            {
               instableLeave = false;
            }

            if (leaveNum < 0)
            {
               // we are not infeasible and have no shift
               if (  shift() <= epsilon()
                  && (  SPxBasis::status() == SPxBasis::REGULAR 
                     || SPxBasis::status() == SPxBasis::DUAL 
                     || SPxBasis::status() == SPxBasis::PRIMAL))
               {
                  Real newpricertol = minpricertol;

                  // is the solution good enough ?
                  // max three times reduced
                  if ((thepricer->epsilon() > minpricertol) && !precisionReached(newpricertol))
                  {  // no
                     // we reduce the pricer tolerance. Note that if the pricer does not find a candiate
                     // with the reduced pricer tolerance, we quit, regardless of the violations.
                     if (newpricertol < minpricertol)
                        newpricertol = minpricertol;

                     thepricer->setEpsilon(newpricertol);

                     MSG_INFO2( spxout << "ISOLVE81 Setting pricer tolerance = "
                                          << thepricer->epsilon()
                                          << std::endl; );
                  }
                  // solution seems good, no check whether we are precise enough
                  else if (lastUpdate() == 0)
                  {
                     priced = true;
                     break;
                  }
                  // We have an iteration limit and everything looks good? Then stop!
                  // 6 is just a number picked.
                  else if (maxIters > 0 && lastUpdate() < 6)
                  {
                     priced = true;
                     break;
                  }
               }
               MSG_INFO3( spxout << "ISOLVE82 solve(leave) triggers refactorization" << std::endl; )

               // if the factorization is not fresh, we better refactorize and call the pricer again; however, this can
               // create cycling, so it is performed only a limited number of times per LEAVE round
               if( lastUpdate() > 0 && leaveFacPivotCount < MAXREFACPIVOTS )
               {
                  factorize();

                  // Inna/Tobi: if the factorization was found out to be singular, we have to quit
                  if (SPxBasis::status() < SPxBasis::REGULAR)
                  {
                     MSG_ERROR( spxout << "ESOLVE10 something wrong with factorization, Basis status: " << SPxBasis::status() << std::endl; )
                     stop = true;
                     break;
                  }

                  // call pricer again
                  leaveNum = thepricer->selectLeave();

                  // count how often the pricer has found something only after refactorizing
                  if( leaveNum >= 0 )
                     leaveFacPivotCount++;
               }

               if (leaveNum < 0)
               {
                  priced = true;
                  break;
               }
            }

            /* check if we have iterations left */
            if (maxIters >= 0 && iterations() >= maxIters)
            {
               MSG_INFO2( spxout << "ISOLVE53l Maximum number of iterations (" << maxIters
                                 << ") reached" << std::endl; )
               m_status = ABORT_ITER;
               stop = true;
               break;
            }

            leave(leaveNum);
            assert((testBounds(), 1));
            thepricer->left4(lastIndex(), lastLeft());
            stop = terminate();
            clearUpdateVecs();
            if( lastIndex() >= 0 ) /* either a successful pivot was performed or a nonbasic variable flipped to the other bound */
            {
               leaveCount++;
               leaveCycleCount = 0;
            }
            else
            {
               leaveCycleCount++;
               if( leaveCycleCount > MAXCYCLES )
               {
                  MSG_INFO2( spxout << "ISOLVE83 Abort solving due to cycling in leaving algorithm" << std::endl; );
                  m_status = ABORT_CYCLING;
                  stop = true;
               }
            }

            /* check every MAXSTALLS iterations whether shift and objective value have not changed */
            if( (iteration() - stallRefIter) % MAXSTALLS == 0 )
            {
               if( fabs(value() - stallRefValue) <= epsilon() && fabs(shift() - stallRefShift) <= epsilon() )
               {
                  if( stallNumRecovers < MAXSTALLRECOVERS )
                  {
                     /* try to recover by switching algorithm up to MAXSTALLRECOVERS times */
                     MSG_INFO3( spxout << "ISOLVE24 Stalling detected - trying to recover by switching to ENTERING algorithm." << std::endl; )

                     ++stallNumRecovers;
                     break;
                  }
                  else
                  {
                     /* giving up */
                     MSG_INFO2( spxout << "ISOLVE25 Abort solving due to stalling in leaving algorithm" << std::endl; );

                     m_status = ABORT_CYCLING;
                     stop = true;
                  }
               }
               else
               {
                  /* merely update reference values */
                  stallRefIter = iteration()-1;
                  stallRefShift = shift();
                  stallRefValue = value();
               }
            }

            //@ assert(isConsistent());
         }
         while (!stop);

         MSG_INFO3(
            spxout << "ISOLVE84 Leave finished. iteration: " << iteration() 
                   << ", value: " << value()
                   << ", shift: " << shift()
                   << ", epsilon: " << epsilon()
                   << ", feastol: " << feastol()
                   << ", opttol: " << opttol()
                   << std::endl
                   << "ISOLVE57 stop: " << stop
                   << ", basis status: " << SPxBasis::status() << " (" << int(SPxBasis::status()) << ")"
                   << ", solver status: " << m_status << " (" << int(m_status) << ")" << std::endl;
         )

         if (!stop)
         {
            if( shift() < minShift )
            {
               minShift = shift();
               cycleCount = 0;
            }
            else
            {
               cycleCount++;
               if( cycleCount > MAXCYCLES )
               {
                  m_status = ABORT_CYCLING;
                  throw SPxStatusException("XSOLVE13 Abort solving due to cycling");
               }
               MSG_INFO3(
                  spxout << "ISOLVE86 maxInfeas: " << maxInfeas()
                         << ", shift: " << shift()
                         << ", leavetol: " << leavetol()
                         << ", cycle count: " << cycleCount << std::endl;
               )
            }

            /**@todo technically it would be ok to finish already when (priced && maxinfeas + shift() <= entertol()) is
             *  satisfied; maybe at least in the case when SoPlex keeps jumping back between ENTER and LEAVE always
             *  shifting (looping), we may relax this condition here;
             *  note also that unShift may increase shift() slightly due to roundoff errors
             */
            if (shift() <= epsilon())
            {
               cycleCount = 0;
               // factorize();
               unShift();

               Real maxinfeas = maxInfeas();

               MSG_INFO3(
                  spxout << "ISOLVE87 maxInfeas: " << maxinfeas
                         << ", shift: " << shift()
                         << ", leavetol: " << leavetol() << std::endl;
               )

               // We stop if we are indeed optimal, or if we have already been
               // two times at this place. In this case it seems futile to
               // continue.
               if (loopCount > 2)
               {
                  m_status = ABORT_CYCLING;
                  throw SPxStatusException("XSOLVE14 Abort solving due to looping");
               }
               else if (priced && maxinfeas + shift() <= leavetol())
               {
                  setBasisStatus(SPxBasis::OPTIMAL);
                  m_status = OPTIMAL;
                  break;
               }
               loopCount++;
            }
            setType(ENTER);
            init();
            thepricer->setType(type());
            theratiotester->setType(type());
         }
      }
      assert(m_status != SINGULAR);

      }
      catch( SPxException E )
      {
         // if we stopped due to a singular basis, we reload the original basis and try again with tighter
         // tolerance (only once)
         if (m_status == SINGULAR && !tightened)
         {
            tightenedtype = type();

            if( tightenedtype == ENTER )
            {
               m_entertol = 0.01 * m_entertol;

               MSG_INFO2( spxout << "ISOLVE26e basis singular: reloading basis and solving with tighter ratio test tolerance " << m_entertol << std::endl; )
            }
            else
            {
               m_leavetol = 0.01 * m_leavetol;

               MSG_INFO2( spxout << "ISOLVE26l basis singular: reloading basis and solving with tighter ratio test tolerance " << m_leavetol << std::endl; )
            }

            // load original basis
            int niters = iterations();
            loadBasis(regulardesc);

            // remember iteration count
            iterCount = niters;

            // try initializing basis (might fail if starting basis was already singular)
            try
            {
               init();
               theratiotester->setType(type());
            }
            catch( SPxException Ex )
            {
               MSG_INFO2( spxout << "ISOLVE27 reloaded basis singular, resetting original tolerances" << std::endl; )

               if( tightenedtype == ENTER )
                  m_entertol = 100.0 * m_entertol;
               else
                  m_leavetol = 100.0 * m_leavetol;

               theratiotester->setType(type());

               throw Ex;
            }

            // reset status and counters
            m_status = RUNNING;
            m_numCycle = 0;
            leaveCount = 0;
            enterCount = 0;
            stallNumRecovers = 0;

            // continue
            stop = false;
            tightened = true;
         }
         // reset tolerance to its original value and pass on the exception
         else if (tightened)
         {
            if( tightenedtype == ENTER )
               m_entertol = 100.0 * m_entertol;
            else
               m_leavetol = 100.0 * m_leavetol;

            theratiotester->setType(type());

            throw E;
         }
         // pass on the exception
         else
            throw E;
      }
   }

   // reset tolerance to its original value
   if (tightened)
   {
      if( tightenedtype == ENTER )
         m_entertol = 100.0 * m_entertol;
      else
         m_leavetol = 100.0 * m_leavetol;

      theratiotester->setType(type());
   }

   if (m_status == RUNNING)
   {
      m_status = ERROR;
      throw SPxStatusException("XSOLVE05 Status is still RUNNING when it shouldn't be");
   }

   MSG_INFO1(
      spxout << "ISOLVE02 Finished solving (status=" << status()
             << ", iters=" << iterCount
             << ", leave=" << leaveCount
             << ", enter=" << enterCount
             << ", flips=" << totalboundflips;
      if( status() == OPTIMAL )
         spxout << ", objValue=" << value();
      spxout << ")" << std::endl;
   )

#ifdef ENABLE_ADDITIONAL_CHECKS
   /* check if solution is really feasible */
   if( status() == OPTIMAL )
   {
      int     c;
      Real    val;
      DVector sol( nCols() );

      getPrimal( sol );

      for(int row = 0; row < nRows(); ++row )
      {
         const SVector& rowvec = rowVector( row );
         val = 0.0;         
         for( c = 0; c < rowvec.size(); ++c )
            val += rowvec.value( c ) * sol[rowvec.index( c )];

         if( LT( val, lhs( row ), feastol() ) ||
             GT( val, rhs( row ), feastol() ) )
         {
            // Minor rhs violations happen frequently, so print these
            // warnings only with verbose level INFO2 and higher.
            MSG_INFO2( spxout << "WSOLVE88 Warning! Constraint " << row
                              << " is violated by solution" << std::endl
                              << "   lhs:" << lhs( row )
                              << " <= val:" << val
                              << " <= rhs:" << rhs( row ) << std::endl; )

            if( type() == LEAVE && isRowBasic( row ) )
            {
               // find basis variable
               for( c = 0; c < nRows(); ++c )
                  if (basis().baseId(c).isSPxRowId()     
                     && (number(basis().baseId(c)) == row))
                     break;

               assert( c < nRows() );

               MSG_WARNING( spxout << "WSOLVE90 basis idx:" << c
                                   << " fVec:" << fVec()[c]
                                   << " fRhs:" << fRhs()[c]
                                   << " fTest:" << fTest()[c] << std::endl; )
            }
         }
      }
      for(int col = 0; col < nCols(); ++col )
      {
         if( LT( sol[col], lower( col ), feastol() ) ||
             GT( sol[col], upper( col ), feastol() ) )
         {
            // Minor bound violations happen frequently, so print these
            // warnings only with verbose level INFO2 and higher.
            MSG_INFO2( spxout << "WSOLVE91 Warning! Bound for column " << col
                                 << " is violated by solution" << std::endl
                                 << "   lower:" << lower( col )
                                 << " <= val:" << sol[col]
                                 << " <= upper:" << upper( col ) << std::endl; )

            if( type() == LEAVE && isColBasic( col ) )
            {
               for( c = 0; c < nRows() ; ++c)
                  if ( basis().baseId( c ).isSPxColId()    
                     && ( number( basis().baseId( c ) ) == col ))
                     break;

               assert( c < nRows() );
               MSG_WARNING( spxout << "WSOLVE92 basis idx:" << c
                                   << " fVec:" << fVec()[c]
                                   << " fRhs:" << fRhs()[c]
                                   << " fTest:" << fTest()[c] << std::endl; )
            }
         }
      }
   }
#endif  // ENABLE_ADDITIONAL_CHECKS

   return status();
}

SPxSolver::Status SPxSolver::solve()
{
   METHOD( "SPxSolver::solve()" );

   theTime.reset();
   theTime.start();

   /* if tolerances are not below irthreshold(), perform a standard floating point solve */
   if( feastol() >= irthreshold() && opttol() >= irthreshold() )
   {
      fpsolve();

      theTime.stop();
      theCumulativeTime += time();

      return status();
   }

   /* perform iterative refinement only if exact arithmetic is available */
   if( !MpqRealIsExact() )
   {
      MSG_WARNING( spxout << "WSOLVE35 Warning: Iterative refinement disabled because of missing GMP support (compile with GMP=true).\n" );

      fpsolve();

      theTime.stop();
      theCumulativeTime += time();

      return status();
   }

   /* otherwise apply iterative refinement */
   SPxBasis::SPxStatus basisstat;
   SPxBasis::Desc basisdesc;
   LPColSet slackcols;
   Real irobjlimit;
   Real irfeastol;
   Real iropttol;
   int iteroffset;
   bool refined;
   bool precisionreached;

   /* refined solution vectors */
   DVector_exact primal_ex;
   DVector_exact slack_ex;
   DVector_exact dual_ex;
   DVector_exact redcost_ex;

   /* store target tolerances */
   irfeastol = feastol();
   iropttol = opttol();

   /* limit floating point tolerance */
   if( feastol() < DEFAULT_BND_VIOL )
   {
      MSG_DEBUG( spxout << "relaxing feastol() to " << DEFAULT_BND_VIOL << "\n" );
      setFeastol(DEFAULT_BND_VIOL);
   }

   if( opttol() < DEFAULT_BND_VIOL )
   {
      MSG_DEBUG( spxout << "relaxing opttol() to " << DEFAULT_BND_VIOL << "\n" );
      setOpttol(DEFAULT_BND_VIOL);
   }

   /* deactivate objective limit; because fpsolve() uses a relaxed tolerance for checking dual feasibility, it might
    * return ABORT_VALUE although dual feasibility is not guaranteed within opttol()
    */
   /**@todo make objLimit consistent in SPxSolver */
   irobjlimit = objLimit;
   if( objLimit < infinity )
   {
      MSG_DEBUG( spxout << "deactivating objective limit\n" );
      objLimit = infinity;
   }

   /* store basis */
   MSG_DEBUG( spxout << "storing basis\n" );

   basisstat = basis().status();
   basisdesc = basis().desc();

   assert(basisstat == SPxBasis::NO_PROBLEM || basis().isDescValid(basisdesc));

   /* add artificial slack variables to convert inequality to equality constraints */
   for( int i = 0; i < nRows(); i++ )
   {
      if( lhs(i) != rhs(i) )
      {
         slackcols.add(0.0, -rhs(i), unitVecs[i], -lhs(i));
         changeLhs(i, 0.0);
         changeRhs(i, 0.0);
      }
   }

   MSG_DEBUG( spxout << "adding slack columns\n" );
   addCols(slackcols);

   /* adjust basis */
   if( basisstat != SPxBasis::NO_PROBLEM )
   {
      basisdesc.reSize(nRows(), nCols());
      for( int i = 0; i < slackcols.num(); i++ )
      {
         int col;
         int row;

         col = nCols() - slackcols.num() + i;
         row = slackcols.colVector(i).index(0);
         assert(row >= 0);
         assert(row < nRows());

         switch( basisdesc.rowStatus(row) )
         {
         case SPxBasis::Desc::P_ON_LOWER:
            basisdesc.colStatus(col) = SPxBasis::Desc::P_ON_UPPER;
            break;
         case SPxBasis::Desc::P_ON_UPPER:
            basisdesc.colStatus(col) = SPxBasis::Desc::P_ON_LOWER;
            break;
         case SPxBasis::Desc::D_ON_UPPER:
            basisdesc.colStatus(col) = SPxBasis::Desc::D_ON_LOWER;
            break;
         case SPxBasis::Desc::D_ON_LOWER:
            basisdesc.colStatus(col) = SPxBasis::Desc::D_ON_UPPER;
            break;
         case SPxBasis::Desc::P_FREE:
         case SPxBasis::Desc::P_FIXED:
         case SPxBasis::Desc::D_FREE:
         case SPxBasis::Desc::D_ON_BOTH:
         case SPxBasis::Desc::D_UNDEFINED:
         default:
            basisdesc.colStatus(col) = basisdesc.rowStatus(row);
            break;
         }

         basisdesc.rowStatus(row) = SPxBasis::Desc::P_FIXED;
      }

      /* load adjusted basis */
      MSG_DEBUG( spxout << "loading adjusted basis\n" );

      loadBasis(basisdesc);
   }

   refined = false;
   precisionreached = false;

   try
   {
      /* perform initial floating point solve */
      MSG_INFO1( spxout << "\nstarting floating-point solve . . .\n\n" );
      fpsolve();

      /* count simplex iterations and update iteration limit */
      iteroffset = iterations();

      if( terminationIter() >= 0 )
      {
         assert(iterations() <= terminationIter());
         setTerminationIter(terminationIter() - iterations());
      }

      /* apply iterative refinement only if initial solve returned optimal */
      refined = (status() == OPTIMAL);
      if( refined )
      {
         /* create refined solution vectors */
         primal_ex = DVector_exact(nCols());
         slack_ex = DVector_exact(nRows());
         dual_ex = DVector_exact(nRows());
         redcost_ex = DVector_exact(nCols());

         precisionreached = refine(irfeastol, iropttol, primal_ex, slack_ex, dual_ex, redcost_ex, 5 * iterations());

         /* count simplex iterations and update iteration limit */
         iteroffset += iterations();

         if( terminationIter() >= 0 )
         {
            assert(iterations() <= terminationIter());
            setTerminationIter(terminationIter() - iterations());
         }
      }
      else
      {
         MSG_INFO1( spxout << "\nskipping refinement because of status " << status() << "\n\n" );
      }
   }
   catch( SPxException E )
   {
      MSG_INFO1( spxout << "iterative refinement threw exception: " << E.what() << "\n" );
   }

   /* adjust basis; restore lhs and rhs */
   assert(basis().status() != SPxBasis::NO_PROBLEM);
   basisdesc = basis().desc();
   for( int i = 0; i < slackcols.num(); i++ )
   {
      int col;
      int row;

      col = nCols() - slackcols.num() + i;
      row = slackcols.colVector(i).index(0);
      assert(row >= 0);
      assert(row < nRows());
      assert(basisdesc.rowStatus(row) == SPxBasis::Desc::P_FIXED || basisdesc.rowStatus(row) == SPxBasis::Desc::D_FREE);

      MSG_DEBUG( spxout << "removing slack column C" << col << ": status=" << basisdesc.colStatus(col) << ", lower=" << lower(col) << ", upper=" << upper(col) << "\n" );

      if( basisdesc.rowStatus(row) == SPxBasis::Desc::P_FIXED )
      {
         switch( basisdesc.colStatus(col) )
         {
         case SPxBasis::Desc::P_ON_LOWER:
            basisdesc.rowStatus(row) = SPxBasis::Desc::P_ON_UPPER;
            break;
         case SPxBasis::Desc::P_ON_UPPER:
            basisdesc.rowStatus(row) = SPxBasis::Desc::P_ON_LOWER;
            break;
         case SPxBasis::Desc::D_ON_UPPER:
            basisdesc.rowStatus(row) = SPxBasis::Desc::D_ON_LOWER;
            break;
         case SPxBasis::Desc::D_ON_LOWER:
            basisdesc.rowStatus(row) = SPxBasis::Desc::D_ON_UPPER;
            break;
         case SPxBasis::Desc::P_FREE:
         case SPxBasis::Desc::P_FIXED:
         case SPxBasis::Desc::D_FREE:
         case SPxBasis::Desc::D_ON_BOTH:
         case SPxBasis::Desc::D_UNDEFINED:
         default:
            basisdesc.rowStatus(row) = basisdesc.colStatus(col);
            break;
         }
      }

      changeLhs(row, -upper(col));
      changeRhs(row, -lower(col));
   }

   if( refined )
   {
      /* adjust slack values */
      for( int i = 0; i < slackcols.num(); i++ )
      {
         int col;
         int row;

         col = nCols() - slackcols.num() + i;
         row = slackcols.colVector(i).index(0);
         assert(row >= 0);
         assert(row < nRows());

         slack_ex[row] -= primal_ex[col];
      }

      /* resize primal and redcost vector */
      primal_ex.reDim(nCols() - slackcols.num());
      redcost_ex.reDim(nCols() - slackcols.num());
   }

   /* remove slack columns */
   MSG_DEBUG( spxout << "removing slack columns\n" );
   removeColRange(nCols() - slackcols.num(), nCols() - 1);

   try
   {
      /* load adjusted basis */
      MSG_DEBUG( spxout << "loading adjusted basis\n" );

      basisdesc.reSize(nRows(), nCols());
      loadBasis(basisdesc);

      /* initialize data structures after removing the slack columns; if refinement has been applied, we are at an
       * optimal basis; otherwise we have to perform a last floating point solve after removing slack columns
       */
      if( refined )
         init();
      else
         fpsolve();
   }
   catch( SPxException E )
   {
      MSG_INFO1( spxout << "iterative refinement threw exception: " << E.what() << "\n" );
      refined = false;
   }

   /* install refined solution */
   if( refined )
   {
      /**@todo use only one DVector as buffer */
      DVector primal_fp(nCols());
      DVector slack_fp(nRows());
      DVector dual_fp(nRows());
      DVector redcost_fp(nCols());

#ifndef NDEBUG
      /**@todo make debug messages */
      DVector tmp;

      tmp = DVector(nCols());
      getPrimal(tmp);
      primal_fp = DVector(primal_ex);
      setPrimal(primal_fp);
      getPrimal(primal_fp);
      tmp -= primal_fp;
      if( tmp.length() > DEFAULT_BND_VIOL )
      {
         MSG_INFO3( spxout << "distance between floating point and refined primal solution is " << tmp.length() << std::endl );
      }

      tmp = DVector(nRows());
      getSlacks(tmp);
      slack_fp = DVector(slack_ex);
      setSlacks(slack_fp);
      getSlacks(slack_fp);
      tmp -= slack_fp;
      if( tmp.length() > DEFAULT_BND_VIOL )
      {
         MSG_INFO3( spxout << "distance between floating point and refined slack vector is " << tmp.length() << std::endl );
      }

      tmp = DVector(nRows());
      getDual(tmp);
      dual_fp = DVector(dual_ex);
      setDual(dual_fp);
      getDual(dual_fp);
      tmp -= dual_fp;
      if( tmp.length() > DEFAULT_BND_VIOL )
      {
         MSG_INFO3( spxout << "distance between floating point and refined dual solution is " << tmp.length() << std::endl );
      }

      tmp = DVector(nCols());
      getRedCost(tmp);
      redcost_fp = DVector(redcost_ex);
      setRedCost(redcost_fp);
      getRedCost(redcost_fp);
      tmp -= redcost_fp;
      if( tmp.length() > DEFAULT_BND_VIOL )
      {
         MSG_INFO3( spxout << "distance between floating point and refined redcost vector is " << tmp.length() << std::endl );
      }
#else
      primal_fp = DVector(primal_ex);
      setPrimal(primal_fp);

      slack_fp = DVector(slack_ex);
      setSlacks(slack_fp);

      dual_fp = DVector(dual_ex);
      setDual(dual_fp);

      redcost_fp = DVector(redcost_ex);
      setRedCost(redcost_fp);
#endif
   }

   /* restore tolerances; this must be done after init(), because the shift procedures can't handle tolerances smaller
    * than epsilon
    */
   MSG_DEBUG( spxout << "restoring tolerances\n" );
   setFeastol(irfeastol);
   setOpttol(iropttol);

   /* restore objective limit */
   if( objLimit != irobjlimit )
   {
      MSG_DEBUG( spxout << "restoring objective limit\n" );
      objLimit = irobjlimit;
   }

   /* restore number of simplex iterations and iteration limit */
   SPxBasis::iterCount = iteroffset;
   if( terminationIter() >= 0 )
      setTerminationIter(terminationIter() + iteroffset);

   /* set refinement has been applied, we are at an optimal basis */
   if( refined )
   {
      m_status = SPxSolver::OPTIMAL;
      setBasisStatus(SPxBasis::OPTIMAL);
   }

   /* print warning if iterative refinement terminated before reaching tolerances */
   if( !precisionreached )
   {
      MSG_INFO1( spxout << "\nwarning: iterative refinement could not reach final precision\n" );
   }

   MSG_DEBUG( spxout << "returning with status " << status() << "\n" );

   theTime.stop();
   theCumulativeTime += time();

   return status();
}

bool SPxSolver::refine(
   Real                  irfeastol,          /**< primal feasibility tolerance */
   Real                  iropttol,           /**< dual feasibility tolerance */
   Vector_exact&         primal_ex,          /**< buffer to return refined primal solution values */
   Vector_exact&         slack_ex,           /**< buffer to return refined slack values */
   Vector_exact&         dual_ex,            /**< buffer to return refined dual solution values */
   Vector_exact&         redcost_ex,         /**< buffer to return refined reduced cost values */
   int                   maxitersround       /**< iteration limit per refinement round */
   )
{
   METHOD( "SPxSolver::refine()" );

   assert(maxRefines == -1 || maxRefines >= 0);
   assert(irfeastol >= 0.0);
   assert(iropttol >= 0.0);
   assert(status() == OPTIMAL);

   /* original problem in floating point precision */
   DVector orilower_fp(nCols());
   DVector oriupper_fp(nCols());
   DVector orirhs_fp(nRows());
   DVector oriobj_fp(nCols());

   /* modified problem in floating point precision */
   DVector modlower_fp(nCols());
   DVector modupper_fp(nCols());
   DVector modrhs_fp(nRows());
   DVector modobj_fp(nCols());

   /* modified problem in exact precision */
   DVector_exact modlower_ex(nCols());
   DVector_exact modupper_ex(nCols());
   DVector_exact modrhs_ex(nRows());
   DVector_exact modobj_ex(nCols());

   /* solution to modified problem in exact precision */
   DVector_exact modprimal_ex(nCols());
   DVector_exact moddual_ex(nRows());

   /* solution in floating point precision needed for cast to exact precision */
   DVector primal_fp(nCols());
   DVector dual_fp(nRows());

   /* maximum violations */
   MpqReal boundsviol_ex;
   MpqReal sidesviol_ex;
   MpqReal redcostviol_ex;

   /* scaling factors */
   MpqReal primalscale_ex;
   MpqReal dualscale_ex;

   /* number of simplex iterations, refinements in total, and refinements with actual simplex iterations being performed */
   int iteroffset;
   int nrefines;
#if defined(DEBUGGING)
   int npivotrefines;
#endif

   /* return status */
   bool precisionreached = false;

   /* flag to mark whether problem has been modified */
   bool changed = false;

   /* store current (regular) basis */
   SPxBasis::Desc basisdesc = basis().desc();
   assert(basis().status() == SPxBasis::OPTIMAL);
   assert(basis().isDescValid(basisdesc));

   /* store original problem */
   orilower_fp = lower();
   oriupper_fp = upper();
   orirhs_fp = rhs();
   getObj(oriobj_fp);

   /* get floating point solution of original problem */
   getPrimal(primal_fp);
   getDual(dual_fp);

   /* store floating point solution in exact precision */
   primal_ex = primal_fp;
   dual_ex = dual_fp;

   /* initialize values */
   nrefines = 0;
   primalscale_ex = 1;
   dualscale_ex = 1;
#if defined(DEBUGGING)
   npivotrefines = 0;
#endif

   /* count iterations performed during refinement */
   iteroffset = 0;

   /* refinement loop */
   do
   {
      MpqReal maxscale_ex;

      assert(status() == OPTIMAL);

      MSG_DEBUG( spxout << "computing violations exactly . . .\n" );

      /* compute bound violations */
      boundsviol_ex = 0;
      for( int c = 0; c < nCols(); c++ )
      {
         modlower_ex[c] = orilower_fp[c];
         modlower_ex[c] -= primal_ex[c];

         modupper_ex[c] = oriupper_fp[c];
         modupper_ex[c] -= primal_ex[c];

         if( modlower_ex[c] > boundsviol_ex )
            boundsviol_ex = modlower_ex[c];

         if( modupper_ex[c] < -boundsviol_ex )
            boundsviol_ex = -modupper_ex[c];
      }

      /* compute sides violation */
      slack_ex = computePrimalActivity(primal_ex);

      sidesviol_ex = 0;
      for( int r = 0; r < nRows(); r++ )
      {
         assert(lhs(r) == rhs(r));

         modrhs_ex[r] = orirhs_fp[r];
         modrhs_ex[r] -= slack_ex[r];

         if( modrhs_ex[r] > sidesviol_ex )
            sidesviol_ex = modrhs_ex[r];
         else if( modrhs_ex[r] < -sidesviol_ex )
            sidesviol_ex = -modrhs_ex[r];
      }
      MSG_DEBUG( spxout << "\n" );

      /* compute reduced costs and reduced cost violation */
      redcost_ex = computeDualActivity(dual_ex);
      redcostviol_ex = 0;

      for( int c = 0; c < nCols(); c++ )
      {
         SPxBasis::Desc::Status basisstat = basisdesc.colStatus(c);

         redcost_ex[c] *= -1;
         redcost_ex[c] += oriobj_fp[c];

         modobj_ex[c] = redcost_ex[c];

         if( (spxSense() == MINIMIZE && basisstat != SPxBasis::Desc::P_ON_UPPER && basisstat != SPxBasis::Desc::P_FIXED)
            || (spxSense() == MAXIMIZE && basisstat != SPxBasis::Desc::P_ON_LOWER && basisstat != SPxBasis::Desc::P_FIXED) )
         {
            if( redcost_ex[c] < -redcostviol_ex )
               redcostviol_ex = -redcost_ex[c];
         }

         if( (spxSense() == MINIMIZE && basisstat != SPxBasis::Desc::P_ON_LOWER && basisstat != SPxBasis::Desc::P_FIXED)
            || (spxSense() == MAXIMIZE && basisstat != SPxBasis::Desc::P_ON_UPPER && basisstat != SPxBasis::Desc::P_FIXED) )
         {
            if( redcost_ex[c] > redcostviol_ex )
               redcostviol_ex = redcost_ex[c];
         }
      }

      /* at this point the vectors primal_ex, slack_ex, dual_ex, redcost_ex are computed and we may exit */

      /* terminate if tolerances are satisfied */
      if( boundsviol_ex <= irfeastol && sidesviol_ex <= irfeastol && redcostviol_ex <= iropttol )
      {
         MSG_INFO1( spxout << "\nrefinement finished: tolerances reached\n\n" );
         assert(status() == OPTIMAL);
         precisionreached = true;
         break;
      }

      /* terminate if maximum number of refinements is reached */
      if( maxRefines >= 0 && nrefines >= maxRefines )
      {
         MSG_INFO1( spxout << "\nrefinement finished: maximum number of refinement rounds reached\n\n" );
         assert(status() == OPTIMAL);
         break;
      }

      /* otherwise continue */
      MSG_INFO1( spxout << "\nstarting refinement round " << nrefines+1 << ": " );

      /* compute primal scaling factor; limit increase in scaling by tolerance used in floating point solve */
      maxscale_ex = primalscale_ex / feastol();

      primalscale_ex = boundsviol_ex > sidesviol_ex ? boundsviol_ex : sidesviol_ex;
      assert(primalscale_ex >= 0);

      if( primalscale_ex > 0 )
      {
         primalscale_ex = 1 / primalscale_ex;
         if( primalscale_ex > maxscale_ex )
            primalscale_ex = maxscale_ex;
      }
      else
         primalscale_ex = maxscale_ex;

      if( primalscale_ex < 1 )
         primalscale_ex = 1;

      MSG_INFO1( spxout << "scaling primal by " << primalscale_ex );

      /* compute dual scaling factor; limit increase in scaling by tolerance used in floating point solve */
      maxscale_ex = dualscale_ex / opttol();

      dualscale_ex = redcostviol_ex;
      assert(dualscale_ex >= 0);

      if( dualscale_ex > 0 )
      {
         dualscale_ex = 1 / dualscale_ex;
         if( dualscale_ex > maxscale_ex )
            dualscale_ex = maxscale_ex;
      }
      else
         dualscale_ex = maxscale_ex;

      if( dualscale_ex < 1 )
         dualscale_ex = 1;

      MSG_INFO1( spxout << ", scaling dual by " << dualscale_ex << " . . .\n\n" );

      /* perform primal scaling */
      modlower_ex *= primalscale_ex;
      modupper_ex *= primalscale_ex;
      modrhs_ex *= primalscale_ex;

      /* perform dual scaling */
      modobj_ex *= dualscale_ex;

      /* change bounds */
      modlower_fp = DVector(modlower_ex);
      modupper_fp = DVector(modupper_ex);
      changeBounds(modlower_fp, modupper_fp);

      /* change sides */
      modrhs_fp = DVector(modrhs_ex);
      changeRange(modrhs_fp, modrhs_fp);

      /* change objective function */
      modobj_fp = DVector(modobj_ex);
      changeObj(modobj_fp);

      /* mark problem as changed */
      changed = true;

      MSG_DEBUG( spxout << "solving modified problem . . .\n" );

      /* load basis */
      loadBasis(basisdesc);

      /* count refinement rounds */
      nrefines++;

      /* solve modified problem */
      try
      {
         int maxiters = terminationIter();

         /* limit number of iterations per refinement round */
         if( maxitersround >= 0 && (terminationIter() < 0 || terminationIter() > maxitersround) )
            setTerminationIter(maxitersround);

         fpsolve();

         setTerminationIter(maxiters);
      }
      catch( SPxException E )
      {
         MSG_DEBUG( spxout << "fpsolve() threw exception: " << E.what() << "\n" );
      }

#if defined(DEBUGGING)
      /* remember whether we moved to a new basis */
      if( iterations() > 0 )
         npivotrefines = nrefines;
#endif

      /* count simplex iterations and update global iteration limit */
      iteroffset += iterations();

      if( terminationIter() >= 0 )
      {
         assert(iterations() <= terminationIter());
         setTerminationIter(terminationIter() - iterations());
      }

      /* if modified problem was not solved to optimality, stop refinement */
      if( status() != OPTIMAL )
      {
         MSG_INFO1( spxout << "\nrefinement finished: modified problem terminated with status " << status() << "\n\n" );

         /* load last optimal basis */
         loadBasis(basisdesc);

         break;
      }
      /* otherwise continue and correct primal and dual solution */
      else
      {
         int nadjusted = 0;

         /* get basis */
         basisdesc = basis().desc();
         assert(basis().isDescValid(basisdesc));

         /* get floating point solution of modified problem */
         getPrimal(primal_fp);
         getDual(dual_fp);

         /* store floating point solution in exact precision */
         modprimal_ex = primal_fp;
         moddual_ex = dual_fp;

         /* correct primal solution */
         MSG_DEBUG( spxout << "correcting primal solution . . ." );

         modprimal_ex *= MpqReal(1 / primalscale_ex);
         primal_ex += modprimal_ex;

         /* force values of nonbasic variables to bounds */
         for( int c = 0; c < nCols(); c++ )
         {
            SPxBasis::Desc::Status basisstat = basisdesc.colStatus(c);

            if( basisstat == SPxBasis::Desc::P_ON_LOWER && primal_ex[c] != orilower_fp[c] )
            {
               primal_ex[c] = orilower_fp[c];
               nadjusted++;
            }
            else if( basisstat == SPxBasis::Desc::P_ON_UPPER && primal_ex[c] != oriupper_fp[c] )
            {
               primal_ex[c] = oriupper_fp[c];
               nadjusted++;
            }
            else if( basisstat == SPxBasis::Desc::P_FIXED )
            {
               assert(orilower_fp[c] == oriupper_fp[c]);

               if( primal_ex[c] != orilower_fp[c] )
               {
                  primal_ex[c] = orilower_fp[c];
                  nadjusted++;
               }
            }
            else if( basisstat == SPxBasis::Desc::P_FREE && primal_ex[c] != 0 )
            {
               primal_ex[c] = 0;
               nadjusted++;
            }
         }

         MSG_DEBUG( spxout << " adjusted " << nadjusted << " nonbasic variables to their bounds\n" );

         /* correct dual solution */
         MSG_DEBUG( spxout << "correcting dual solution . . .\n" );

         moddual_ex *= MpqReal(1 / dualscale_ex);
         dual_ex += moddual_ex;
      }
   }
   while( true );

   /* output violations; the reduced cost violations for artificially introduced slack columns are actually violations of the dual multipliers */
   MSG_INFO3( spxout << "maximum violations: bounds=" << boundsviol_ex << ", sides=" << sidesviol_ex << ", duals/redcosts=" << redcostviol_ex << "\n" );

   MSG_DEBUG(
      MpqReal viol_ex = boundsviol_ex;
      if( sidesviol_ex > viol_ex )
         viol_ex = sidesviol_ex;
      if( redcostviol_ex > viol_ex )
         viol_ex = redcostviol_ex;

      spxout << "maxboundsviol: " << boundsviol_ex
      << " maxsidesviol: " << sidesviol_ex
      << " maxredcostviol: " << redcostviol_ex
      << " maxviol: " << viol_ex
      << " after ncorrections: " << nrefines
      << " --"
      << " maxcorrections: " << maxRefines
      << " primaltol: " << irfeastol
      << " dualtol: " << iropttol
      << " pivotiters: " << npivotrefines
      << "\n"
      );

   /* restore original problem and load last regular basis */
   if( changed )
   {
      MSG_DEBUG( spxout << "restoring original problem . . .\n" );

      changeBounds(orilower_fp, oriupper_fp);
      changeRange(orirhs_fp, orirhs_fp);
      changeObj(oriobj_fp);

      loadBasis(basisdesc);
   }

   /* remember total number of simplex iterations during refinement */
   SPxBasis::iterCount = iteroffset;

   /* restore iteration limit */
   if( terminationIter() >= 0 )
      setTerminationIter(terminationIter() + iteroffset);

   return precisionreached;
}

void SPxSolver::testVecs()
{
   METHOD( "SPxSolver::testVecs()" );

   assert(SPxBasis::status() > SPxBasis::SINGULAR);

   DVector tmp(dim());

   tmp = *theCoPvec;
   multWithBase(tmp);
   tmp -= *theCoPrhs;
   if (tmp.length() > leavetol())
   {
      MSG_INFO3( spxout << "ISOLVE93 " << iteration() << ":\tcoP error = \t"
                        << tmp.length() << std::endl; )

      tmp.clear();
      SPxBasis::coSolve(tmp, *theCoPrhs);
      multWithBase(tmp);
      tmp -= *theCoPrhs;
      MSG_INFO3( spxout << "ISOLVE94\t\t" << tmp.length() << std::endl; )

      tmp.clear();
      SPxBasis::coSolve(tmp, *theCoPrhs);
      tmp -= *theCoPvec;
      MSG_INFO3( spxout << "ISOLVE95\t\t" << tmp.length() << std::endl; )
   }

   tmp = *theFvec;
   multBaseWith(tmp);
   tmp -= *theFrhs;
   if (tmp.length() > entertol())
   {
      MSG_INFO3( spxout << "ISOLVE96 " << iteration() << ":\t  F error = \t"
                           << tmp.length() << std::endl; )

      tmp.clear();
      SPxBasis::solve(tmp, *theFrhs);
      tmp -= *theFvec;
      MSG_INFO3( spxout << "ISOLVE97\t\t" << tmp.length() << std::endl; )
   }

   if (type() == ENTER)
   {
      for (int i = 0; i < dim(); ++i)
      {
         if (theCoTest[i] < -leavetol() && isCoBasic(i))
         {
            /// @todo Error message "this shalt not be": shalt this be an assert (also below)?
            MSG_ERROR( spxout << "ESOLVE98 testVecs: theCoTest: this shalt not be!"
                              << std::endl
                              << "  i=" << i 
                              << ", theCoTest[i]=" << theCoTest[i]
                              << ", leavetol()=" << leavetol() << std::endl; )
         }
      }

      for (int i = 0; i < coDim(); ++i)
      {
         if (theTest[i] < -leavetol() && isBasic(i))
         {
            MSG_ERROR( spxout << "ESOLVE99 testVecs: theTest: this shalt not be!"
                              << std::endl
                              << "  i=" << i 
                              << ", theTest[i]=" << theTest[i]
                              << ", leavetol()=" << leavetol() << std::endl; )
         }
      }
   }
}

bool SPxSolver::terminate()
{
   METHOD( "SPxSolver::terminate()" );
#ifdef ENABLE_ADDITIONAL_CHECKS
   if (SPxBasis::status() > SPxBasis::SINGULAR)
      testVecs();
#endif

   int redo = dim();

   if (redo < 1000)
      redo = 1000;

   if (iteration() > 10 && iteration() % redo == 0)
   {
#ifdef ENABLE_ADDITIONAL_CHECKS
      DVector cr(*theCoPrhs);
      DVector fr(*theFrhs);
#endif 

      if (type() == ENTER)
         computeEnterCoPrhs();
      else
         computeLeaveCoPrhs();

      computeFrhs();

#ifdef ENABLE_ADDITIONAL_CHECKS
      cr -= *theCoPrhs;
      fr -= *theFrhs;
      if (cr.length() > leavetol())
         MSG_WARNING( spxout << "WSOLVE50 unexpected change of coPrhs " 
                             << cr.length() << std::endl; )
      if (fr.length() > entertol())
         MSG_WARNING( spxout << "WSOLVE51 unexpected change of   Frhs " 
                             << fr.length() << std::endl; )
#endif

      if (updateCount > 1)
      {
         MSG_INFO3( spxout << "ISOLVE52 terminate triggers refactorization" 
                           << std::endl; )
         factorize();
      }
      SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
      SPxBasis::solve (*theFvec, *theFrhs);

      if (pricing() == FULL)
      {
         computePvec();
         if (type() == ENTER)
            computeTest();
      }

      if (shift() > 0.0)
         unShift();
   }

   if ( maxTime >= 0 && maxTime < infinity && time() >= maxTime )
   {
      MSG_INFO2( spxout << "ISOLVE54 Timelimit (" << maxTime
                        << ") reached" << std::endl; )
      m_status = ABORT_TIME;
      return true;   
   }

   // objLimit is set and we are running DUAL:
   // - objLimit is set if objLimit < infinity
   // - DUAL is running if rep() * type() > 0 == DUAL (-1 == PRIMAL)
   //
   // In this case we have given a objective value limit, e.g, through a
   // MIP solver, and we want stop solving the LP if we figure out that the
   // optimal value of the current LP can not be better then this objective
   // limit. More precisely:
   // - MINIMIZATION Problem
   //   We want stop the solving process if
   //   objLimit <= current objective value of the DUAL LP
   // - MAXIMIZATION Problem
   //   We want stop the solving process if 
   //   objLimit >= current objective value of the DUAL LP
   if (objLimit < infinity && type() * rep() > 0)
   {
      // We have no bound shifts; therefore, we can trust the current
      // objective value.
      // It might be even possible to use this termination value in case of
      // bound violations (shifting) but in this case it is quite difficult
      // to determine if we already reached the limit.
      if( shift() < epsilon() && maxInfeas() + shift() <= opttol() )
      {
         // SPxSense::MINIMIZE == -1, so we have sign = 1 on minimizing
         if( spxSense() * value() <= spxSense() * objLimit ) 
         {
            MSG_INFO2( spxout << "ISOLVE55 Objective value limit (" << objLimit
               << ") reached" << std::endl; )
            MSG_DEBUG(
               spxout << "DSOLVE56 Objective value limit reached" << std::endl
                      << " (value: " << value()
                      << ", limit: " << objLimit << ")" << std::endl
                      << " (spxSense: " << int(spxSense())
                      << ", rep: " << int(rep())
                      << ", type: " << int(type()) << ")" << std::endl;
            )
            
            m_status = ABORT_VALUE;
            return true;
         }
      }
   }

   if( SPxBasis::status() >= SPxBasis::OPTIMAL  ||
       SPxBasis::status() <= SPxBasis::SINGULAR )
   {
      m_status = UNKNOWN;
      return true;
   }
   return false;
}

SPxSolver::Status SPxSolver::getPrimal (Vector& p_vector) const
{
   METHOD( "SPxSolver::getPrimal()" );

   if (!isInitialized())
   {
      /* exit if presolving/simplifier cleared the problem */
      if (status() == NO_PROBLEM)
         return status();
      throw SPxStatusException("XSOLVE06 Not Initialized");
   }
   if (rep() == ROW)
      p_vector = coPvec();
   else
   {
      const SPxBasis::Desc& ds = desc();

      for (int i = 0; i < nCols(); ++i)
      {
         switch (ds.colStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER :
            p_vector[i] = SPxLP::lower(i);
            break;
         case SPxBasis::Desc::P_ON_UPPER :
         case SPxBasis::Desc::P_FIXED :
            p_vector[i] = SPxLP::upper(i);
            break;
         case SPxBasis::Desc::P_FREE :
            p_vector[i] = 0;
            break;
         case SPxBasis::Desc::D_FREE :
         case SPxBasis::Desc::D_ON_UPPER :
         case SPxBasis::Desc::D_ON_LOWER :
         case SPxBasis::Desc::D_ON_BOTH :
         case SPxBasis::Desc::D_UNDEFINED :
            break;
         default:
            throw SPxInternalCodeException("XSOLVE07 This should never happen.");
         }
      }
      for (int j = 0; j < dim(); ++j)
      {
         if (baseId(j).isSPxColId())
            p_vector[ number(SPxColId(baseId(j))) ] = fVec()[j];
      }
   }
   return status();
}

SPxSolver::Status SPxSolver::getDual (Vector& p_vector) const
{
   METHOD( "SPxSolver::getDual()" );

   assert(isInitialized());

   if (!isInitialized()) 
   {
      /* exit if presolving/simplifier cleared the problem */
      if (status() == NO_PROBLEM)
         return status();
      throw SPxStatusException("XSOLVE08 No Problem loaded");
   }

   if (rep() == ROW)
   {
      int i;
      p_vector.clear ();
      for (i = nCols() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
            p_vector[ number(SPxRowId(baseId(i))) ] = fVec()[i];
      }
   }
   else
      p_vector = coPvec();

   p_vector *= Real(spxSense());

   return status();
}

SPxSolver::Status SPxSolver::getRedCost (Vector& p_vector) const
{
   METHOD( "SPxSolver::getRedCost()" );

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE09 No Problem loaded");    
      // return NOT_INIT;
   }

   if (rep() == ROW)
   {
      int i;
      p_vector.clear();
      if (spxSense() == SPxLP::MINIMIZE)
      {
         for (i = dim() - 1; i >= 0; --i)
         {
            if (baseId(i).isSPxColId())
               p_vector[ number(SPxColId(baseId(i))) ] = -fVec()[i];
         }
      }
      else
      {
         for (i = dim() - 1; i >= 0; --i)
         {
            if (baseId(i).isSPxColId())
               p_vector[ number(SPxColId(baseId(i))) ] = fVec()[i];
         }
      }
   }
   else
   {
      p_vector = maxObj();
      p_vector -= pVec();
      if (spxSense() == SPxLP::MINIMIZE)
         p_vector *= -1.0;
   }

   return status();
}

SPxSolver::Status SPxSolver::getPrimalray (Vector& p_vector) const
{
   METHOD( "SPxSolver::getPrimalray()" );

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE10 No Problem loaded");
      // return NOT_INIT;
   }

   assert(SPxBasis::status() == SPxBasis::UNBOUNDED);
   p_vector.clear();
   p_vector = primalRay;

   return status();
}

SPxSolver::Status SPxSolver::getDualfarkas (Vector& p_vector) const
{
   METHOD( "SPxSolver::getDualfarkas()" );

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE10 No Problem loaded");
      // return NOT_INIT;
   }

   assert(SPxBasis::status() == SPxBasis::INFEASIBLE);
   p_vector.clear();
   p_vector = dualFarkas;

   return status();
}

SPxSolver::Status SPxSolver::getSlacks (Vector& p_vector) const
{
   METHOD( "SPxSolver::getSlacks()" );

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE11 No Problem loaded");
      // return NOT_INIT;
   }

   if (rep() == COLUMN)
   {
      int i;
      const SPxBasis::Desc& ds = desc();
      for (i = nRows() - 1; i >= 0; --i)
      {
         switch (ds.rowStatus(i))
         {
         case SPxBasis::Desc::P_ON_LOWER :
            p_vector[i] = lhs(i);
            break;
         case SPxBasis::Desc::P_ON_UPPER :
         case SPxBasis::Desc::P_FIXED :
            p_vector[i] = rhs(i);
            break;
         case SPxBasis::Desc::P_FREE :
            p_vector[i] = 0;
            break;
         case SPxBasis::Desc::D_FREE :
         case SPxBasis::Desc::D_ON_UPPER :
         case SPxBasis::Desc::D_ON_LOWER :
         case SPxBasis::Desc::D_ON_BOTH :
         case SPxBasis::Desc::D_UNDEFINED :
            break;
         default:
            throw SPxInternalCodeException("XSOLVE12 This should never happen.");
         }
      }
      for (i = dim() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
            p_vector[ number(SPxRowId(baseId(i))) ] = -(*theFvec)[i];
      }
   }
   else
      p_vector = pVec();

   return status();
}

void SPxSolver::setPrimal(Vector& p_vector)
{
   METHOD( "SPxSolver::setPrimal()" );

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE20 Not Initialized");
   }

   if (rep() == ROW)
      coPvec() = p_vector;
   else
   {
      for (int j = 0; j < dim(); ++j)
      {
         if (baseId(j).isSPxColId())
            fVec()[j] = p_vector[ number(SPxColId(baseId(j))) ];
      }
   }
}

void SPxSolver::setDual(Vector& p_vector)
{
   METHOD( "SPxSolver::setDual()" );

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE21 Not Initialized");
   }

   if (rep() == ROW)
   {
      for (int i = nCols() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
         {
            if (spxSense() == SPxLP::MAXIMIZE)
               fVec()[i] = p_vector[ number(SPxRowId(baseId(i))) ];
            else
               fVec()[i] = -p_vector[ number(SPxRowId(baseId(i))) ];
         }
      }
   }
   else
   {
      coPvec() = p_vector;
      if (spxSense() == SPxLP::MINIMIZE)
         coPvec() *= -1.0;
   }
}

void SPxSolver::setRedCost(Vector& p_vector)
{
   METHOD( "SPxSolver::setRedCost()" );

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE22 Not Initialized");
   }

   if (rep() == ROW)
   {
      for (int i = dim() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxColId())
         {
            if (spxSense() == SPxLP::MINIMIZE)
               fVec()[i] = -p_vector[ number(SPxColId(baseId(i))) ];
            else
               fVec()[i] = p_vector[ number(SPxColId(baseId(i))) ];
         }
      }
   }
   else
   {
      pVec() = maxObj();

      if (spxSense() == SPxLP::MINIMIZE)
         pVec() += p_vector;
      else
         pVec() -= p_vector;
   }
}

void SPxSolver::setSlacks(Vector& p_vector)
{
   METHOD( "SPxSolver::getSlacks()" );

   assert(isInitialized());

   if (!isInitialized())
   {
      throw SPxStatusException("XSOLVE23 Not Initialized");
   }

   if (rep() == COLUMN)
   {
      for (int i = dim() - 1; i >= 0; --i)
      {
         if (baseId(i).isSPxRowId())
            (*theFvec)[i] = -p_vector[ number(SPxRowId(baseId(i))) ];
      }
   }
   else
      pVec() = p_vector;
}

SPxSolver::Status SPxSolver::status() const
{
   METHOD( "SPxSolver::status()" );
   switch( m_status )
   {
   case UNKNOWN :      
      switch (SPxBasis::status())
      {
      case SPxBasis::NO_PROBLEM :
         return NO_PROBLEM;
      case SPxBasis::SINGULAR :
         return SINGULAR;
      case SPxBasis::REGULAR :
      case SPxBasis::DUAL :
      case SPxBasis::PRIMAL :
         return UNKNOWN;
      case SPxBasis::OPTIMAL :
         return OPTIMAL;
      case SPxBasis::UNBOUNDED :
         return UNBOUNDED;
      case SPxBasis::INFEASIBLE :
         return INFEASIBLE;
      default:
         return ERROR;
      }
   case SINGULAR : 
      return m_status;
   case OPTIMAL :
      assert( SPxBasis::status() == SPxBasis::OPTIMAL );
      /*lint -fallthrough*/
   case ABORT_CYCLING :
   case ABORT_TIME :
   case ABORT_ITER :
   case ABORT_VALUE :
   case RUNNING :
   case REGULAR :
   case NOT_INIT :
   case NO_SOLVER :
   case NO_PRICER :
   case NO_RATIOTESTER :
   case ERROR:
      return m_status;
   default:
      return ERROR;
   }
}

SPxSolver::Status SPxSolver::getResult(
   Real* p_value,
   Vector* p_primal,
   Vector* p_slacks,
   Vector* p_dual,
   Vector* reduCosts) const
{
   METHOD( "SPxSolver::getResult()" );
   if (p_value)
      *p_value = this->value();
   if (p_primal)
      getPrimal(*p_primal);
   if (p_slacks)
      getSlacks(*p_slacks);
   if (p_dual)
      getDual(*p_dual);
   if (reduCosts)
      getRedCost(*reduCosts);
   return status();
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
