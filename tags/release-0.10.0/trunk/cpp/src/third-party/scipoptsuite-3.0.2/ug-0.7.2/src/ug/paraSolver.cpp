/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*             This file is part of the program and software framework       */
/*                  UG --- Ubquity Generator Framework                       */
/*                                                                           */
/*    Copyright (C) 2010-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  UG is distributed under the terms of the ZIB Academic Licence.           */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with UG; see the file COPYING. If not email to scip@zib.de.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    paraSolver.cpp
 * @brief   Base class for solver: Generic parallelized solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cfloat>
#include <climits>
#include <cassert>
#include <algorithm>
#include "paraComm.h"
#include "paraNode.h"
#include "paraInstance.h"
#include "paraSolver.h"
#include "paraSolution.h"
#include "paraSolverTerminationState.h"
#include "paraCalculationState.h"
#include "paraSolverState.h"
#ifdef _COMM_PTH
#include "paraCommPth.h"
#endif

using namespace UG;

ParaSolver::ParaSolver(
      int argc,
      char **argv,
      ParaComm     *comm,
      ParaParamSet *inParaParamSet,
      ParaInstance *inParaInstance,
      ParaDeterministicTimer *inParaDetTimer
      ) : notificationIdGenerator(0),
      paraComm(comm),
      paraParams(inParaParamSet),
      racingParams(0),
      winnerRacingParams(0),
      paraDetTimer(inParaDetTimer),
      globalBestDualBoundValueAtWarmStart(-DBL_MAX),
      globalBestIncumbentValue(DBL_MAX),
      lcBestDualBoundValue(-DBL_MAX),
      globalBestIncumbentSolution(0),
      localIncumbentSolution(0),
      pendingSolution(0),
      pendingIncumbentValue(DBL_MAX),
      paraInstance(inParaInstance),
      currentNode(0),
      newNode(0),
      collectingMode(false),
      aggressiveCollecting(false),
      nSendInCollectingMode(0),
      nCollectOnce(0),
      collectingManyNodes(false),
      terminationMode(NoTerminationMode),
      warmStarted(false),
      rampUp(false),
      racingIsInterrupted(false),
      racingWinner(false),
      anotherNodeIsRequested(false),
      waitingSpecificMessage(false),
      lightWeightRootNodeComputation(false),
      onceBreak(false),
      rootNodeTime(0.0),
      totalRootNodeTime(0.0),
      minRootNodeTime(DBL_MAX),
      maxRootNodeTime(-DBL_MAX),
      previousNotificationTime(0.0),
      paraNodeStartTime(0.0),
      previousStopTime(-DBL_MAX),
      idleTimeToFirstParaNode(0.0),
      idleTimeBetweenParaNodes(0.0),
      idleTimeAfterLastParaNode(0.0),
      idleTimeToWaitNotificationId(0.0),
      idelTimeToWaitAckCompletion(0.0), idleTimeToWaitToken(0.0), previousIdleTimeToWaitToken(0.0), offsetTimeToWaitToken(0.0),
      nSolved(0),
      nSent(0),
      nImprovedIncumbent(0),
      nSolvedWithNoPreprocesses(0),
      totalNSolved(0),
      minNSolved(INT_MAX),
      maxNSolved(INT_MIN),
      totalNSent(0),
      totalNImprovedIncumbent(0),
      nParaNodesReceived(0),
      nParaNodesSolved(0),
      nParaNodesSolvedAtRoot(0),
      nParaNodesSolvedAtPreCheck(0),
      nSimplexIterRoot(0),
      minIisum(DBL_MAX),
      maxIisum(0.0),
      minNii(INT_MAX),
      maxNii(0),
      globalIncumbnetValueUpdateFlag(false),
      notificationProcessed(false),
      eps(0.0),
      targetBound(-DBL_MAX),
      nTransferLimit(-1),
      nTransferredNodes(-1),
      previousCommTime(0.0)
{
   /** create timer for this ParaSolver */
   paraTimer = paraComm->createParaTimer();
   paraTimer->init(paraComm);
   if( paraParams->getBoolParamValue(Deterministic) )
   {
      assert(paraDetTimer);
   }

   /** register message handlers */
   for( int i = 0; i < N_TAGS; i++ )
   {
      messageHandler[i] = 0;
   }
   messageHandler[TagNode] = &UG::ParaSolver::processTagNode;
   messageHandler[TagNodeReceived] = &UG::ParaSolver::processTagNodeReceived;
   messageHandler[TagRampUp] = &UG::ParaSolver::processTagRampUp;
   messageHandler[TagRetryRampUp] = &UG::ParaSolver::processTagRetryRampUp;
   messageHandler[TagSolution] = &UG::ParaSolver::processTagSolution;
   messageHandler[TagIncumbentValue] = &UG::ParaSolver::processTagIncumbentValue;
   messageHandler[TagGlobalBestDualBoundValueAtWarmStart] = &UG::ParaSolver::processTagGlobalBestDualBoundValueAtWarmStart;
   messageHandler[TagNoNodes] = &UG::ParaSolver::processTagNoNodes;
   messageHandler[TagInCollectingMode] = &UG::ParaSolver::processTagInCollectingMode;
   messageHandler[TagCollectAllNodes] = &UG::ParaSolver::processTagCollectAllNodes;
   messageHandler[TagOutCollectingMode] = &UG::ParaSolver::processTagOutCollectingMode;
   messageHandler[TagLCBestBoundValue] = &UG::ParaSolver::processTagLCBestBoundValue;
   messageHandler[TagNotificationId] = &UG::ParaSolver::processTagNotificationId;
   messageHandler[TagTerminateRequest] = &UG::ParaSolver::processTagTerminateRequest;
   messageHandler[TagInterruptRequest] = &UG::ParaSolver::processTagInterruptRequest;
   messageHandler[TagRacingRampUpParamSet] = &UG::ParaSolver::processTagWinnerRacingRampUpParamSet;
   messageHandler[TagWinner] = &UG::ParaSolver::processTagWinner;
   messageHandler[TagLightWeightRootNodeProcess] = &UG::ParaSolver::processTagLightWeightRootNodeProcess;
   messageHandler[TagBreaking] = &UG::ParaSolver::processTagBreaking;
   messageHandler[TagToken] = &UG::ParaSolver::processTagToken;

   offsetTimeToWaitToken = ( paraParams->getRealParamValue(NotificationInterval) / (paraComm->getSize() - 1) )
         * (paraComm->getRank() - 1);
   pthread_mutex_init(&tokenAccessLock, NULL);
}

ParaSolver::~ParaSolver(
      )
{
   iReceiveMessages();
   if(racingParams) delete racingParams;
   if(winnerRacingParams) delete winnerRacingParams;
   if( globalBestIncumbentSolution ) delete globalBestIncumbentSolution;
   if(localIncumbentSolution) delete localIncumbentSolution;
   if(currentNode) delete currentNode;
   if(newNode) delete newNode;

   if( paraParams->getBoolParamValue(Deterministic) )
   {
      for(;;)
      {
         while( !waitToken(paraComm->getRank()) )
         {
            iReceiveMessages();
         }
         // paraDetTimer->update(1.0);
         // previousCommTime = paraDetTimer->getElapsedTime();
         if( paraComm->passTermToken(paraComm->getRank()) ) break;
      }
   }

   double stopTime = paraTimer->getElapsedTime();
   idleTimeAfterLastParaNode = stopTime - previousStopTime - ( idleTimeToWaitToken - previousIdleTimeToWaitToken );
   int interrupted = terminationMode == InterruptedTerminationMode ? 1 : 0;

   double detTime = -1.0;
   /*
   if( paraDetTimer )
   {
      detTime = paraDetTimer->getElapsedTime();
   } */
   ParaSolverTerminationState *paraSolverTerminationState = paraComm->createParaSolverTerminationState(
          interrupted,
          paraComm->getRank(),
          totalNSolved, minNSolved, maxNSolved, totalNSent, totalNImprovedIncumbent,
          nParaNodesReceived, nParaNodesSolved, nParaNodesSolvedAtRoot, nParaNodesSolvedAtPreCheck,
          stopTime, idleTimeToFirstParaNode, idleTimeBetweenParaNodes, idleTimeAfterLastParaNode,
          idleTimeToWaitNotificationId, idelTimeToWaitAckCompletion, idleTimeToWaitToken,
          totalRootNodeTime, minRootNodeTime, maxRootNodeTime, detTime);
   paraSolverTerminationState->send(paraComm, 0, TagTerminated);
   delete paraSolverTerminationState;
   delete paraTimer;

   if( paraInstance ) delete paraInstance;
}

int
ParaSolver::processTagNode(
      int source,
      int tag
      )
{
   if( currentNode )
   {
      newNode = paraComm->createParaNode();
      newNode->receive(paraComm, source);
      if( newNode->getMergingStatus() == 3 )  // This means, the received node is separated
      {
         newNode->setMergingStatus(-1);
      }
   }
   else
   {
      currentNode = paraComm->createParaNode();
      currentNode->receive(paraComm, source);
      if( currentNode->getMergingStatus() == 3 )  // This means, the received node is separated
      {
         currentNode->setMergingStatus(-1);
      }
   }
   anotherNodeIsRequested = false;
   nParaNodesReceived++;
   return 0;
}

int
ParaSolver::processTagNodeReceived(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagNodeReceived )
         );
   return 0;
}

int
ParaSolver::processTagRampUp(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagRampUp)
         );
   if( isRacingStage() )
   {
      assert(!racingWinner);
   }
   rampUp = true;
#ifdef _DEBUG_CHECK_RECEIVE
   std::cout << paraTimer->getElapsedTime() << " Solver" << paraComm->getRank() << " received TagRampUp" << std::endl;
#endif
   return 0;
}

int
ParaSolver::processTagRetryRampUp(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagRetryRampUp)
         );
   rampUp = false;
   return 0;
}

int
ParaSolver::processTagSolution(
      int source,
      int tag
      )
{
   ParaSolution *sol = paraComm->createParaSolution();
   sol->receive(paraComm, source);
   if( sol->getObjectiveFuntionValue()  < globalBestIncumbentValue )
      //   updateGlobalBestIncumbentValue( sol->getObjectiveFuntionValue() ) ) //  DO NOT UPDATE!!
      //   The timing of the update depends on solver upsed
   {
      if( pendingSolution )
      {
         delete pendingSolution;
      }
      pendingSolution = sol;
      pendingIncumbentValue = sol->getObjectiveFuntionValue();
   }
   else
   {
      delete sol;
   }
   return 0;
}

int
ParaSolver::processTagIncumbentValue(
      int source,
      int tag
      )
{
   double incumbent;
   PARA_COMM_CALL(
         paraComm->receive( &incumbent, 1, ParaDOUBLE, source, TagIncumbentValue)
         );
   if( paraParams->getBoolParamValue(Deterministic) )
   {
      if( incumbent  < globalBestIncumbentValue && incumbent < pendingIncumbentValue )
      {
         pendingIncumbentValue = incumbent;
      }
   }
   else
   {
      updateGlobalBestIncumbentValue(incumbent);
   }
   return 0;
}

int
ParaSolver::processTagGlobalBestDualBoundValueAtWarmStart(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( &globalBestDualBoundValueAtWarmStart, 1, ParaDOUBLE, source, TagGlobalBestDualBoundValueAtWarmStart)
         );
   return 0;
}

int
ParaSolver::processTagNoNodes(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagNoNodes)
         );
   anotherNodeIsRequested = false;
   return 0;
}

int
ParaSolver::processTagInCollectingMode(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( &nSendInCollectingMode, 1, ParaINT, source, TagInCollectingMode)
         );
   if( nSendInCollectingMode < 0 )
   {
      nSendInCollectingMode = ( 0 - nSendInCollectingMode ) - 1;
      aggressiveCollecting = true;
   }
   collectingMode = true;
   return 0;
}

int
ParaSolver::processTagCollectAllNodes(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( &nCollectOnce, 1, ParaINT, source, TagCollectAllNodes)
         );
   collectingManyNodes = true;
   return 0;
}

int
ParaSolver::processTagOutCollectingMode(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagOutCollectingMode)
         );
   collectingMode = false;
   aggressiveCollecting = false;
   nSendInCollectingMode = 0;
   return 0;
}

int
ParaSolver::processTagLCBestBoundValue(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( &lcBestDualBoundValue, 1, ParaDOUBLE, source, TagLCBestBoundValue)
         );
   return 0;
}

int
ParaSolver::processTagNotificationId(
      int source,
      int tag
      )
{
   assert(notificationProcessed);
   unsigned int notificationId;
   PARA_COMM_CALL(
         paraComm->receive( &notificationId, 1, ParaUNSIGNED, source, TagNotificationId)
         );
   if( notificationId == notificationIdGenerator) notificationProcessed = false;
   else {
      THROW_LOGICAL_ERROR4("notificationId received is ", notificationId, ", but generator value is ", notificationIdGenerator);
   }
   return 0;
}

int
ParaSolver::processTagTerminateRequest(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagTerminateRequest)
         );
   terminationMode = NormalTerminationMode;
   return 0;
}

int
ParaSolver::processTagInterruptRequest(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagInterruptRequest)
         );
   terminationMode = InterruptedTerminationMode;
   return 0;
}

int
ParaSolver::processTagWinnerRacingRampUpParamSet(
      int source,
      int tag
      )
{
   /* received racing parameter set is always winner one,
    * because the initail racing parameter set is broadcasted. */
   winnerRacingParams = paraComm->createParaRacingRampUpParamSet();
   PARA_COMM_CALL(
         winnerRacingParams->receive(paraComm, 0)
         );
   if( isRacingStage() )
   {
      assert(!racingWinner);
      racingIsInterrupted = true;
   }
#ifdef _DEBUG_CHECK_RECEIVE
   std::cout << paraTimer->getElapsedTime() << " Solver" << paraComm->getRank() << " received racing winner params" << std::endl;
#endif
   return 0;
}

int
ParaSolver::processTagWinner(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagWinner)
         );
   racingWinner = true;
   delete racingParams;
   racingParams = 0;      // No racing stage now
   return 0;
}

int
ParaSolver::processTagLightWeightRootNodeProcess(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagLightWeightRootNodeProcess)
         );
   lightWeightRootNodeComputation = true;
   setLightWeightRootNodeProcess();
   return 0;
}

int
ParaSolver::processTagBreaking(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( &targetBound, 1, ParaDOUBLE, source, TagBreaking )
         );
   PARA_COMM_CALL(
         paraComm->receive( &nTransferLimit, 1, ParaINT, source, TagBreaking )
         );
   nTransferredNodes = 0;
   collectingManyNodes = true;
   return 0;
}

int
ParaSolver::processTagToken(
      int source,
      int tag
      )
{
   assert(paraParams->getBoolParamValue(Deterministic));
   int token[2];
   PARA_COMM_CALL(
         paraComm->receive( token, 2, ParaINT, source, TagToken )
         );
   paraComm->setToken(paraComm->getRank(), token);
   return 0;
}

void
ParaSolver::run(){
   for(;;)
   {
      /***************************************************
       *  Wait a new ParaNode from ParaLoadCoordinator   *
       *  If Termination message is received, then break *
       ***************************************************/
      if( !currentNode )
      {
         if( receiveNewNodeAndReactivate() == false ) break;
      }

      if( winnerRacingParams ) // winner racing parameter set is received
      {
         setWinnerRacingParams(winnerRacingParams);
         delete winnerRacingParams;
         winnerRacingParams = 0;
      }

      if( paraParams->getBoolParamValue(SetAllDefaultsAfterRacing) &&
            racingWinner &&
            nParaNodesReceived == 2 )   // after racing, parameters should be set to default values
      {
         setWinnerRacingParams(0);
      }
      racingWinner = false;

      /** set collecting mode */
      collectingMode = false;  /* begin with out-collecting mode: NOTE: LC clear collecting mode for new solver */
      nSendInCollectingMode = 0;
      aggressiveCollecting = false;
      collectingManyNodes = false;
      nCollectOnce = 0;
      resetBreakingInfo();     // set false on collectingManyNodes in the  resetBreakingInfo
      onceBreak = false;

      /** set start time and ilde times */
      paraNodeStartTime = paraTimer->getElapsedTime();
      if( previousStopTime < 0.0 )
      {
         idleTimeToFirstParaNode = paraNodeStartTime - (idleTimeToWaitToken - previousIdleTimeToWaitToken);
      }
      else
      {
         idleTimeBetweenParaNodes += ( paraNodeStartTime - previousStopTime - ( idleTimeToWaitToken - previousIdleTimeToWaitToken ) );
      }

      /****************************************************
       * create subproblem into target solver environment *
       ***************************************************/
      createSubproblem();
      updatePendingSolution();
      if( globalBestIncumbentSolution )
      {
         tryNewSolution(globalBestIncumbentSolution);
      }
      globalIncumbnetValueUpdateFlag = false;

      /******************
       * start solving  *
       ******************/
      assert(!newNode);
      solve();

      /*****************************************************
       * notify completion of a calculation of a ParaNode  *
       *****************************************************/
      previousStopTime = paraTimer->getElapsedTime();
      double compTime = previousStopTime - paraNodeStartTime;
      previousIdleTimeToWaitToken = idleTimeToWaitToken;
      nSolved = getNNodesSolved();
      if( paraParams->getBoolParamValue(CheckEffectOfRootNodePreprocesses) && nSolved == 1)
      {
         solveToCheckEffectOfRootNodePreprocesses();
      }

      /****************************************************************************
      * send completion of calculation and update counters and accumulation time  *
      ****************************************************************************/
      if( paraParams->getBoolParamValue(Deterministic) )
      {
         do
         {
            iReceiveMessages();
         } while ( !waitToken(paraComm->getRank()) );
      }
      iReceiveMessages();     /** Before sending completion state, receiving message should be checked.
                               *   When subproblem terminated with no branch, solver lost a timing for receiving new node */
      sendCompletionOfCalculation(compTime);

      if( paraParams->getBoolParamValue(Deterministic) )
      {
         waitAckCompletion();
         // if( hasToken() ) passToken();
         paraDetTimer->update(1.0);
         previousCommTime = paraDetTimer->getElapsedTime();
#ifdef _DEBUG_DET
         std::cout << previousCommTime << " run2 R." << paraComm->getRank() << ": token passed" << std::endl;
#endif
         passToken(paraComm->getRank());
      }

      /*******************************************
      * free solving environment for subproblem  *
      ********************************************/
      freeSubproblem();

      /* if light wait root node computation is applied, rest it */
      if( lightWeightRootNodeComputation )
      {
         setOriginalRootNodeProcess();
         lightWeightRootNodeComputation = false;
      }

      /**************************
      * update current ParaNode *
      ***************************/
      assert( currentNode );
      delete currentNode;
      if( newNode )
      {
         currentNode = newNode;
         newNode = 0;
      }
      else
      {
         currentNode = 0;
      }
      if( terminationMode )
      {
         break;
      }
   }
   return;
}

bool
ParaSolver::receiveNewNodeAndReactivate()
{
   for(;;)
   {
      int source;
      int tag;
      int status;
      /*******************************************
       *  waiting for any message form anywhere  *
       *******************************************/
      if( paraParams->getBoolParamValue(Deterministic) )
      {
         do
         {
            iReceiveMessages();
         } while( !waitToken(paraComm->getRank()) );
         iReceiveMessages();
         paraDetTimer->update(1.0);
         previousCommTime = paraDetTimer->getElapsedTime();
#ifdef _DEBUG_DET
         std::cout << previousCommTime << " receiveNewNodeAndReactivate R." << paraComm->getRank() << ": token passed" << std::endl;
#endif
         passToken(paraComm->getRank());
      }
      else
      {
         paraComm->probe(&source, &tag);
         if( messageHandler[tag] )
         {
            status = (this->*messageHandler[tag])(source, tag);
            if( status )
            {
               std::ostringstream s;
               s << "[ERROR RETURN form Message Hander]:" <<  __FILE__ <<  "] func = "
                 << __func__ << ", line = " << __LINE__ << " - "
                 << "process tag = " << tag << std::endl;
               abort();
            }
         }
         else
         {
            THROW_LOGICAL_ERROR3( "No message hander for ", tag, " is not registered" );
         }
      }
      if( currentNode )
      {
         if( paraParams->getBoolParamValue(Deterministic)  )
         {
            previousNotificationTime = paraDetTimer->getElapsedTime();
         }
         else
         {
            previousNotificationTime = paraTimer->getElapsedTime();
         }
         return true;
      }
      if( terminationMode ) break;
   }
   return false;
}

void
ParaSolver::iReceiveMessages(
      )
{
#ifdef _DEBUG_CHECK_RECEIVE
   static double previousRreceiveCheckTime = DBL_MAX;
   double currentTime = paraTimer->getElapsedTime();
   if( ( currentTime - previousRreceiveCheckTime ) < -1.0 )
   {
      std::cout << currentTime << " Solver" << paraComm->getRank() << " No check receiving message over 500 (sec.) is logging." << std::endl;
   }
   if( ( currentTime - previousRreceiveCheckTime ) > 500.0 )
   {
      std::cout << currentTime << " Solver" << paraComm->getRank() << " did not check receiving message over 500 (sec.)" << std::endl;
      writeSubproblem();
   }
   previousRreceiveCheckTime = currentTime;
#endif
   /************************************************************************
    * This fucntion is called from a callback routine of the target solver *
    * **********************************************************************/
   int source;
   int tag;
   int status;
   /************************************
    * check if there are some messages *
    ************************************/
   while( paraComm->iProbe(&source, &tag) )
   {
      if( messageHandler[tag] )
      {
         status = (this->*messageHandler[tag])(source, tag);
         if( status )
         {
            std::ostringstream s;
            s << "[ERROR RETURN form Message Hander]:" <<  __FILE__ <<  "] func = "
              << __func__ << ", line = " << __LINE__ << " - "
              << "process tag = " << tag << std::endl;
            abort();
         }
      }
      else
      {
         THROW_LOGICAL_ERROR3( "No message hander for ", tag, " is not registered" );
      }

   }
}

void
ParaSolver::setRootNodeTime(
      )
{
   rootNodeTime = paraTimer->getElapsedTime() - paraNodeStartTime;
   if( rootNodeTime < minRootNodeTime )
   {
      minRootNodeTime = rootNodeTime;
   }
   if( rootNodeTime > maxRootNodeTime )
   {
      maxRootNodeTime = rootNodeTime;
   }
}

void
ParaSolver::sendCompletionOfCalculation(
      double compTime
      )
{
   int terminationState = CompTerminatedNormally;
   if( terminationMode == InterruptedTerminationMode )
   {
      terminationState = CompTerminatedByInterruptRequest;
   }
   else
   {
      if ( newNode )
      {
         terminationState = CompTerminatedByAnotherNode;
      }
      else
      {
         if( anotherNodeIsRequested )
         {
            for(;;)
            {
               int source;
               int tag;
               int status;
               /*******************************************
                *  waiting for any message form anywhere  *
                *******************************************/
               paraComm->probe(&source, &tag);
                if( messageHandler[tag] )
                {
                   status = (this->*messageHandler[tag])(source, tag);
                   if( status )
                   {
                      std::ostringstream s;
                      s << "[ERROR RETURN form Message Hander]:" <<  __FILE__ <<  "] func = "
                        << __func__ << ", line = " << __LINE__ << " - "
                        << "process tag = " << tag << std::endl;
                      abort();
                   }
                }
                else
                {
                   THROW_LOGICAL_ERROR3( "No message hander for ", tag, " is not registered" );
                }
               if( newNode )
               {
                  terminationState = CompTerminatedByAnotherNode;
                  break;
               }
               if( !anotherNodeIsRequested ) break;
            }
         }
         if( currentNode->getMergingStatus() == 3 )
         {
            terminationState = CompInterruptedInMerging;
         }
      }

      if( isRacingStage() &&
            !newNode && !racingWinner)
      {
         if( !racingIsInterrupted && wasTerminatedNormally() )
         {
            // assert(getNNodesLeft() == 0);   // hard time limit case, this happned.
            /* CompTerminatedInRacingStage means computation is finished, so terminates all solvers */
            terminationState = CompTerminatedInRacingStage;
         }
         else
         {
            /* CompInterruptedInRacingStage means computation is interrupted by termination of the other solvers */
            terminationState = CompInterruptedInRacingStage;   // even if a racing solver terminated badly, just keep running.
            racingIsInterrupted = false;
         }

      }
      else
      {
         if( !wasTerminatedNormally() )
         {
            THROW_LOGICAL_ERROR3( "ParaSolver", paraComm->getRank(), " was terminated abnormally." );
         }
      }
   }

   double averageSimplexIter = 0.0;

   if( nSolved <= 1 )
   {
      /** nSolved > 1 is set within callback routine */
      setRootNodeTime();
      setRootNodeSimplexIter(getSimplexIter());
   }
   else
   {
      averageSimplexIter = static_cast<double>( ( getSimplexIter() - nSimplexIterRoot )/(nSolved - 1) );
   }
   ParaCalculationState *paraCalculationState = paraComm->createParaCalculationState(
         compTime, rootNodeTime, nSolved,nSent, nImprovedIncumbent, terminationState, nSolvedWithNoPreprocesses,
         nSimplexIterRoot, averageSimplexIter, getNRestarts(), minIisum, maxIisum, minNii, maxNii);
   paraCalculationState->send(paraComm, 0, TagCompletionOfCalculation);
   delete paraCalculationState;

   /*******************
    * update counters *
    *******************/
   if( nSolved < minNSolved )
   {
      minNSolved = nSolved;
   }
   if( nSolved > maxNSolved )
   {
      maxNSolved = nSolved;
   }
   totalNSolved += nSolved;
   totalNSent += nSent;
   totalNImprovedIncumbent += nImprovedIncumbent;
   nParaNodesSolved++;
   if( nSolved == 1)
   {
      nParaNodesSolvedAtRoot++;
   }

   nSolved = 0;
   nSent = 0;
   nImprovedIncumbent = 0;
   nSolvedWithNoPreprocesses = 0;
   /**********************************
   * accumulate total root node time *
   ***********************************/
   totalRootNodeTime += rootNodeTime;
   rootNodeTime = 0.0;

   minIisum = DBL_MAX;
   maxIisum = 0.0;
   minNii = INT_MAX;
   maxNii = 0;

   double detTime = -1.0;
   if( paraParams->getBoolParamValue(Deterministic)  )
   {
      detTime = paraDetTimer->getElapsedTime();
   }
   double stopTime = paraTimer->getElapsedTime();
   if( isRacingStage() )
   {
      if( newNode )
      {
         nParaNodesReceived--;
      }
      /** Transfer SolverTermination state during racing ramp-up */
      ParaSolverTerminationState *paraSolverTerminationState = paraComm->createParaSolverTerminationState(
    		3,    /** interupted flag == 3 means the information for racing ramp-up */
            paraComm->getRank(),
            totalNSolved, minNSolved, maxNSolved, totalNSent, totalNImprovedIncumbent,
            nParaNodesReceived, nParaNodesSolved, nParaNodesSolvedAtRoot, nParaNodesSolvedAtPreCheck,
            stopTime, idleTimeToFirstParaNode, idleTimeBetweenParaNodes, 0.0,
            idleTimeToWaitNotificationId, idelTimeToWaitAckCompletion, idleTimeToWaitToken,
            totalRootNodeTime, minRootNodeTime, maxRootNodeTime, detTime);
      paraSolverTerminationState->send(paraComm, 0, TagTerminated);
      delete paraSolverTerminationState;
      assert(! racingWinner);
      /** re-initialize all counters, winner counts on current counters */
      minNSolved = INT_MAX;
      maxNSolved = INT_MIN;
      totalNSolved = 0;
      totalNSent = 0;
      totalNImprovedIncumbent = 0;
      if( newNode )
         nParaNodesReceived = 1;
      else
         nParaNodesReceived = 0;
      nParaNodesSolved = 0;
      nParaNodesSolvedAtRoot = 0;
      nParaNodesSolvedAtPreCheck = 0;
      totalRootNodeTime = 0.0;
      minRootNodeTime = DBL_MAX;
      maxRootNodeTime = -DBL_MAX;
      terminateRacing();
   }
   else
   {
      /** Transfer SolverTermination state to save statistic information for checkpoint */
      ParaSolverTerminationState *paraSolverTerminationState = paraComm->createParaSolverTerminationState(
    		2,    /** interupted flag == 2 means the information for checkpoint */
            paraComm->getRank(),
            totalNSolved, minNSolved, maxNSolved, totalNSent, totalNImprovedIncumbent,
            nParaNodesReceived, nParaNodesSolved, nParaNodesSolvedAtRoot, nParaNodesSolvedAtPreCheck,
            stopTime, idleTimeToFirstParaNode, idleTimeBetweenParaNodes, 0.0,
            idleTimeToWaitNotificationId, idelTimeToWaitAckCompletion, idleTimeToWaitToken,
            totalRootNodeTime, minRootNodeTime, maxRootNodeTime, detTime);
      paraSolverTerminationState->send(paraComm, 0, TagTerminated);
      delete paraSolverTerminationState;
   }
   nSendInCollectingMode = 0;
}

/*** sendIfImprovedSolutionWasFound routine should be removed in the future **/
void
ParaSolver::sendIfImprovedSolutionWasFound(
      ParaSolution *sol
      )
{
   if( EPSLT( sol->getObjectiveFuntionValue(), globalBestIncumbentValue, eps ) )
   {
      globalBestIncumbentValue = sol->getObjectiveFuntionValue();
      globalIncumbnetValueUpdateFlag = true;
      assert( localIncumbentSolution == 0 );
      sol->send(paraComm, 0);
      delete sol;
      nImprovedIncumbent++;
   }
   else
   {
      delete sol;
   }
}

void
ParaSolver::saveIfImprovedSolutionWasFound(
      ParaSolution *sol
      )
{
   if( EPSLT( sol->getObjectiveFuntionValue(), globalBestIncumbentValue, eps ) )  // compare to globalBestIncumbentValue
   {                                                                              // no solution sending is accepted!!!
      globalBestIncumbentValue = sol->getObjectiveFuntionValue();
      globalIncumbnetValueUpdateFlag = true;
      if( localIncumbentSolution  )
      {
         if( EPSLT( sol->getObjectiveFuntionValue(), localIncumbentSolution->getObjectiveFuntionValue(), eps ) )
         {
            // there is a possibility to be updated several times
            delete localIncumbentSolution;
            localIncumbentSolution = sol;
            nImprovedIncumbent++;
         }
         else
         {
            delete localIncumbentSolution;
            localIncumbentSolution = 0;
         }
      }
      else
      {
         localIncumbentSolution = sol;
         nImprovedIncumbent++;
      }
   }
   else
   {
      delete sol;
   }
}

void
ParaSolver::sendLocalSolution(
      )
{
   if( localIncumbentSolution && (!notificationProcessed) ) // if solution is sent in notification is processed,
                                                            //  dead lock may be occurred depending of MPI system buffer size
   {
      if( !globalBestIncumbentSolution )
      {
         localIncumbentSolution->send(paraComm, 0);
         globalBestIncumbentSolution = localIncumbentSolution;
      }
      else
      {
         if( EPSLT(localIncumbentSolution->getObjectiveFuntionValue(), globalBestIncumbentSolution->getObjectiveFuntionValue(), eps) )
            // NOTE: globalBestIncumbnetValue may be an objective function value of localIncumbentSolution
         {
            localIncumbentSolution->send(paraComm, 0);
            delete globalBestIncumbentSolution;
            globalBestIncumbentSolution = localIncumbentSolution;
         }
         else
         {
            delete localIncumbentSolution;
         }
      }
      localIncumbentSolution = 0;
   }
}

bool
ParaSolver::notificationIsNecessary(
      )
{
   // if( paraParams->getBoolParamValue(Deterministic) )
   // {
   //  return true;  // always true, in case of deterministic run
   // }
   if( !rampUp )
   {
      if( isRacingStage() )
      {
         if( lcBestDualBoundValue < getDualBoundValue() )
         {
            return true;
         }
         else
         {
            if( paraParams->getBoolParamValue(Deterministic) )
            {
               if( ( paraDetTimer->getElapsedTime() -  previousNotificationTime )
                  > paraParams->getRealParamValue(NotificationInterval) )
               {
                  return true;
               }
               else
               {
                  return false;
               }
            }
            else
            {
               if( ( paraTimer->getElapsedTime() -  previousNotificationTime )
                  > paraParams->getRealParamValue(NotificationInterval) )
               {
                  return true;
               }
               else
               {
                  return false;
               }
            }
         }
      }
      else
      {  // normal ramp-up phase
         return true;
      }
   }
   else
   {
      if( collectingMode || collectingManyNodes )
      {
         return true;
      }
      if( paraParams->getBoolParamValue(Deterministic) )
      {
         if(  ( paraDetTimer->getElapsedTime() -  previousNotificationTime )
               > paraParams->getRealParamValue(NotificationInterval) )
         {
            return true;
         }
         else
         {
            return false;
         }
      }
      else
      {
         if( ( paraTimer->getElapsedTime() -  previousNotificationTime )
               > paraParams->getRealParamValue(NotificationInterval) )
         {
            return true;
         }
         else
         {
            return false;
         }
      }
   }
}

void
ParaSolver::sendSolverState(
      long long nNodesSolved,
      int nNodesLeft,
      double bestDualBoundValue,
      double detTime
      )
{
   if(!notificationProcessed)
   {
      int racingStage = 0;   /** assume not racing stage */
      if( isRacingStage() )
      {
         racingStage = 1;
      }

      double tempGlobalBestPrimalBound = DBL_MAX;
      if( globalBestIncumbentSolution )
      {
         tempGlobalBestPrimalBound = globalBestIncumbentSolution->getObjectiveFuntionValue();
      }

      ParaSolverState *solverState = paraComm->createParaSolverState(
            racingStage,
            ++notificationIdGenerator,
            currentNode->getLcId(), currentNode->getGlobalSubtreeIdInLc(),
            nNodesSolved, nNodesLeft, std::max( bestDualBoundValue,currentNode->getDualBoundValue() ),
            tempGlobalBestPrimalBound,
            detTime
      );
      solverState->send(paraComm, 0, TagSolverState);
      delete solverState;
      if( paraParams->getBoolParamValue(Deterministic) )
      {
         previousNotificationTime = paraDetTimer->getElapsedTime();
      }
      else
      {
         previousNotificationTime = paraTimer->getElapsedTime();
      }
      notificationProcessed = true;
   }
}

int
ParaSolver::getThresholdValue(
      int nNodes      /**< number of processed nodes, including the focus node */
      )
{
   if( paraParams->getBoolParamValue(Deterministic) )
   {
      return 10;
   }

   double compTime = paraTimer->getElapsedTime() - (idleTimeToFirstParaNode + idleTimeBetweenParaNodes);
   double meanRootNodeTime = ( totalRootNodeTime + rootNodeTime )/(nParaNodesSolved + 1);  // +1 is this ParaNode
   double meanNodeTime = -1.0;             // unrealiable
   if( ( compTime - (totalRootNodeTime + rootNodeTime ) ) > 10.0  // less equal than 10.0 sec. treats as unreliable
         && meanRootNodeTime > 0.0001 ) // meanRootNode time less equal than 0.0001 sec. treats as unreliable
   {
      meanNodeTime = ( compTime - (totalRootNodeTime + rootNodeTime ) )  / (totalNSolved + nNodes - nParaNodesSolved);
                                                                      // compute mean node comp. time except root nodes(nParaNodesSolved+1) and the focus node(-1)
      if( meanNodeTime < 0.0 ) meanNodeTime = 0.0;
   }

   int n;
   if( meanNodeTime < 0.000001 )
   {   // heuristic initial value setting
      n = ( rootNodeTime < 0.01 ) ? 100 : ( ( rootNodeTime < 0.1 ) ? 50 : ( (rootNodeTime < 1.0) ? 10 : 5 ) );
   }
   else
   {
      if( ( meanRootNodeTime / meanNodeTime ) > 5.0 )
      {  // base value at least 5.0
         n = (int) (
               paraParams->getRealParamValue(MultiplierToDetermineThresholdValue)
               * ( meanRootNodeTime / meanNodeTime ) );
         if( n > 100 )
         {
             if( meanNodeTime > 1.0 )
             {
                n = 3;
             }
             else
             {
                n = 100;
             }
         }
      }
      else
      {   // heuristic value setting
         n = ( rootNodeTime < 0.01 ) ? 100 : ( ( rootNodeTime < 0.1 ) ? 50 : ( (rootNodeTime < 1.0) ? 10 : 5 ) );
      }
   }
   if( n < 2 ) n = 2;
   return n;
}

void
ParaSolver::sendParaNode(
      long long n,
      int depth,
      double dualBound,
      double estimateValue,
      ParaDiffSubproblem *diffSubproblem
      )
{
   ParaNode *node = paraComm->createParaNode(
         NodeId( SubtreeId(currentNode->nodeId.subtreeId.lcId,
                           currentNode->nodeId.subtreeId.globalSubtreeIdInLc,
                           paraComm->getRank() ),
                           0),
         NodeId( currentNode->nodeId.subtreeId, n),
         currentNode->getDepth() + depth,
         std::max( dualBound, currentNode->getDualBoundValue() ),
         dualBound,
         estimateValue,
         diffSubproblem
         );
   double inIdleTime = paraTimer->getElapsedTime();
   node->send(paraComm, 0);
   delete node;
   int tag;
   int status;
   waitingSpecificMessage = true;
   while( paraComm->waitSpecTagFromSpecSource(0, ParaBYTE, TagNodeReceived, &tag) != 0 )
   {
      idleTimeToWaitNotificationId += ( paraTimer->getElapsedTime() - inIdleTime );
      status = (this->*messageHandler[tag])(0, tag);
      if( status )
      {
         std::ostringstream s;
         s << "[ERROR RETURN form Message Hander]:" <<  __FILE__ <<  "] func = "
           << __func__ << ", line = " << __LINE__ << " - "
           << "process tag = " << tag << std::endl;
         abort();
      }
      inIdleTime = paraTimer->getElapsedTime();
   }
   idleTimeToWaitNotificationId += ( paraTimer->getElapsedTime() - inIdleTime );
   //
   // receive the NULL message for TagNodeReceived
   //
   status = (this->*messageHandler[tag])(0, tag);
   if( status )
   {
      std::ostringstream s;
      s << "[ERROR RETURN form Message Hander]:" <<  __FILE__ <<  "] func = "
        << __func__ << ", line = " << __LINE__ << " - "
        << "process tag = " << tag << std::endl;
      abort();
   }
   waitingSpecificMessage = false;
   nSent++;
   if( nSendInCollectingMode > 0 ) nSendInCollectingMode--;
   if( nCollectOnce > 0 ) nCollectOnce--;
   if( aggressiveCollecting && nSendInCollectingMode == 0 ) aggressiveCollecting = false;
   if( collectingManyNodes && nCollectOnce == 0 ) collectingManyNodes = false;
   if( isBreaking() ) nTransferredNodes++;
}

void
ParaSolver::sendAnotherNodeRequest(
      double bestDualBoundValue
      )
{
    if( anotherNodeIsRequested ) return;
    PARA_COMM_CALL(
          paraComm->send( &bestDualBoundValue, 1, ParaDOUBLE, 0, TagAnotherNodeRequest)
          );
    anotherNodeIsRequested = true;

}

bool
ParaSolver::updateGlobalBestIncumbentValue(
      double newValue
      )
{
    if( newValue < globalBestIncumbentValue ){
        globalBestIncumbentValue = newValue;
        globalIncumbnetValueUpdateFlag = true;
        return true;
    } else {
        return false;
    }
}

bool
ParaSolver::updateGlobalBestIncumbentSolution(
      ParaSolution *sol
      )
{
   if( globalBestIncumbentSolution  )
   {
      if( sol->getObjectiveFuntionValue() < globalBestIncumbentSolution->getObjectiveFuntionValue() )
      {
         delete globalBestIncumbentSolution;
         globalBestIncumbentSolution = sol;
         if( sol->getObjectiveFuntionValue() < globalBestIncumbentValue )
         {
            globalBestIncumbentValue = sol->getObjectiveFuntionValue();
            globalIncumbnetValueUpdateFlag = true;
         }
         return true;
      }
      else
      {
         return false;
      }
   }
   else
   {
      globalBestIncumbentSolution = sol;
      if( sol->getObjectiveFuntionValue() < globalBestIncumbentValue )
      {
         globalBestIncumbentValue = sol->getObjectiveFuntionValue();
         globalIncumbnetValueUpdateFlag = true;
      }
      return true;
   }
}

void
ParaSolver::waitMessageIfNecessary(
      )
{
   if( paraParams->getIntParamValue(NotificationSynchronization) == 0 ||
          !rampUp ||
         ( paraParams->getIntParamValue(NotificationSynchronization) == 1 && ( collectingMode || collectingManyNodes ) ) )
   {
      waitNotificationIdMessage();
   }
}

void
ParaSolver::waitNotificationIdMessage(
      )
{
   double inIdleTime = paraTimer->getElapsedTime();
   int tag;
   int status;
   waitingSpecificMessage = true;
   while( paraComm->waitSpecTagFromSpecSource(0, ParaUNSIGNED, TagNotificationId, &tag) != 0 )
   {
      idleTimeToWaitNotificationId += ( paraTimer->getElapsedTime() - inIdleTime );
      status = (this->*messageHandler[tag])(0, tag);
      if( status )
      {
         std::ostringstream s;
         s << "[ERROR RETURN form Message Hander]:" <<  __FILE__ <<  "] func = "
           << __func__ << ", line = " << __LINE__ << " - "
           << "process tag = " << tag << std::endl;
         abort();
      }
      inIdleTime = paraTimer->getElapsedTime();
   }
   idleTimeToWaitNotificationId += ( paraTimer->getElapsedTime() - inIdleTime );
   //
   // receive the notification Id message
   //
   status = (this->*messageHandler[tag])(0, tag);
   if( status )
   {
      std::ostringstream s;
      s << "[ERROR RETURN form Message Hander]:" <<  __FILE__ <<  "] func = "
        << __func__ << ", line = " << __LINE__ << " - "
        << "process tag = " << tag << std::endl;
      abort();
   }
   waitingSpecificMessage = false;
}

void
ParaSolver::waitAckCompletion(
      )
{
   double inIdleTime = paraTimer->getElapsedTime();
   int tag;
   int hstatus;
   waitingSpecificMessage = true;
   while( paraComm->waitSpecTagFromSpecSource(0, ParaBYTE, TagAckCompletion, &tag) != 0 )
   {
      idelTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - inIdleTime );
      hstatus = (this->*messageHandler[tag])(0, tag);
      if( hstatus )
      {
         std::ostringstream s;
         s << "[ERROR RETURN form Message Hander]:" <<  __FILE__ <<  "] func = "
           << __func__ << ", line = " << __LINE__ << " - "
           << "process tag = " << tag << std::endl;
         abort();
      }
      inIdleTime = paraTimer->getElapsedTime();
   }
   idelTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - inIdleTime );
   //
   // receive the notification Id message
   //
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, 0, TagAckCompletion)
         );
   waitingSpecificMessage = false;
}
