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

/**@file    paraLoadCoordinator.cpp
 * @brief   Load coordinator.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifdef _COMM_PTH
#include <unistd.h>
#endif
#include <unistd.h>
#include <math.h>
#include <ctime>
#include <cfloat>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <climits>
#include <algorithm>
#include <iomanip>
#include "gzstream.h"
#include "paraLoadCoordinator.h"
#include "paraInitialStat.h"

using namespace UG;

ParaLoadCoordinator::ParaLoadCoordinator(
      ParaComm *inComm,
      ParaParamSet *inParaParamSet,
      ParaInitiator *inParaInitiator,
      bool *inRacingSolversExist,
      ParaTimer *inParaTimer,
      ParaDeterministicTimer *inParaDetTimer
      ) :
      globalSubtreeIdGen(0),
      paraParamSet(inParaParamSet),
      paraInitiator(inParaInitiator),
      racingSolversExist(inRacingSolversExist),
      restarted(false),
      initialNodesGenerated(false),
      runningPhase(RampUpPhase),
      firstCollectingModeState(-1),
      isCollectingModeRestarted(false),
      breakingSolverId(-1),
      nReplaceToBetterNode(0),
      computationIsInterrupted(false),
      interruptedFromControlTerminal(false),
      hardTimeLimitIsReached(false),
      winnerSolverNodesCollected(false),
      interruptIsRequested(false),
      paraRacingSolverPool(0),
      // bestDualBoundValueInInterruptedRacingSolvers(-DBL_MAX),
      nSolvedInInterruptedRacingSolvers(-1),
      nNodesLeftInInterruptedRacingSolvers(-1),
      previousCheckpointTime(0.0),
      statEmptyNodePoolTime(DBL_MAX),
      eps(0.0),
      racingWinner(-1),
      minDepthInWinnerSolverNodes(INT_MAX),
      maxDepthInWinnerSolverNodes(-1),
      racingWinnerParams(0),
      racingTermination(false),
      nSolvedRacingTermination(0),
      merging(false),
      varIndexTable(0),
      mergeInfoHead(0),
      mergeInfoTail(0),
      nTerminated(0),
      paraTimer(inParaTimer),
      paraDetTimer(inParaDetTimer),
      previousTabularOutputTime(0.0),
      osLogSolvingStatus(0),
      osLogNodesTransfer(0),
      osTabularSolvingStatus(0),
      osLogSubtreeInfo(0),
      osStatisticsFinalRun(0),
      osStatisticsRacingRampUp(0)
{

   paraComm = inComm;
   lcts.rank = paraComm->getRank();

   /** register message handlers */
   for( int i = 0; i < N_TAGS; i++ )
   {
      messageHandler[i] = 0;
   }
   messageHandler[TagNode] = &UG::ParaLoadCoordinator::processTagNode;
   messageHandler[TagSolution] = &UG::ParaLoadCoordinator::processTagSolution;
   messageHandler[TagSolverState] = &UG::ParaLoadCoordinator::processTagSolverState;
   messageHandler[TagCompletionOfCalculation] = &UG::ParaLoadCoordinator::processTagCompletionOfCalculation;
   messageHandler[TagAnotherNodeRequest] = &UG::ParaLoadCoordinator::processTagAnotherNodeRequest;
   messageHandler[TagTerminated] = &UG::ParaLoadCoordinator::processTagTerminated;
   messageHandler[TagHardTimeLimit] = &UG::ParaLoadCoordinator::processTagHardTimeLimit;
   messageHandler[TagInitialStat] = &UG::ParaLoadCoordinator::processTagInitialStat;
   if( paraParamSet->getBoolParamValue(Deterministic) )
   {
      messageHandler[TagToken] = &UG::ParaLoadCoordinator::processTagToken;
   }

   /* set up status log and transfer log */
   logSolvingStatusFlag = paraParamSet->getBoolParamValue(LogSolvingStatus);
   if( logSolvingStatusFlag )
   {
      std::ostringstream s;
      s << paraParamSet->getStringParamValue(LogSolvingStatusFilePath)
      << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".status";
      ofsLogSolvingStatus.open(s.str().c_str(), std::ios::app );
      if( !ofsLogSolvingStatus )
      {
         std::cout << "Solving status log file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osLogSolvingStatus = &ofsLogSolvingStatus;

      /* if initial solution is given, output the primal value */
      if( paraInitiator->getGlobalBestIncumbentSolution() )
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " LC" << " INITIAL_PRIMAL_VALUE "
         << paraInitiator->convertToExternalValue(
               paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue()
               ) << std::endl;
      }
   }

   logSubtreeInfoFlag = paraParamSet->getBoolParamValue(LogSubtreeInfo);
   if( logSubtreeInfoFlag )
   {
      std::ostringstream s;
      s << paraParamSet->getStringParamValue(LogSolvingStatusFilePath)
      << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".treelog";
      ofsLogSubtreeInfo.open(s.str().c_str(), std::ios::app );
      if( !ofsLogSubtreeInfo )
      {
         std::cout << "Sub tree info. log file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osLogSubtreeInfo = &ofsLogSubtreeInfo;
   }

   logNodesTransferFlag = paraParamSet->getBoolParamValue(LogNodesTransfer);
   if( logNodesTransferFlag )
   {
      std::ostringstream s;
      s << paraParamSet->getStringParamValue(LogNodesTransferFilePath)
      << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".transfer";
      ofsLogNodesTransfer.open(s.str().c_str(), std::ios::app);
      if( !ofsLogNodesTransfer )
      {
         std::cout << "Node transfer log file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osLogNodesTransfer = &ofsLogNodesTransfer;
   }

   outputTabularSolvingStatusFlag = paraParamSet->getBoolParamValue(OutputTabularSolvingStatus);
   if( outputTabularSolvingStatusFlag || paraParamSet->getBoolParamValue(Quiet) )
   {
      if( paraParamSet->getBoolParamValue(Quiet) )
      {
         osTabularSolvingStatus = &std::cout;
      }
      else
      {
         std::ostringstream s;
         s << paraParamSet->getStringParamValue(LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << "_T.status";
         ofsTabularSolvingStatus.open(s.str().c_str(), std::ios::app );
         if( !ofsTabularSolvingStatus )
         {
            std::cout << "Tabular solving status file cannot open : file name = " << s.str() << std::endl;
            exit(1);
         }
         osTabularSolvingStatus = &ofsTabularSolvingStatus;
      }
      if( outputTabularSolvingStatusFlag )
      {
         // output title line 1
         *osTabularSolvingStatus << std::setw(1) << " ";
         *osTabularSolvingStatus << std::setw(8) << std::right << " ";
         *osTabularSolvingStatus << std::setw(15) << std::right << " ";
         *osTabularSolvingStatus << std::setw(12) << std::right << "Nodes";
         *osTabularSolvingStatus << std::setw(10) << std::right << "Active";
         *osTabularSolvingStatus << std::setw(17) << std::right << " ";
         *osTabularSolvingStatus << std::setw(17) << std::right << " ";
         *osTabularSolvingStatus << std::setw(10) << std::right << " ";
         *osTabularSolvingStatus << std::endl;
         // output title line 2
         *osTabularSolvingStatus << std::setw(1) << " ";
         *osTabularSolvingStatus << std::setw(8) << std::right << "Time";
         *osTabularSolvingStatus << std::setw(15) << std::right << "Nodes";
         *osTabularSolvingStatus << std::setw(12) << std::right << "Left";
         *osTabularSolvingStatus << std::setw(10) << std::right << "Solvers";
         *osTabularSolvingStatus << std::setw(17) << std::right << "Best Integer";
         *osTabularSolvingStatus << std::setw(17) << std::right << "Best Node";
         *osTabularSolvingStatus << std::setw(10) << std::right << "Gap";
         *osTabularSolvingStatus << std::setw(17) << std::right << "Best Node(S)";
         *osTabularSolvingStatus << std::setw(10) << std::right << "Gap(S)";
         *osTabularSolvingStatus << std::endl;
      }
   }

   if( !paraParamSet->getBoolParamValue(Quiet) )
   {
      //
      // open statistic files
      //
      std::ostringstream ssfr;
      ssfr << paraParamSet->getStringParamValue(LogSolvingStatusFilePath)
      << paraInitiator->getParaInstance()->getProbName() << "_statistics_final_LC" << paraComm->getRank();
      ofsStatisticsFinalRun.open(ssfr.str().c_str(), std::ios::app);
      if( !ofsStatisticsFinalRun )
      {
         std::cout << "Statistics file for final run cannot open : file name = " << ssfr.str() << std::endl;
         exit(1);
      }
      osStatisticsFinalRun = &ofsStatisticsFinalRun;

      if( paraParamSet->getIntParamValue(RampUpPhaseProcess) == 1 ||
            paraParamSet->getIntParamValue(RampUpPhaseProcess) == 2
            )  /** racing ramp-up */
      {
         std::ostringstream ssrru;
         ssrru << paraParamSet->getStringParamValue(LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_statistics_racing_LC" << paraComm->getRank();
         ofsStatisticsRacingRampUp.open(ssrru.str().c_str(), std::ios::app);
         if( !ofsStatisticsRacingRampUp )
         {
            std::cout << "Statistics file for racing ramp-up cannot open : file name = " << ssrru.str() << std::endl;
            exit(1);
         }
         osStatisticsRacingRampUp = &ofsStatisticsRacingRampUp;
      }
   }

   paraSolverPool = new ParaSolverPoolForMinimization(    // always minimization problem
                    paraParamSet->getRealParamValue(MultiplierForCollectingMode),
                    paraParamSet->getRealParamValue(BgapCollectingMode),
                    paraParamSet->getRealParamValue(MultiplierForBgapCollectingMode),
                    1,                // paraSolver origin rank
                    paraComm, paraParamSet);
   // always minimization problem
   paraNodePool = new ParaNodePoolForMinimization(paraParamSet->getRealParamValue(BgapCollectingMode));

   eps = paraInitiator->getEpsilon();

   lastCheckpointTimeStr[0] = ' ';
   lastCheckpointTimeStr[1] = '\0';

   if( paraParamSet->getIntParamValue(NSolverNodesStartBreaking) == 0 ||
         paraParamSet->getIntParamValue(NStopBreaking) == 0 )
   {
      isBreakingFinised = true;
   }
   else
   {
      isBreakingFinised = false;
   }

   if( paraParamSet->getBoolParamValue(NChangeIntoCollectingModeNSolvers) )
   {
      paraParamSet->setIntParamValue(NChangeIntoCollectingMode, std::max((paraSolverPool->getNSolvers()/2), 1) );
   }

   if( paraParamSet->getBoolParamValue(Deterministic) )
   {
      assert(paraDetTimer);
   }

}

ParaLoadCoordinator::~ParaLoadCoordinator(
      )
{
   if( nTerminated == 0 &&
         ( ( paraRacingSolverPool && paraRacingSolverPool->getNumActiveSolvers() > 0 ) ||
               paraSolverPool->getNumActiveSolvers() > 0 ||
               interruptIsRequested
               )
     )
   {
      terminateAllSolvers();
      runningPhase = TerminationPhase;
      int source;
      int tag;

      for(;;)
      {
         /*******************************************
          *  waiting for any message form anywhere  *
          *******************************************/
         double inIdleTime = paraTimer->getElapsedTime();
         paraComm->probe(&source, &tag);
         lcts.idleTime += ( paraTimer->getElapsedTime() - inIdleTime );
         if( messageHandler[tag] )
         {
            int status = (this->*messageHandler[tag])(source, tag);
            if( status )
            {
               std::ostringstream s;
               s << "[ERROR RETURN form termination message handler]:" <<  __FILE__ <<  "] func = "
                 << __func__ << ", line = " << __LINE__ << " - "
                 << "process tag = " << tag << std::endl;
               abort();
            }
         }
         else
         {
            THROW_LOGICAL_ERROR3( "No message handler for ", tag, " is not registered" );
         }
         if( nTerminated == paraSolverPool->getNSolvers() ) break;
      }
   }

   /** write final solution */
   paraInitiator->writeSolution("Final Solution");

   if( paraNodePool )
   {
      if( logNodesTransferFlag )
      {
         *osLogNodesTransfer << std::endl << "ParaLoadCoordinator: # received = " << lcts.nReceived
         << ", # sent = " << lcts.nSent << ", # sent immediately = " << lcts.nSentBackImmediately << ", # deleted = " << lcts.nDeletedInLc
         << ", # failed to send back = " << lcts.nFailedToSendBack
         << ", Max usage of node pool = " << paraNodePool->getMaxUsageOfPool() << std::endl;
         *osLogNodesTransfer << "# sent immediately ( another node ) = " << lcts.nSentBackImmediatelyAnotherNode
         << ", # # failed to send back ( another node ) = " << lcts.nFailedToSendBackAnotherNode << std::endl;
         if( !paraNodePool->isEmpty() ){
            *osLogNodesTransfer << "LoadCoodibator: NodePool in LoadCoordinatator is not empty: "
            << paraNodePool->getNumOfNodes()  << " nodes remained" << std::endl;
         }
      }
      lcts.nMaxUsageOfNodePool = paraNodePool->getMaxUsageOfPool();
      lcts.nNodesInNodePool = paraNodePool->getNumOfNodes();
   }

   // set final solving status
   if( initialNodesGenerated )
   {
      paraInitiator->setFinalSolverStatus(InitialNodesGenerated);
   }
   else
   {
      if( hardTimeLimitIsReached )
      {
         paraInitiator->setFinalSolverStatus(HardTimeLimitIsReached);
      }
      else
      {
         if( interruptedFromControlTerminal ||
               (!racingTermination && computationIsInterrupted ) )
         {
            paraInitiator->setFinalSolverStatus(ComputingWasInterrupted);
         }
         else
         {
            paraInitiator->setFinalSolverStatus(ProblemWasSolved);
            if( outputTabularSolvingStatusFlag )
            {
               outputTabularSolvingStatus(' ');
            }
         }
      }
   }


   // set number of nodes solved and final dual bound value
   if( initialNodesGenerated )
   {
      paraInitiator->setNumberOfNodesSolved( std::max(1LL, paraSolverPool->getTotalNodesSolved()) );
      paraInitiator->setDualBound(paraNodePool->getBestDualBoundValue());
   }
   else
   {
      if( racingTermination )
      {
         if( ( interruptedFromControlTerminal || hardTimeLimitIsReached ) && paraRacingSolverPool )
         {
            paraInitiator->setNumberOfNodesSolved(paraRacingSolverPool->getNnodesSolvedInBestSolver());
         }
         else
         {
            paraInitiator->setNumberOfNodesSolved(nSolvedRacingTermination);
         }
         paraInitiator->setDualBound(lcts.globalBestDualBoundValue);
      }
      else
      {
         if( !computationIsInterrupted && !hardTimeLimitIsReached )
         {
            paraInitiator->setNumberOfNodesSolved( std::max(1LL, paraSolverPool->getTotalNodesSolved()) );
            paraInitiator->setDualBound(lcts.globalBestDualBoundValue);
         }
         else
         {
            if( isRacingStage() )
            {
               if( paraRacingSolverPool )
               {
                  paraInitiator->setDualBound(lcts.globalBestDualBoundValue);
                  paraInitiator->setNumberOfNodesSolved(  paraRacingSolverPool->getNnodesSolvedInBestSolver() );
               }
               else
               {
                  THROW_LOGICAL_ERROR1("Comutation is interrupted in racing stage, but no paraRacingSolverPool");
               }
            }
            else
            {
               paraInitiator->setNumberOfNodesSolved( std::max(1LL, paraSolverPool->getTotalNodesSolved()) );
               paraInitiator->setDualBound(lcts.globalBestDualBoundValue);
            }
         }
      }
   }

   paraInitiator->outputFinalSolverStatistics( osTabularSolvingStatus, paraTimer->getElapsedTime() );

   if( racingWinnerParams )
   {
      delete racingWinnerParams;
   }

   if( paraSolverPool && osStatisticsFinalRun )
   {
      lcts.nNodesLeftInAllSolvers = paraSolverPool->getNnodesInSolvers();
      *osStatisticsFinalRun << "######### The number of nodes solved in all solvers: "
            << paraSolverPool->getTotalNodesSolved();
      if( paraSolverPool->getTotalNodesSolved() != paraSolverPool->getNnodesSolvedInSolvers() )
      {
         *osStatisticsFinalRun << " : "
               << paraSolverPool->getNnodesSolvedInSolvers();
      }
      *osStatisticsFinalRun << " #########" << std::endl;
   }
   if( paraSolverPool )  delete paraSolverPool;

   lcts.runningTime = paraTimer->getElapsedTime();
   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus << lcts.runningTime << " ParaLoadCoordinator_TERMINATED" << std::endl;
      if( paraRacingSolverPool )
      {
         *osLogSolvingStatus << lcts.runningTime << " "
               << paraRacingSolverPool->getNumActiveSolvers() << " active racing ramp-up solvers exist." << std::endl;
      }
   }
#ifdef _DEBUG_LB
   std::cout << lcts.runningTime << " ParaLoadCoordinator_TERMINATED" << std::endl;
   if( paraRacingSolverPool )
   {
      std::cout << lcts.runningTime << " "
            << paraRacingSolverPool->getNumActiveSolvers() << " active racing ramp-up solvers exist." << std::endl;
   }
#endif


   lcts.isCheckpointState = false;
   if( paraParamSet->getBoolParamValue(StatisticsToStdout) )
   {
      std::cout << lcts.toString() << std::endl;
   }

   if( osStatisticsFinalRun )
   {
      *osStatisticsFinalRun << lcts.toString();
   }

   if( paraRacingSolverPool && osStatisticsFinalRun )
   {
      *osStatisticsFinalRun << "***** <<< NOTE >>> "
            << paraRacingSolverPool->getNumActiveSolvers()
            << " active racing ramp-up solvers exist." << std::endl;
      *osStatisticsFinalRun << paraRacingSolverPool->getStrActiveSolerNumbers() << std::endl;
   }

   if( paraRacingSolverPool )
   {
      delete paraRacingSolverPool;
      *racingSolversExist = true;
   }

   if( osStatisticsFinalRun )
   {
      osStatisticsFinalRun->flush();
   }

   if( paraNodePool ) delete paraNodePool;
}

int
ParaLoadCoordinator::processTagNode(
      int source,
      int tag
      )
{
   ParaNode *paraNode = paraComm->createParaNode();
   paraNode->receive(paraComm, source);
   lcts.nReceived++;
   if( ( paraInitiator->getGlobalBestIncumbentSolution() &&
         paraNode->getDualBoundValue() < paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue() ) ||
         !( paraInitiator->getGlobalBestIncumbentSolution() ) )
   {
      paraNode->setGlobalSubtreeId(paraComm->getRank(), createNewGlobalSubtreeId());
      /** in the case that ParaNode received from LoadCoordinator is not implemented yet */
      ParaNode *ancestorNode = paraSolverPool->getCurrentNode(source);
      paraNode->setAncestor(
            new ParaNodeGenealogicalLocalPtr( ancestorNode->getNodeId(), ancestorNode ));
      ancestorNode->addDescendant(
            new ParaNodeGenealogicalLocalPtr( paraNode->getNodeId(), paraNode) );
      paraNodePool->insert(paraNode);
      if( merging )
      {
         addNodeToMergeNodeStructs(paraNode);
      }

      if( paraParamSet->getBoolParamValue(RacingStatBranching) &&
            source == racingWinner &&
            !winnerSolverNodesCollected )
      {
         if( minDepthInWinnerSolverNodes > paraNode->getDepth() )
         {
            minDepthInWinnerSolverNodes = paraNode->getDepth();
         }
      }

      if( !( paraSolverPool->getNumInactiveSolvers() > 0 && sendParaNodesToIdleSolvers() ) )
                                              // In racing stage, some solver could take time for solving a node
                                              // Therefore, some solver could stay idle so long time in paraSolverPool
                                              // However, node collecting has to terminated to sending switch out message
      {
         // paraNodePool->insert(paraNode);
         /** switch-out request might be delayed to reach the target solver.
          * This behavior affects load balancing sometimes.
          *  Therefore, if additional nodes are received, sends switch-out request again.
         //if( runningPhase != RampUpPhase && !(paraSolverPool->isInCollectingMode()) )
         //{
            // paraSolverPool->enforcedSwitchOutCollectingMode(source);
         //}
         ** I don't want to do this. So, I decided to wait message after each notification when LC in collecting mode if necessary */
         double globalBestDualBoundValue =
            std::max (
                  std::min( paraSolverPool->getGlobalBestDualBoundValue(), paraNodePool->getBestDualBoundValue() ),
                  lcts.globalBestDualBoundValue );
         if( runningPhase != RampUpPhase && paraSolverPool->isInCollectingMode() &&
              ( paraNodePool->getNumOfGoodNodes( globalBestDualBoundValue )
               > ( paraParamSet->getRealParamValue(MultiplierForCollectingMode) *
                     paraParamSet->getIntParamValue(NChangeIntoCollectingMode) )
                     )
                     )
         {
            paraSolverPool->switchOutCollectingMode();
            firstCollectingModeState = 1;
            isCollectingModeRestarted = false;
         }

      }
   }
   else
   {
      assert( !( EPSLT( paraNode->getDualBoundValue(), paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue(), eps) ) );
      delete paraNode;
      lcts.nDeletedInLc++;
      if( runningPhase != RampUpPhase && !(paraSolverPool->isInCollectingMode()) &&
            paraSolverPool->getNumInactiveSolvers() > ( paraSolverPool->getNSolvers() * 0.1 )  &&
            paraNodePool->isEmpty() )
      { // inactive solver exists but cannot send a ParaNode to it
         paraSolverPool->switchInCollectingMode(paraNodePool);
         if( firstCollectingModeState == -1 ) firstCollectingModeState = 0;
      }
   }

   PARA_COMM_CALL(
          paraComm->send( NULL, 0, ParaBYTE, source, TagNodeReceived);
          );

   return 0;
}

int
ParaLoadCoordinator::processTagSolution(
      int source,
      int tag
      )
{
   ParaSolution *sol = paraComm->createParaSolution();
   sol->receive(paraComm, source);

#ifdef _DEBUG_DET
   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus << paraTimer->getElapsedTime()
      << " S." << source << " R.SOL "
      << paraInitiator->convertToExternalValue(
            sol->getObjectiveFuntionValue()
            ) << std::endl;
   }
#endif

   if( updateSolution(sol) ) 
   {
      delete sol;
      sendIncumbentValue(source);
      if( logSolvingStatusFlag )
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " I.SOL "
         << paraInitiator->convertToExternalValue(
               paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue()
               ) << std::endl;
      }
#ifdef _DEBUG_LB
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " I.SOL "
      << paraInitiator->convertToExternalValue(
            paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue()
            ) << std::endl;
#endif
      /** output tabular solving status */
      if( outputTabularSolvingStatusFlag )
      {
         outputTabularSolvingStatus('*');
      }
      /* Do not have to remove ParaNodes from NodePool. It is checked and removed before sending them */
      /** save incumbent solution */
      char solutionFileNameTemp[256];
      char solutionFileName[256];
      if( paraParamSet->getBoolParamValue(Checkpoint) && paraComm->getRank() == 0  )
      {
         sprintf(solutionFileNameTemp,"%s%s_after_checkpointing_solution_t.gz",
               paraParamSet->getStringParamValue(CheckpointFilePath),
               paraInitiator->getParaInstance()->getProbName() );
         paraInitiator->writeCheckpointSolution(std::string(solutionFileNameTemp));
         sprintf(solutionFileName,"%s%s_after_checkpointing_solution.gz",
               paraParamSet->getStringParamValue(CheckpointFilePath),
               paraInitiator->getParaInstance()->getProbName() );
         if ( rename(solutionFileNameTemp, solutionFileName) )
         {
            std::cout << "after checkpointing solution file cannot be renamed: errno = " << strerror(errno) << std::endl;
            exit(1);
         }
      }
   }
   else
   {
      delete sol;
   }
   return 0;
}

int
ParaLoadCoordinator::processTagSolverState(
      int source,
      int tag
      )
{

   double globalBestDualBoundValue = -DBL_MAX;

   ParaSolverState *solverState = paraComm->createParaSolverState();
   solverState->receive(paraComm, source, tag);

   if( paraDetTimer
         && paraDetTimer->getElapsedTime() < solverState->getDeterministicTime() )

   {
      paraDetTimer->update( solverState->getDeterministicTime() - paraDetTimer->getElapsedTime() );
   }

   if( solverState->isRacingStage() )
   {
      globalBestDualBoundValue =
               std::min( paraSolverPool->getGlobalBestDualBoundValue(), paraNodePool->getBestDualBoundValue() );
      /** not update paraSolverPool. The solver is inactive in paraSolverPool */
      if( logSolvingStatusFlag )
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " | "
         << paraInitiator->convertToExternalValue(
               solverState->getSolverLocalBestDualBoundValue()
               );
         if( !paraInitiator->getGlobalBestIncumbentSolution() ||
               paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) > displayInfOverThisValue )
         {
            *osLogSolvingStatus << " ( Inf )";
         }
         else
         {
            *osLogSolvingStatus << " ( " << paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) * 100 << "% )";
         }
         *osLogSolvingStatus << " [ " << solverState->getNNodesLeft() << " ]";
         *osLogSolvingStatus << " ** G.B.: " << paraInitiator->convertToExternalValue(globalBestDualBoundValue);
         if( !paraInitiator->getGlobalBestIncumbentSolution() ||
               paraInitiator->getGap(globalBestDualBoundValue) > displayInfOverThisValue )
         {
            *osLogSolvingStatus << " ( Inf ) ";
         }
         else
         {
            *osLogSolvingStatus << " ( " << paraInitiator->getGap(globalBestDualBoundValue) * 100 << "% ) ";
         }
         *osLogSolvingStatus << "[ " << paraSolverPool->getNnodesInSolvers() << ", " << paraNodePool->getNumOfNodes()
         << " ( " <<  paraNodePool->getNumOfGoodNodes( globalBestDualBoundValue )
//          <<" ) ] ** RR" << std::endl;
         <<" ) ] ** RR " << solverState->getDeterministicTime() << std::endl;  // for debug

      }
#ifdef _DEBUG_LB
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " | "
      << paraInitiator->convertToExternalValue(
            solverState->getSolverLocalBestDualBoundValue()
            );
      if( !paraInitiator->getGlobalBestIncumbentSolution() ||
            paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) > displayInfOverThisValue )
      {
         std::cout << " ( Inf )";
      }
      else
      {
         std::cout << " ( " << paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) * 100 << "% )";
      }
      std::cout << " [ " << solverState->getNNodesLeft() << " ]";
      globalBestDualBoundValue =
               std::min( paraSolverPool->getGlobalBestDualBoundValue(), paraNodePool->getBestDualBoundValue() );
      std::cout << " ** G.B.: " << paraInitiator->convertToExternalValue(globalBestDualBoundValue);
      if( !paraInitiator->getGlobalBestIncumbentSolution() ||
            paraInitiator->getGap(globalBestDualBoundValue) > displayInfOverThisValue )
      {
         std::cout << " ( Inf ) ";
      }
      else
      {
         std::cout << " ( " << paraInitiator->getGap(globalBestDualBoundValue) * 100 << "% ) ";
      }
      std::cout << "[ " << paraSolverPool->getNnodesInSolvers() << ", " << paraNodePool->getNumOfNodes()
      << " ( " <<  paraNodePool->getNumOfGoodNodes( globalBestDualBoundValue )
      <<" ) ] ** RR" << std::endl;
#endif
   }
   else
   {
      paraSolverPool->updateSolverStatus(source,
                                      solverState->getNNodesSolved(),
                                      solverState->getNNodesLeft(),
                                      solverState->getSolverLocalBestDualBoundValue(),
                                      paraNodePool
                                      );
      globalBestDualBoundValue =
         std::max (
               std::min( paraSolverPool->getGlobalBestDualBoundValue(), paraNodePool->getBestDualBoundValue() ),
               lcts.globalBestDualBoundValue );
      if( logSolvingStatusFlag )
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " | "
         << paraInitiator->convertToExternalValue(
               solverState->getSolverLocalBestDualBoundValue()
               );
         if( !paraInitiator->getGlobalBestIncumbentSolution() ||
               paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) > displayInfOverThisValue )
         {
            *osLogSolvingStatus << " ( Inf )";
         }
         else
         {
            *osLogSolvingStatus << " ( " << paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) * 100 << "% )";
         }
         *osLogSolvingStatus << " [ " << solverState->getNNodesLeft() << " ]";
         *osLogSolvingStatus << " ** G.B.: " << paraInitiator->convertToExternalValue(globalBestDualBoundValue);
         if( !paraInitiator->getGlobalBestIncumbentSolution() ||
               paraInitiator->getGap(globalBestDualBoundValue) > displayInfOverThisValue )
         {
            *osLogSolvingStatus << " ( Inf ) ";
         }
         else
         {
            *osLogSolvingStatus << " ( " << paraInitiator->getGap(globalBestDualBoundValue) * 100 << "% ) ";
         }
         *osLogSolvingStatus << "[ " << paraSolverPool->getNnodesInSolvers() << ", " << paraNodePool->getNumOfNodes()
         << " ( " <<  paraNodePool->getNumOfGoodNodes( globalBestDualBoundValue )
         <<" ) ] **";
         if( runningPhase == RampUpPhase ) *osLogSolvingStatus << " R";
         if( paraSolverPool->isInCollectingMode() )
         {
            *osLogSolvingStatus << " C";
            if( paraSolverPool->isSolverInCollectingMode(source) )
            {
               *osLogSolvingStatus << " 1";
            }
            else
            {
               *osLogSolvingStatus << " 0";
            }

         }
         *osLogSolvingStatus << " " << solverState->getDeterministicTime();   // for debug
         *osLogSolvingStatus << std::endl;
      }

      ParaNode *node = paraSolverPool->getCurrentNode(source);
      if( node->getMergeNodeInfo() != 0 && solverState->getNNodesSolved() > 2 )  // I stand on the safety side. we can write "> 1"
      {
         mergeNodes(node);
      }

#ifdef _DEBUG_LB
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " | "
      << paraInitiator->convertToExternalValue(
            solverState->getSolverLocalBestDualBoundValue()
            );
      if( !paraInitiator->getGlobalBestIncumbentSolution() ||
            paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) > displayInfOverThisValue )
      {
         std::cout << " ( Inf )";
      }
      else
      {
         std::cout << " ( " << paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) * 100 << "% )";
      }
      std::cout << " [ " << solverState->getNNodesLeft() << " ]";
      globalBestDualBoundValue =
         std::max (
               std::min( paraSolverPool->getGlobalBestDualBoundValue(), paraNodePool->getBestDualBoundValue() ),
               lcts.globalBestDualBoundValue );
      std::cout << " ** G.B.: " << paraInitiator->convertToExternalValue(globalBestDualBoundValue);
      if( !paraInitiator->getGlobalBestIncumbentSolution() ||
            paraInitiator->getGap(globalBestDualBoundValue) > displayInfOverThisValue )
      {
         std::cout << " ( Inf ) ";
      }
      else
      {
         std::cout << " ( " << paraInitiator->getGap(globalBestDualBoundValue) * 100 << "% ) ";
      }
      std::cout << "[ " << paraSolverPool->getNnodesInSolvers() << ", " << paraNodePool->getNumOfNodes()
      << " ( " <<  paraNodePool->getNumOfGoodNodes( globalBestDualBoundValue )
      <<" ) ] **";
      if( runningPhase == RampUpPhase ) std::cout << " R";
      if( paraSolverPool->isInCollectingMode() )
      {
         std::cout << " C";
         if( paraSolverPool->isSolverInCollectingMode(source) )
         {
            std::cout << " 1";
         }
         else
         {
            std::cout << " 0";
         }

      }
      std::cout << std::endl;
#endif

   }

   /** the following should be before noticationId back to the source solver */
   if( paraParamSet->getBoolParamValue(DistributeBestPrimalSolution) )
   {
      if( paraInitiator->getGlobalBestIncumbentSolution() &&
            paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue()
            < solverState->getGlobalBestPrimalBoundValue() )
      {
         paraInitiator->getGlobalBestIncumbentSolution()->send(paraComm, source);
      }
   }

   double lcBestDualBoundValue = paraNodePool->getBestDualBoundValue();
   PARA_COMM_CALL(
         paraComm->send( &lcBestDualBoundValue, 1, ParaDOUBLE, source, TagLCBestBoundValue)
         );
   unsigned int notificationId = solverState->getNotificaionId();
   PARA_COMM_CALL(
         paraComm->send( &notificationId, 1, ParaUNSIGNED, source, TagNotificationId)
         );

   if( paraParamSet->getBoolParamValue(InitialNodesGeneration) &&
         ( paraSolverPool->getNnodesInSolvers() + paraNodePool->getNumOfNodes() ) >= paraParamSet->getIntParamValue(NumberOfInitialNodes) )
   {
      for(int i = 1; i <= paraSolverPool->getNSolvers(); i++ )
      {
         int nCollect = -1;
         if( paraSolverPool->isSolverActive(i) )
         {
            PARA_COMM_CALL(
                  paraComm->send( &nCollect, 1, ParaINT, i, TagCollectAllNodes )
            );
         }
      }
      initialNodesGenerated = true;
   }
   else
   {
      if( runningPhase != RampUpPhase  )
      {
         if( paraNodePool->getNumOfGoodNodes(
               globalBestDualBoundValue
               ) > 0 )
         {
            statEmptyNodePoolTime = DBL_MAX;
         }
         else
         {
            if( paraParamSet->getBoolParamValue(Deterministic) )
            {
               if( ( paraDetTimer->getElapsedTime() - statEmptyNodePoolTime ) < 0 )
               {
                  statEmptyNodePoolTime = paraDetTimer->getElapsedTime();
               }
            }
            else
            {
               if( ( paraTimer->getElapsedTime() - statEmptyNodePoolTime ) < 0 )
               {
                  statEmptyNodePoolTime = paraTimer->getElapsedTime();
               }
            }

         }

         if( ( paraParamSet->getBoolParamValue(Deterministic) &&
               ( paraDetTimer->getElapsedTime() - statEmptyNodePoolTime ) > paraParamSet->getRealParamValue(TimeToIncreaseCMS) ) ||
               ( !paraParamSet->getBoolParamValue(Deterministic) &&
               ( paraTimer->getElapsedTime() - statEmptyNodePoolTime ) > paraParamSet->getRealParamValue(TimeToIncreaseCMS) ) )
         {
            if( !isCollectingModeRestarted )
            {
               paraSolverPool->switchOutCollectingMode();
               paraSolverPool->switchInCollectingMode(paraNodePool);
               if( logSolvingStatusFlag )
               {
                  *osLogSolvingStatus << paraTimer->getElapsedTime()
                  << " Collecting mode is restarted with " << paraSolverPool->getNLimitCollectingModeSolvers()
                  << std::endl;
               }
               isCollectingModeRestarted = true;
            }
            else
            {
               if( paraNodePool->getNumOfNodes() == 0 && paraSolverPool->canIncreaseLimitNLimitCollectingModeSolvers() )
                  // ramp-up may collect nodes having not so good nodes. As long as nodes exist, the limit number should not be increased.
               {
                  paraSolverPool->switchOutCollectingMode();
                  paraSolverPool->incNLimitCollectingModeSolvers();
                  paraSolverPool->switchInCollectingMode(paraNodePool);
                  if( logSolvingStatusFlag )
                  {
                     *osLogSolvingStatus << paraTimer->getElapsedTime()
                     << " Limit number of collecting mode solvers extends to " << paraSolverPool->getNLimitCollectingModeSolvers()
                     << std::endl;
                     if( outputTabularSolvingStatusFlag )
                     {
                        *osTabularSolvingStatus <<
                              "Limit number of collecting mode solvers extends to " <<
                              paraSolverPool->getNLimitCollectingModeSolvers() <<
                              " after " << paraTimer->getElapsedTime() << " seconds." << std::endl;
                     }
   #ifdef _DEBUG_LB
                     std::cout << paraTimer->getElapsedTime()
                           << " Limit number of collecting mode solvers extends to " << paraSolverPool->getNLimitCollectingModeSolvers()
                           << std::endl;
   #endif
                  }
                  isCollectingModeRestarted = false;
               }
               else   // cannot increase the number of collecting mode solvers
               {
                  paraSolverPool->switchOutCollectingMode();
                  paraSolverPool->switchInCollectingMode(paraNodePool);
                  if( logSolvingStatusFlag )
                  {
                     *osLogSolvingStatus << paraTimer->getElapsedTime()
                     << " Collecting mode is restarted with " << paraSolverPool->getNLimitCollectingModeSolvers()
                     << std::endl;
                  }
                  isCollectingModeRestarted = true;
               }
            }
            statEmptyNodePoolTime = DBL_MAX;
         }
      }

      if( !( solverState->isRacingStage() ) &&
            runningPhase != RampUpPhase && !(paraSolverPool->isInCollectingMode()) &&
            ( paraNodePool->getNumOfGoodNodes( globalBestDualBoundValue )
            < paraParamSet->getIntParamValue(NChangeIntoCollectingMode) ) )
      {
         paraSolverPool->switchInCollectingMode(paraNodePool);
         if( firstCollectingModeState == -1 ) firstCollectingModeState = 0;
      }

      if( !isBreakingFinised )
      {
         if( (!solverState->isRacingStage()) && runningPhase == NormalRunningPhase )
         {
            if( paraNodePool->getNumOfNodes() > paraParamSet->getIntParamValue(NStopBreaking) ||
                  paraSolverPool->getNnodesInSolvers() < paraParamSet->getIntParamValue(NStopBreaking) )
            {
               isBreakingFinised = true;
            }
            else
            {
               if( breakingSolverId == -1 )
               {
                  if( paraSolverPool->getNumOfNodesLeftInBestSolver()
                        > paraParamSet->getIntParamValue(NSolverNodesStartBreaking) )
                  {
                     breakingSolverId = paraSolverPool->getBestSolver();
                     assert( breakingSolverId != -1 );
                     double targetBound = ( paraSolverPool->getGlobalBestDualBoundValue()*
                           paraParamSet->getRealParamValue(MultiplierForBreakingTargetBound) );
                     int nLimitTransfer = paraParamSet->getIntParamValue(NTransferLimitForBreaking);
                     PARA_COMM_CALL(
                           paraComm->send( &targetBound, 1, ParaDOUBLE, breakingSolverId, TagBreaking )
                     );
                     PARA_COMM_CALL(
                           paraComm->send( &nLimitTransfer, 1, ParaINT, breakingSolverId, TagBreaking )
                     );
                  }
                  else
                  {
                     if( ( ( paraSolverPool->getGlobalBestDualBoundValue()
                           + paraParamSet->getRealParamValue(ABgapForSwitchingToBestSolver)*3 ) >
                           solverState->getSolverLocalBestDualBoundValue() ) &&
                                 paraSolverPool->getNnodesInSolvers() >
                                 std::max(paraParamSet->getIntParamValue(NStopBreaking)*2, paraSolverPool->getNSolvers()*2 ) &&
                         solverState->getNNodesLeft() >  ( paraSolverPool->getNnodesInSolvers()*0.5 ) )
                     {
                        breakingSolverId = source;
                        double targetBound = ( solverState->getSolverLocalBestDualBoundValue()*
                              paraParamSet->getRealParamValue(MultiplierForBreakingTargetBound) );
                        int nLimitTransfer = paraParamSet->getIntParamValue(NTransferLimitForBreaking);
                        PARA_COMM_CALL(
                              paraComm->send( &targetBound, 1, ParaDOUBLE, breakingSolverId, TagBreaking )
                        );
                        PARA_COMM_CALL(
                              paraComm->send( &nLimitTransfer, 1, ParaINT, breakingSolverId, TagBreaking )
                        );
                     }
                  }
               }
            }
         }
      }
      else   // isBootstrapFinised
      {
         if( runningPhase == NormalRunningPhase &&
               paraNodePool->getNumOfNodes() < paraParamSet->getIntParamValue(NStopBreaking) &&
               paraSolverPool->getNnodesInSolvers()
               > std::max(paraParamSet->getIntParamValue(NStopBreaking)*2, paraSolverPool->getNSolvers()*2 ) )
         {
            // break again. several solvers can be in breaking situation. That is, braking message can be sent to breaking solver
            isBreakingFinised = false;
            breakingSolverId = -1;
         }
      }
   }

   if( lcts.globalBestDualBoundValue < globalBestDualBoundValue )
   {
      lcts.globalBestDualBoundValue = globalBestDualBoundValue;
      lcts.externalGlobalBestDualBoundValue = paraInitiator->convertToExternalValue(globalBestDualBoundValue);
   }

   delete solverState;
   return 0;
}

int
ParaLoadCoordinator::processTagCompletionOfCalculation(
      int source,
      int tag
      )
{
   ParaCalculationState *calcState = paraComm->createParaCalculationState();
   calcState->receive(paraComm, source, tag);
   if( paraRacingSolverPool && paraRacingSolverPool->isActive(source) ) // racing root node termination
   {
      writeTransferLogInRacing(source, calcState);
   }
   else
   {
      writeTransferLog(source, calcState);
   }

   if( !winnerSolverNodesCollected &&
         racingWinner == source )
   {
      winnerSolverNodesCollected = true;
      if( merging )
      {
         generateMergeNodesCandidates();   // Anyway,merge nodes candidates have to be generated,
                                           // even if running with InitialNodesGeneration
         merging = false;
      }
      if( paraParamSet->getBoolParamValue(InitialNodesGeneration) &&
           paraNodePool->getNumOfNodes() >= paraParamSet->getIntParamValue(NumberOfInitialNodes) )
      {
         initialNodesGenerated = true;
      }
   }

   int calcTerminationState = calcState->getTerminationState();

   if( logSolvingStatusFlag )
   {
      switch ( calcTerminationState )
      {
      case CompTerminatedNormally:
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " >";
         break;
      }
      case CompTerminatedByAnotherNode:
      {
         /** When starts with racing ramp-up, solvers except winner should be terminated in this state */
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " >(INTERRUPTED_BY_ANOTHER_NODE)";
         break;
      }
      case CompTerminatedByInterruptRequest:
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " >(INTERRUPTED)";
         break;
      }
      case CompTerminatedInRacingStage:
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " >(TERMINATED_IN_RACING_STAGE)";
         break;
      }
      case CompInterruptedInRacingStage:
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " >(INTERRUPTED_IN_RACING_STAGE)";
         break;
      }
      case CompInterruptedInMerging:
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " >(INTERRUPTED_IN_MERGING)";
         break;
      }
      default:
         THROW_LOGICAL_ERROR2("Invalid termination: termination state = ", calcState->getTerminationState() )
      }

      if( paraParamSet->getBoolParamValue(CheckEffectOfRootNodePreprocesses) &&
            calcState->getNSolvedWithNoPreprocesses() > 0 )
      {
         *osLogSolvingStatus << " SOLVED_AT_ROOT ( DEPTH = " << paraSolverPool->getCurrentNode(source)->getDepth()
         << ", Gap = " << paraInitiator->getGap(paraSolverPool->getCurrentNode(source)->getDualBoundValue()) * 100  << "%, TrueGap = "
         << paraInitiator->getGap(paraSolverPool->getCurrentNode(source)->getInitialDualBoundValue()) * 100  << "% ) [ "
         << calcState->getNSolvedWithNoPreprocesses() << " ]";
      }
      *osLogSolvingStatus << std::endl;
   }

#ifdef _DEBUG_LB
   switch ( calcTerminationState )
   {
   case CompTerminatedNormally:
   {
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " >";
      break;
   }
   case CompTerminatedByAnotherNode:
   {
      /** When starts with racing ramp-up, solvers except winner should be terminated in this state */
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " >(INTERRUPTED_BY_ANOTHER_NODE)";
      break;
   }
   case CompTerminatedByInterruptRequest:
   {
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " >(INTERRUPTED)";
      break;
   }
   case CompTerminatedInRacingStage:
   {
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " >(TERMINATED_IN_RACING_STAGE)";
      break;
   }
   case CompInterruptedInRacingStage:
   {
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " >(INTERRUPTED_IN_RACING_STAGE)";
      break;
   }
   default:
      THROW_LOGICAL_ERROR2("Invalid termination: termination state = ", calcState->getTerminationState() )
   }

   if( paraParamSet->getBoolParamValue(CheckEffectOfRootNodePreprocesses) &&
         calcState->getNSolvedWithNoPreprocesses() > 0 )
   {
      std::cout << " SOLVED_AT_ROOT ( DEPTH = " << paraSolverPool->getCurrentNode(source)->getDepth()
      << ", Gap = " << paraInitiator->getGap(paraSolverPool->getCurrentNode(source)->getDualBoundValue()) * 100  << "%, TrueGap = "
      << paraInitiator->getGap(paraSolverPool->getCurrentNode(source)->getInitialDualBoundValue()) * 100  << "% ) [ "
      << calcState->getNSolvedWithNoPreprocesses() << " ]";
   }
   std::cout << std::endl;
#endif

   switch ( calcTerminationState )
   {
   case CompTerminatedNormally:
   {
      writeSubtreeInfo(source, calcState);
      ParaNode *node = paraSolverPool->getCurrentNode(source);
      if( node->getMergeNodeInfo() )
      {
         mergeNodes(node);
      }
      if( (!paraRacingSolverPool) || (!paraRacingSolverPool->isActive(source) ) )  // RacingSolverPool is inactivated below
      {
         paraSolverPool->inactivateSolver(source, calcState->getNSolved(),paraNodePool);
         paraSolverPool->addTotalNodesSolved(calcState->getNSolved());
      }
      // DO NOT send ParaNode here!
      // Send ParaNode after solver termination state is received.
      break;
   }
   case CompInterruptedInRacingStage:
   {
      // DO NOT send ParaNode here!
      // Send ParaNode after solver termination state is received and RacingSolver is inactivated.
      // Do not have to update counters of ParaSolverPool.
      break;
   }
   case CompTerminatedByAnotherNode:
   {
      /** in this case the following two numbers should be different */
      /** # Total > # Solved */
      paraSolverPool->addNumNodesSolved(calcState->getNSolved());
      paraSolverPool->addTotalNodesSolved(calcState->getNSolved());
      break;
   }
   case CompTerminatedByInterruptRequest:
   {
      if( (!paraRacingSolverPool) || (!paraRacingSolverPool->isActive(source) ) )  // RacingSolverPool is inactivated below
      {
         // paraRacingSolverPool entry is inactivated, when it receives ParaSolverTerminationState message in below.
         ParaNode *solvingNode = paraSolverPool->extractCurrentNodeAndInactivate(source, paraNodePool);
         paraNodePool->insert(solvingNode);
      }
      break;
   }
   case CompTerminatedInRacingStage:
   {
      racingTermination = true; // even if interruptIsRequested, 
                                // solver should have been terminated before receiveing it
      if( osStatisticsRacingRampUp )
      {
         *osStatisticsRacingRampUp << "######### Solver Rank = " <<
               source << " is terminated in racing stage #########" << std::endl;
      }
      nSolvedRacingTermination = calcState->getNSolved();
      if( !interruptedFromControlTerminal && !hardTimeLimitIsReached &&
            !computationIsInterrupted
            && paraInitiator->getGlobalBestIncumbentSolution() )
      {
         lcts.globalBestDualBoundValue = paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue();
         lcts.externalGlobalBestDualBoundValue = paraInitiator->convertToExternalValue(lcts.globalBestDualBoundValue);
      }
      break;
   }
   case CompInterruptedInMerging:
   {
      ParaNode *solvingNode = paraSolverPool->extractCurrentNodeAndInactivate(source, paraNodePool);
      assert(solvingNode->getMergeNodeInfo());
      regenerateMergeNodesCandidates(solvingNode);
      paraNodePool->insert(solvingNode);
      break;
   }
   default:
      THROW_LOGICAL_ERROR2("Invalid termination: termination state = ", calcState->getTerminationState() )
   }

   delete calcState;

   ParaSolverTerminationState *termState = paraComm->createParaSolverTerminationState();
   termState->receive(paraComm, source, TagTerminated);

   if( paraDetTimer )
   {
      if( paraDetTimer->getElapsedTime() < termState->getDeterministicTime() )
      {
         paraDetTimer->update( termState->getDeterministicTime() - paraDetTimer->getElapsedTime() );
      }
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, source, TagAckCompletion )
      );
   }

   switch( termState->getInterruptedMode() )
   {
   case 2: /** checkpoint; This is normal termination */
   {
      /** in order to save termination status to check point file, keep this information to solver pool */
      paraSolverPool->setTermState(source, termState);
      // don't delete termState! it is saved in paraSolverPool
      if( runningPhase != TerminationPhase )
      {
         if( paraNodePool->isEmpty() )
         {
            lcts.nFailedToSendBack++;
            if( runningPhase != RampUpPhase && !(paraSolverPool->isInCollectingMode()) )
            {
               paraSolverPool->switchInCollectingMode(paraNodePool);
               if( firstCollectingModeState == -1 ) firstCollectingModeState = 0;
            }
         }
         else
         {
            if( sendParaNodesToIdleSolvers() )
            {
               lcts.nSentBackImmediately++;
            }
         }
      }
      break;
   }
   case 3: /** racing ramp-up */
   {
      if( osStatisticsRacingRampUp )
      {
         *osStatisticsRacingRampUp << termState->toString();
         osStatisticsRacingRampUp->flush();
      }
      // nTerminated++;      We should not count this, We should always send Term from LC!
      inactivateRacingSolverPool(source);
      /*  Anyway, incumbent value is sent to all Solvers except its generator. In such case, this is not necessary.
      double globalBestIncumbentValue = paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue();
      PARA_COMM_CALL(
            paraComm->send( &globalBestIncumbentValue, 1, ParaDOUBLE, source, TagIncumbentValue )
      );
      */
      if( runningPhase != TerminationPhase )
      {
         if( paraNodePool->isEmpty() )
         {
            lcts.nFailedToSendBack++;
            if( runningPhase != RampUpPhase && !(paraSolverPool->isInCollectingMode()) )
            {
               paraSolverPool->switchInCollectingMode(paraNodePool);
               if( firstCollectingModeState == -1 ) firstCollectingModeState = 0;
            }
         }
         else
         {
            if( sendParaNodesToIdleSolvers() )
            {
               lcts.nSentBackImmediately++;
            }
         }
      }
      delete termState;
      break;
   }
   default:  /** unexpected mode */
      THROW_LOGICAL_ERROR4("Unexpected termination state received from rank = ", source,
            ", interrupted mode = ", termState->getInterruptedMode());
   }

   if( paraParamSet->getBoolParamValue(Quiet) && racingTermination )
   {
      /** in this case, do not have to wait statistical information from the other solvers */
      nTerminated = 1;
      delete this;
#ifdef _COMM_PTH
      _exit(0);
#else
      exit(0);
#endif
   }

   if( source == breakingSolverId )
   {
      breakingSolverId = -1;
      isBreakingFinised = false;
   }
   return 0;
}

int
ParaLoadCoordinator::processTagAnotherNodeRequest(
      int source,
      int tag
      )
{
   double bestDualBoundValue;
   PARA_COMM_CALL(
         paraComm->receive( &bestDualBoundValue, 1, ParaDOUBLE, source, TagAnotherNodeRequest)
         );
   if( paraNodePool->isEmpty() || paraSolverPool->currentSolvingNodehaeDescendant(source) )
   {
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, source, TagNoNodes)
            );
      lcts.nFailedToSendBackAnotherNode++;
   }
   else
   {
      ParaNode *paraNode = 0;
      while( !paraNodePool->isEmpty() )
      {
         paraNode = paraNodePool->extractNode();
         if( !paraNode ) break;
         if( ( paraInitiator->getGlobalBestIncumbentSolution() &&
               paraNode->getDualBoundValue() < paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue() ) ||
               !( paraInitiator->getGlobalBestIncumbentSolution() ) )
         {
            break;
         }
         else
         {
            delete paraNode;
            paraNode = 0;
            lcts.nDeletedInLc++;
         }
      }
      if( paraNode )
      {
         if( ( paraSolverPool->getDualBoundValue(source) - paraNode->getDualBoundValue()) > 0.0 &&
               ( REALABS( ( paraSolverPool->getDualBoundValue(source) - paraNode->getDualBoundValue() )
                     / paraNode->getDualBoundValue() ) > paraParamSet->getRealParamValue(BgapStopSolvingMode) ) )
         {
            paraSolverPool->sendSwitchOutCollectingModeIfNecessary(source);
            paraNode->send(paraComm, source);
            lcts.nSent++;
            lcts.nSentBackImmediatelyAnotherNode++;
            ParaNode *solvingNode = paraSolverPool->extractCurrentNodeAndInactivate(source, paraNodePool);
            solvingNode->setDualBoundValue(bestDualBoundValue);
            solvingNode->setInitialDualBoundValue(bestDualBoundValue);
            paraNodePool->insert(solvingNode);
            paraSolverPool->activateSolver(source, paraNode);
            writeTransferLog(source);
            if( logSolvingStatusFlag )
            {
               *osLogSolvingStatus << paraTimer->getElapsedTime()
               << " S." << source << " <(ANOTHER_NODE) "
               << paraInitiator->convertToExternalValue(
                     paraNode->getDualBoundValue() );
               if( paraInitiator->getGlobalBestIncumbentSolution() )
               {
                  if( paraInitiator->getGap(paraNode->getDualBoundValue()) > displayInfOverThisValue )
                  {
                     *osLogSolvingStatus << " ( Inf )";
                  }
                  else
                  {
                     *osLogSolvingStatus << " ( " << paraInitiator->getGap(paraNode->getDualBoundValue()) * 100 << "% )";
                  }
               }
               *osLogSolvingStatus << std::endl;
            }
#ifdef _DEBUG_LB
            std::cout << paraTimer->getElapsedTime()
            << " S." << source << " <(ANOTHER_NODE) "
            << paraInitiator->convertToExternalValue(
                  paraNode->getDualBoundValue() );
            if( paraInitiator->getGlobalBestIncumbentSolution() )
            {
               if( paraInitiator->getGap(paraNode->getDualBoundValue()) > displayInfOverThisValue )
               {
                  std::cout << " ( Inf )";
               }
               else
               {
                  std::cout << " ( " << paraInitiator->getGap(paraNode->getDualBoundValue()) * 100 << "% )";
               }
            }
            std::cout << std::endl;
#endif
         }
         else
         {
            paraNodePool->insert(paraNode);
            PARA_COMM_CALL(
                  paraComm->send( NULL, 0, ParaBYTE, source, TagNoNodes)
                  );
            lcts.nFailedToSendBackAnotherNode++;
         }
      }
      else
      {
         PARA_COMM_CALL(
               paraComm->send( NULL, 0, ParaBYTE, source, TagNoNodes)
               );
         lcts.nFailedToSendBackAnotherNode++;
      }
   }
   return 0;
}

int
ParaLoadCoordinator::processTagTerminated(
      int source,
      int tag
      )
{
   ParaSolverTerminationState *paraSolverTerminationState = paraComm->createParaSolverTerminationState();
   paraSolverTerminationState->receive(paraComm, source, tag);

   if( paraDetTimer )
   {
      if( paraDetTimer->getElapsedTime() < paraSolverTerminationState->getDeterministicTime() )
      {
         paraDetTimer->update( paraSolverTerminationState->getDeterministicTime() - paraDetTimer->getElapsedTime() );
      }
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, source, TagAckCompletion )
      );
   }

   if( osStatisticsFinalRun )
   {
      *osStatisticsFinalRun << paraSolverTerminationState->toString();
      osStatisticsFinalRun->flush();
   }
   if( paraParamSet->getBoolParamValue(StatisticsToStdout) )
   {
      std::cout << paraSolverTerminationState->toString() << std::endl;
   }
   /* do not have to do this. All ParaSolvers should be inactivated.
   if( paraSolverPool->isSolverActive(source) )   // if TimeLimit is specified, solver can be active
   {
      paraSolverPool->inactivateSolver(source, 0, paraNodePool );   // 0: the number of solved should be counted in completion of calculation process
   }
   */
   if( (!racingTermination) && paraSolverTerminationState->getInterruptedMode() == 1 )
   {
      computationIsInterrupted = true;
   }
   nTerminated++;

   delete paraSolverTerminationState;

   return 0;
}

int
ParaLoadCoordinator::processTagHardTimeLimit(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, ParaBYTE, source, TagHardTimeLimit)
         );
   hardTimeLimitIsReached = true;
   return 0;
}

int
ParaLoadCoordinator::processTagInitialStat(
      int source,
      int tag
      )
{
   ParaInitialStat *initialStat = paraComm->createParaInitialStat();
   initialStat->receive(paraComm, source);
   if( maxDepthInWinnerSolverNodes < initialStat->getMaxDepth() )
   {
      maxDepthInWinnerSolverNodes = initialStat->getMaxDepth();
   }
   paraInitiator->accumulateInitialStat(initialStat);
   delete initialStat;
   return 0;
}

int
ParaLoadCoordinator::processTagToken(
      int source,
      int tag
      )
{
   int token[2];
   PARA_COMM_CALL(
         paraComm->receive( token, 2, ParaINT, source, TagToken)
         );
   PARA_COMM_CALL(
         paraComm->send( token, 2, ParaINT, token[0], TagToken )
   );

   paraComm->setToken(0, token);    // for debug

   return 0;
}

void
ParaLoadCoordinator::outputTabularSolvingStatus(
      char incumbent
      )
{
   *osTabularSolvingStatus << std::setw(1) << incumbent;
   *osTabularSolvingStatus << std::setw(8) << std::right << std::setprecision(0) << std::fixed << paraTimer->getElapsedTime();
   if( !restarted &&
         ( paraParamSet->getIntParamValue(RampUpPhaseProcess) == 1 ||
               paraParamSet->getIntParamValue(RampUpPhaseProcess) == 2 )
               && !racingWinnerParams )
   {
      /** racing ramp-up stage now */
      if( paraRacingSolverPool )
      {
         *osTabularSolvingStatus << std::setw(15) << std::right << paraRacingSolverPool->getNnodesSolvedInBestSolver();
         *osTabularSolvingStatus << std::setw(12) << std::right << paraRacingSolverPool->getNnodesLeftInBestSolver();
         *osTabularSolvingStatus << std::setw(10) << std::right << paraRacingSolverPool->getNumActiveSolvers();
         if( paraInitiator->getGlobalBestIncumbentSolution() )
         {
            *osTabularSolvingStatus << std::setw(17) << std::right << std::setprecision(4) <<
                  paraInitiator->convertToExternalValue(
                  paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue());
         }
         else
         {
            *osTabularSolvingStatus << std::setw(17) << std::right << "-";
         }
         if( paraRacingSolverPool->getNnodesSolvedInBestSolver() == 0 )
         {
            *osTabularSolvingStatus << std::setw(17) << std::right << "-";
         }
         else
         {
            if( EPSEQ( lcts.globalBestDualBoundValue,-DBL_MAX, paraInitiator->getEpsilon() ))
            {
               *osTabularSolvingStatus << std::setw(17) << std::right << "-";
            }
            else
            {
               *osTabularSolvingStatus << std::setw(17) << std::right << std::setprecision(4) <<
                     lcts.externalGlobalBestDualBoundValue;
            }
         }
      }
      else  // One of ParaSolvers terminates in racing stage
      {
         if( nSolvedRacingTermination > 0 )
         {
            *osTabularSolvingStatus << std::setw(15) << std::right << nSolvedRacingTermination;
            *osTabularSolvingStatus << std::setw(12) << std::right << 0;
            *osTabularSolvingStatus << std::setw(10) << std::right << 0;
            if( paraInitiator->getGlobalBestIncumbentSolution() )
            {
               *osTabularSolvingStatus << std::setw(17) << std::right << std::setprecision(4) <<
                     paraInitiator->convertToExternalValue(
                     paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue());
            }
            else
            {
               *osTabularSolvingStatus << std::setw(17) << std::right << "-";
            }
            *osTabularSolvingStatus << std::setw(17) << std::right << "-";
         }
         else   // should be interrupted
         {
            *osTabularSolvingStatus << std::setw(15) << std::right << nSolvedInInterruptedRacingSolvers;
            *osTabularSolvingStatus << std::setw(12) << std::right << nNodesLeftInInterruptedRacingSolvers;
            *osTabularSolvingStatus << std::setw(10) << std::right << 0;
            if( paraInitiator->getGlobalBestIncumbentSolution() )
            {
               *osTabularSolvingStatus << std::setw(17) << std::right << std::setprecision(4) <<
                     paraInitiator->convertToExternalValue(
                           paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue());
            }
            else
            {
               *osTabularSolvingStatus << std::setw(17) << std::right << "-";
            }
            *osTabularSolvingStatus << std::setw(17) << std::right << std::setprecision(4) <<
                  lcts.externalGlobalBestDualBoundValue;
         }
      }
      if( !paraInitiator->getGlobalBestIncumbentSolution() ||
            paraInitiator->getGap(lcts.globalBestDualBoundValue) > displayInfOverThisValue ||
            ( paraSolverPool->isActive() && paraSolverPool->getNumActiveSolvers() == 0 && paraNodePool->getNumOfNodes() == 0 )
            )
      {
         *osTabularSolvingStatus << std::setw(10) << std::right << "-";
      }
      else
      {
         *osTabularSolvingStatus << std::setw(9) << std::right << std::setprecision(2) <<
               paraInitiator->getGap(lcts.globalBestDualBoundValue) * 100 << "%";
      }
   }
   else
   {
      *osTabularSolvingStatus << std::setw(15) << std::right << paraSolverPool->getNnodesSolvedInSolvers();
      *osTabularSolvingStatus << std::setw(12) << std::right << ( paraSolverPool->getNnodesInSolvers() + paraNodePool->getNumOfNodes() );
      *osTabularSolvingStatus << std::setw(10) << std::right << paraSolverPool->getNumActiveSolvers();
      if( paraInitiator->getGlobalBestIncumbentSolution() )
      {
         *osTabularSolvingStatus << std::setw(17) << std::right << std::setprecision(4) <<
               paraInitiator->convertToExternalValue(
               paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue());
      }
      else
      {
         *osTabularSolvingStatus << std::setw(17) << std::right << "-";
      }

      if( paraSolverPool->getNnodesSolvedInSolvers() == 0 ||
            ( paraSolverPool->getNumActiveSolvers() == 0 && paraNodePool->getNumOfNodes() == 0 )
            )
      {
         *osTabularSolvingStatus << std::setw(17) << std::right << "-";
      }
      else
      {
         *osTabularSolvingStatus << std::setw(17) << std::right << std::setprecision(4) <<
               lcts.externalGlobalBestDualBoundValue;
      }
      if( !paraInitiator->getGlobalBestIncumbentSolution() ||
            paraInitiator->getGap(lcts.globalBestDualBoundValue) > displayInfOverThisValue ||
            ( paraSolverPool->isActive() && paraSolverPool->getNumActiveSolvers() == 0 && paraNodePool->getNumOfNodes() == 0 )
            )
      {
         *osTabularSolvingStatus << std::setw(10) << std::right << "-";
      }
      else
      {
         *osTabularSolvingStatus << std::setw(9) << std::right << std::setprecision(2) <<
               paraInitiator->getGap(lcts.globalBestDualBoundValue) * 100 << "%";
      }
      if( paraSolverPool->getNumActiveSolvers() > 0 )
      {
         if( paraSolverPool->getGlobalBestDualBoundValue() >= -1e+10 )
         {
            *osTabularSolvingStatus << std::setw(17) << std::right << std::setprecision(4) <<
                  paraInitiator->convertToExternalValue( paraSolverPool->getGlobalBestDualBoundValue() );
         }
         else
         {
            *osTabularSolvingStatus << std::setw(17) << std::right << "-";
         }
         if(  paraInitiator->getGap( paraSolverPool->getGlobalBestDualBoundValue() ) > displayInfOverThisValue )
         {
            *osTabularSolvingStatus << std::setw(10) << std::right << "-";
         }
         else
         {
            *osTabularSolvingStatus << std::setw(10) << std::right << std::setprecision(2) <<
                  paraInitiator->getGap( paraSolverPool->getGlobalBestDualBoundValue() ) * 100 << "%";
         }
      }
      else
      {
         *osTabularSolvingStatus << std::setw(17) << std::right << "-";
      }
   }
   *osTabularSolvingStatus << std::endl;
}

void
ParaLoadCoordinator::run(
      )
{

   int source;
   int tag;

   for(;;)
   {
      if( paraSolverPool->getNumActiveSolvers() == 0 )
      {
         if( paraNodePool->isEmpty() )                             // paraNodePool has to be checked
                                                                   // because node cannot send in a parameter settings
         {
            if( runningPhase != TerminationPhase )
            {
               if( !interruptedFromControlTerminal
                     && !computationIsInterrupted
                     && !hardTimeLimitIsReached
                     && paraInitiator->getGlobalBestIncumbentSolution() )
               {
                  lcts.globalBestDualBoundValue = paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue();
                  lcts.externalGlobalBestDualBoundValue = paraInitiator->convertToExternalValue(lcts.globalBestDualBoundValue);
               }
               /* No active solver exists */
               terminateAllSolvers();
               runningPhase = TerminationPhase;
            }
            else // runningPhase == TerminationPhase
            {
               if( ( paraRacingSolverPool &&
                     paraSolverPool->getNumInactiveSolvers() == (paraRacingSolverPool->getNumActiveSolvers() + nTerminated ) ) ||
                     paraSolverPool->getNumInactiveSolvers() == nTerminated  )
               {
                  break;
               }
            }
         }
         else
         {
            if( initialNodesGenerated  )
            {
               if( runningPhase != TerminationPhase )
               {
                  lcts.globalBestDualBoundValue = paraNodePool->getBestDualBoundValue();
                  lcts.externalGlobalBestDualBoundValue = paraInitiator->convertToExternalValue(lcts.globalBestDualBoundValue);
                  updateCheckpointFiles();
                  /* No active solver exists */
                  terminateAllSolvers();
                  runningPhase = TerminationPhase;
               }
               else  // runningPhase == TerminationPhase
               {
                  if( paraSolverPool->getNumInactiveSolvers() == nTerminated  )
                  {
                     break;
                  }
               }
            }
         }
      }

      if( !paraRacingSolverPool && paraSolverPool->getNumActiveSolvers() == 0 &&
               paraParamSet->getRealParamValue(TimeLimit) > 0.0 &&
               paraTimer->getElapsedTime() > paraParamSet->getRealParamValue(TimeLimit) )
      {
         hardTimeLimitIsReached = true;
         break;
      }

      if( racingTermination && !paraRacingSolverPool && paraSolverPool->getNumActiveSolvers() == 0 && paraNodePool->getNumOfNodes() == 1 )
      {
         /*
          * special timining problem
          *
          * 1113.58 S.4 I.SOL 0
          * 1113.58 S.3 is the racing winner! Selected strategy 2.
          * 1113.58 S.4 >(TERMINATED_IN_RACING_STAGE)
          *
          */
         break;
      }

      /*******************************************
       *  waiting for any message form anywhere  *
       *******************************************/
      double inIdleTime = paraTimer->getElapsedTime();
      paraComm->probe(&source, &tag);
      lcts.idleTime += ( paraTimer->getElapsedTime() - inIdleTime );
      if( messageHandler[tag] )
      {
         int status = (this->*messageHandler[tag])(source, tag);
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

      /** completion message may delay */
      if( paraRacingSolverPool && racingTermination &&  paraRacingSolverPool->getNumActiveSolvers() == 0 )
      {
         delete paraRacingSolverPool;
         paraRacingSolverPool = 0;
         break;
      }

      /** output tabular solving status */
      if( outputTabularSolvingStatusFlag &&
            paraSolverPool->getNumActiveSolvers() != 0 &&
            ( ( ( paraParamSet->getBoolParamValue(Deterministic) &&
                  paraParamSet->getBoolParamValue(DeterministicTabularSolvingStatus) ) &&
                  ( paraDetTimer->getElapsedTime() - previousTabularOutputTime ) >
                                 paraParamSet->getRealParamValue(TabularSolvingStatusInterval) ) ||
            ( ( !paraParamSet->getBoolParamValue(Deterministic) ||
                  !paraParamSet->getBoolParamValue(DeterministicTabularSolvingStatus) ) &&
                  ( paraTimer->getElapsedTime() - previousTabularOutputTime ) >
               paraParamSet->getRealParamValue(TabularSolvingStatusInterval) ) ) )
      {
         outputTabularSolvingStatus(' ');
         if( paraParamSet->getBoolParamValue(Deterministic) )
         {
            if( paraParamSet->getBoolParamValue(DeterministicTabularSolvingStatus) )
            {
               previousTabularOutputTime = paraDetTimer->getElapsedTime();
            }
            else
            {
               previousTabularOutputTime = paraTimer->getElapsedTime();
            }
         }
         else
         {
            previousTabularOutputTime = paraTimer->getElapsedTime();
         }
      }

      switch ( runningPhase )
      {
      case RampUpPhase:
      {
         if( racingTermination )
         {
            sendInterruptRequest();
            runningPhase = TerminationPhase;
         }
         else
         {
            double globalBestDualBoundValue =
               std::max (
                     std::min( paraSolverPool->getGlobalBestDualBoundValue(), paraNodePool->getBestDualBoundValue() ),
                     lcts.globalBestDualBoundValue );
            if( paraSolverPool->isActive() &&
                 ( paraNodePool->getNumOfGoodNodes( globalBestDualBoundValue )
                                 > paraParamSet->getIntParamValue(NChangeIntoCollectingMode) ||
                  ( paraNodePool->getNumOfNodes()
                                 > paraParamSet->getIntParamValue(NChangeIntoCollectingMode) * 2  &&
                                 paraNodePool->getNumOfGoodNodes( globalBestDualBoundValue ) > 0  ) ) )
            {
               sendRampUpToAllSolvers();
               runningPhase = NormalRunningPhase;
            }
            (void) sendParaNodesToIdleSolvers();
         }
         break;
      }
      case NormalRunningPhase:
      {
         if( racingTermination ||
             (  paraParamSet->getRealParamValue(TimeLimit) > 0.0 &&
               paraTimer->getElapsedTime() > paraParamSet->getRealParamValue(TimeLimit) ) )
         {
            sendInterruptRequest();
            runningPhase = TerminationPhase;
         }
         else
         {
            (void) sendParaNodesToIdleSolvers();
         }
         if( isCollectingModeRestarted && paraNodePool->isEmpty() &&
               ( (!paraRacingSolverPool) || ( paraRacingSolverPool && paraRacingSolverPool->getWinner() > 0 ) ) )
         {
            sendRetryRampUpToAllSolvers();
            runningPhase = RampUpPhase;
         }
         break;
      }
      case TerminationPhase:
      {
         break;
      }
      default:
      {
         THROW_LOGICAL_ERROR2( "Undefined running phase: ", static_cast<int>(runningPhase) );
      }
      }
      if( paraParamSet->getBoolParamValue(Checkpoint) &&
            ( paraTimer->getElapsedTime() - previousCheckpointTime )
            > paraParamSet->getRealParamValue(CheckpointInterval) )
      {
          updateCheckpointFiles();
          previousCheckpointTime = paraTimer->getElapsedTime();
      }
  }
}

bool
ParaLoadCoordinator::sendParaNodesToIdleSolvers(
      )
{
   if( merging || initialNodesGenerated ||
       runningPhase == TerminationPhase ||
         ( !restarted &&
            paraParamSet->getBoolParamValue(RacingStatBranching) &&
            ( !winnerSolverNodesCollected ||
                  ( paraRacingSolverPool &&
                        paraRacingSolverPool->getNumInactiveSolvers() < paraRacingSolverPool->getNumActiveSolvers() )
            )
         )
       )
   {
      return false;
   }

   bool sentNode = false;
   while( paraSolverPool->getNumInactiveSolvers() > 0 && !paraNodePool->isEmpty() )
   {
      ParaNode *paraNode = 0;
      while( !paraNodePool->isEmpty() )
      {
         paraNode = paraNodePool->extractNode();
         if( !paraNode ) break;
         assert( !paraNode->getMergeNodeInfo() ||
               ( paraNode->getMergeNodeInfo() &&
               paraNode->getMergeNodeInfo()->status == ParaMergeNodeInfo::PARA_MERGED_RPRESENTATIVE &&
               paraNode->getMergeNodeInfo()->mergedTo == 0 ) );
         if( ( paraInitiator->getGlobalBestIncumbentSolution() &&
               paraNode->getDualBoundValue() < paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue() ) ||
               !( paraInitiator->getGlobalBestIncumbentSolution() ) )
         {
            break;
         }
         else
         {
            delete paraNode;
            paraNode = 0;
            lcts.nDeletedInLc++;
         }
      }
      if( paraNode )
      {
         if( paraNode->getMergeNodeInfo() && paraNode->getMergeNodeInfo()->nMergedNodes == 0 )
         {
            ParaMergeNodeInfo *mNode = paraNode->getMergeNodeInfo();
            paraNode->setDiffSubproblem(mNode->origDiffSubproblem);
            paraNode->setMergeNodeInfo(0);
            paraNode->setMergingStatus(-1);
            delete mNode->mergedDiffSubproblem;
            deleteMergeNodeInfo(mNode);
         }
         if( paraParamSet->getBoolParamValue(RacingStatBranching) &&
               paraNode->isSameParetntNodeSubtreeIdAs( NodeId() ) &&  //  if parent is the root node
               paraNode->getDiffSubproblem()                          //  paraNode deos not root
               )
         {
            paraInitiator->setInitialStatOnDiffSubproblem(
                  minDepthInWinnerSolverNodes, maxDepthInWinnerSolverNodes,
                  paraNode->getDiffSubproblem());
         }
         int destination = paraSolverPool->activateSolver(paraNode, paraRacingSolverPool, (runningPhase==RampUpPhase) );
         if( destination < 0 )
         {
            /** cannot activate */
            paraNodePool->insert(paraNode);
            return sentNode;
         }
         else
         {
            lcts.nSent++;
            writeTransferLog(destination);
            sentNode = true;
            if( logSolvingStatusFlag )
            {
               *osLogSolvingStatus << paraTimer->getElapsedTime()
               << " S." << destination << " < "
               << paraInitiator->convertToExternalValue(
                     paraNode->getDualBoundValue() );
               if( paraInitiator->getGlobalBestIncumbentSolution() )
               {
                  if( paraInitiator->getGap(paraNode->getDualBoundValue()) > displayInfOverThisValue )
                  {
                     *osLogSolvingStatus << " ( Inf )";
                  }
                  else
                  {
                     *osLogSolvingStatus << " ( " << paraInitiator->getGap(paraNode->getDualBoundValue()) * 100 << "% )";
                  }
               }
               if( paraParamSet->getBoolParamValue(LightWeightRootNodeProcess) &&
                     runningPhase != RampUpPhase && (!paraRacingSolverPool) &&
                     paraSolverPool->getNumInactiveSolvers() >
                         ( paraSolverPool->getNSolvers() * paraParamSet->getRealParamValue(RatioToApplyLightWeightRootProcess) ) )
               {
                  *osLogSolvingStatus << " L";
               }
               if( paraNode->getMergeNodeInfo() )
               {
                  *osLogSolvingStatus << " M(" << paraNode->getMergeNodeInfo()->nMergedNodes + 1 << ")";
               }
               if( paraNode->getDiffSubproblem() )
               {
                  *osLogSolvingStatus << " " << paraNode->getDiffSubproblem()->getNBoundChanges();
               }
               *osLogSolvingStatus << std::endl;
            }
#ifdef _DEBUG_LB
            std::cout << paraTimer->getElapsedTime()
            << " S." << destination << " < "
            << paraInitiator->convertToExternalValue(
                  paraNode->getDualBoundValue() );
            if( paraInitiator->getGlobalBestIncumbentSolution() )
            {
               if( paraInitiator->getGap(paraNode->getDualBoundValue()) > displayInfOverThisValue )
               {
                  std::cout << " ( Inf )";
               }
               else
               {
                  std::cout << " ( " << paraInitiator->getGap(paraNode->getDualBoundValue()) * 100 << "% )";
               }
            }
            if( paraParamSet->getBoolParamValue(LightWeightRootNodeProcess) &&
                  runningPhase != RampUpPhase && (!paraRacingSolverPool) &&
                  paraSolverPool->getNumInactiveSolvers() >
                     ( paraSolverPool->getNSolvers() * paraParams->getRealParamValue(RatioToApplyLightWeightRootProcess) ) )
            {
               std::cout << " L";
            }
            std::cout << std::endl;
#endif
         }
      }
      else
      {
         break;
      }
   }
   return sentNode;
}

void
ParaLoadCoordinator::updateCheckpointFiles(
      )
{
   time_t timer;
   char timeStr[26];

   /** get checkpoint time */
   time(&timer);
   /** make checkpoint time string */
   ctime_r(&timer, timeStr);
   for( int i = 0; timeStr[i] != '\0' && i < 26; i++ )
   {
      if( timeStr[i] == ' ') timeStr[i] = '_';
      if( timeStr[i] == '\n' ) timeStr[i] = '\0';
   }
   char *newCheckpointTimeStr = &timeStr[4];    // remove a day of the week

   /** save nodes information */
   char nodesFileName[256];
   sprintf(nodesFileName,"%s%s_%s_nodes_LC%d.gz",
         paraParamSet->getStringParamValue(CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr, paraComm->getRank());
   ogzstream checkpointNodesStream;
   checkpointNodesStream.open(nodesFileName, std::ios::out | std::ios::binary);
   if( !checkpointNodesStream )
   {
      std::cout << "Checkpoint file for ParaNodes cannot open. file name = " << nodesFileName << std::endl;
      exit(1);
   }
   paraSolverPool->updateDualBoundsForSavingNodes();
   paraNodePool->updateDualBoundsForSavingNodes();
   int n;
   n = paraSolverPool->writeParaNodesToCheckpointFile(checkpointNodesStream);
   n += paraNodePool->writeParaNodesToCheckpointFile(checkpointNodesStream);
   checkpointNodesStream.close();
   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus << paraTimer->getElapsedTime()
      << " Checkpoint: " << n << " ParaNodes were saved" <<  std::endl;
   }
#ifdef _DEBUG_LB
   std::cout << paraTimer->getElapsedTime()
   << " Checkpoint: " << n << " ParaNodes were saved" <<  std::endl;
#endif

   if( outputTabularSolvingStatusFlag )
   {
      *osTabularSolvingStatus <<
            "Storing check-point data after " <<
            paraTimer->getElapsedTime() << " seconds. " <<
            n << " nodes were saved." << std::endl;
   }

   /** save incumbent solution */
   char solutionFileName[256];
   if( paraComm->getRank() == 0 )
   {
      sprintf(solutionFileName,"%s%s_%s_solution.gz",
            paraParamSet->getStringParamValue(CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr);
      paraInitiator->writeCheckpointSolution(std::string(solutionFileName));
   }

   /** save Solver statistics */
   char solverStatisticsFileName[256];
   sprintf(solverStatisticsFileName,"%s%s_%s_solverStatistics_LC%d.gz",
         paraParamSet->getStringParamValue(CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr, paraComm->getRank());
   ogzstream checkpointSolverStatisticsStream;
   checkpointSolverStatisticsStream.open(solverStatisticsFileName, std::ios::out | std::ios::binary);
   if( !checkpointSolverStatisticsStream )
   {
      std::cout << "Checkpoint file for SolverStatistics cannot open. file name = " << solverStatisticsFileName << std::endl;
      exit(1);
   }
   int nSolverInfo = paraSolverPool->writeSolverStatisticsToCheckpointFile(checkpointSolverStatisticsStream);
   checkpointSolverStatisticsStream.close();

   /** save LoadCoordinator statistics */
   char loadCoordinatorStatisticsFileName[256];
   sprintf(loadCoordinatorStatisticsFileName,"%s%s_%s_loadCoordinatorStatistics_LC%d.gz",
         paraParamSet->getStringParamValue(CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr, paraComm->getRank());
   ogzstream loadCoordinatorStatisticsStream;
   loadCoordinatorStatisticsStream.open(loadCoordinatorStatisticsFileName, std::ios::out | std::ios::binary);
   if( !loadCoordinatorStatisticsStream )
   {
      std::cout << "Checkpoint file for SolverStatistics cannot open. file name = " << loadCoordinatorStatisticsFileName << std::endl;
      exit(1);
   }
   // double globalBestDualBoundValue =
   //   std::max (
   //      std::min( paraSolverPool->getGlobalBestDualBoundValue(), paraNodePool->getBestDualBoundValue() ),
   //      lcts.globalBestDualBoundValue );
   // double externalGlobalBestDualBoundValue = paraInitiator->convertToExternalValue(globalBestDualBoundValue);
   writeLoadCoordinatorStatisticsToCheckpointFile(loadCoordinatorStatisticsStream, nSolverInfo,
         lcts.globalBestDualBoundValue, lcts.externalGlobalBestDualBoundValue );
   loadCoordinatorStatisticsStream.close();

   if( lastCheckpointTimeStr[0] == ' ' )
   {
      /** the first time for checkpointing */
      if( racingWinnerParams )
      {
         /** save racing winner params */
         char racingWinnerParamsName[256];
         sprintf(racingWinnerParamsName,"%s%s_racing_winner_params.gz",
               paraParamSet->getStringParamValue(CheckpointFilePath),
               paraInitiator->getParaInstance()->getProbName());
         ogzstream racingWinnerParamsStream;
         racingWinnerParamsStream.open(racingWinnerParamsName, std::ios::out | std::ios::binary);
         if( !racingWinnerParamsStream )
         {
            std::cout << "Racing winner parameter file cannot open. file name = " << racingWinnerParamsName << std::endl;
            exit(1);
         }
         racingWinnerParams->write(racingWinnerParamsStream);
         racingWinnerParamsStream.close();
      }
   }
   else
   {
      /** remove old check point files */
      sprintf(nodesFileName,"%s%s_%s_nodes_LC%d.gz",
            paraParamSet->getStringParamValue(CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr, paraComm->getRank());
      if( paraComm->getRank() == 0 )
      {
         sprintf(solutionFileName,"%s%s_%s_solution.gz",
               paraParamSet->getStringParamValue(CheckpointFilePath),
               paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr);
      }
      sprintf(solverStatisticsFileName,"%s%s_%s_solverStatistics_LC%d.gz",
            paraParamSet->getStringParamValue(CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr, paraComm->getRank());
      sprintf(loadCoordinatorStatisticsFileName,"%s%s_%s_loadCoordinatorStatistics_LC%d.gz",
            paraParamSet->getStringParamValue(CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr, paraComm->getRank());
      if( remove(nodesFileName) )
      {
         std::cout << "checkpoint nodes file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(solutionFileName) )
      {
         std::cout << "checkpoint solution file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(solverStatisticsFileName) )
      {
         std::cout << "checkpoint SolverStatistics file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(loadCoordinatorStatisticsFileName) )
      {
         std::cout << "checkpoint LoadCoordinatorStatistics file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      char afterCheckpointingSolutionFileName[256];
      sprintf(afterCheckpointingSolutionFileName,"%s%s_after_checkpointing_solution.gz",
            paraParamSet->getStringParamValue(CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName() );
      igzstream afterCheckpointingSolutionStream;
      afterCheckpointingSolutionStream.open(afterCheckpointingSolutionFileName, std::ios::in | std::ios::binary);
      if( afterCheckpointingSolutionStream  )
      {
         /** afater checkpointing solution file exists */
         afterCheckpointingSolutionStream.close();
         if ( remove(afterCheckpointingSolutionFileName) )
         {
            std::cout << "after checkpointing solution file cannot be removed: errno = " << strerror(errno) << std::endl;
            exit(1);
         }
      }
   }

   /** update last checkpoint time string */
   strcpy(lastCheckpointTimeStr,newCheckpointTimeStr);
}

void
ParaLoadCoordinator::writeLoadCoordinatorStatisticsToCheckpointFile(
      ogzstream &loadCoordinatorStatisticsStream,
      int nSolverInfo,
      double globalBestDualBoundValue,
      double externalGlobalBestDualBoundValue
      )
{
   loadCoordinatorStatisticsStream.write((char *)&nSolverInfo, sizeof(int));
   lcts.isCheckpointState = true;
   lcts.nMaxUsageOfNodePool = paraNodePool->getMaxUsageOfPool();
   lcts.nNodesInNodePool = paraNodePool->getNumOfNodes();
   lcts.nNodesLeftInAllSolvers = paraSolverPool->getNnodesInSolvers();
   lcts.globalBestDualBoundValue = globalBestDualBoundValue;
   lcts.externalGlobalBestDualBoundValue = externalGlobalBestDualBoundValue;
   lcts.runningTime = paraTimer->getElapsedTime();
   lcts.write(loadCoordinatorStatisticsStream);
}

void
ParaLoadCoordinator::warmStart(
      )
{
   restarted = true;
   /** write previous statistics information */
   writePreviousStatisticsInformation();

   /** try to read racing winner params */
   char racingWinnerParamsName[256];
   sprintf(racingWinnerParamsName,"%s%s_racing_winner_params.gz",
         paraParamSet->getStringParamValue(CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName());
   igzstream racingWinnerParamsStream;
   racingWinnerParamsStream.open(racingWinnerParamsName, std::ios::in | std::ios::binary);
   if( racingWinnerParamsStream )
   {
      assert(!racingWinnerParams);
      racingWinnerParams = paraComm->createParaRacingRampUpParamSet();
      racingWinnerParams->read(paraComm, racingWinnerParamsStream);
      racingWinnerParamsStream.close();
      for( int i = 1; i < paraComm->getSize(); i++ )
      {
         /** send racing winner params: NOTE: should not broadcast. if we do it, solver routine need to recognize staring process */
         PARA_COMM_CALL(
               racingWinnerParams->send(paraComm, i)
         );
      }
   }

   /** set solution and get internal incumbent value */
   char afterCheckpointingSolutionFileName[256];
   sprintf(afterCheckpointingSolutionFileName,"%s%s_after_checkpointing_solution.gz",
         paraParamSet->getStringParamValue(CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName() );
   double incumbentValue = paraInitiator->readSolutionFromCheckpointFile(afterCheckpointingSolutionFileName);
   for( int i = 1; i < paraComm->getSize(); i++ )
   {
      /** send internal incumbent value */
      PARA_COMM_CALL(
            paraComm->send( &incumbentValue, 1, ParaDOUBLE, i, TagIncumbentValue )
      );
      if( paraParamSet->getBoolParamValue(DistributeBestPrimalSolution) && paraInitiator->getGlobalBestIncumbentSolution() )
      {
         paraInitiator->getGlobalBestIncumbentSolution()->send(paraComm, i);
      }
      /** send internal global dual bound value */
      PARA_COMM_CALL(
            paraComm->send( &lcts.globalBestDualBoundValue, 1, ParaDOUBLE, i, TagGlobalBestDualBoundValueAtWarmStart )
      );
   }
   if( !paraParamSet->getBoolParamValue(Quiet) )
   {
      paraInitiator->writeSolution("[Warm started from "+std::string(paraInitiator->getPrefixWarm())+" : the solution from the checkpoint file]");
   }

   int n = 0;
   ParaNode *paraNode;
   bool onlyBoundChanges = false;
   if( !paraParamSet->getBoolParamValue(TransferLocalCuts) && !paraParamSet->getBoolParamValue(TransferConflicts) )
   {
      onlyBoundChanges = true;
   }
   if( paraParamSet->getBoolParamValue(MergeNodesAtRestart) )
   {
      initMergeNodesStructs();
   }

   ParaNodePoolForMinimization tempParaNodePool(paraParamSet->getRealParamValue(BgapCollectingMode));

   if( paraParamSet->getIntParamValue(AddDualBoundCons) == 3 )
   {
      /** NOT implemented */
      for(;;)
      {
         paraNode = paraInitiator->readParaNodeFromCheckpointFile(onlyBoundChanges);
         if( paraNode == 0 )
            break;
         n++;
         paraNode->setDualBoundValue(paraNode->getInitialDualBoundValue());
         paraNodePool->insert(paraNode);   /** in order to sort ParaNodes, insert paraNodePool once */
         if( paraParamSet->getBoolParamValue(MergeNodesAtRestart) )
         {
            addNodeToMergeNodeStructs(paraNode);
         }
      }
   }
   else
   {
      for(;;)
      {
         paraNode = paraInitiator->readParaNodeFromCheckpointFile(onlyBoundChanges);
         if( paraNode == 0 )
            break;
         n++;
         paraNode->setDualBoundValue(paraNode->getInitialDualBoundValue());
         // paraNodePool->insert(paraNode);   /** in order to sort ParaNodes, insert paraNodePool once */
         if( paraParamSet->getBoolParamValue(MergeNodesAtRestart) )
         {
            tempParaNodePool.insert(paraNode);  /** in order to sort ParaNodes with new dual value, insert it to tempParaNodePool once */
            // addNodeToMergeNodeStructs(paraNode);
         }
         else
         {
            paraNodePool->insert(paraNode);   /** in order to sort ParaNodes, insert paraNodePool once */
         }
      }
   }

   if( paraParamSet->getBoolParamValue(MergeNodesAtRestart) )
   {
      ParaNode *tempNode = 0;
      while( ( tempNode = tempParaNodePool.extractNode() ) )
      {
         paraNodePool->insert(tempNode);
         addNodeToMergeNodeStructs(tempNode);
      }
      generateMergeNodesCandidates();
   }

   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus << paraTimer->getElapsedTime()
      << " Warm started from "
      << paraInitiator->getPrefixWarm()
      << " : " << n << " ParaNodes read. Current incumbent value = "
      <<  paraInitiator->convertToExternalValue(
            paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue() )
            << std::endl;
   }
#ifdef _DEBUG_LB
   std::cout << paraTimer->getElapsedTime()
   << " Warm started from "
   << paraInitiator->getPrefixWarm()
   << " : " << n << " ParaNodes read. Current incumbent value = "
   <<  paraInitiator->convertToExternalValue(
         paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue() )
         << std::endl;
#endif
   (void) sendParaNodesToIdleSolvers();
   if( paraSolverPool->getNumInactiveSolvers() > 0 )
   {
      runningPhase = RampUpPhase;
   }
   else
   {
      sendRampUpToAllSolvers();
      runningPhase = NormalRunningPhase;
   }
   run();
}

void
ParaLoadCoordinator::writePreviousStatisticsInformation(
      )
{
   /* read previous LoadCoordinator statistics */
   char loadCoordinatorStatisticsFileName[256];
   sprintf(loadCoordinatorStatisticsFileName,"%s_loadCoordinatorStatistics_LC0.gz", paraInitiator->getPrefixWarm() );
   igzstream  loadCoordinatorStatisticsStream;
   loadCoordinatorStatisticsStream.open(loadCoordinatorStatisticsFileName, std::ios::in | std::ios::binary);
   if( !loadCoordinatorStatisticsStream )
   {
      std::cout << "checkpoint LoadCoordinatorStatistics file cannot open: file name = " <<  loadCoordinatorStatisticsFileName << std::endl;
      exit(1);
   }
   int nSolverStatistics;
   loadCoordinatorStatisticsStream.read((char *)&nSolverStatistics, sizeof(int));
   ParaLoadCoordinatorTerminationState *prevLcts = new ParaLoadCoordinatorTerminationState();
   if( !prevLcts->read(paraComm, loadCoordinatorStatisticsStream) )
   {
      std::cout << "checkpoint LoadCoordinatorStatistics file cannot read: file name = " <<  loadCoordinatorStatisticsFileName << std::endl;
      exit(1);
   }
   loadCoordinatorStatisticsStream.close();

   /* open Solver statistics file */
   char solverStatisticsFileName[256];
   sprintf(solverStatisticsFileName,"%s_solverStatistics_LC0.gz", paraInitiator->getPrefixWarm() );
   igzstream  solverStatisticsStream;
   solverStatisticsStream.open(solverStatisticsFileName, std::ios::in | std::ios::binary);
   if( !solverStatisticsStream )
   {
      std::cout << "checkpoint SolverStatistics file cannot open: file name = " <<  solverStatisticsFileName << std::endl;
      exit(1);
   }

   /* opne output statistics file */
   char previousStatisticsFileName[256];
   sprintf(previousStatisticsFileName, "%s_statistics_w%02d_LC0",
         paraInitiator->getPrefixWarm(),
         prevLcts->nWarmStart);
   std::ofstream ofsStatistics;
   ofsStatistics.open(previousStatisticsFileName);
   if( !ofsStatistics )
   {
      std::cout << "previous statistics file cannot open : file name = " << previousStatisticsFileName << std::endl;
      exit(1);
   }

   /* read and write solver statistics */
   for( int i = 0; i < nSolverStatistics; i++ )
   {
      ParaSolverTerminationState *psts = paraComm->createParaSolverTerminationState();
      if( !psts->read(paraComm, solverStatisticsStream) )
      {
         std::cout << "checkpoint SolverStatistics file cannot read: file name = " <<  solverStatisticsFileName << std::endl;
         exit(1);
      }
      ofsStatistics << psts->toString();
      delete psts;
   }

   /* write LoadCoordinator statistics */
   ofsStatistics << prevLcts->toString();

   /* update warm start counter */
   lcts.nWarmStart = prevLcts->nWarmStart + 1;
   lcts.globalBestDualBoundValue = prevLcts->globalBestDualBoundValue;
   lcts.externalGlobalBestDualBoundValue = prevLcts->externalGlobalBestDualBoundValue;
   delete prevLcts;

   /* close solver statistics file and output file */
   solverStatisticsStream.close();
   ofsStatistics.close();
}

void
ParaLoadCoordinator::run(
      ParaNode *paraNode
      )
{
   assert(!paraRacingSolverPool);
   int destination = paraSolverPool->activateSolver(paraNode, paraRacingSolverPool, ( runningPhase == RampUpPhase ) );
   lcts.nSent++;
   if( paraParamSet->getBoolParamValue(Deterministic) )
   {
      int token[2];
      token[0] = 1;
      token[1] = -1;
      PARA_COMM_CALL(
            paraComm->send( token, 2, ParaINT, token[0], TagToken )
      );
   }
   writeTransferLog(destination);
   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus << paraTimer->getElapsedTime()
      << " S." << destination << " < "
      << paraInitiator->convertToExternalValue(
            paraNode->getDualBoundValue() );
      if( paraInitiator->getGlobalBestIncumbentSolution() )
      {
         if( paraInitiator->getGap(paraNode->getDualBoundValue()) > displayInfOverThisValue )
         {
            *osLogSolvingStatus << " ( Inf )";
         }
         else
         {
            *osLogSolvingStatus << " ( " << paraInitiator->getGap(paraNode->getDualBoundValue()) * 100 << "% )";
         }
      }
      if( paraParamSet->getBoolParamValue(LightWeightRootNodeProcess) &&
            runningPhase != RampUpPhase && (!paraRacingSolverPool) &&
            paraSolverPool->getNumInactiveSolvers() >
               ( paraSolverPool->getNSolvers() * paraParamSet->getRealParamValue(RatioToApplyLightWeightRootProcess) ) )
      {
         *osLogSolvingStatus << " L";
      }
      *osLogSolvingStatus << std::endl;
   }
#ifdef DEBUG_LB
   std::cout << paraTimer->getElapsedTime()
   << " S." << destination << " > "
   << paraInitiator->convertToExternalValue(
         paraNode->getDualBoundValue() );
   if( paraInitiator->getGlobalBestIncumbentSolution() )
   {
      if( paraInitiator->getGap(paraNode->getDualBoundValue()) > displayInfOverThisValue )
      {
         std::cout << " ( Inf )";
      }
      else
      {
         std::cout << " ( " << paraInitiator->getGap(paraNode->getDualBoundValue()) * 100 << "% )";
      }
   }
   if( paraParamSet->getBoolParamValue(LightWeightRootNodeProcess) &&
         runningPhase != RampUpPhase && (!paraRacingSolverPool) &&
         paraSolverPool->getNumInactiveSolvers() >
            ( paraSolverPool->getNSolvers() * paraParams->getRealParamValue(RatioToApplyLightWeightRootProcess) ) )
   {
      std::cout << " L";
   }
   std::cout << std::endl;
#endif
   run();
}

int
ParaLoadCoordinator::processRacingRampUpTagSolverState(
      int source,
      int tag
      )
{

   ParaSolverState *solverState = paraComm->createParaSolverState();
   solverState->receive(paraComm, source, tag);

#ifdef _DEBUG_DET
   if( paraDetTimer )
   {
      std::cout << "Rank " << source << ": ET = " << paraDetTimer->getElapsedTime() << ", Det time = " << solverState->getDeterministicTime() << std::endl;
   }
#endif

   if( paraDetTimer
         && paraDetTimer->getElapsedTime() < solverState->getDeterministicTime() )

   {
      paraDetTimer->update( solverState->getDeterministicTime() - paraDetTimer->getElapsedTime() );
   }

   assert(solverState->isRacingStage());
   paraRacingSolverPool->updateSolverStatus(source,
                                   solverState->getNNodesSolved(),
                                   solverState->getNNodesLeft(),
                                   solverState->getSolverLocalBestDualBoundValue());
   assert( paraRacingSolverPool->getNumNodesLeft(source) == solverState->getNNodesLeft() );
   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus << paraTimer->getElapsedTime()
      << " S." << source << " | "
      << paraInitiator->convertToExternalValue(
            solverState->getSolverLocalBestDualBoundValue()
            );
      if( !paraInitiator->getGlobalBestIncumbentSolution() ||
            paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) > displayInfOverThisValue
            || solverState->getNNodesLeft() == 0 )
      {
         *osLogSolvingStatus << " ( Inf )";
      }
      else
      {
         *osLogSolvingStatus << " ( " << paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) * 100 << "% )";
      }
      *osLogSolvingStatus << " [ " << solverState->getNNodesLeft() << " ]";
      double globalBestDualBoundValue = paraRacingSolverPool->getGlobalBestDualBoundValue();
      *osLogSolvingStatus << " ** G.B.: " << paraInitiator->convertToExternalValue(globalBestDualBoundValue);
      if( !paraInitiator->getGlobalBestIncumbentSolution() ||
            paraInitiator->getGap(globalBestDualBoundValue) > displayInfOverThisValue )
      {
         *osLogSolvingStatus << " ( Inf ) ";
      }
      else
      {
         *osLogSolvingStatus << " ( " << paraInitiator->getGap(globalBestDualBoundValue) * 100 << "% ) ";
      }
      *osLogSolvingStatus << "[ " << paraRacingSolverPool->getNnodesLeftInBestSolver()
 //     <<" ] ** RR" << std::endl;
      <<" ] ** RR " << solverState->getDeterministicTime() << std::endl;   // for debug
   }
#ifdef _DEBUG_LB
   std::cout << paraTimer->getElapsedTime()
   << " S." << source << " | "
   << paraInitiator->convertToExternalValue(
         solverState->getSolverLocalBestDualBoundValue()
         );
   if( !paraInitiator->getGlobalBestIncumbentSolution() ||
         paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) > displayInfOverThisValue
         || solverState->getNNodesLeft() == 0 )
   {
      std::cout << " ( Inf )";
   }
   else
   {
      std::cout << " ( " << paraInitiator->getGap(solverState->getSolverLocalBestDualBoundValue()) * 100 << "% )";
   }
   std::cout << " [ " << solverState->getNNodesLeft() << " ]";
   double globalBestDualBoundValue = paraRacingSolverPool->getGlobalBestDualBoundValue();
   std::cout << " ** G.B.: " << paraInitiator->convertToExternalValue(globalBestDualBoundValue);
   if( !paraInitiator->getGlobalBestIncumbentSolution() ||
         paraInitiator->getGap(globalBestDualBoundValue) > displayInfOverThisValue )
   {
      std::cout << " ( Inf ) ";
   }
   else
   {
      std::cout << " ( " << paraInitiator->getGap(globalBestDualBoundValue) * 100 << "% ) ";
   }
   std::cout << "[ " << paraRacingSolverPool->getNnodesLeftInBestSolver()
   <<" ] ** RR" << std::endl;
#endif

   if( !paraParamSet->getBoolParamValue(NoUpperBoundTransferInRacing) )
   {
      /** the following should be before noticationId back to the source solver */
      if( paraParamSet->getBoolParamValue(DistributeBestPrimalSolution) )
      {
         if( paraInitiator->getGlobalBestIncumbentSolution() &&
               paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue()
               < solverState->getGlobalBestPrimalBoundValue() )
         {
            paraInitiator->getGlobalBestIncumbentSolution()->send(paraComm, source);
         }
      }
   }

   double lcBestDualBoundValue = paraRacingSolverPool->getGlobalBestDualBoundValue();
   PARA_COMM_CALL(
         paraComm->send( &lcBestDualBoundValue, 1, ParaDOUBLE, source, TagLCBestBoundValue)
         );
   unsigned int notificationId = solverState->getNotificaionId();
   PARA_COMM_CALL(
         paraComm->send( &notificationId, 1, ParaUNSIGNED, source, TagNotificationId)
         );

   if( lcts.globalBestDualBoundValue < paraRacingSolverPool->getGlobalBestDualBoundValue() )
   {
      lcts.globalBestDualBoundValue = paraRacingSolverPool->getGlobalBestDualBoundValue();
      lcts.externalGlobalBestDualBoundValue = paraInitiator->convertToExternalValue(lcts.globalBestDualBoundValue);
   }
   delete solverState;
   return 0;
}

int
ParaLoadCoordinator::processRacingRampUpTagCompletionOfCalculation(
      int source,
      int tag
      )
{
   ParaCalculationState *calcState = paraComm->createParaCalculationState();
   calcState->receive(paraComm, source, tag);
   writeTransferLogInRacing(source, calcState);
   if( logSolvingStatusFlag )
   {
      switch ( calcState->getTerminationState() )
      {
      case CompTerminatedInRacingStage:
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " >(TERMINATED_IN_RACING_STAGE)";
         break;
      }
      case CompTerminatedByInterruptRequest:
      case CompInterruptedInRacingStage:
      {
         *osLogSolvingStatus << paraTimer->getElapsedTime()
         << " S." << source << " >(INTERRUPTED_BY_TIME_LIMIT or INTERRUPTED_BY_SOME_SOLVER_TERMINATED_IN_RACING_STAGE)";
         break;
      }
      default:
         THROW_LOGICAL_ERROR2("Invalid termination: termination state = ", calcState->getTerminationState() )
      }
      *osLogSolvingStatus << std::endl;
   }
#ifdef _DEBUG_LB
   switch ( calcState->getTerminationState() )
   {
   case CompTerminatedInRacingStage:
   {
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " >(TERMINATED_IN_RACING_STAGE)";
      break;
   }
   case CompTerminatedByInterruptRequest:
   case CompInterruptedInRacingStage:
   {
      std::cout << paraTimer->getElapsedTime()
      << " S." << source << " >(INTERRUPTED_BY_TIME_LIMIT or INTERRUPTED_BY_SOME_SOLVER_TERMINATED_IN_RACING_STAGE)";
      break;
   }
   default:
      THROW_LOGICAL_ERROR2("Invalid termination: termination state = ", calcState->getTerminationState() )
   }
   std::cout << std::endl;
#endif

   if( calcState->getTerminationState() == CompTerminatedInRacingStage )
   {
      racingTermination = true; // even if interruptIsRequested, 
                                // solver should have been terminated before receiveing it
      if( osStatisticsRacingRampUp )
      {
         *osStatisticsRacingRampUp << "######### Solver Rank = " <<
               source << " is terminated in racing stage #########" << std::endl;
      }
      nSolvedRacingTermination = calcState->getNSolved();
      if(  !interruptedFromControlTerminal && !hardTimeLimitIsReached &&
            !computationIsInterrupted && paraInitiator->getGlobalBestIncumbentSolution() )
      {
         lcts.globalBestDualBoundValue = paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue();
         lcts.externalGlobalBestDualBoundValue = paraInitiator->convertToExternalValue(lcts.globalBestDualBoundValue);
      }
   }

   delete calcState;

   ParaSolverTerminationState *termState = paraComm->createParaSolverTerminationState();
   termState->receive(paraComm, source, TagTerminated);

   if( paraDetTimer )
   {
      if( paraDetTimer->getElapsedTime() < termState->getDeterministicTime() )
      {
         paraDetTimer->update( termState->getDeterministicTime() - paraDetTimer->getElapsedTime() );
      }
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, source, TagAckCompletion )
      );
   }

   if( osStatisticsRacingRampUp )
   {
      *osStatisticsRacingRampUp << termState->toString();
   }

   if( paraParamSet->getBoolParamValue(StatisticsToStdout) )
   {
      std::cout << termState->toString() << std::endl;
   }

   delete termState;
   // nTerminated++;      We should not count this, We should always send Term from LC!
   inactivateRacingSolverPool(source);

   if( paraParamSet->getBoolParamValue(Quiet) && racingTermination )
   {
      /** in this case, do not have to wait statistical information from the other solvers */
      nTerminated = 1;
      delete this;
#ifdef _COMM_PTH
      _exit(0);
#else
      exit(0);
#endif
   }

   return 0;
}


void
ParaLoadCoordinator::run(
      ParaNode *paraNode,
      int nRacingSolvers,
      ParaRacingRampUpParamSet **racingRampUpParams
      )
{
   /** register message handlers */
   for( int i = 0; i < N_TAGS; i++ )
   {
      racingRampUpMessageHandler[i] = 0;
   }
   racingRampUpMessageHandler[TagSolution] = &UG::ParaLoadCoordinator::processTagSolution;
   racingRampUpMessageHandler[TagSolverState] = &UG::ParaLoadCoordinator::processRacingRampUpTagSolverState;
   racingRampUpMessageHandler[TagCompletionOfCalculation] = &UG::ParaLoadCoordinator::processRacingRampUpTagCompletionOfCalculation;
   racingRampUpMessageHandler[TagAnotherNodeRequest] = &UG::ParaLoadCoordinator::processTagAnotherNodeRequest;
   racingRampUpMessageHandler[TagTerminated] = &UG::ParaLoadCoordinator::processTagTerminated;
   racingRampUpMessageHandler[TagHardTimeLimit] = &UG::ParaLoadCoordinator::processTagHardTimeLimit;
   racingRampUpMessageHandler[TagInitialStat] = &UG::ParaLoadCoordinator::processTagInitialStat;
   if( paraParamSet->getBoolParamValue(Deterministic) )
   {
      racingRampUpMessageHandler[TagToken] = &UG::ParaLoadCoordinator::processTagToken;
   }

   /** creates racing solver pool */
   paraRacingSolverPool = new ParaRacingSolverPool(
         1,                // paraSolver origin rank
         paraComm, paraParamSet, paraTimer, paraDetTimer);

   /** activate racing solver with root node */
   PARA_COMM_CALL(
         paraNode->bcast(paraComm, 0)
         );
   paraRacingSolverPool->activate(paraNode);
   lcts.nSent++;
   if( paraParamSet->getBoolParamValue(Deterministic) )
   {
      int token[2];
      token[0] = 1;
      token[1] = -1;
      PARA_COMM_CALL(
            paraComm->send( token, 2, ParaINT, token[0], TagToken )
      );
   }
   if( logNodesTransferFlag )
   {
      for(int i = 1; i < paraComm->getSize(); i++ )
      {
         writeTransferLogInRacing(i);
      }
   }

   /** output start racing to log file, if it is necessary */
   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus << paraTimer->getElapsedTime()
      << " All Solvers starts racing "
      << paraInitiator->convertToExternalValue(
            paraNode->getDualBoundValue() ) << std::endl;
   }
#ifdef _DEBUG_LB
   std::cout << paraTimer->getElapsedTime()
   << " All Solvers starts racing "
   << paraInitiator->convertToExternalValue(
         paraNode->getDualBoundValue() ) << std::endl;
#endif

   int source;
   int tag;

   for(;;)
   {
      /*******************************************
       *  waiting for any message form anywhere  *
       *******************************************/
      double inIdleTime = paraTimer->getElapsedTime();
      paraComm->probe(&source, &tag);
      lcts.idleTime += ( paraTimer->getElapsedTime() - inIdleTime );
      if( racingRampUpMessageHandler[tag] )
      {
         int status = (this->*racingRampUpMessageHandler[tag])(source, tag);
         if( status )
         {
            std::ostringstream s;
            s << "[ERROR RETURN form Racing Ramp-up Message Hander]:" <<  __FILE__ <<  "] func = "
              << __func__ << ", line = " << __LINE__ << " - "
              << "process tag = " << tag << std::endl;
            abort();
         }
      }
      else
      {
         THROW_LOGICAL_ERROR3( "No racing ramp-up message hander for ", tag, " is not registered" );
      }

      /** output tabular solving status */
      if( outputTabularSolvingStatusFlag && (!racingTermination) &&
            ( ( ( paraParamSet->getBoolParamValue(Deterministic) &&
                  paraParamSet->getBoolParamValue(DeterministicTabularSolvingStatus) ) &&
                  ( ( paraDetTimer->getElapsedTime() - previousTabularOutputTime ) >
               paraParamSet->getRealParamValue(TabularSolvingStatusInterval) ) ) ||
               ( (!paraParamSet->getBoolParamValue(Deterministic) ||
                     !paraParamSet->getBoolParamValue(DeterministicTabularSolvingStatus) )  &&
                  ( ( paraTimer->getElapsedTime() - previousTabularOutputTime ) >
               paraParamSet->getRealParamValue(TabularSolvingStatusInterval) ) ) ) )
      {
         outputTabularSolvingStatus(' ');
         if( paraParamSet->getBoolParamValue(Deterministic) )
         {
            if( paraParamSet->getBoolParamValue(DeterministicTabularSolvingStatus) )
            {
               previousTabularOutputTime = paraDetTimer->getElapsedTime();
            }
            else
            {
               previousTabularOutputTime = paraTimer->getElapsedTime();
            }
         }
         else
         {
            previousTabularOutputTime = paraTimer->getElapsedTime();
         }
      }

      if( hardTimeLimitIsReached )
          break;

      switch ( runningPhase )
      {
      case RampUpPhase:
      {
         if( !paraRacingSolverPool )
         {
            paraSolverPool->activate();
            run();
            return;
         }
         if( racingTermination )
         {
            if( paraRacingSolverPool->getNumActiveSolvers() == 0 )
            {
               delete paraRacingSolverPool;
               paraRacingSolverPool = 0;
               run();
               return;
            }
            sendInterruptRequest();
            break;
         }

         if( paraParamSet->getRealParamValue(TimeLimit) > 0.0 &&
               paraTimer->getElapsedTime() > paraParamSet->getRealParamValue(TimeLimit) )
         {
            hardTimeLimitIsReached = true;
            sendInterruptRequest();
            // runningPhase = TerminationPhase;     waits until paraRacingSolverPool becomes empty
            break;
         }

         if( paraRacingSolverPool->isWinnerDecided(
               ( paraInitiator->getGlobalBestIncumbentSolution() &&
                 EPSLT(paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue(),DBL_MAX, DEFAULT_NUM_EPSILON ) ) ) )
         {
            racingWinner = paraRacingSolverPool->getWinner();
            assert( racingWinner >0 );
            int numNodesLeft = paraRacingSolverPool->getNumNodesLeft(racingWinner);
            paraSolverPool->activateSolver(
                  racingWinner, paraRacingSolverPool->extractNode(), numNodesLeft );
            paraRacingSolverPool->inactivateSolver(racingWinner);
            assert(paraRacingSolverPool->getNumActiveSolvers() >= 0);
            PARA_COMM_CALL(
                  paraComm->send( NULL, 0, ParaBYTE, racingWinner, TagWinner )
            );

            if( numNodesLeft >
                 2.0 * paraParamSet->getIntParamValue(StopRacingNumberOfNodesLeft)*paraParamSet->getRealParamValue(StopRacingNumberOfNodesLeftMultiplier)
                 ||
                 numNodesLeft <= ( paraSolverPool->getNSolvers() * paraParamSet->getRealParamValue(ProhibitCollectOnceMultiplier) )
                 )
            {
               paraParamSet->setBoolParamValue(CollectOnce,false);
               paraParamSet->setBoolParamValue(MergeNodesAtRestart,false);
               paraParamSet->setBoolParamValue(RacingStatBranching,false);
               paraParamSet->setIntParamValue(RampUpPhaseProcess, 1);
               std::cout << "Warning: Ramp-Up Phase Process is switched to 1. CollectOnce, MergeNodesAtRestart and RacingStatBranching are switched to FALSE." << std::endl;
               std::cout << "You should check the following parameter values: StopRacingNumberOfNodesLeft, StopRacingNumberOfNodesLeftMultiplier, ProhibitCollectOnceMultiplier" << std::endl;
            }

            if( numNodesLeft > paraSolverPool->getNSolvers()
                  ||
                  numNodesLeft > ( paraSolverPool->getNSolvers() * paraParamSet->getRealParamValue(ProhibitCollectOnceMultiplier) )
                  )
            {
               if( paraParamSet->getBoolParamValue(CollectOnce) )
               {
                  int nCollect = paraParamSet->getIntParamValue(NCollectOnce);
                  if( nCollect == 0 )
                  {
                     nCollect = ( paraSolverPool->getNSolvers() * 5 );
                  }
                  PARA_COMM_CALL(
                        paraComm->send( &nCollect, 1, ParaINT, racingWinner, TagCollectAllNodes )
                  );
               }
               if( paraParamSet->getIntParamValue(RampUpPhaseProcess) == 2 )
               {
                  merging = true;
                  initMergeNodesStructs();
               }
            }
            else
            {
               winnerSolverNodesCollected = true;   // do not wait until all nodes are collected
            }
            racingWinnerParams = racingRampUpParams[racingWinner - 1];
            racingRampUpParams[racingWinner - 1] = 0;
            for(int i = 1; i < paraComm->getSize(); i++)
            {
               if( racingRampUpParams[i - 1] )
               {
                  PARA_COMM_CALL(
                        racingWinnerParams->send(paraComm, i)
                        );
               }
            }
            /** output winner to log file, if it is necessary */
            if( logSolvingStatusFlag )
            {
               *osLogSolvingStatus << paraTimer->getElapsedTime()
               << " S." << racingWinner << " is the racing winner!"
               << " Selected strategy " << racingWinnerParams->getStrategy()
               << "." << std::endl;
            }
#ifdef _DEBUG_LB
            std::cout << paraTimer->getElapsedTime()
            << " S." << racingWinner << " is the racing winner!"
            << " Selected strategy " << racingWinnerParams->getStrategy()
            << "." << std::endl;
#endif
            if( outputTabularSolvingStatusFlag )
            {
               *osTabularSolvingStatus <<
                     "Racing ramp-up finished after " <<
                     paraTimer->getElapsedTime() << " seconds." <<
                     " Selected strategy " << racingWinnerParams->getStrategy() <<
                     "." << std::endl;
            }
            // runningPhase = NormalRunningPhase;
            // Keep running as RampUpPhase, but in the run() switching into RampUpPhase in normal running mode
            // delete paraRacingSolverPool;
            run();
            return;
         }
         break;
      }
      default:
      {
         THROW_LOGICAL_ERROR2( "Undefined running phase: ", static_cast<int>(runningPhase) );
      }
      }
   }
   return;
}

void
ParaLoadCoordinator::sendRampUpToAllSolvers(
      )
{
   for( int i = 1; i < paraComm->getSize(); i++ )
   {
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, i, TagRampUp )
      );
   }
}

void
ParaLoadCoordinator::sendRetryRampUpToAllSolvers(
      )
{
   for( int i = 1; i < paraComm->getSize(); i++ )
   {
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, i, TagRetryRampUp )
      );
   }
}

void
ParaLoadCoordinator::sendInterruptRequest(
      )
{
   if( interruptIsRequested ) return;
   for( int i = 1; i < paraComm->getSize(); i++ )
   {
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, i, TagInterruptRequest )
      );
   }
   interruptIsRequested = true;
}

void
ParaLoadCoordinator::terminateAllSolvers(
      )
{
   for( int i = 1; i < paraComm->getSize(); i++ )
   {
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, i, TagTerminateRequest )
      );
   }
}

void
ParaLoadCoordinator::writeTransferLog(
      int rank
      )
{
   // output comp infomation to tree log file
   if( logNodesTransferFlag )
   {
      *osLogNodesTransfer << "[Solver-ID: " << rank
      << "] ParaNode was sent " << (paraSolverPool->getCurrentNode(rank))->toString() << std::endl;
   }
}

void
ParaLoadCoordinator::writeTransferLog(
      int rank,
      ParaCalculationState *state
      )
{
   // output comp infomation to tree log file
   if( logNodesTransferFlag )
   {
      *osLogNodesTransfer << "[Solver-ID: " << rank
      << "] Solved " << (paraSolverPool->getCurrentNode(rank))->toString() << std::endl;
      *osLogNodesTransfer << "[Solver-ID: " << rank
      << "] " << state->toString() << std::endl;
   }
}

void
ParaLoadCoordinator::writeTransferLogInRacing(
      int rank
      )
{
   // output comp infomation to tree log file
   if( logNodesTransferFlag )
   {
      *osLogNodesTransfer << "[Solver-ID: " << rank
      << "] ParaNode was sent " << (paraRacingSolverPool->getCurrentNode(rank))->toString() << std::endl;
   }
}

void
ParaLoadCoordinator::writeTransferLogInRacing(
      int rank,
      ParaCalculationState *state
      )
{
   // output comp infomation to tree log file
   if( logNodesTransferFlag )
   {
      *osLogNodesTransfer << "[Solver-ID: " << rank
      << "] Solved " << (paraRacingSolverPool->getCurrentNode(rank))->toString() << std::endl;
      *osLogNodesTransfer << "[Solver-ID: " << rank
      << "] " << state->toString() << std::endl;
   }
}

bool
ParaLoadCoordinator::updateSolution(
      ParaSolution *sol
      )
{
   if( !paraInitiator->getGlobalBestIncumbentSolution() )
      return paraInitiator->tryToSetIncumbentSolution(sol->clone(paraComm));
   if( sol->getObjectiveFuntionValue()
         < paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue() )
      return paraInitiator->tryToSetIncumbentSolution(sol->clone(paraComm));
   else
      return false;
}

void
ParaLoadCoordinator::sendIncumbentValue(
      int receivedRank
      )
{
   double globalBestIncumbentValue = paraInitiator->getGlobalBestIncumbentSolution()->getObjectiveFuntionValue();
   if( !paraParamSet->getBoolParamValue(NoUpperBoundTransferInRacing) || !isRacingStage() )
   {
      for( int i = 1; i < paraComm->getSize(); i++ )
      {
         if( i !=  receivedRank )
         {
            PARA_COMM_CALL(
                  paraComm->send( &globalBestIncumbentValue, 1, ParaDOUBLE, i, TagIncumbentValue )
            );
         }

      }
   }
   lcts.nDeletedInLc += paraNodePool->removeBoundedNodes(globalBestIncumbentValue);
}

void
ParaLoadCoordinator::inactivateRacingSolverPool(
      int rank
      )
{
   nSolvedInInterruptedRacingSolvers = paraRacingSolverPool->getNnodesSolvedInBestSolver();
   nNodesLeftInInterruptedRacingSolvers = paraRacingSolverPool->getNnodesLeftInBestSolver();
   if( paraRacingSolverPool->isActive(rank) )   // if rank is the winner, it should be inactive
   {
      paraRacingSolverPool->inactivateSolver(rank);
   }
   if( paraSolverPool->isSolverActive(rank) )
   {
      /*
       * special timining problem
       *
       * 1113.58 S.4 I.SOL 0
       * 1113.58 S.3 is the racing winner! Selected strategy 2.
       * 1113.58 S.4 >(TERMINATED_IN_RACING_STAGE)
       *
       */
      paraSolverPool->inactivateSolver(rank, 0, paraNodePool);
   }

   if ( (!interruptIsRequested) && paraRacingSolverPool->getNumActiveSolvers() == 0 )
   {
      /** if the computation is interrupted,
       *  paraRacingSolverPool is needed to output statistics */
      delete paraRacingSolverPool;
      paraRacingSolverPool = 0;
   }

}

void
ParaLoadCoordinator::sendTagToAllSolvers(
      const int tag
      )
{
   for( int i = 1; i < paraComm->getSize(); i++ )
   {
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, i, tag )
      );
   }
}

void
ParaLoadCoordinator::writeSubtreeInfo(
      int source,
      ParaCalculationState *calcState
      )
{
   if( logSubtreeInfoFlag )
   {
     ParaNode *node = paraSolverPool->getCurrentNode(source);
     *osLogSubtreeInfo  << paraTimer->getElapsedTime()
           << ", "
           << source
           << ", "
           << node->toSimpleString()
           << ", "
           << calcState->toSimpleString()
           << std::endl;
   }
}

void
ParaLoadCoordinator::initMergeNodesStructs(
      )
{
   varIndexTable = new ParaFixedValuePtr[paraInitiator->getParaInstance()->getNVars()];
   for( int i = 0; i < paraInitiator->getParaInstance()->getNVars(); i++ )
   {
      varIndexTable[i] = 0;
   }
   mergeInfoHead = 0;
   mergeInfoTail = 0;
}

void
ParaLoadCoordinator::addNodeToMergeNodeStructs(
      ParaNode *node
      )
{
   //
   // create mergeNodeInfo and linked to paraNode
   //
   ParaMergeNodeInfo *mNode = new ParaMergeNodeInfo();
   mNode->status = ParaMergeNodeInfo::PARA_MERGING;
   mNode->nSameValueVariables = -1;
   mNode->nMergedNodes = -1;
   mNode->keyIndex = -1;
   mNode->nFixedVariables = 0;
   mNode->fixedVariables = 0;
   mNode->mergedTo = 0;
   mNode->paraNode = node;
   mNode->origDiffSubproblem = 0;
   mNode->mergedDiffSubproblem = 0;
   mNode->next = 0;
   node->setMergingStatus(0);      // checking
   node->setMergeNodeInfo(mNode);  // set merge node info.

   //
   // make the same value fixed variables links
   //
   /* get fixed variables array */
   if( node->getDiffSubproblem() )
   {
      mNode->nFixedVariables = node->getDiffSubproblem()->getFixedVariables(
            paraInitiator->getParaInstance(),
            &(mNode->fixedVariables));
   }
   else
   {
      mNode->nFixedVariables = 0;
   }
   if( mNode->nFixedVariables == 0 )  // cannot merge!
   {
      delete mNode;
      node->setMergingStatus(3);     // cannot be merged
      node->setMergeNodeInfo(0);
      return;
   }

   //
   // add mergeNode to mergeNodeInfo list
   //
   if( mergeInfoTail == 0 )
   {
      mergeInfoTail = mNode;
      mergeInfoHead = mNode;
   }
   else
   {
      mergeInfoTail->next = mNode;
      mergeInfoTail = mNode;
   }

   for( int i = 0; i < mNode->nFixedVariables; i++ )
   {
      mNode->fixedVariables[i].mnode = mNode;
      ParaFixedValue *fixedValue = 0;
      if( varIndexTable[mNode->fixedVariables[i].index] == 0 )
      {
         fixedValue = new ParaFixedValue();
         fixedValue->value = mNode->fixedVariables[i].value;
         fixedValue->head = 0;
         fixedValue->tail = 0;
         fixedValue->next = 0;
         varIndexTable[mNode->fixedVariables[i].index] = fixedValue;
      }
      else
      {
         ParaFixedValue *prev = varIndexTable[mNode->fixedVariables[i].index];
         for( fixedValue = varIndexTable[mNode->fixedVariables[i].index];
               fixedValue != 0 &&  !EPSEQ( fixedValue->value, mNode->fixedVariables[i].value, DEFAULT_NUM_EPSILON );
               fixedValue = fixedValue->next )
         {
            prev = fixedValue;
         }
         if( fixedValue == 0 )
         {
            fixedValue = new ParaFixedValue();
            fixedValue->value = mNode->fixedVariables[i].value;
            fixedValue->head = 0;
            fixedValue->tail = 0;
            fixedValue->next = 0;
            prev->next = fixedValue;
         }
      }
      assert( fixedValue );
      if( fixedValue->tail == 0 )
      {
         fixedValue->head = &(mNode->fixedVariables[i]);
         fixedValue->tail = &(mNode->fixedVariables[i]);
      }
      else
      {
         fixedValue->tail->next = &(mNode->fixedVariables[i]);
         fixedValue->tail->next->prev = fixedValue->tail;
         fixedValue->tail = &(mNode->fixedVariables[i]);
      }
      for( ParaFixedVariable *p = fixedValue->head; p != fixedValue->tail; p = p->next )
      {
         (p->nSameValue)++;
      }
   }
}

void
ParaLoadCoordinator::generateMergeNodesCandidates(
      )
{
   ParaMergeNodeInfo *mPre = mergeInfoHead;
   ParaMergeNodeInfo *mNode = mergeInfoHead;
   mNode = mergeInfoHead;
   while( mNode )
   {
      assert( mNode->paraNode->getMergeNodeInfo() == mNode );
      if( mNode->status == ParaMergeNodeInfo::PARA_MERGING && mNode->nMergedNodes < 0 )
      {
         // make sorted variables list
         std::multimap<int, ParaSortedVariable, std::greater<int> > descendent;
         for( int i = 0; i < mNode->nFixedVariables; i++ )
         {
            ParaSortedVariable sortedVar;
            sortedVar.idxInFixedVariabes = i;
            sortedVar.fixedVariable = &(mNode->fixedVariables[i]);
            descendent.insert(std::make_pair(mNode->fixedVariables[i].nSameValue, sortedVar));

         }
         //
         //  try to make merge candidates
         //
         std::multimap<int, ParaSortedVariable, std::greater<int> >::iterator pos;
         pos = descendent.begin();
         mNode->keyIndex = pos->second.idxInFixedVariabes;
         mNode->nSameValueVariables = 1;
         ParaFixedVariable *traverse = mNode->fixedVariables[mNode->keyIndex].next;
         int nmNodes = 0;
         for( ;
               traverse;
               traverse = traverse->next )
         {
            if( traverse->mnode->status == ParaMergeNodeInfo::PARA_MERGING && traverse->mnode->nMergedNodes < 0 )
            {
               assert( traverse->mnode != mNode );
               traverse->mnode->mergedTo = mNode;
               traverse->mnode->nMergedNodes = 0;
               traverse->mnode->nSameValueVariables = 1;
               nmNodes++;
            }
         }
         ++pos;
         for( ; pos != descendent.end(); ++pos )
         {
            // check if there are merged nodes in case adding one more variable
            for( traverse = mNode->fixedVariables[pos->second.idxInFixedVariabes].next;
                  traverse;
                  traverse = traverse->next )
            {
               if( traverse->mnode->nMergedNodes == 0 && traverse->mnode->mergedTo == mNode )
               {
                  if( traverse->mnode->nSameValueVariables == mNode->nSameValueVariables )
                  {
                     break;   // at least one node can be merged
                  }
               }
            }
            if( traverse == 0 )  // cannot merge any nodes
            {
               break;
            }

            // merge nodes
            mNode->nSameValueVariables++;
            for( traverse = mNode->fixedVariables[pos->second.idxInFixedVariabes].next;
                  traverse;
                  traverse = traverse->next )
            {
               if( traverse->mnode->nMergedNodes == 0 && traverse->mnode->mergedTo == mNode )
               {
                  if( traverse->mnode->nSameValueVariables == (mNode->nSameValueVariables - 1) )
                  {
                     traverse->mnode->nSameValueVariables++;
                  }
               }
            }
         }
         // if the number of fixed variables is too small, then the merged node is not created
         if( nmNodes < 2 ||     // no merging nodes
               static_cast<int>((mNode->nFixedVariables)*paraParamSet->getRealParamValue(FixedVariablesRatioInMerging)) < 1 || //  0 same value variables are not allowed
               mNode->nSameValueVariables < (mNode->nFixedVariables)*paraParamSet->getRealParamValue(FixedVariablesRatioInMerging) )
         {
            for( ParaFixedVariable *cleanup = mNode->fixedVariables[mNode->keyIndex].next;
                  cleanup;
                  cleanup = cleanup->next )
            {
               if( cleanup->mnode->mergedTo == mNode )
               {
                  cleanup->mnode->nSameValueVariables = -1;
                  cleanup->mnode->nMergedNodes = -1;
                  cleanup->mnode->keyIndex = -1;
                  cleanup->mnode->mergedTo = 0;
                  assert( cleanup->mnode->status == ParaMergeNodeInfo::PARA_MERGING );
               }
            }
            assert( !(mNode->origDiffSubproblem) );
            assert( !(mNode->mergedDiffSubproblem) );
            mNode->paraNode->setMergeNodeInfo(0);
            mNode->paraNode->setMergingStatus(3);   // cannot merged
            ParaMergeNodeInfo *doomed = mNode;
            if( mNode == mergeInfoHead )
            {
               mergeInfoHead = mNode->next;
               mPre = mergeInfoHead;
               mNode = mergeInfoHead;
            }
            else
            {
               mPre->next = mNode->next;
               mNode = mNode->next;
            }
            if( mNode == mergeInfoTail )
            {
               mergeInfoTail = mPre;
            }
            deleteMergeNodeInfo(doomed);
         }
         else  // cleanup and merge nodes
         {
            int nMergedNodes = 0;
            for( ParaFixedVariable *cleanup = mNode->fixedVariables[mNode->keyIndex].next;
                  cleanup;
                  cleanup = cleanup->next )
            {
               if( cleanup->mnode->mergedTo == mNode )
               {
                  if( mNode->nSameValueVariables == cleanup->mnode->nSameValueVariables )
                  {
                     nMergedNodes++;
                     cleanup->mnode->status = ParaMergeNodeInfo::PARA_MERGE_CHECKING_TO_OTHER_NODE;
                  }
                  else
                  {
                     assert( cleanup->mnode->status == ParaMergeNodeInfo::PARA_MERGING );
                     cleanup->mnode->nSameValueVariables = -1;
                     cleanup->mnode->nMergedNodes = -1;
                     cleanup->mnode->keyIndex = -1;
                     cleanup->mnode->mergedTo = 0;
                  }
               }
            }
            mNode->nMergedNodes = nMergedNodes;
            assert(nMergedNodes > 0);
            int n = 0;
            ParaFixedVariable *fixedVariables = new ParaFixedVariable[mNode->nSameValueVariables];
            for( pos = descendent.begin(); pos != descendent.end(); ++pos )
            {
               fixedVariables[n] = *(pos->second.fixedVariable);
               n++;
               if( n == mNode->nSameValueVariables ) break;
            }
            mNode->origDiffSubproblem = mNode->paraNode->getDiffSubproblem();
            mNode->mergedDiffSubproblem = mNode->origDiffSubproblem->createDiffSubproblem(paraComm, n, fixedVariables );
            delete [] fixedVariables;
            mNode->paraNode->setDiffSubproblem(mNode->mergedDiffSubproblem);
            mNode->status = ParaMergeNodeInfo::PARA_MERGED_RPRESENTATIVE;
            assert( mNode->mergedTo == 0 );
            mPre = mNode;
            mNode = mNode->next;
         }
      }
      else
      {
         mPre = mNode;
         mNode = mNode->next;
      }
   }

   // remove data which are not used anymore.
   if( varIndexTable )
   {
      for( int i = 0; i < paraInitiator->getParaInstance()->getNVars(); i++ )
      {
         if( varIndexTable[i] )
         {
            while ( varIndexTable[i] )
            {
               ParaFixedValue *del = varIndexTable[i];
               varIndexTable[i] = varIndexTable[i]->next;
               delete del;
            }
         }
      }
      delete [] varIndexTable;
      varIndexTable = 0;
   }
}

void
ParaLoadCoordinator::regenerateMergeNodesCandidates(
      ParaNode *node
      )
{
   ParaMergeNodeInfo *mNode = node->getMergeNodeInfo();
   assert(mNode);
   assert(mNode->paraNode == node);
   assert(mNode->status == ParaMergeNodeInfo::PARA_MERGED_RPRESENTATIVE);
   assert(mNode->mergedTo == 0);
   node->setMergeNodeInfo(0);
   node->resetDualBoundValue();
   node->setDiffSubproblem(mNode->origDiffSubproblem);
   delete mNode->mergedDiffSubproblem;
   mNode->mergedDiffSubproblem = 0;
   assert( mNode->status != ParaMergeNodeInfo::PARA_MERGING );
   // set new range
   mergeInfoHead =  0;
   mergeInfoTail = 0;
   ParaMergeNodeInfo *mPrev = 0;
   for( ParaFixedVariable *traverse = mNode->fixedVariables[mNode->keyIndex].next;
         traverse;
         traverse = traverse->next )
   {
      if( mergeInfoTail )
      {
         mPrev->next = traverse->mnode;
         mergeInfoTail = traverse->mnode;
         mPrev = traverse->mnode;
      }
      if( mNode == traverse->mnode->mergedTo )
      {
         if( !mergeInfoHead )
         {
            mergeInfoHead = traverse->mnode;
            mergeInfoTail = traverse->mnode;
            mPrev = traverse->mnode;
         }
         assert( traverse->mnode->status == ParaMergeNodeInfo::PARA_MERGE_CHECKING_TO_OTHER_NODE );
      }
   }
   mergeInfoTail->next = 0;
   // remove mnode
   mNode->paraNode->setMergingStatus(-1);  // no merging node
   deleteMergeNodeInfo(mNode);
   generateMergeNodesCandidates();
}

void
ParaLoadCoordinator::deleteMergeNodeInfo(
      ParaMergeNodeInfo *mNode
      )
{

   if( mNode->nMergedNodes == 0 && mNode->mergedTo )
   {
      assert(mNode->status == ParaMergeNodeInfo::PARA_MERGE_CHECKING_TO_OTHER_NODE);
      assert(mNode->mergedTo->status ==  ParaMergeNodeInfo::PARA_MERGED_RPRESENTATIVE);
      assert(mNode->mergedTo->mergedTo == 0);
      mNode->mergedTo->nMergedNodes--;
      if( mNode->mergedTo->nMergedNodes == 0 && mNode->mergedTo->paraNode->getMergeNodeInfo() )
      {
         mNode->mergedTo->paraNode->setDiffSubproblem(mNode->origDiffSubproblem);
         mNode->mergedTo->paraNode->setMergeNodeInfo(0);
         mNode->mergedTo->paraNode->setMergingStatus(-1);
         delete mNode->mergedDiffSubproblem;
         deleteMergeNodeInfo(mNode->mergedTo);
      }
      mNode->mergedTo = 0;

   }

   if( mNode->status == ParaMergeNodeInfo::PARA_MERGED_RPRESENTATIVE)
   {
      assert( mNode->mergedTo == 0 );
      if( mNode->paraNode->getMergingStatus() == -1 )  // merging failed
      {
         for( ParaFixedVariable *traverse = mNode->fixedVariables[mNode->keyIndex].next;
               traverse;
               traverse = traverse->next )
         {
            if( traverse->mnode->nMergedNodes == 0 && mNode == traverse->mnode->mergedTo )
            {
               traverse->mnode->mergedTo->nMergedNodes--;
               if( traverse->mnode->mergedTo->nMergedNodes == 0 && traverse->mnode->mergedTo->paraNode->getMergeNodeInfo() )
               {
                  traverse->mnode->mergedTo->paraNode->setDiffSubproblem(mNode->origDiffSubproblem);
                  traverse->mnode->mergedTo->paraNode->setMergeNodeInfo(0);
                  traverse->mnode->mergedTo->paraNode->setMergingStatus(-1);
                  delete traverse->mnode->mergedDiffSubproblem;
                  deleteMergeNodeInfo(traverse->mnode->mergedTo);
               }
               traverse->mnode->mergedTo = 0;
               traverse->mnode->paraNode->setMergingStatus(0);
               traverse->mnode->status = ParaMergeNodeInfo::PARA_MERGING;
               traverse->mnode->nMergedNodes = -1;
               traverse->mnode->nSameValueVariables = -1;
               traverse->mnode->keyIndex = -1;
            }
         }
      }
      else
      {
         for( ParaFixedVariable *traverse = mNode->fixedVariables[mNode->keyIndex].next;
               traverse;
               traverse = traverse->next )
         {
            if( traverse->mnode->nMergedNodes == 0 && mNode == traverse->mnode->mergedTo )
            {
               traverse->mnode->mergedTo->nMergedNodes--;
               if( traverse->mnode->mergedTo->nMergedNodes == 0 && traverse->mnode->mergedTo->paraNode->getMergeNodeInfo() )
               {
                  traverse->mnode->mergedTo->paraNode->setDiffSubproblem(mNode->origDiffSubproblem);
                  traverse->mnode->mergedTo->paraNode->setMergeNodeInfo(0);
                  traverse->mnode->mergedTo->paraNode->setMergingStatus(-1);
                  delete traverse->mnode->mergedDiffSubproblem;
                  deleteMergeNodeInfo(traverse->mnode->mergedTo);
               }
               traverse->mnode->mergedTo = 0;
               if( traverse->mnode->paraNode->getDualBoundValue() < mNode->paraNode->getDualBoundValue() )
               {
                  traverse->mnode->paraNode->setMergingStatus(0);
                  traverse->mnode->status = ParaMergeNodeInfo::PARA_MERGING;
                  traverse->mnode->nMergedNodes = -1;
                  traverse->mnode->nSameValueVariables = -1;
                  traverse->mnode->keyIndex = -1;
               }
               else
               {
                  traverse->mnode->paraNode->setMergingStatus(4);  // merging representative was deleted -> this node should be deleted
                  traverse->mnode->status = ParaMergeNodeInfo::PARA_DELETED;
                  traverse->mnode->nMergedNodes = -1;
                  traverse->mnode->nSameValueVariables = -1;
                  traverse->mnode->keyIndex = -1;
               }
            }
         }
      }
   }

   if( mNode->fixedVariables )
   {
      for( int i = 0; i < mNode->nFixedVariables; i++ )
      {
         for( ParaFixedVariable *traverse = mNode->fixedVariables[i].prev;
               traverse;
               traverse = traverse->prev
               )
         {
            traverse->nSameValue--;
         }
         if( mNode->fixedVariables[i].prev )
         {
            mNode->fixedVariables[i].prev->next = mNode->fixedVariables[i].next;
            if( mNode->fixedVariables[i].next )
            {
               mNode->fixedVariables[i].next->prev = mNode->fixedVariables[i].prev;
            }
            else
            {
               mNode->fixedVariables[i].prev->next = 0;
            }
         }
         else  // prev == 0
         {
            if( mNode->fixedVariables[i].next )
            {
               mNode->fixedVariables[i].next->prev = 0;
            }
         }

      }
      delete [] mNode->fixedVariables;
   }
   delete mNode;
}

void
ParaLoadCoordinator::mergeNodes(
      ParaNode *node
      )
{
   ParaMergeNodeInfo *mNode = node->getMergeNodeInfo();
   assert( mNode->status == ParaMergeNodeInfo::PARA_MERGED_RPRESENTATIVE );
   assert( mNode->mergedTo == 0 );
   ParaMergedNodeListElement *head = 0;
   ParaMergedNodeListElement *cur = 0;
   int nMerged = 0;
   for( ParaFixedVariable *traverse = mNode->fixedVariables[mNode->keyIndex].next;
         traverse;
         traverse = traverse->next )
   {
      if( traverse->mnode->nMergedNodes == 0 && mNode == traverse->mnode->mergedTo )
      {
         if( head == 0 )
         {
            head = cur = new ParaMergedNodeListElement();
         }
         else
         {
            cur->next = new ParaMergedNodeListElement();
            cur = cur->next;
         }
         cur->node = traverse->mnode->paraNode;
         cur->node->setMergingStatus(2);
         cur->next = 0;
         nMerged++;
      }
   }
   assert( mNode->nMergedNodes == nMerged);
   int delNodes = 0;
   if ( head )
   {
      delNodes = paraNodePool->removeMergedNodes(head);
   }
   assert( delNodes == nMerged );
   if( delNodes != nMerged )
   {
	   THROW_LOGICAL_ERROR4("delNodes != nMerged, delNodes = ", delNodes, ", nMerged = ", nMerged );
   }
   node->setDiffSubproblem(mNode->mergedDiffSubproblem);
   delete mNode->origDiffSubproblem;
   node->setMergeNodeInfo(0);
   node->setMergingStatus(1);
   deleteMergeNodeInfo(mNode);
   lcts.nDeletedByMerging += nMerged;
   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus << (nMerged + 1) <<
            " nodes are merged at " <<
            paraTimer->getElapsedTime() << " seconds." << std::endl;
   }
}
