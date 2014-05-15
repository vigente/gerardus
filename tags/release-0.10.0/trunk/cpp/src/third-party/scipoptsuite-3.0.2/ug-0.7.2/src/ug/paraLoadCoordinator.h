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

/**@file    paraLoadCoordinator.h
 * @brief   Load Coordinator.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_LOADCOORDINATOR_H__
#define __PARA_LOADCOORDINATOR_H__

#include <list>
#include <queue>
#include "paraDef.h"
#include "paraComm.h"
#include "paraCalculationState.h"
#include "paraLoadCoordinatorTerminationState.h"
#include "paraNode.h"
#include "paraNodePool.h"
#include "paraSolverState.h"
#include "paraSolverPool.h"
#include "paraInitiator.h"
#include "paraTimer.h"
#include "paraDeterministicTimer.h"

namespace UG
{

enum RunningPhase
{
   RampUpPhase,
   NormalRunningPhase,
   TerminationPhase               // includes interrupting phase
};

static const double displayInfOverThisValue = 5.0;

class ParaLoadCoordinator
{
   int            (ParaLoadCoordinator::*messageHandler[N_TAGS])(int source, int tag);
   int            (ParaLoadCoordinator::*racingRampUpMessageHandler[N_TAGS])(int source, int tag);

   int            globalSubtreeIdGen;
   ParaComm       *paraComm;
   ParaParamSet   *paraParamSet;
   ParaInitiator  *paraInitiator;
   bool           *racingSolversExist;
   bool           restarted;                    /**< indicates that this run is restarted from checkpoint files */
   bool           initialNodesGenerated;        /**< indicates that initial nodes have been generated */
   /** Status of LoadCoordinator */
   RunningPhase   runningPhase;
   int            firstCollectingModeState;     /**< status of first collecting mode
                                                            -1 : have not been in collecting mode
                                                             0 : once in collecting mode
                                                             1 : collecting mode is terminated once. */
   bool           isCollectingModeRestarted;     /**< this flag indicate if a collecting mode is restarted or not */
   bool           isBreakingFinised;             /**< indicate that breaking is finished or not
                                                           if bootstrap ramp-up is not specified, this flag should be always true */
   int            breakingSolverId;             /**< all nodes collecting solver Id: -1: no collecting */
   int            nReplaceToBetterNode;         /**< the number of replacing to a better nodes */

   bool           computationIsInterrupted;     /**< indicate that current computation is interrupted or not */
   bool           interruptedFromControlTerminal; /**< interrupted from control terminal        */
   bool           hardTimeLimitIsReached;       /**< indicate that hard time limit is reached or not */
   bool           winnerSolverNodesCollected;   /**< indicate that all winner solver nodes has been collected */
   bool           interruptIsRequested;         /**< indicate that all solver interrupt message is requested or not */

   /** Pools in LoadCorrdinator  */
   ParaSolverPool   *paraSolverPool;
   ParaNodePool     *paraNodePool;
   ParaRacingSolverPool *paraRacingSolverPool;
   long long        nSolvedInInterruptedRacingSolvers;
   long long        nNodesLeftInInterruptedRacingSolvers;

   /** For checkpoint */
   double           previousCheckpointTime;       /**< previous checkpoint time */
   char             lastCheckpointTimeStr[26];    /**< lastCheckpointTimeStr[0] == ' ' means no checkpoint */

   /** To measure how long does node pool stay in empty situation */
   double          statEmptyNodePoolTime;         /** start time that node pool become empty. initialized by max double */

   /** epsilon **/
   double            eps;  /**< absolute values smaller than this are considered zero  */

   /** LoadCoordinatorTerminationState: counters and time */
   ParaLoadCoordinatorTerminationState lcts;

   /** racing winner */
   int                racingWinner;
   int                minDepthInWinnerSolverNodes;
   int                maxDepthInWinnerSolverNodes;
   ParaRacingRampUpParamSet *racingWinnerParams;

   /** racing termination flag */
   bool               racingTermination;
   int                nSolvedRacingTermination;   /**< number of nodes solved at racing termination solver */

   /***********************/
   /** for merging nodes  */
   /***********************/
   bool              merging;                     /**< indicate that merging is processing */
   /** The followings are used temporary to generate merge nodes info */
   ParaFixedValue    **varIndexTable;             /**< variable indices table. */
   ParaMergeNodeInfo *mergeInfoHead;              /**< head of ParaMergeNodeInfo list */
   ParaMergeNodeInfo *mergeInfoTail;              /**< tail of ParaMergeNodeInfo list */

   /** counter to check all solver termination */
   int                nTerminated;

   /** Timer for LoadCoordinator */
   ParaTimer          *paraTimer;
   ParaDeterministicTimer *paraDetTimer;    /**< deterministic timer for this ParaSolver */

   /** Keep tabular solving status output time */
   double             previousTabularOutputTime;

   bool               logSolvingStatusFlag;
   std::ofstream      ofsLogSolvingStatus;
   std::ostream       *osLogSolvingStatus;
   bool               logNodesTransferFlag;
   std::ofstream      ofsLogNodesTransfer;
   std::ostream       *osLogNodesTransfer;
   bool               outputTabularSolvingStatusFlag;
   std::ofstream      ofsTabularSolvingStatus;
   std::ostream       *osTabularSolvingStatus;
   bool               logSubtreeInfoFlag;
   std::ofstream      ofsLogSubtreeInfo;
   std::ostream       *osLogSubtreeInfo;
   std::ofstream      ofsStatisticsFinalRun;
   std::ostream       *osStatisticsFinalRun;
   std::ofstream      ofsStatisticsRacingRampUp;
   std::ostream       *osStatisticsRacingRampUp;

   void            writeTransferLog(int rank, ParaCalculationState *state);
   void            writeTransferLog(int rank);
   void            writeTransferLogInRacing(int rank, ParaCalculationState *state);
   void            writeTransferLogInRacing(int rank);
   void            sendRampUpToAllSolvers();
   void            sendRetryRampUpToAllSolvers();
   void            sendInterruptRequest();
   void            terminateAllSolvers();
   bool            updateSolution(ParaSolution *);
   void            sendIncumbentValue(int receivedRank);
   int             createNewGlobalSubtreeId(){ return ++globalSubtreeIdGen; }

   /** message handlers */
   int              processTagNode(int source, int tag);
   int              processTagSolution(int source, int tag);
   int              processTagSolverState(int source, int tag);
   int              processTagCompletionOfCalculation(int source, int tag);
   int              processTagAnotherNodeRequest(int source, int tag);
   int              processTagTerminated(int source, int tag);
   int              processTagHardTimeLimit(int source, int tag);
   int              processTagInitialStat(int source, int tag);
   int              processTagToken(int source, int tag);

   /** message handlers specialized for racing ramp-up */
   int              processRacingRampUpTagSolverState(int source, int tag);
   int              processRacingRampUpTagCompletionOfCalculation(int source, int tag);

   void            run();

   void            updateCheckpointFiles();
   void            writeLoadCoordinatorStatisticsToCheckpointFile(
                       ogzstream &loadCoordinatorStatisticsStream, int nSolverInfo,
                       double globalBestDualBoundValue, double externalGlobalBestDualBoundValue );
   bool            sendParaNodesToIdleSolvers();
   bool            isRacingStage()
                   {
                      if( paraParamSet->getIntParamValue(RampUpPhaseProcess) > 0 &&  racingWinner < 0 )
                         return true;
                      else
                         return false;
                   }
   void            writePreviousStatisticsInformation();
   void            inactivateRacingSolverPool(int rank);
   void            outputTabularSolvingStatus(char incumbent);
   void            sendTagToAllSolvers(const int tag);
   void            writeSubtreeInfo(int source,  ParaCalculationState *calcState);

   void            initMergeNodesStructs();
   void            addNodeToMergeNodeStructs(ParaNode *node);
   void            generateMergeNodesCandidates();
   void            regenerateMergeNodesCandidates(ParaNode *node);
   void            mergeNodes(ParaNode *node);
   void            deleteMergeNodeInfo(ParaMergeNodeInfo *mNode);

public:
   ParaLoadCoordinator(ParaComm *inComm, ParaParamSet *inParaParamSet,
         ParaInitiator *paraInitiator, bool *racingSolversExist, ParaTimer *paraTimer, ParaDeterministicTimer *detTimer);
   ~ParaLoadCoordinator();
   void interrupt(){ interruptedFromControlTerminal = true; }
   void warmStart();
   void run(ParaNode *paraNode);
   void run(ParaNode *paraNode, int nRacingSolvers, ParaRacingRampUpParamSet **racingRampUpParams);
   void setGlobalBestIncumbentSolution(ParaSolution *sol);
};

}

#endif // __PARA_LOADCOORDINATOR_H__

