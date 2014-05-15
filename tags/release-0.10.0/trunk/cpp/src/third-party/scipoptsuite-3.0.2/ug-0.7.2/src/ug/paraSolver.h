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

/**@file    paraSolver.h
 * @brief   Base class for solver: Generic parallelized solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_SOLVER_H__
#define __PARA_SOLVER_H__

#include "paraDef.h"
#include "paraComm.h"
#include "paraParamSet.h"
#include "paraRacingRampUpParamSet.h"
#include "paraNode.h"
#include "paraTimer.h"
#include "paraDeterministicTimer.h"
#include "paraSolution.h"

#define ENFORCED_THRESHOLD 5

namespace UG
{

static const int NoTerminationMode          = 0;
static const int NormalTerminationMode      = 1;
static const int InterruptedTerminationMode = 2;

class ParaSolver
{
protected:

   int     (ParaSolver::*messageHandler[N_TAGS])(int source, int tag);
   unsigned int notificationIdGenerator;

   ParaComm     *paraComm;                  /**< ParaCommunicator object */
   ParaParamSet *paraParams;                /**< ParaParamSet object */
   ParaRacingRampUpParamSet *racingParams;  /**< ParaRacingRampUpParamSet object.
                                                 This is also a flag to indicate running with racing ramp-up */
   ParaRacingRampUpParamSet *winnerRacingParams; /**< Winner ParaRacingRampUpParamSet object */

   ParaTimer    *paraTimer;                 /**< timer for this ParaSolver */
   ParaDeterministicTimer *paraDetTimer;    /**< deterministic timer for this ParaSolver */

   double       globalBestDualBoundValueAtWarmStart;   /**< global best dual bound value which is set when system warm starts */
   double       globalBestIncumbentValue;   /**< global best incumbent value */
   double       lcBestDualBoundValue;       /**< LoadCoordinator best dual bound value */
   ParaSolution *globalBestIncumbentSolution;   /**< global best solution. However, this is not always feasible for the current sub-MIP */
   ParaSolution *localIncumbentSolution;    /**< incumbent solution generated in local solver */
   ParaSolution *pendingSolution;           /**< solution which is pending to update in case of deterministic runs */
   double        pendingIncumbentValue;     /**< incumbent value which is pending to update in case of deterministic runs */
   ParaInstance *paraInstance;              /**< root problem instance */
   ParaNode     *currentNode;               /**< solving node */
   ParaNode     *newNode;                   /**< new node to solve */
   bool         collectingMode;             /**< indicate whether if this solver is in collecting mode or not */
   bool         aggressiveCollecting;       /**< indicate that if this solver has two nodes, this solver sends one to LC */
   int          nSendInCollectingMode;      /**< number of nodes need to send in collecting mode */
   int          nCollectOnce;               /**< number of nodes need to collect once */
   bool         collectingManyNodes;        /**< indicate that many nodes collecting is requested by LC */
   int          terminationMode;            /**< indicate that termination mode */
                                            /**<     0: no termination mode    */
                                            /**<     1: normal termination mode */
                                            /**<     2: interrupted termination */
   bool         warmStarted;                /**< indicate whether if system is warm started or not */
   bool         rampUp;                     /**< indicate whether if ramp-up phase is finished or not: true - finish */
   bool         racingIsInterrupted;        /**< indicate whether if racing phases is interrupted or not: true - interrupted */
   bool         racingWinner;               /**< indicate racing ramp-up winner or not: true - winner */
   bool         anotherNodeIsRequested;     /**< indicate that another node is requested or not */
   bool         waitingSpecificMessage;     /**< indicate that this solver is waiting for a specific message */
   bool         lightWeightRootNodeComputation; /**< indicate that fast root node computation is required */

   bool         onceBreak;                  /**< indicate that the sub-MIP is broken down once */

   /** Times */
   double       rootNodeTime;               /**< root node process time of current ParaNode */
   double       totalRootNodeTime;          /**< accumulated root node process time solved by this solver so far */
   double       minRootNodeTime;            /**< minimum time consumed by root node process */
   double       maxRootNodeTime;            /**< maximum time consumed by root node process */
   double       previousNotificationTime;   /**< previous notification time */
   double       paraNodeStartTime;          /**< start time of current ParaNode */
   /** Idle Times */
   double       previousStopTime;           /**< previous stop solving time of this Solver: For measurement */
   double       idleTimeToFirstParaNode;    /**< idle time to start solving the first ParaNode */
   double       idleTimeBetweenParaNodes;   /**< idle time between ParaNodes processing */
   double       idleTimeAfterLastParaNode;  /**< idle time after the last ParaNode was solved */
   double       idleTimeToWaitNotificationId;   /**< idle time to wait a message within collecting mode */
   double       idelTimeToWaitAckCompletion;    /**< idle time to wait acknowledgment of completion */
   double       idleTimeToWaitToken;         /**< idle time to wait token */
   double       previousIdleTimeToWaitToken; /**< previous idle time to wait token */
   double       offsetTimeToWaitToken;       /**< offset time to wait token */

   /** Counters related to the current ParaNode */
   int          nSolved;                    /**< number of nodes solved, that is, number of subtree nodes rooted from ParaNode */
   int          nSent;                      /**< number of ParaNodes sent from this subtree rooted from the current ParaNode */
   int          nImprovedIncumbent;         /**< number of improvements of incumbent value */
   int          nSolvedWithNoPreprocesses;  /**< number of nodes solved when it is solved with no preprocesses */
   /** Counters related to this ParaSolver */
   int          totalNSolved;               /**< accumulated number of nodes solved in this ParaSolver */
   int          minNSolved;                 /**< minimum number of subtree nodes rooted from ParaNode */
   int          maxNSolved;                 /**< maximum number of subtree nodes rooted from ParaNode */
   int          totalNSent;                 /**< accumulated number of nodes sent from this ParaSolver */
   int          totalNImprovedIncumbent;    /**< accumulated number of improvements of incumbent value in this ParaSolver */
   int          nParaNodesReceived;         /**< number of ParaNodes received in this ParaSolver */
   int          nParaNodesSolved;           /**< number of ParaNodes solved ( received ) in this ParaSolver */
   int          nParaNodesSolvedAtRoot;     /**< number of ParaNodes solved at root node */
   int          nParaNodesSolvedAtPreCheck; /**< number of ParaNodes solved at pre-checking of root node solvability */

   int          nSimplexIterRoot;           /**< number of simplex iteration at root node */
   // long long    totalSimplexIter;           /**< total number of simplex iteration */
   double       minIisum;                   /**< minimum sum of integer infeasibility */
   double       maxIisum;                   /**< maximum sum of integer infeasibility */
   int          minNii;                     /**< minimum number of integer infeasibility */
   int          maxNii;                     /**< maximum number of integer infeasibility */

   bool         globalIncumbnetValueUpdateFlag; /**< indicate that global incumbent value is updated in iReceiveMessages() routine */
   bool         notificationProcessed;          /**< if true, notification is issued but not receive the corresponding LCB */

   double       eps;                            /**< absolute values smaller than this are considered zero  */
                                                /**< esp should be set in the constructor of the derived class of ParaSolver */

   double       targetBound;               /**< target bound value for breaking */
   int          nTransferLimit;             /**< limit number of transferring nodes for breaking */
   int          nTransferredNodes;          /**< keep track number of transferred nodes for breaking  */
   /* for deterministic execution */
   double       previousCommTime;           /**< previous communication time for deterministic execution */

   pthread_mutex_t   tokenAccessLock;       /**< if a solver is multi-threaded, several threads manipulate token at the same time */

   int processTagNode(int source, int tag);
   int processTagNodeReceived(int source, int tag);
   int processTagRampUp(int source, int tag);
   int processTagRetryRampUp(int source, int tag);
   int processTagSolution(int source, int tag);
   int processTagIncumbentValue(int source, int tag);
   int processTagGlobalBestDualBoundValueAtWarmStart(int source, int tag);
   int processTagNoNodes(int source, int tag);
   int processTagInCollectingMode(int source, int tag);
   int processTagInEnforcedCollectingMode(int source, int tag);
   int processTagCollectAllNodes(int source, int tag);
   int processTagOutCollectingMode(int source, int tag);
   int processTagLCBestBoundValue(int source, int tag);
   int processTagNotificationId(int source, int tag);
   int processTagTerminateRequest(int source, int tag);
   int processTagInterruptRequest(int source, int tag);
   int processTagWinnerRacingRampUpParamSet(int source, int tag);
   int processTagWinner(int source, int tag);
   int processTagLightWeightRootNodeProcess(int source, int tag);
   int processTagBreaking(int source, int tag);
   int processTagToken(int source, int tag);

   bool receiveNewNodeAndReactivate();
   void waitNotificationIdMessage();
   void waitAckCompletion();

   void sendCompletionOfCalculation(double stopTime);
   bool updateGlobalBestIncumbentSolution(ParaSolution *sol);
   bool updateGlobalBestIncumbentValue(double newValue);

   virtual void setRacingParams(ParaRacingRampUpParamSet *racingParms, bool winnerParam) = 0;
   virtual void setWinnerRacingParams(ParaRacingRampUpParamSet *racingParms) = 0;
   virtual void createSubproblem() = 0;
   virtual void freeSubproblem() = 0;
   virtual void solve() = 0;
   virtual long long getNNodesSolved() = 0;
   virtual int getNNodesLeft() = 0;
   virtual double getDualBoundValue() = 0;

   virtual void solveToCheckEffectOfRootNodePreprocesses(){ /* set nSolvedWithNoPreprocesses */ }

public:
   ParaSolver(){ THROW_LOGICAL_ERROR1("Default constructor of ParaSolver is called"); }
   ParaSolver( int argc, char **argv,
         ParaComm *comm, ParaParamSet *inParaParamSet, ParaInstance *paraInstance, ParaDeterministicTimer *detTimer );
   virtual ~ParaSolver();

   bool isWarmStarted(){ return warmStarted; }
   void run();

   void run(
         ParaNode *paraNode
         )
   {
      currentNode = paraNode;
      run();
   }

   void run(
         ParaRacingRampUpParamSet *inRacingRampUpParamSet
         )
   {
      ParaNode *rootNode = paraComm->createParaNode();
      PARA_COMM_CALL(
            rootNode->bcast(paraComm, 0)
            );
      nParaNodesReceived++;
      racingParams = inRacingRampUpParamSet;
      setRacingParams(racingParams, false);
      if( paraParams->getBoolParamValue(Deterministic) )
      {
         do
         {
            iReceiveMessages();
         } while( !waitToken(paraComm->getRank()) );
      }
      iReceiveMessages();   // Feasible solution may be received.
      if( paraParams->getBoolParamValue(Deterministic) )
      {
         passToken(paraComm->getRank());
      }
      run( rootNode );
   }

   /** the following functions may be called from callback routines of the target Solver */
   double getElapsedTimeOfNodeSolving(){ return (paraTimer->getElapsedTime() - paraNodeStartTime); }
   double getGlobalBestDualBoundValueAtWarmStart(){ return globalBestDualBoundValueAtWarmStart; }
   double getLcBestDualBoundValue(){ return lcBestDualBoundValue; }
   int getNStopSolvingMode(){ return paraParams->getIntParamValue(NStopSolvingMode);}
   double getRootNodeTime(){ return rootNodeTime; }
   double getBoundGapForStopSolving(){ return paraParams->getRealParamValue(BgapStopSolvingMode); }
   double getBoundGapForCollectingMode(){ return paraParams->getRealParamValue(BgapCollectingMode); }
   void iReceiveMessages();
   bool isGlobalIncumbentUpdated(){ return globalIncumbnetValueUpdateFlag; }
   void globalIncumbnetValueIsReflected() { globalIncumbnetValueUpdateFlag = false; }
   bool isRampUp(){ return rampUp; }
   bool isRacingWinnerParamReceived(){ return ( winnerRacingParams != 0 ); }
   bool isRacingRampUp()
   {
      return ( ( paraParams->getIntParamValue(RampUpPhaseProcess) == 1 ) ||
            ( paraParams->getIntParamValue(RampUpPhaseProcess) == 2 ) );
   }
   bool isRacingWinner(){ return racingWinner; }
   void setRootNodeTime();
   void sendIfImprovedSolutionWasFound(ParaSolution *sol);
   void saveIfImprovedSolutionWasFound(ParaSolution *sol);
   void sendLocalSolution();
   bool notificationIsNecessary();
   void sendSolverState(long long nNodesSolved, int nNodesLeft, double bestDualBoundValue, double detTime);
   bool isInterrupting(){ return ( terminationMode == InterruptedTerminationMode ); }
   bool isTerminationRequested(){ return  ( terminationMode == NormalTerminationMode ); }
   bool newParaNodeExists(){ return (newNode != 0); }
   bool isInCollectingMode(){ return ( collectingMode || collectingManyNodes ); }
   bool isAggressiveCollecting(){ return aggressiveCollecting; }
   bool isManyNodesCollectionRequested(){ return collectingManyNodes; }
   int  getThresholdValue(
         int nNodes      /**< number of processed nodes, including the focus node */
         );
   void sendParaNode(long long n, int depth, double dualBound, double estimateValue, ParaDiffSubproblem *diffSubproblem);
   void sendAnotherNodeRequest(double bestDualBoundValue);
   bool getNotificaionProcessed(){ return notificationProcessed; }
   double getGlobalBestIncumbentValue(){ return globalBestIncumbentValue; }
   ParaNode *getCurrentNode(){ return currentNode; }
   ParaInstance *getParaInstance(){ return paraInstance; }
   ParaParamSet *getParaParamSet(){ return paraParams; }
   int getRank(){ return paraComm->getRank(); }
   void countInPrecheckSolvedParaNodes(){ nParaNodesSolvedAtPreCheck++; }
   void waitMessageIfNecessary();
   int getNSendInCollectingMode(){ return nSendInCollectingMode; }
   bool isRacingStage()
   {
      return (racingParams &&
            (paraParams->getIntParamValue(RampUpPhaseProcess) == 1 ||
             paraParams->getIntParamValue(RampUpPhaseProcess) == 2 ) );
   }
   void terminateRacing()
   {
      assert(racingParams);
      delete racingParams;
      racingParams = 0;
      racingIsInterrupted = false;    // rampUp message might have been received before terminate racing
                                      // Then, this flag should be set false
   }

   ParaSolution *getGlobalBestIncumbentSolution(
         )
   {
      return globalBestIncumbentSolution;
   }

   bool isWaitingForSpecificMessage(){ return waitingSpecificMessage; }
   bool isBreaking(){ return (nTransferLimit > 0); }
   double getTargetBound(){ return targetBound; }
   bool isTransferLimitReached(){ return (nTransferredNodes >= nTransferLimit); }
   void resetBreakingInfo(
         )
   {
      targetBound = -DBL_MAX;
      nTransferLimit = -1;
      nTransferredNodes = -1;
      collectingManyNodes = false;
   }

   bool isOnceBreak(){ return onceBreak; }
   void setOnceBreak()
   {
      nCollectOnce = -1;
      collectingManyNodes = true;
      onceBreak = true;
   }
   bool isAggressivePresolvingSpecified()
   {
      return ( paraParams->getIntParamValue(AggressivePresolveDepth) >= 0 );
   }
   int getAggresivePresolvingDepth()
   {
      return paraParams->getIntParamValue(AggressivePresolveDepth);
   }
   int getAggresivePresolvingStopDepth()
   {
      return paraParams->getIntParamValue(AggressivePresolveStopDepth);
   }
   int getSubMipDepth()
   {
      return currentNode->getDepth();
   }
   void setSendBackAllNodes()
   {
      nCollectOnce = -1;          // collect all
      collectingManyNodes = true;
   }
   bool isCollectingAllNodes()
   {
      return( collectingManyNodes && (nCollectOnce < 0) );
   }
   int getBigDualGapSubtreeHandlingStrategy()
   {
      return paraParams->getIntParamValue(BigDualGapSubtreeHandling);
   }

   bool isIterativeBreakDownApplied()
   {
      return paraParams->getBoolParamValue(IterativeBreakDown);
   }

   void setII(
         double sum,    /**< sum of integer infeasibility */
         int count      /**< number of integer infeasibility */
         )
   {
      if( minIisum > sum ) minIisum = sum;
      if( maxIisum < sum ) maxIisum = sum;
      if( minNii > count ) minNii = count;
      if( maxNii < count ) maxNii = count;
   }

   void setRootNodeSimplexIter(
         int iter
         )
   {
      nSimplexIterRoot = iter;
   }

   bool waitToken(int rank){
      bool result;
      double startTimeToWaitToken = paraTimer->getElapsedTime();
      result = paraComm->waitToken(rank);
      idleTimeToWaitToken += (paraTimer->getElapsedTime() - startTimeToWaitToken);
      return result;
   }
   void passToken(
         int rank
         )
   {
      paraComm->passToken(rank);
   }

   int getCurrentSolivingNodeMergingStatus(){ return currentNode->getMergingStatus(); }
   double getCurrentSolvingNodeInitialDualBound(){ return currentNode->getInitialDualBoundValue(); }

   ParaDeterministicTimer *getDeterministicTimer(){ return paraDetTimer; }

   double getOffsetTimeToWaitToken(){ return offsetTimeToWaitToken; }

   void updatePendingSolution()
   {
      if( !pendingSolution )
      {
         updateGlobalBestIncumbentValue(pendingIncumbentValue);
         pendingIncumbentValue = DBL_MAX;
         return;
      }
      if( updateGlobalBestIncumbentSolution(pendingSolution) )
      {
         tryNewSolution(pendingSolution);
      }
      else
      {
         delete pendingSolution;
      }
      pendingIncumbentValue = DBL_MAX;
      pendingSolution = 0;
   }

   virtual bool wasTerminatedNormally() = 0;

   /** this method is always useful for debugging, so we should implement this method */
   virtual void writeCurrentNodeProblem(const std::string& filename) = 0;
   virtual void tryNewSolution(ParaSolution *sol) = 0;
   virtual void setLightWeightRootNodeProcess() = 0;
   virtual void setOriginalRootNodeProcess() = 0;
   virtual void writeSubproblem() = 0;
   virtual long long getSimplexIter() = 0;
   virtual int getNRestarts() { return 0; } // Derived class for SCIP should override this function

   /* for deterministic execution */
   void setPreviousCommTime(double detTime){ previousCommTime = detTime; }
   double getPreviousCommTime(){ return previousCommTime; }

};

}

#endif // __PARA_SOLVER_H__
