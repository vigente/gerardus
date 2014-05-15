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

/**@file    paraComm.h
 * @brief   Base class of communicator for UG Framework
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_COMM_H__
#define __PARA_COMM_H__
#include <sstream>
#include "paraTagDef.h"
#include "paraParamSet.h"

namespace UG
{

#define PARA_COMM_CALL( paracommcall ) \
  { \
          int status = paracommcall; \
          if( status )  \
          { \
             std::ostringstream s_; \
             s_ << "[PARA_COMM_CALL ERROR: " << __FILE__  << "] func = " \
             << __func__ << ", line = " << __LINE__ << ": " \
             << "error_code = " << status << std::endl; \
             throw std::logic_error( s_.str() ); \
          }\
  }


static const int TYPE_FIRST             =               0;
static const int ParaCHAR               = TYPE_FIRST +  0;
static const int ParaSHORT              = TYPE_FIRST +  1;
static const int ParaINT                = TYPE_FIRST +  2;
static const int ParaLONG               = TYPE_FIRST +  3;
static const int ParaLONG_LONG          = TYPE_FIRST +  4;
static const int ParaSIGNED_CHAR        = TYPE_FIRST +  5;
static const int ParaUNSIGNED_CHAR      = TYPE_FIRST +  6;
static const int ParaUNSIGNED_SHORT     = TYPE_FIRST +  7;
static const int ParaUNSIGNED           = TYPE_FIRST +  8;
static const int ParaUNSIGNED_LONG      = TYPE_FIRST +  9;
static const int ParaUNSIGNED_LONG_LONG = TYPE_FIRST + 10;
static const int ParaFLOAT              = TYPE_FIRST + 11;
static const int ParaDOUBLE             = TYPE_FIRST + 12;
static const int ParaLONG_DOUBLE        = TYPE_FIRST + 13;
static const int ParaBOOL               = TYPE_FIRST + 14;
static const int ParaBYTE               = TYPE_FIRST + 15;
static const int TYPE_LAST              = TYPE_FIRST + 15;
static const int TYPE_LIST_SIZE         = TYPE_LAST - TYPE_FIRST + 1;

static const int NumMaxWorkers          = 20000;

class ParaCalculationState;
class ParaParamSet;
class ParaSolverState;
class ParaSolverTerminationState;
class ParaInstance;
class ParaDiffSubproblem;
class ParaSolution;
class ParaInitialStat;
class ParaRacingRampUpParamSet;
class ParaNode;
class ParaTimer;
class NodeId;

class ParaComm
{
public:
   ParaComm(){}
   virtual void init( int argc, char **argv ) = 0;
   virtual ~ParaComm(){}
   virtual int getRank() = 0;
   virtual int getSize() = 0;
   virtual void lcInit(ParaParamSet *paraParamSet) = 0;
   virtual void solverInit(ParaParamSet *paraParamSet) = 0;
   virtual void abort() = 0;
   virtual bool waitTerminatedMessage() = 0;
   virtual bool waitToken(int rank){ return true; }
   virtual void passToken(int rank){}
   virtual bool passTermToken(int rank){ return true; }
   virtual void setToken(int rank, int *token){}
   virtual void lockApp(){}
   virtual void unlockApp(){}
   /*******************************************************************************
   * transfer object factory
   *******************************************************************************/
   virtual ParaCalculationState *createParaCalculationState() = 0;
   virtual ParaCalculationState *createParaCalculationState(
               double compTime,                   /**< computation time of this ParaNode */
               double rootTime,                   /**< computation time of the root node */
               int    nSolved,                    /**< the number of nodes solved   */
               int    nSent,                      /**< the number of ParaNodes sent */
               int    nImprovedIncumbent,         /**< the number of improved solution generated in this ParaSolver */
               int    terminationState,           /**< indicate whether if this computation is terminationState or not. 0: no, 1: terminationState */
               int    nSolvedWithNoPreprocesses,  /**< number of nodes solved when it is solved with no preprocesses */
               int    nSimplexIterRoot,           /**< number of simplex iteration at root node */
               double averageSimplexIter,         /**< average number of simplex iteration except root node */
               int    nRestarts,                  /**< number of restarts */
               double minIisum,                   /**< minimum sum of integer infeasibility */
               double maxIisum,                   /**< maximum sum of integer infeasibility */
               int    minNii,                     /**< minimum number of integer infeasibility */
               int    maxNii                      /**< maximum number of integer infeasibility */
           ) = 0;
   virtual ParaInitialStat* createParaInitialStat() = 0;
   virtual ParaRacingRampUpParamSet* createParaRacingRampUpParamSet() = 0;
   virtual ParaNode *createParaNode() = 0;
   virtual ParaNode *createParaNode(
               NodeId inNodeId,
               NodeId inGeneratorNodeId,
               int inDepth,
               double inDualBoundValue,
               double inOriginalDualBoundValue,
               double inEstimatedValue,
               ParaDiffSubproblem *inDiffSubproblem
           ) = 0;
   virtual ParaParamSet *createParaParamSet() = 0;
   virtual ParaSolverState *createParaSolverState() = 0;
   virtual ParaSolverState *createParaSolverState(
               int racingStage,
               unsigned int notificationId,
               int lcId,
               int globalSubtreeId,
               long long nodesSolved,
               int nodesLeft,
               double bestDualBoundValue,
               double globalBestPrimalBoundValue,
               double detTime
           ) = 0;
   virtual ParaSolverTerminationState *createParaSolverTerminationState() = 0;
   virtual ParaSolverTerminationState *createParaSolverTerminationState(
               int    interrupted,                /**< indicate that this solver is interrupted or not. 0: not interrupted, 1: interrputed
                                                                                                       2: checkpoint, 3: racing-ramp up */
               int    rank,                       /**< rank of this solver */
               int    totalNSolved,               /**< accumulated number of nodes solved in this ParaSolver */
               int    minNSolved,                 /**< minimum number of subtree nodes rooted from ParaNode */
               int    maxNSolved,                 /**< maximum number of subtree nodes rooted from ParaNode */
               int    totalNSent,                 /**< accumulated number of nodes sent from this ParaSolver */
               int    totalNImprovedIncumbent,    /**< accumulated number of improvements of incumbent value in this ParaSolver */
               int    nParaNodesReceived,         /**< number of ParaNodes received in this ParaSolver */
               int    nParaNodesSolved,           /**< number of ParaNodes solved ( received ) in this ParaSolver */
               int    nParaNodesSolvedAtRoot,     /**< number of ParaNodes solved at root node before sending  */
               int    nParaNodesSolvedAtPreCheck, /**< number of ParaNodes solved at pre-checking of root node solvability */
               double runningTime,                /**< this solver running time */
               double idleTimeToFirstParaNode,    /**< idle time to start solving the first ParaNode */
               double idleTimeBetweenParaNodes,   /**< idle time between ParaNodes processing */
               double iddleTimeAfterLastParaNode, /**< idle time after the last ParaNode was solved */
               double idleTimeToWaitNotificationId,   /**< idle time to wait notification Id messages  */
               double idleTimeToWaitAckCompletion,    /**< idle time to wait ack completion message  */
               double idleTimeToWaitToken,        /**< idle time to wait token  */
               double totalRootNodeTime,          /**< total time consumed by root node processes */
               double minRootNodeTime,            /**< minimum time consumed by root node processes */
               double maxRootNodeTime,            /**< maximum time consumed by root node processes */
               double detTime                     /**< deterministic time, -1: should be non-deterministic */
           ) = 0;
   virtual ParaTimer *createParaTimer() = 0;
   virtual ParaDiffSubproblem *createParaDiffSubproblem() = 0;
   virtual ParaInstance *createParaInstance() = 0;
   virtual ParaSolution *createParaSolution() = 0;
   /*******************************************************************************
   *  Some action need to be taken for fault tolerant, when the functions return. *
   *  So, they rerun status value                                                 *
   *******************************************************************************/
   virtual int bcast( void* buffer, int count, const int datatypeId, int root ) = 0;
   virtual int send( void* bufer, int count, const int datatypeId, int dest, const int tag ) = 0;
   virtual int receive( void* bufer, int count, const int datatypeId, int source, const int tag ) = 0;
   virtual int waitSpecTagFromSpecSource(const int source, const int datatypeId, const int tag, int *receivedTag) = 0;
   /*************************************************************************
   /  No need to take action for fault tolerant, when the functions return. *
   *  So, they do not rerun status value                                    *
   *************************************************************************/
   virtual void probe(int *source, int *tag) = 0;
   virtual bool iProbe(int *source, int *tag) = 0;  /** retrun value: true => there is message, false => no message */
};

}

#endif // __PARA_COMM_H__
