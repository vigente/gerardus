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

/**@file    paraCommMpiWorld.cpp
 * @brief   ParaComm extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraCommMpiWorld.h"

using namespace UG;

MPI::Datatype
ParaCommMpiWorld::datatypes[TYPE_LIST_SIZE];

const char *
ParaCommMpiWorld::tagStringTable[] = {
  TAG_STR(TagNode),
  TAG_STR(TagNodeReceived),
  TAG_STR(TagDiffSubproblem),
  TAG_STR(TagRampUp),
  TAG_STR(TagRetryRampUp),
  TAG_STR(TagSolution),
  TAG_STR(TagIncumbentValue),
  TAG_STR(TagGlobalBestDualBoundValue),
  TAG_STR(TagSolverState),
  TAG_STR(TagCompletionOfCalculation),
  TAG_STR(TagAnotherNodeRequest),
  TAG_STR(TagNoNodes),
  TAG_STR(TagInCollectingMode),
  TAG_STR(TagCollectAllNodes),
  TAG_STR(TagOutCollectingMode),
  TAG_STR(TagLCBestBoundValue),
  TAG_STR(TagNotificationId),
  TAG_STR(TagTerminateRequest),
  TAG_STR(TagInterruptRequest),
  TAG_STR(TagTerminated),
  TAG_STR(TagRacingRampUpParamSets),
  TAG_STR(TagWinner),
  TAG_STR(TagLightWeightRootNodeProcess),
  TAG_STR(TagBreaking),
  TAG_STR(TagHardTimeLimit),
  TAG_STR(TagInitialStat),
  TAG_STR(TagAckCompletion),
  TAG_STR(TagToken),
  TAG_STR(TagSolverDiffParamSet1),
  TAG_STR(TagSolverDiffParamSet2),
  TAG_STR(TagDiffSubproblem1),
  TAG_STR(TagDiffSubproblem2),
  TAG_STR(TagSolution1),
  TAG_STR(TagInitialStat1)
};

void
ParaCommMpiWorld::lcInit(
      ParaParamSet *paraParamSet
      )
{
   tagTraceFlag = paraParamSet->getBoolParamValue(TagTrace);
   if( tagTraceFlag )
   {
      if( paraParamSet->isStringParamDefaultValue(TagTraceFileName) )
      {
         tos = &std::cout;
      }
      else
      {
         std::ostringstream s;
         s << paraParamSet->getStringParamValue(TagTraceFileName) << myRank;
         ofs.open(s.str().c_str());
         tos = &ofs;
      }
   }
   if( paraParamSet->getBoolParamValue(Deterministic) )
   {
      token[0] = 0;
      token[1] = -1;
   }
}

void
ParaCommMpiWorld::solverInit(
      ParaParamSet *paraParamSet
      )
{
   tagTraceFlag = paraParamSet->getBoolParamValue(TagTrace);
   if( tagTraceFlag )
   {
      if( paraParamSet->isStringParamDefaultValue(TagTraceFileName) )
      {
         tos = &std::cout;
      }
      else
      {
         std::ostringstream s;
         s << paraParamSet->getStringParamValue(TagTraceFileName) << myRank;
         ofs.open(s.str().c_str());
         tos = &ofs;
      }
   }
}

void
ParaCommMpiWorld::abort(
      )
{
   MPI::COMM_WORLD.Abort(0);
}

bool
ParaCommMpiWorld::waitToken(
      int tempRank
      )
{
   pthread_mutex_lock(&tokenAccessLock);
   if( token[0] == myRank )
   {
      pthread_mutex_unlock(&tokenAccessLock);
      return true;
   }
   else
   {
      int previousRank = myRank - 1;
      if( previousRank == 0 )
      {
         if( token[0] != -1 )
         {
            previousRank = comSize - 1;
         }
      }
      int receivedTag;
      MPI::Status mpiStatus;
      MPI::COMM_WORLD.Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, mpiStatus);
      receivedTag = mpiStatus.Get_tag();
      TAG_TRACE (Probe, From, mpiStatus.Get_source(), receivedTag);
      if( receivedTag == TagToken )
      {
         receive(token, 2, ParaINT, 0, TagToken);
         assert(token[0] == myRank);
         pthread_mutex_unlock(&tokenAccessLock);
         return true;
      }
      else
      {
         pthread_mutex_unlock(&tokenAccessLock);
         return false;
      }
   }
}

void
ParaCommMpiWorld::passToken(
      int tempRank
      )
{
   pthread_mutex_lock(&tokenAccessLock);
   assert( token[0] == myRank );
   token[0] = ( token[0]  % (comSize - 1) ) + 1;
   token[1] = -1;
   send(token, 2, ParaINT, 0, TagToken);
   pthread_mutex_unlock(&tokenAccessLock);
}

bool
ParaCommMpiWorld::passTermToken(
      int tempRank
      )
{
   pthread_mutex_lock(&tokenAccessLock);
   if( myRank == token[0] )
   {
      if( token[1] == token[0] ) token[1] = -2;
      else if( token[1] == -1 ) token[1] = token[0];
      token[0] = ( token[0]  % (comSize - 1) ) + 1;
   }
   else
   {
      THROW_LOGICAL_ERROR4("Invalid token update. Rank = ", getRank(), ", token = ", token[0] );
   }
   send(token, 2, ParaINT, 0, TagToken);
   if( token[1] == -2 )
   {
      pthread_mutex_unlock(&tokenAccessLock);
      return true;
   }
   else
   {
      pthread_mutex_unlock(&tokenAccessLock);
      return false;
   }
}

/** MPI call wrappers */
void
ParaCommMpiWorld::init( int argc, char **argv )
{
   MPI::Init( argc, argv );
   startTime = MPI::Wtime();
   MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
   MPI_CALL_WITH_RET_VAL(
      myRank, MPI::COMM_WORLD.Get_rank()
   );
   MPI_CALL_WITH_RET_VAL(
      comSize, MPI::COMM_WORLD.Get_size()
   );
   char *pprocName = procName;
   MPI_CALL(
      MPI::Get_processor_name(pprocName,namelen)
   );

   /** if you add tag, you should add tagStringTale too */
   assert( sizeof(tagStringTable)/sizeof(char*) == N_MPI_TAGS );

   /** Data Types */
   datatypes[ParaCHAR] = MPI::CHAR;
   datatypes[ParaSHORT] = MPI::SHORT;
   datatypes[ParaINT] = MPI::INT;
   datatypes[ParaLONG] = MPI::LONG;
   datatypes[ParaUNSIGNED_CHAR] = MPI::UNSIGNED_CHAR;
   datatypes[ParaUNSIGNED_SHORT] = MPI::UNSIGNED_SHORT;
   datatypes[ParaUNSIGNED] = MPI::UNSIGNED;
   datatypes[ParaUNSIGNED_LONG] = MPI::UNSIGNED_LONG;
   datatypes[ParaFLOAT] = MPI::FLOAT;
   datatypes[ParaDOUBLE] = MPI::DOUBLE;
   datatypes[ParaLONG_DOUBLE] = MPI::LONG_DOUBLE;
   datatypes[ParaBYTE] = MPI::BYTE;

#ifdef _ALIBABA
   datatypes[ParaSIGNED_CHAR] = MPI::CHAR;
   datatypes[ParaLONG_LONG] = MPI::LONG;
   datatypes[ParaUNSIGNED_LONG_LONG] = MPI::UNSIGNED_LONG;
   datatypes[ParaBOOL] = MPI::INT;
#else
   datatypes[ParaSIGNED_CHAR] = MPI::SIGNED_CHAR;
   datatypes[ParaLONG_LONG] = MPI::LONG_LONG;
   datatypes[ParaUNSIGNED_LONG_LONG] = MPI::UNSIGNED_LONG_LONG;
   datatypes[ParaBOOL] = MPI::BOOL;
#endif

}

ParaCommMpiWorld::~ParaCommMpiWorld()
{
   MPI::Finalize();
}

ParaCalculationState *
ParaCommMpiWorld::createParaCalculationState(
      )
{
   return new ParaCalculationStateMpi();
}

ParaCalculationState *
ParaCommMpiWorld::createParaCalculationState(
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
           )
{
   return new ParaCalculationStateMpi(
                  compTime,
                  rootTime,
                  nSolved,
                  nSent,
                  nImprovedIncumbent,
                  terminationState,
                  nSolvedWithNoPreprocesses,
                  nSimplexIterRoot,
                  averageSimplexIter,
                  nRestarts,
                  minIisum,
                  maxIisum,
                  minNii,
                  maxNii
              );
}

ParaNode *
ParaCommMpiWorld::createParaNode(
      )
{
   return new ParaNodeMpi();
}

ParaNode *
ParaCommMpiWorld::createParaNode(
               NodeId inNodeId,
               NodeId inGeneratorNodeId,
               int inDepth,
               double inDualBoundValue,
               double inOriginalDualBoundValue,
               double inEstimatedValue,
               ParaDiffSubproblem *inDiffSubproblem
            )
{
    return new ParaNodeMpi(
                  inNodeId,
                  inGeneratorNodeId,
                  inDepth,
                  inDualBoundValue,
                  inOriginalDualBoundValue,
                  inEstimatedValue,
                  inDiffSubproblem
              );
}

ParaParamSet *
ParaCommMpiWorld::createParaParamSet(
      )
{
   return new ParaParamSetMpi();
}

ParaSolverState *
ParaCommMpiWorld::createParaSolverState(
      )
{
   return new ParaSolverStateMpi();
}

ParaSolverState *
ParaCommMpiWorld::createParaSolverState(
               int racingStage,
               unsigned int notificationId,
               int lcId,
               int globalSubtreeId,
               long long nodesSolved,
               int nodesLeft,
               double bestDualBoundValue,
               double globalBestPrimalBoundValue,
               double detTime
           )
{
   return new ParaSolverStateMpi(
                  racingStage,
                  notificationId,
                  lcId,
                  globalSubtreeId,
                  nodesSolved,
                  nodesLeft,
                  bestDualBoundValue,
                  globalBestPrimalBoundValue,
                  detTime
              );
}


ParaSolverTerminationState *
ParaCommMpiWorld::createParaSolverTerminationState(
      )
{
   return new ParaSolverTerminationStateMpi();
}

ParaSolverTerminationState *
ParaCommMpiWorld::createParaSolverTerminationState(
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
           )
{
   return new ParaSolverTerminationStateMpi(
                 interrupted,
                 rank,
                 totalNSolved,
                 minNSolved,
                 maxNSolved,
                 totalNSent,
                 totalNImprovedIncumbent,
                 nParaNodesReceived,
                 nParaNodesSolved,
                 nParaNodesSolvedAtRoot,
                 nParaNodesSolvedAtPreCheck,
                 runningTime,
                 idleTimeToFirstParaNode,
                 idleTimeBetweenParaNodes,
                 iddleTimeAfterLastParaNode,
                 idleTimeToWaitNotificationId,
                 idleTimeToWaitAckCompletion,
                 idleTimeToWaitToken,
                 totalRootNodeTime,
                 minRootNodeTime,
                 maxRootNodeTime,
                 detTime
              );
}

int
ParaCommMpiWorld::bcast(
   void* buffer,
   int count,
   const int datatypeId,
   int root
   )
{
   MPI_CALL(
      MPI::COMM_WORLD.Bcast( buffer, count, datatypes[datatypeId], root )
   );
   return 0;
}

int
ParaCommMpiWorld::send(
   void* buffer,
   int count,
   const int datatypeId,
   int dest,
   const int tag
   )
{
   MPI::COMM_WORLD.Send( buffer, count, datatypes[datatypeId], dest, tag );
   TAG_TRACE (Send, To, dest, tag);
   return 0;
}

int
ParaCommMpiWorld::receive(
   void* buffer,
   int count,
   const int datatypeId,
   int source,
   const int tag
   )
{
   MPI_CALL (
      MPI::COMM_WORLD.Recv( buffer, count, datatypes[datatypeId], source, tag )
   );
   TAG_TRACE (Recv, From, source, tag);
   return 0;
}

int
ParaCommMpiWorld::waitSpecTagFromSpecSource(
      const int source,
      const int datatypeId,
      const int tag,
      int *receivedTag
      )
{
   MPI::Status mpiStatus;
   MPI::COMM_WORLD.Probe(source, MPI::ANY_TAG, mpiStatus);
   (*receivedTag) = mpiStatus.Get_tag();
   TAG_TRACE (Probe, From, source, (*receivedTag));
   if( tag == (*receivedTag) )
   {
      return 0;
   }
   else
   {
      return 1;
   }
}

void
ParaCommMpiWorld::probe(
   int* source,
   int* tag
   )
{
   MPI::Status mpiStatus;
   MPI::COMM_WORLD.Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, mpiStatus);
   *source = mpiStatus.Get_source();
   *tag = mpiStatus.Get_tag();
   TAG_TRACE (Probe, From, *source, *tag);
}

bool
ParaCommMpiWorld::iProbe(
   int* source,
   int* tag
   )
{
   bool flag;
   MPI::Status mpiStatus;
   flag = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, MPI::ANY_TAG, mpiStatus);
   if( flag )
   {
      *source = mpiStatus.Get_source();
      *tag = mpiStatus.Get_tag();
      TAG_TRACE (Iprobe, From, *source, *tag);
   }
   return flag;
}

int
ParaCommMpiWorld::bcast(
   void* buffer,
   int count,
   MPI::Datatype datatype,
   int root
   )
{
   MPI_CALL(
      MPI::COMM_WORLD.Bcast( buffer, count, datatype, root )
   );
   return 0;
}

int
ParaCommMpiWorld::send(
   void* buffer,
   int count,
   MPI::Datatype datatype,
   int dest,
   const int tag
   )
{
   MPI::COMM_WORLD.Send( buffer, count, datatype, dest, tag );
   TAG_TRACE (Send, To, dest, tag);
   return 0;
}

int
ParaCommMpiWorld::receive(
   void* buffer,
   int count,
   MPI::Datatype datatype,
   int source,
   const int tag
   )
{
   MPI::COMM_WORLD.Recv( buffer, count, datatype, source, tag );
   TAG_TRACE (Recv, From, source, tag);
   return 0;
}
