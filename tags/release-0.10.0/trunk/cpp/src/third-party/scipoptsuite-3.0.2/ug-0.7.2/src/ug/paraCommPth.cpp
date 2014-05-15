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

/**@file    paraCommPth.cpp
 * @brief   ParaComm extension for Pthreads communication
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cstring>
#include "paraCommPth.h"
#include "paraNodePth.h"
#include "paraSolution.h"
#include "paraSolverTerminationStatePth.h"
#include "paraCalculationStatePth.h"
#include "paraSolverStatePth.h"
#include "paraRacingRampUpParamSet.h"
#include "paraInitialStat.h"

using namespace UG;

pthread_mutex_t hashtableLock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t hashtableUpdated = PTHREAD_COND_INITIALIZER;

HashTableElement *
ParaCommPth::hashtable[HashTableSize];

const char *
ParaCommPth::tagStringTable[] = {
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
  TAG_STR(TagParaInstance),
  TAG_STR(TagSolverDiffParamSet)
};

void
ParaCommPth::init( int argc, char **argv )
{
   // don't have to take any lock, because only LoadCoordinator call this function

   timer.start();
   comSize = 0;

   for( int i = 1; i < argc; i++ )
   {
      if( strcmp(argv[i], "-sth") == 0 )
      {
         i++;
         if( i < argc )
            comSize = atoi(argv[i]);   // if -sth 0, then it is considered as use the number of cores system has
         else
         {
            std::cerr << "missing the number of solver threads after parameter '-sth'" << std::endl;
            exit(1);
         }
      }
   }

   if( comSize > 0 )
   {
      comSize++;
   }
   else
   {
      comSize = sysconf(_SC_NPROCESSORS_CONF) + 1;
   }

   tokenAccessLock = new pthread_mutex_t[comSize];
   token = new int*[comSize];
   for( int i = 0; i < comSize; i++ )
   {
      pthread_mutex_init(&tokenAccessLock[i], NULL);
      token[i] = new int[2];
      token[i][0] = 0;
      token[i][1] = -1;
   }

   /** if you add tag, you should add tagStringTale too */
   assert( sizeof(tagStringTable)/sizeof(char*) == N_PTH_TAGS );

   /** initialize hashtable */
   for(int i = 0; i < HashTableSize; i++ )
   {
      hashtable[i] = 0;
   }

   messageQueueTable = new MessageQueueTableElement *[comSize + 1];  // +1 for TimeLimitMonitor
   for( int i = 0; i < ( comSize + 1 ); i++ )
   {
      messageQueueTable[i] = new MessageQueueTableElement;
   }

}

void
ParaCommPth::lcInit(
      ParaParamSet *paraParamSet
      )
{
   // don't have to take any lock, because only LoadCoordinator call this function
   pthread_t tid = pthread_self();
 
   pthread_mutex_lock(&hashtableLock);
   hashtable[HashEntry(tid)] = new HashTableElement(tid, 0, paraParamSet);
   pthread_mutex_unlock(&hashtableLock);

   tagTraceFlag = paraParamSet->getBoolParamValue(TagTrace);
}

void
ParaCommPth::solverInitPth(
      pthread_t tid,
      int rank,
      ParaParamSet *paraParamSet
      )
{
   // don't have to take any lock, because only LoadCoordinator call this function
   // CHANGED in multi-threaded solver case
   pthread_mutex_lock(&hashtableLock);
   int index = HashEntry(tid);
   if( hashtable[index] == 0 )
   {
      hashtable[index] = new HashTableElement(tid, rank, paraParamSet);
   }
   else
   {
      HashTableElement *elem = hashtable[index];
      while( elem->getNext() != 0 )
      {
         elem = elem->getNext();
      }
      elem->link(new HashTableElement(tid, rank, paraParamSet));
   }
   pthread_mutex_unlock(&hashtableLock);
}

void
ParaCommPth::waitUntilRegistered(
      )
{
   pthread_t tid = pthread_self();
   bool registered = false;
   pthread_mutex_lock(&hashtableLock);
   while( !registered )
   {
      HashTableElement *elem = hashtable[HashEntry(tid)];
      while( elem && !pthread_equal( elem->getTid(), tid ) )
      {
         elem = elem->getNext();
      }
      if( elem )
      {
         assert( pthread_equal( elem->getTid(), tid ) );
         registered = true;
         break;
      }
      pthread_cond_wait(&hashtableUpdated, &hashtableLock);
   }
   pthread_mutex_unlock(&hashtableLock);
}

void
ParaCommPth::registedAllSolvers(
      )
{
   pthread_mutex_lock(&hashtableLock);
   pthread_cond_broadcast(&hashtableUpdated);
   pthread_mutex_unlock(&hashtableLock);

}

bool
ParaCommPth::waitToken(
      int rank
      )
{

   // int rank = getRank();   // multi-thread solver may change rank here
   pthread_mutex_lock(&tokenAccessLock[rank]);
   if( token[rank][0] == rank )
   {
      pthread_mutex_unlock(&tokenAccessLock[rank]);
      return true;
   }
   else
   {
      int receivedTag;
      int source;
      probe(&source, &receivedTag);
      TAG_TRACE (Probe, From, source, receivedTag);
      if( source == 0 && receivedTag == TagToken )
      {
         receive(token[rank], 2, ParaINT, 0, TagToken);
         assert( token[rank][0] == rank );
         pthread_mutex_unlock(&tokenAccessLock[rank]);
         return true;
      }
      else
      {
         pthread_mutex_unlock(&tokenAccessLock[rank]);
         return false;
      }
   }
}

void
ParaCommPth::passToken(
      int rank
      )
{
   // int rank = getRank();   // multi-thread solver may change rank here
   pthread_mutex_lock(&tokenAccessLock[rank]);
   assert( token[rank][0] == rank && rank != 0 );
   token[rank][0] = ( token[rank][0]  % (comSize - 1) ) + 1;
   token[rank][1] = -1;
   send(token[rank], 2, ParaINT, 0, TagToken);
   pthread_mutex_unlock(&tokenAccessLock[rank]);
}

bool
ParaCommPth::passTermToken(
      int rank
      )
{
   // int rank = getRank();   // multi-thread solver may change rank here
   pthread_mutex_lock(&tokenAccessLock[rank]);
   if( rank == token[rank][0] )
   {
      if( token[rank][1] == token[rank][0] ) token[rank][1] = -2;
      else if( token[rank][1] == -1 ) token[rank][1] = token[rank][0];
      token[rank][0] = ( token[rank][0]  % (comSize - 1) ) + 1;
   }
   else
   {
      THROW_LOGICAL_ERROR4("Invalid token update. Rank = ", getRank(), ", token = ", token[0] );
   }
   send(token[rank], 2, ParaINT, 0, TagToken);
   if( token[rank][1] == -2 )
   {
      pthread_mutex_unlock(&tokenAccessLock[rank]);
      return true;
   }
   else
   {
      pthread_mutex_unlock(&tokenAccessLock[rank]);
      return false;
   }
}

void
ParaCommPth::setToken(
      int rank,
      int *inToken
      )
{
   // int rank = getRank();
   pthread_mutex_lock(&tokenAccessLock[rank]);
   assert( rank == 0 || ( rank != 0 && inToken[0] == rank ) );
   token[rank][0] = inToken[0];
   token[rank][1] = inToken[1];
   pthread_mutex_unlock(&tokenAccessLock[rank]);
}

ParaCommPth::~ParaCommPth()
{

   for(int i = 0; i < HashTableSize; i++ )
   {
      if( hashtable[i] )
      {
         while( hashtable[i] )
         {
            HashTableElement  *next = hashtable[i]->getNext();
            delete hashtable[i];
            hashtable[i] = next;
         }
      }
   }

   for( int i = 0; i < comSize; i++ )
   {
      delete [] token[i];
   }
   delete [] token;
   delete [] tokenAccessLock;

   for(int i = 0; i < (comSize + 1); i++)
   {
      MessageQueueElement *elem = messageQueueTable[i]->extarctElement();
      while( elem )
      {
         if( elem->getData() )
         {
            if( elem->getDataTypeId() < USER_TYPE_FIRST )
            {
               freeMem(elem->getData(), elem->getCount(), elem->getDataTypeId() );
            }
            else
            {
               switch( elem->getDataTypeId())
               {
               case ParaInstanceType:
               {
                  delete reinterpret_cast<ParaInstance *>(elem->getData());
                  break;
               }
               case ParaSolutionType:
               {
                  delete reinterpret_cast<ParaSolution *>(elem->getData());
                  break;
               }
               case ParaParamSetType:
               {
                  delete reinterpret_cast<ParaParamSet *>(elem->getData());
                  break;
               }
               case ParaNodeType:
               {
                  delete reinterpret_cast<ParaNode *>(elem->getData());
                  break;
               }
               case ParaSolverStateType:
               {
                  delete reinterpret_cast<ParaSolverState *>(elem->getData());
                  break;
               }
               case ParaCalculationStateType:
               {
                  delete reinterpret_cast<ParaCalculationState *>(elem->getData());
                  break;
               }
               case ParaSolverTerminationStateType:
               {
                  delete reinterpret_cast<ParaSolverTerminationState *>(elem->getData());
                  break;
               }
               case ParaRacingRampUpParamType:
               {
                  delete reinterpret_cast<ParaRacingRampUpParamSet *>(elem->getData());
                  break;
               }
               case ParaSolverDiffParamType:
               {
                  // Not supported now.
                  break;
               }
               case ParaInitialStatType:
               {
                  delete reinterpret_cast<ParaInitialStat *>(elem->getData());
                  break;
               }
               default:
               {
                  THROW_LOGICAL_ERROR2("Requested type is not implemented. Type = ", elem->getDataTypeId() );
               }
               }
            }
         }
         delete elem;
         elem = messageQueueTable[i]->extarctElement();
      }
      delete messageQueueTable[i];
   }
   delete [] messageQueueTable;
}

unsigned int
ParaCommPth::hashCode(
      pthread_t tid
      )
{
   union {
      pthread_t tid;
      unsigned char      cTid[sizeof(pthread_t)];
   } reinterpret;

   reinterpret.tid = tid;
   unsigned int h = 0;
   for (unsigned int i = 0; i < sizeof(pthread_t); i++) {
       h = 31*h + reinterpret.cTid[i];
   }
   return h;
}

int
ParaCommPth::getRank(
      )
{
   pthread_mutex_lock(&hashtableLock);          // it is necessary in case multi-threaded solver
   pthread_t tid = pthread_self();
   // std::cout << "gerRank tid = " << tid << std::endl;
   HashTableElement *elem = hashtable[HashEntry(tid)];
   // assert(elem);
   while( elem && !pthread_equal( elem->getTid(), tid ) )
   {
      elem = elem->getNext();
      // assert(elem);
   }
   if( elem )
   {
      pthread_mutex_unlock(&hashtableLock);
      return elem->getRank();
   }
   else
   {
      pthread_mutex_unlock(&hashtableLock);
      return -1; // No ug threads
   }
}

int
ParaCommPth::getRank(
      pthread_t tid
      )
{
   pthread_mutex_lock(&hashtableLock);          // it is necessary in case multi-threaded solver
   // std::cout << "gerRank tid = " << tid << std::endl;
   HashTableElement *elem = hashtable[HashEntry(tid)];
   // assert(elem);
   while( elem && !pthread_equal( elem->getTid(), tid ) )
   {
      elem = elem->getNext();
      // assert(elem);
   }
   if( elem )
   {
      pthread_mutex_unlock(&hashtableLock);
      return elem->getRank();
   }
   else
   {
      pthread_mutex_unlock(&hashtableLock);
      return -1; // No ug threads
   }
}

std::ostream *
ParaCommPth::getOstream(
      )
{
   pthread_t tid = pthread_self();
   pthread_mutex_lock(&hashtableLock);
   HashTableElement *elem = hashtable[HashEntry(tid)];
   assert(elem);
   while( !pthread_equal( elem->getTid(), tid ) )
   {
      elem = elem->getNext();
      assert(elem);
   }
   pthread_mutex_unlock(&hashtableLock);
   return elem->getOstream();
}

ParaCalculationState *
ParaCommPth::createParaCalculationState(
      )
{
   return new ParaCalculationStatePth();
}

ParaCalculationState *
ParaCommPth::createParaCalculationState(
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
   return new ParaCalculationStatePth(
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
ParaCommPth::createParaNode(
      )
{
   return new ParaNodePth();
}

ParaNode *
ParaCommPth::createParaNode(
               NodeId inNodeId,
               NodeId inGeneratorNodeId,
               int inDepth,
               double inDualBoundValue,
               double inOriginalDualBoundValue,
               double inEstimatedValue,
               ParaDiffSubproblem *inDiffSubproblem
            )
{
    return new ParaNodePth(
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
ParaCommPth::createParaParamSet(
      )
{
   return new ParaParamSetPth();
}

ParaSolverState *
ParaCommPth::createParaSolverState(
      )
{
   return new ParaSolverStatePth();
}

ParaSolverState *
ParaCommPth::createParaSolverState(
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
   return new ParaSolverStatePth(
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
ParaCommPth::createParaSolverTerminationState(
      )
{
   return new ParaSolverTerminationStatePth();
}

ParaSolverTerminationState *
ParaCommPth::createParaSolverTerminationState(
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
   return new ParaSolverTerminationStatePth(
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

void *
ParaCommPth::allocateMemAndCopy(
      void* buffer,
      int count,
      const int datatypeId
      )
{
   void *newBuf = 0;
   if( count == 0 ) return newBuf;

   switch(datatypeId)
   {
   case ParaCHAR :
   {
      newBuf = new char[count];
      memcpy(newBuf, buffer, sizeof(char)*count);
      break;
   }
   case ParaSHORT :
   {
      newBuf = new short[count];
      memcpy(newBuf, buffer, sizeof(short)*count);
      break;
   }
   case ParaINT :
   {
      newBuf = new int[count];
      memcpy(newBuf, buffer, sizeof(int)*count);
      break;
   }
   case ParaLONG :
   {
      newBuf = new long[count];
      memcpy(newBuf, buffer, sizeof(long)*count);
      break;
   }
   case ParaUNSIGNED_CHAR :
   {
      newBuf = new unsigned char[count];
      memcpy(newBuf, buffer, sizeof(unsigned char)*count);
      break;
   }
   case ParaUNSIGNED_SHORT :
   {
      newBuf = new unsigned short[count];
      memcpy(newBuf, buffer, sizeof(unsigned short)*count);
      break;
   }
   case ParaUNSIGNED :
   {
      newBuf = new unsigned int[count];
      memcpy(newBuf, buffer, sizeof(unsigned int)*count);
      break;
   }
   case ParaUNSIGNED_LONG :
   {
      newBuf = new unsigned long[count];
      memcpy(newBuf, buffer, sizeof(unsigned long)*count);
      break;
   }
   case ParaFLOAT :
   {
      newBuf = new float[count];
      memcpy(newBuf, buffer, sizeof(float)*count);
      break;
   }
   case ParaDOUBLE :
   {
      newBuf = new double[count];
      memcpy(newBuf, buffer, sizeof(double)*count);
      break;
   }
   case ParaLONG_DOUBLE :
   {
      newBuf = new long double[count];
      memcpy(newBuf, buffer, sizeof(long double)*count);
      break;
   }
   case ParaBYTE :
   {
      newBuf = new char[count];
      memcpy(newBuf, buffer, sizeof(char)*count);
      break;
   }
   case ParaSIGNED_CHAR :
   {
      newBuf = new char[count];
      memcpy(newBuf, buffer, sizeof(char)*count);
      break;
   }
   case ParaLONG_LONG :
   {
      newBuf = new long long[count];
      memcpy(newBuf, buffer, sizeof(long long)*count);
      break;
   }
   case ParaUNSIGNED_LONG_LONG :
   {
      newBuf = new unsigned long long[count];
      memcpy(newBuf, buffer, sizeof(unsigned long long)*count);
      break;
   }
   case ParaBOOL :
   {
      newBuf = new bool[count];
      memcpy(newBuf, buffer, sizeof(bool)*count);
      break;
   }
   default :
      THROW_LOGICAL_ERROR2("This type is not implemented. Type = ", datatypeId);
   }

   return newBuf;
}

void
ParaCommPth::copy(
      void *dest, void *src, int count, int datatypeId
      )
{

   if( count == 0 ) return;

   switch(datatypeId)
   {
   case ParaCHAR :
   {
      memcpy(dest, src, sizeof(char)*count);
      break;
   }
   case ParaSHORT :
   {
      memcpy(dest, src, sizeof(short)*count);
      break;
   }
   case ParaINT :
   {
      memcpy(dest, src, sizeof(int)*count);
      break;
   }
   case ParaLONG :
   {
      memcpy(dest, src, sizeof(long)*count);
      break;
   }
   case ParaUNSIGNED_CHAR :
   {
      memcpy(dest, src, sizeof(unsigned char)*count);
      break;
   }
   case ParaUNSIGNED_SHORT :
   {
      memcpy(dest, src, sizeof(unsigned short)*count);
      break;
   }
   case ParaUNSIGNED :
   {
      memcpy(dest, src, sizeof(unsigned int)*count);
      break;
   }
   case ParaUNSIGNED_LONG :
   {
      memcpy(dest, src, sizeof(unsigned long)*count);
      break;
   }
   case ParaFLOAT :
   {
      memcpy(dest, src, sizeof(float)*count);
      break;
   }
   case ParaDOUBLE :
   {
      memcpy(dest, src, sizeof(double)*count);
      break;
   }
   case ParaLONG_DOUBLE :
   {
      memcpy(dest, src, sizeof(long double)*count);
      break;
   }
   case ParaBYTE :
   {
      memcpy(dest, src, sizeof(char)*count);
      break;
   }
   case ParaSIGNED_CHAR :
   {
      memcpy(dest, src, sizeof(char)*count);
      break;
   }
   case ParaLONG_LONG :
   {
      memcpy(dest, src, sizeof(long long)*count);
      break;
   }
   case ParaUNSIGNED_LONG_LONG :
   {
      memcpy(dest, src, sizeof(unsigned long long)*count);
      break;
   }
   case ParaBOOL :
   {
      memcpy(dest, src, sizeof(bool)*count);
      break;
   }
   default :
      THROW_LOGICAL_ERROR2("This type is not implemented. Type = ", datatypeId);
   }

}

void
ParaCommPth::freeMem(
      void* buffer,
      int count,
      const int datatypeId
      )
{

   if( count == 0 ) return;

   switch(datatypeId)
   {
   case ParaCHAR :
   {
      delete [] static_cast<char *>(buffer);
      break;
   }
   case ParaSHORT :
   {
      delete [] static_cast<short *>(buffer);
      break;
   }
   case ParaINT :
   {
      delete [] static_cast<int *>(buffer);
      break;
   }
   case ParaLONG :
   {
      delete [] static_cast<long *>(buffer);
      break;
   }
   case ParaUNSIGNED_CHAR :
   {
      delete [] static_cast<unsigned char *>(buffer);
      break;
   }
   case ParaUNSIGNED_SHORT :
   {
      delete [] static_cast<unsigned short *>(buffer);
      break;
   }
   case ParaUNSIGNED :
   {
      delete [] static_cast<unsigned int *>(buffer);
      break;
   }
   case ParaUNSIGNED_LONG :
   {
      delete [] static_cast<unsigned long *>(buffer);
      break;
   }
   case ParaFLOAT :
   {
      delete [] static_cast<float *>(buffer);
      break;
   }
   case ParaDOUBLE :
   {
      delete [] static_cast<double *>(buffer);
      break;
   }
   case ParaLONG_DOUBLE :
   {
      delete [] static_cast<long double *>(buffer);
      break;
   }
   case ParaBYTE :
   {
      delete [] static_cast<char *>(buffer);
      break;
   }
   case ParaSIGNED_CHAR :
   {
      delete [] static_cast<char *>(buffer);
      break;
   }
   case ParaLONG_LONG :
   {
      delete [] static_cast<long long *>(buffer);
      break;
   }
   case ParaUNSIGNED_LONG_LONG :
   {
      delete [] static_cast<unsigned long long *>(buffer);
      break;
   }
   case ParaBOOL :
   {
      delete [] static_cast<bool *>(buffer);;
      break;
   }
   default :
      THROW_LOGICAL_ERROR2("This type is not implemented. Type = ", datatypeId);
   }

}

int
ParaCommPth::bcast(
   void* buffer,
   int count,
   const int datatypeId,
   int root
   )
{
   if( getRank() == root )
   {
      for(int i=0; i < comSize; i++)
      {
         if( i != root )
         {
            send(buffer, count, datatypeId, i, -1);
         }
      }
   }
   else
   {
      receive(buffer, count, datatypeId, root, -1);
   }
   return 0;
}

int
ParaCommPth::send(
   void* buffer,
   int count,
   const int datatypeId,
   int dest,
   const int tag
   )
{
   messageQueueTable[dest]->enqueue(
         new MessageQueueElement(getRank(), count, datatypeId, tag,
               allocateMemAndCopy(buffer, count, datatypeId) ) );
   TAG_TRACE (Send, To, dest, tag);
   return 0;
}

int
ParaCommPth::receive(
   void* buffer,
   int count,
   const int datatypeId,
   int source,
   const int tag
   )
{
   if( !messageQueueTable[getRank()]->checkElement(source, datatypeId, tag) )
   {
      messageQueueTable[getRank()]->watiMessage(source, datatypeId, tag);
   }
   MessageQueueElement *elem = messageQueueTable[getRank()]->extarctElement(source, datatypeId, tag);
   assert(elem);
   copy( buffer, elem->getData(), count, datatypeId );
   freeMem(elem->getData(), count, datatypeId );
   delete elem;
   TAG_TRACE (Recv, From, source, tag);
   return 0;
}

int
ParaCommPth::waitSpecTagFromSpecSource(
      const int source,
      const int datatypeId,
      const int tag,
      int *receivedTag
      )
{
   /*
   // Just wait, iProbe and receive will be performed after this call
   messageQueueTable[getRank()]->watiMessage(source, datatypeId, tag);
   TAG_TRACE (Probe, From, source, tag);
   return 0;
   */

   messageQueueTable[getRank()]->watiMessage(source, receivedTag);
   TAG_TRACE (Probe, From, source, *receivedTag);
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
ParaCommPth::probe(
   int* source,
   int* tag
   )
{
   messageQueueTable[getRank()]->waitMessage();
   MessageQueueElement *elem = messageQueueTable[getRank()]->getHead();
   *source = elem->getSource();
   *tag = elem->getTag();
   TAG_TRACE (Probe, From, *source, *tag);
}

bool
ParaCommPth::iProbe(
   int* source,
   int* tag
   )
{
   bool flag;
   flag = !(messageQueueTable[getRank()]->isEmpty());
   if( flag )
   {
      MessageQueueElement *elem = messageQueueTable[getRank()]->getHead();
      *source = elem->getSource();
      *tag = elem->getTag();
      TAG_TRACE (Iprobe, From, *source, *tag);
   }
   return flag;
}

int
ParaCommPth::uTypeSend(
   void* buffer,
   const int datatypeId,
   int dest,
   const int tag
   )
{
   messageQueueTable[dest]->enqueue(
         new MessageQueueElement(getRank(), 1, datatypeId, tag, buffer ) );
   TAG_TRACE (Send, To, dest, tag);
   return 0;
}

int
ParaCommPth::uTypeReceive(
   void** buffer,
   const int datatypeId,
   int source,
   const int tag
   )
{
   if( !messageQueueTable[getRank()]->checkElement(source, datatypeId, tag) )
   {
      messageQueueTable[getRank()]->watiMessage(source, datatypeId, tag);
   }
   MessageQueueElement *elem = messageQueueTable[getRank()]->extarctElement(source, datatypeId, tag);
   assert(elem);
   *buffer = elem->getData();
   delete elem;
   TAG_TRACE (Recv, From, source, tag);
   return 0;
}
