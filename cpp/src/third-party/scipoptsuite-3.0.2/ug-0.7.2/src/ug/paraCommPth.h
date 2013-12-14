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

/**@file    paraCommPth.h
 * @brief   ParaComm extension for Pthreads communication
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_COMM_PTH_H__
#define __PARA_COMM_PTH_H__

#include <pthread.h>
#include <stdexcept>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cassert>
#include "paraDef.h"
#include "paraComm.h"
#include "paraSysTimer.h"
#include "paraParamSetPth.h"
#include "paraTimerPth.h"

#define HashEntry( tid ) ( hashCode(tid) % UG::HashTableSize )

#define TAG_TRACE( call, fromTo, sourceDest, tag ) \
   if( tagTraceFlag )  \
   {  \
      if( tag >= 0 ) \
      { \
         /* std::cout << " call = " << #call << ", Rank = " << getRank() << ", tag = " << tag << ", " << tagStringTable[tag] << std::endl; */ \
         *getOstream() << timer.getRTimeInterval() << " [Rank = " << getRank() << "] " << #call << " " << #fromTo  \
         << " " << sourceDest << " with Tag = " << tagStringTable[tag] << std::endl; \
      } \
      else \
      { \
         /* std::cout << " call = " << #call << ", Rank = " << getRank() << ", tag = " << tag << std::endl; */  \
         *getOstream() << timer.getRTimeInterval() << " [Rank = " << getRank() << "] " << #call << " " << #fromTo  \
         << " " << sourceDest << " as broadcast" << std::endl; \
      } \
}

namespace UG
{

static const int USER_TYPE_FIRST                = TYPE_LAST + 1;
static const int ParaInstanceType               = USER_TYPE_FIRST +  0;
static const int ParaSolutionType               = USER_TYPE_FIRST +  1;
static const int ParaParamSetType               = USER_TYPE_FIRST +  2;
static const int ParaNodeType                   = USER_TYPE_FIRST +  3;
static const int ParaSolverStateType            = USER_TYPE_FIRST +  4;
static const int ParaCalculationStateType       = USER_TYPE_FIRST +  5;
static const int ParaSolverTerminationStateType = USER_TYPE_FIRST +  6;
static const int ParaRacingRampUpParamType      = USER_TYPE_FIRST +  7;
static const int ParaSolverDiffParamType        = USER_TYPE_FIRST +  8;
static const int ParaInitialStatType            = USER_TYPE_FIRST +  9;

static const int HashTableSize = 751;

class MessageQueueElement
{
   int source;                 // source thread rank of this message
   int count;                  // number of elements of the data type
   int dataTypeId;             // data type id
   int tag;                    // -1 : in case of broadcast message
   void *data;                 // NOTE : For basic data types, this is copy of sender side memory.
                               //        When the memory is copied at receive function, the memory
                               //        have to be freed.
                               //        For user defined data type, this is the receiver side memory,
                               //        because it is better to allocate memory in the sender side for
                               //        mutex locking. Sender side functions have to allocate memory.
                               //        In this case, memory do not hvae to freed. The memory is for
                               //        receiver side.
   MessageQueueElement *next;
public:
   MessageQueueElement(
         ) : source(-1), count(0), dataTypeId(-1), tag(-1), data(0), next(0)
   {
   }
   MessageQueueElement(
         int inSource,
         int inCount,
         int inDataTypeId,
         int inTag,
         void *inData
         ) : source(inSource), count(inCount), dataTypeId(inDataTypeId), tag(inTag), data(inData), next(0)
   {
   }
   ~MessageQueueElement(
         )
   {
   }

   int getSource(
         )
   {
      return source;
   }

   int getCount(
         )
   {
      return count;
   }

   int getDataTypeId(
         )
   {
      return dataTypeId;
   }

   int getTag(
         )
   {
      return tag;
   }

   void *getData(
         )
   {
      return data;
   }

   MessageQueueElement *getNext(
         )
   {
      return next;
   }

   void link(
         MessageQueueElement *nextElement
         )
   {
      next = nextElement;
   }
};

class MessageQueueTableElement
{
   bool                sentMessage;
   pthread_mutex_t     queueLock;
   pthread_cond_t      sentMsg;
   MessageQueueElement *head;
   MessageQueueElement *tail;
public:
   MessageQueueTableElement(
         ) : sentMessage(false), head(0), tail(0)
   {
      pthread_mutex_init(&queueLock, NULL);
      pthread_cond_init(&sentMsg, NULL);
   }

   ~MessageQueueTableElement(
         )
   {
      pthread_mutex_lock(&queueLock);
      while( head )
      {
         MessageQueueElement *next = head->getNext();
         delete head;
         head = next;
      }
      pthread_mutex_unlock(&queueLock);
   }

   MessageQueueElement *checkElement(
         int source,
         int datatypeId,
         int tag
         )
   {
      MessageQueueElement *ret = 0;

      pthread_mutex_lock(&queueLock);
      for( MessageQueueElement *current = head; current; current = current->getNext() )
      {
         if( current->getSource() == source
               && current->getDataTypeId() == datatypeId
               && current->getTag() == tag )
         {
            ret = current;
            break;
         }
      }
      pthread_mutex_unlock(&queueLock);
      return ret;
   }

   MessageQueueElement *extarctElement(
         int source,
         int datatypeId,
         int tag
         )
   {
      MessageQueueElement *ret = 0;

      pthread_mutex_lock(&queueLock);
      MessageQueueElement *prev = head;
      MessageQueueElement *current = head;
      while( current )
      {
         MessageQueueElement *next = current->getNext();
         if( current->getSource() == source
               && current->getDataTypeId() == datatypeId
               && current->getTag() == tag )
         {
            ret = current;
            ret->link(0);
            if( current == head )
            {
               if( current == tail )
               {
                  head = 0;
                  tail = 0;
                  sentMessage = false;
               }
               else
               {
                  head=next;
               }
            }
            else
            {
               if( current == tail )
               {
                  tail = prev;
               }
               prev->link(next);
            }
            break;
         }
         else
         {
            prev = current;
            current = next;
         }
      }
      pthread_mutex_unlock(&queueLock);
      return ret;
   }

   // This method is only for desctructor or ParaCommPth. No lock is necessary.
   MessageQueueElement *extarctElement(
         )
   {
      if( head == 0 ) return 0;

      MessageQueueElement *current = head;
      MessageQueueElement *next = current->getNext();
      current->link(0);
      if( current == tail )
      {
         head = 0;
         tail = 0;
         sentMessage = false;
      }
      else
      {
         head=next;
      }
      return current;
   }

   void enqueue(
         MessageQueueElement *newElement
         )
   {
      pthread_mutex_lock(&queueLock);
      if( tail == 0 )
      {
         head = tail = newElement;
      }
      else
      {
         tail->link(newElement);
         tail = newElement;
      }
      sentMessage = true;
      pthread_cond_signal(&sentMsg);
      pthread_mutex_unlock(&queueLock);
   }

   MessageQueueElement *getHead(
         )
   {
      return head;
   }

   bool isEmpty(
         )
   {
      bool empty = true;
      pthread_mutex_lock(&queueLock);
      if( head ) empty = false;
      pthread_mutex_unlock(&queueLock);
      return empty;
   }

   void waitMessage(
         )
   {
      pthread_mutex_lock(&queueLock);

      while( sentMessage == false )
      {
         pthread_cond_wait(&sentMsg, &queueLock);
      }

      pthread_mutex_unlock(&queueLock);
   }

   void watiMessage(
         int source,
         int datatypeId,
         int tag
         )
   {
      pthread_mutex_lock(&queueLock);

      for(;;)
      {
         while( sentMessage == false )
         {
            pthread_cond_wait(&sentMsg, &queueLock);
         }
         MessageQueueElement *current = head;
         while( current )
         {
            MessageQueueElement *next = current->getNext();
            if( current->getSource() == source
                  && current->getDataTypeId() == datatypeId
                  && current->getTag() == tag )
            {
               break;
            }
            current = next;
         }
         if( current ) break;
         sentMessage = false;
      }

      pthread_mutex_unlock(&queueLock);
   }

   void watiMessage(
         int source,
         int *tag
         )
   {
      pthread_mutex_lock(&queueLock);

      for(;;)
      {
         while( sentMessage == false )
         {
            pthread_cond_wait(&sentMsg, &queueLock);
         }
         MessageQueueElement *current = head;
         while( current )
         {
            MessageQueueElement *next = current->getNext();
            if( current->getSource() == source )
            {
               *tag = current->getTag();
               break;
            }
            current = next;
         }
         if( current ) break;
         sentMessage = false;
      }

      pthread_mutex_unlock(&queueLock);
   }

};

class HashTableElement
{
   pthread_t          tid;          /*< thread id of this thread */
   int                rank;         /*< rank of this thread */
   std::ostream       *tos;         /*< tag trace stream for this thread */
   HashTableElement   *next;        /*< next HashTableElement pointer */
public:


   HashTableElement(
         ) : rank(0), tos(0), next(0)
   {
   }

   HashTableElement(
         pthread_t inTid,
         int       inRank,
         ParaParamSet *paraParamSet
         ) : tid(inTid), rank(inRank), tos(0), next(0)
   {
      if( paraParamSet->getBoolParamValue(TagTrace) )
      {
         if( paraParamSet->isStringParamDefaultValue(TagTraceFileName) )
         {
            tos = &std::cout;
         }
         else
         {
            std::ostringstream s;
            s << paraParamSet->getStringParamValue(TagTraceFileName) << inRank;
            std::ofstream  *ofs = new std::ofstream();
            ofs->open(s.str().c_str());
            tos = ofs;
         }
      }
   }

   ~HashTableElement(
         )
   {
      if( tos )
      {
         std::ofstream *ofs = dynamic_cast<std::ofstream *>(tos);
         ofs->close();
         delete ofs;
      }
   }

   pthread_t getTid(
         )
   {
      return tid;
   }

   int getRank(
         )
   {
      return rank;
   }

   std::ostream *getOstream(
         )
   {
      return tos;
   }

   HashTableElement *getNext(
         )
   {
      return next;
   }

   void link(
         HashTableElement *nextElement
         )
   {
      next = nextElement;
   }
};

class ParaCalculationState;
class ParaNode;
class ParaSolverState;
class ParaSolverTerminationState;
class ParaDiffSubproblem;
class ParaInstance;
class ParaSolution;

class ParaCommPth : public ParaComm
{
protected:
   int                      comSize;
   bool                     tagTraceFlag;
   int                      **token;    // index 0: token
                                        // index 1: token color
                                        //           -1: green
                                        //           > 0: yellow ( termination origin solver number )
                                        //           -2: red ( means the solver can terminate )
   ParaSysTimer             timer;
   static const char        *tagStringTable[];
   static HashTableElement  *hashtable[HashTableSize];
   MessageQueueTableElement **messageQueueTable;

   pthread_mutex_t           *tokenAccessLock;
   pthread_mutex_t           applicationLock;
   pthread_mutex_t           rankLock;

   unsigned int hashCode(pthread_t tid);
   void *allocateMemAndCopy(void* buffer, int count, const int datatypeId);
   void copy(void *dest, void *src, int count, int datatypeId);
   void freeMem(void* buffer, int count, const int datatypeId);

public:
   ParaCommPth() : comSize(-1),
         tagTraceFlag(false), token(0), messageQueueTable(0)
   {
      pthread_mutex_init(&applicationLock, NULL);
      pthread_mutex_init(&rankLock, NULL);
   }

   void init( int argc, char **argv );
   virtual ~ParaCommPth();
   double getStartTime(){ // should not used
      return timer.getStartTime();
   }
   int getRank();
   int getRank( pthread_t tid);
   int getSize(){ return comSize; }
   void lcInit(ParaParamSet *paraParamSet);
   void solverInit(ParaParamSet *paraParamSet){}
   void solverInitPth(pthread_t tid, int rank, ParaParamSet *paraParamSet);
   void abort(){abort();}
   bool waitTerminatedMessage(){ return false;}
   bool waitToken(int rank);
   void passToken(int rank);
   bool passTermToken(int rank);
   void setToken(int rank, int *inToken);
   void waitUntilRegistered();
   void registedAllSolvers();
   std::ostream *getOstream();

   void lockApp()
   {
      pthread_mutex_lock(&applicationLock);
   }
   void unlockApp()
   {
      pthread_mutex_unlock(&applicationLock);
   }

   void lockRank()
   {
      pthread_mutex_lock(&rankLock);
   }
   void unlockRank()
   {
      pthread_mutex_unlock(&rankLock);
   }

   ParaCalculationState *createParaCalculationState();
   ParaCalculationState *createParaCalculationState(
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
           );
   ParaNode *createParaNode();
   ParaNode *createParaNode(
               NodeId inNodeId,
               NodeId inGeneratorNodeId,
               int inDepth,
               double inDualBoundValue,
               double inOriginalDualBoundValue,
               double inEstimatedValue,
               ParaDiffSubproblem *inDiffSubproblem 
            );
   ParaParamSet *createParaParamSet();
   ParaSolverState *createParaSolverState();
   ParaSolverState *createParaSolverState(
               int racingStage,
               unsigned int notificationId,
               int lcId,
               int globalSubtreeId,
               long long nodesSolved,
               int nodesLeft,
               double bestDualBoundValue,
               double globalBestPrimalBoundValue,
               double detTime
           );
   ParaSolverTerminationState *createParaSolverTerminationState();
   ParaSolverTerminationState *createParaSolverTerminationState(
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
           );
   ParaTimer *createParaTimer(){ return new ParaTimerPth(); }
   int bcast( void* buffer, int count, const int datatypeId, int root );
   int send( void* bufer, int count, const int datatypeId, int dest, const int tag );
   int receive( void* bufer, int count, const int datatypeId, int source, const int tag );
   int waitSpecTagFromSpecSource(const int source, const int datatypeId, const int tag, int *receivedTag);
   void probe(int *source, int *tag);
   bool iProbe(int *source, int *tag);  /** return value: true => there is message, false => no message */

   /** For Created Datatypes */
   int uTypeSend( void* bufer, const int datatypeId, int dest, int tag );
   int uTypeReceive( void** bufer, const int datatypeId, int source, int tag );
};

#define DEF_PARA_COMM( para_comm, comm ) ParaCommPth *para_comm = dynamic_cast< ParaCommPth* >(comm)

}

#endif  // __PARA_COMM_PTH_H__
