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

/**@file    paraSolverTerminationStatePth.cpp
 * @brief   ParaSolverTerminationState extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraCommPth.h"
#include "paraSolverTerminationStatePth.h"

using namespace UG;

ParaSolverTerminationStatePth *
ParaSolverTerminationStatePth::createDatatype(
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
         idleTimeAfterLastParaNode,
         idleTimeToWaitNotificationId,
         idleTimeToWaitAckCompletion,
         idleTimeToWaitToken,
         totalRootNodeTime,
         minRootNodeTime,
         maxRootNodeTime,
         detTime
         );
}

void
ParaSolverTerminationStatePth::send(
      ParaComm *comm,
      int destination,
      int tag
      )
{
   DEF_PARA_COMM( commPth, comm);

   PARA_COMM_CALL(
      commPth->uTypeSend((void *)createDatatype(), ParaSolverTerminationStateType, destination, tag)
   );
}

void
ParaSolverTerminationStatePth::receive(
      ParaComm *comm,
      int source,
      int tag
      )
{
   DEF_PARA_COMM( commPth, comm);

   ParaSolverTerminationStatePth *received;
   PARA_COMM_CALL(
      commPth->uTypeReceive((void **)&received, ParaSolverTerminationStateType, source, tag)
   );
   interrupted = received->interrupted;
   rank = received->rank;
   totalNSolved = received->totalNSolved;
   minNSolved = received->minNSolved;
   maxNSolved = received->maxNSolved;
   totalNSent = received->totalNSent;
   totalNImprovedIncumbent = received->totalNImprovedIncumbent;
   nParaNodesReceived = received->nParaNodesReceived;
   nParaNodesSolved = received->nParaNodesSolved;
   nParaNodesSolvedAtRoot = received->nParaNodesSolvedAtRoot;
   nParaNodesSolvedAtPreCheck = received->nParaNodesSolvedAtPreCheck;
   runningTime = received->runningTime;
   idleTimeToFirstParaNode = received->idleTimeToFirstParaNode;
   idleTimeBetweenParaNodes = received->idleTimeBetweenParaNodes;
   idleTimeAfterLastParaNode = received->idleTimeAfterLastParaNode;
   idleTimeToWaitNotificationId = received->idleTimeToWaitNotificationId;
   idleTimeToWaitAckCompletion = received->idleTimeToWaitAckCompletion;
   idleTimeToWaitToken = received->idleTimeToWaitToken;
   totalRootNodeTime = received->totalRootNodeTime;
   minRootNodeTime = received->minRootNodeTime;
   maxRootNodeTime = received->maxRootNodeTime;
   detTime = received->detTime;

   delete received;
}
