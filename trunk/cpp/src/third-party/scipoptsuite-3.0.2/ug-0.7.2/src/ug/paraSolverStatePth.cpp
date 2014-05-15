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

/**@file    paraSolverStatePth.cpp
 * @brief   ParaSolverState extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraCommPth.h"
#include "paraSolverStatePth.h"

using namespace UG;

ParaSolverStatePth *
ParaSolverStatePth::createDatatype(
      )
{
   return new ParaSolverStatePth(
         racingStage,
         notificationId,
         lcId,
         globalSubtreeIdInLc,
         nNodesSolved,
         nNodesLeft,
         bestDualBoundValue,
         globalBestPrimalBoundValue,
         detTime
         );
}

void
ParaSolverStatePth::send(
      ParaComm *comm,
      int destination,
      int tag
      )
{
   assert(nNodesLeft >= 0);
   assert(bestDualBoundValue >= -1e+10);
   DEF_PARA_COMM( commPth, comm);

   PARA_COMM_CALL(
      commPth->uTypeSend((void *)createDatatype(), ParaSolverStateType, destination, tag)
   );
}

void
ParaSolverStatePth::receive(
      ParaComm *comm,
      int source,
      int tag
      )
{
   DEF_PARA_COMM( commPth, comm);

   ParaSolverStatePth *received;
   PARA_COMM_CALL(
      commPth->uTypeReceive((void **)&received, ParaSolverStateType, source, tag)
   );

   racingStage = received->racingStage;
   notificationId = received->notificationId;
   lcId = received->lcId;
   globalSubtreeIdInLc = received->globalSubtreeIdInLc;
   nNodesSolved = received->nNodesSolved;
   nNodesLeft = received->nNodesLeft;
   bestDualBoundValue = received->bestDualBoundValue;
   globalBestPrimalBoundValue = received->globalBestPrimalBoundValue;
   detTime = received->detTime;

   delete received;

}
