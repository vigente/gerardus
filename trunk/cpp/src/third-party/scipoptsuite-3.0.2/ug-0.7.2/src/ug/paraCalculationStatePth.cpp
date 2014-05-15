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

/**@file    paraCalculationStatePth.cpp
 * @brief   CalcutationStte object extension for Pthreads communication
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraCommPth.h"
#include "paraCalculationStatePth.h"

using namespace UG;

ParaCalculationStatePth*
ParaCalculationStatePth::createDatatype(
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

void
ParaCalculationStatePth::send(
      ParaComm *comm,
      int destination,
      int tag
      )
{
   DEF_PARA_COMM( commPth, comm);

   PARA_COMM_CALL(
      commPth->uTypeSend(createDatatype(), ParaCalculationStateType, destination, tag)
   );
}

void
ParaCalculationStatePth::receive(
      ParaComm *comm,
      int source,
      int tag
      )
{
   DEF_PARA_COMM( commPth, comm);

   ParaCalculationStatePth *received;

   PARA_COMM_CALL(
      commPth->uTypeReceive((void **)&received, ParaCalculationStateType, source, tag)
   );

   compTime = received->compTime;
   rootTime = received->rootTime;
   nSolved = received->nSolved;
   nSent = received->nSent;
   nImprovedIncumbent = received->nImprovedIncumbent;
   terminationState = received->terminationState;
   nSolvedWithNoPreprocesses = received->nSolvedWithNoPreprocesses;
   nSimplexIterRoot = received->nSimplexIterRoot;
   averageSimplexIter = received->averageSimplexIter;
   nRestarts = received->nRestarts;
   minIisum = received->minIisum;
   maxIisum = received->maxIisum;
   minNii = received->minNii;
   maxNii = received->maxNii;

   delete received;

}
