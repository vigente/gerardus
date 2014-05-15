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

/**@file    scipParaRacingRampUpParamSetPth.cpp
 * @brief   ScipParaRacingRampUpParamSet extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scipParaCommPth.h"
#include "scipDiffParamSetPth.h"
#include "scipParaRacingRampUpParamSetPth.h"
#include <cstring>

using namespace ParaSCIP;

/** create Datatype */
ScipParaRacingRampUpParamSetPth *
ScipParaRacingRampUpParamSetPth::createDatatype(
      )
{
   ScipDiffParamSetPth *scipDiffParamSetPth = dynamic_cast<ScipDiffParamSetPth *>(scipDiffParamSet);

   if( scipDiffParamSetInfo )
   {
      return new ScipParaRacingRampUpParamSetPth(
            terminationCriteria,
            nNodesLeft,
            timeLimit,
            scipRacingParamSeed,
            permuteProbSeed,
            generateBranchOrderSeed,
            scipDiffParamSetPth->clone()
            );
   }
   else
   {
      return new ScipParaRacingRampUpParamSetPth(
            terminationCriteria,
            nNodesLeft,
            timeLimit,
            scipRacingParamSeed,
            permuteProbSeed,
            generateBranchOrderSeed,
            0
            );
   }
}

int
ScipParaRacingRampUpParamSetPth::send(
      UG::ParaComm *comm,
      int dest)
{

   DEF_SCIP_PARA_COMM( commPth, comm);

   PARA_COMM_CALL(
      commPth->uTypeSend((void *)createDatatype(), UG::ParaRacingRampUpParamType, dest, UG::TagRacingRampUpParamSet)
   );

   return 0;

}

int
ScipParaRacingRampUpParamSetPth::receive(
      UG::ParaComm *comm,
      int source)
{

   DEF_SCIP_PARA_COMM( commPth, comm);

   ScipParaRacingRampUpParamSetPth *received;
   PARA_COMM_CALL(
      commPth->uTypeReceive((void **)&received, UG::ParaRacingRampUpParamType, source, UG::TagRacingRampUpParamSet)
   );

   terminationCriteria = received->terminationCriteria;
   nNodesLeft = received->nNodesLeft;
   timeLimit = received->timeLimit;
   scipRacingParamSeed = received->scipRacingParamSeed;
   permuteProbSeed = received->permuteProbSeed;
   generateBranchOrderSeed = received->generateBranchOrderSeed;
   scipDiffParamSetInfo = received->scipDiffParamSetInfo;
   if( scipDiffParamSetInfo )
   {
      ScipDiffParamSetPth *scipDiffParamSetPth = dynamic_cast<ScipDiffParamSetPth *>(received->scipDiffParamSet);
      scipDiffParamSet = scipDiffParamSetPth->clone();
   }

   delete received;

   return 0;

}
