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

/**@file    scipParaRacingRampUpParamSetMpi.h
 * @brief   ScipParaRacingRampUpParamSet extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_RACING_RAMP_UP_PARAM_SET_MPI_H__
#define __SCIP_PARA_RACING_RAMP_UP_PARAM_SET_MPI_H__

#include <mpi.h>
#include "ug/paraTagDefMpi.h"
#include "ug/paraCommMpiWorld.h"
#include "scipParaRacingRampUpParamSet.h"

namespace ParaSCIP
{

/** The difference between instance and subproblem: this is base class */
class ScipParaRacingRampUpParamSetMpi : public ScipParaRacingRampUpParamSet
{

   /** create ScipParaRacingRampUpParamSet datatype */
   MPI::Datatype createDatatype();

public:
   /** default constructor */
   ScipParaRacingRampUpParamSetMpi(
         )
   {
   }

   /** Constructor */
   ScipParaRacingRampUpParamSetMpi(
         int inTerminationCriteria,
         int inNNodesLeft,
         double inTimeLimit,
         int inScipRacingParamSeed,
         int inPermuteProbSeed,
         int inGenerateBranchOrderSeed,
         ScipDiffParamSet *inScipDiffParamSet
         ) : ScipParaRacingRampUpParamSet(inTerminationCriteria, inNNodesLeft, inTimeLimit,
         inScipRacingParamSeed,inPermuteProbSeed, inGenerateBranchOrderSeed, inScipDiffParamSet)
   {
   }

   /** destructor */
   ~ScipParaRacingRampUpParamSetMpi(
         )
   {
   }

   int send(UG::ParaComm *comm, int dest);
   int receive(UG::ParaComm *comm, int source);
};

}

#endif    // __SCIP_PARA_RACING_RAMP_UP_PARAM_SET_MPI_H__

