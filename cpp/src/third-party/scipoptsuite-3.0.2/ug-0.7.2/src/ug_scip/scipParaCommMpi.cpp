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

/**@file    scipParaCommMpi.cpp
 * @brief   SCIP ParaComm extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scipParaCommMpi.h"
#include "scipParaInstanceMpi.h"
#include "scipParaDiffSubproblemMpi.h"
#include "scipParaSolutionMpi.h"
#include "scipParaInitialStatMpi.h"
#include "scipParaRacingRampUpParamSetMpi.h"
#include "scipDiffParamSetMpi.h"

using namespace ParaSCIP;

/*******************************************************************************
* transfer object factory
*******************************************************************************/
UG::ParaDiffSubproblem *
ScipParaCommMpi::createParaDiffSubproblem(
    )
{ 
    return new ScipParaDiffSubproblemMpi(); 
}

UG::ParaInitialStat* 
ScipParaCommMpi::createParaInitialStat(
    )
{ 
    return new ScipParaInitialStatMpi(); 
}

UG::ParaRacingRampUpParamSet* 
ScipParaCommMpi::createParaRacingRampUpParamSet(
    )
{ 
    return new ScipParaRacingRampUpParamSetMpi(); 
}

UG::ParaInstance*
ScipParaCommMpi::createParaInstance(
    )
{ 
    return new ScipParaInstanceMpi(); 
}

UG::ParaSolution*
ScipParaCommMpi::createParaSolution(
    )
{ 
    return new ScipParaSolutionMpi(); 
}

ScipParaInstance*
ScipParaCommMpi::createScipParaInstance(
    SCIP *scip, 
    int method
    )
{
    return new ScipParaInstanceMpi(scip, method);
}

ScipParaSolution*
ScipParaCommMpi::createScipParaSolution(
    SCIP_Real objval, 
    int inNvars, 
    SCIP_VAR ** vars, 
    SCIP_Real *vals
    )
{
    return new ScipParaSolutionMpi(objval, inNvars, vars, vals);
}

ScipParaDiffSubproblem*
ScipParaCommMpi::createScipParaDiffSubproblem(
         SCIP *scip,
         ScipParaSolver *scipParaSolver,
         int nNewBranchVars,
         SCIP_VAR **newBranchVars,
         SCIP_Real *newBranchBounds,
         SCIP_BOUNDTYPE *newBoundTypes
         )
{
    return new ScipParaDiffSubproblemMpi(
         scip,
         scipParaSolver,
         nNewBranchVars,
         newBranchVars,
         newBranchBounds,
         newBoundTypes
         );
}

ScipParaInitialStat*
ScipParaCommMpi::createScipParaInitialStat(
         SCIP *scip
         )
{
    return new ScipParaInitialStatMpi(
         scip
         );
}

ScipParaInitialStat*
ScipParaCommMpi::createScipParaInitialStat(
            int inMaxDepth,
            int inMaxTotalDepth,
            int inNVarBranchStatsDown,
            int inNVarBranchStatsUp,
            int *inIdxLBranchStatsVarsDown,
            int *inNVarBranchingDown,
            int *inIdxLBranchStatsVarsUp,
            int *inNVarBranchingUp,
            SCIP_Real *inDownpscost,
            SCIP_Real *inDownvsids,
            SCIP_Real *inDownconflen,
            SCIP_Real *inDowninfer,
            SCIP_Real *inDowncutoff,
            SCIP_Real *inUppscost,
            SCIP_Real *inUpvsids,
            SCIP_Real *inUpconflen,
            SCIP_Real *inUpinfer,
            SCIP_Real *inUpcutoff
         )
{
    return new ScipParaInitialStatMpi(
            inMaxDepth,
            inMaxTotalDepth,
            inNVarBranchStatsDown,
            inNVarBranchStatsUp,
            inIdxLBranchStatsVarsDown,
            inNVarBranchingDown,
            inIdxLBranchStatsVarsUp,
            inNVarBranchingUp,
            inDownpscost,
            inDownvsids,
            inDownconflen,
            inDowninfer,
            inDowncutoff,
            inUppscost,
            inUpvsids,
            inUpconflen,
            inUpinfer,
            inUpcutoff
         );
}

ScipParaRacingRampUpParamSet *
ScipParaCommMpi::createScipParaRacingRampUpParamSet(
         int inTerminationCriteria,
         int inNNodesLeft,
         double inTimeLimit,
         int inScipRacingParamSeed,
         int inPermuteProbSeed,
         int inGenerateBranchOrderSeed,
         ScipDiffParamSet *inScipDiffParamSet
         )
{
    return new ScipParaRacingRampUpParamSetMpi(
               inTerminationCriteria,
               inNNodesLeft,
               inTimeLimit,
               inScipRacingParamSeed,
               inPermuteProbSeed,
               inGenerateBranchOrderSeed,
               inScipDiffParamSet
               );
}

ScipDiffParamSet *
ScipParaCommMpi::createScipDiffParamSet()
{
    return new ScipDiffParamSetMpi();
}

ScipDiffParamSet *
ScipParaCommMpi::createScipDiffParamSet(
        SCIP *scip
        )
{
    return new ScipDiffParamSetMpi(scip);
}

