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

/**@file    scipParaCommPth.cpp
 * @brief   SCIP ParaComm extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scipParaCommPth.h"
#include "scipParaInstancePth.h"
#include "scipParaDiffSubproblemPth.h"
#include "scipParaSolutionPth.h"
#include "scipParaInitialStatPth.h"
#include "scipParaRacingRampUpParamSetPth.h"
#include "scipDiffParamSetPth.h"

using namespace ParaSCIP;

/*******************************************************************************
* transfer object factory
*******************************************************************************/
UG::ParaDiffSubproblem *
ScipParaCommPth::createParaDiffSubproblem(
    )
{ 
    return new ScipParaDiffSubproblemPth(); 
}

UG::ParaInitialStat* 
ScipParaCommPth::createParaInitialStat(
    )
{ 
    return new ScipParaInitialStatPth(); 
}

UG::ParaRacingRampUpParamSet* 
ScipParaCommPth::createParaRacingRampUpParamSet(
    )
{ 
    return new ScipParaRacingRampUpParamSetPth(); 
}

UG::ParaInstance*
ScipParaCommPth::createParaInstance(
    )
{ 
    return new ScipParaInstancePth(); 
}

UG::ParaSolution*
ScipParaCommPth::createParaSolution(
    )
{ 
    return new ScipParaSolutionPth(); 
}

ScipParaInstance*
ScipParaCommPth::createScipParaInstance(
    SCIP *scip, 
    int method
    )
{
    return new ScipParaInstancePth(scip, method);
}

ScipParaSolution*
ScipParaCommPth::createScipParaSolution(
    SCIP_Real objval, 
    int inNvars, 
    SCIP_VAR ** vars, 
    SCIP_Real *vals
    )
{
    return new ScipParaSolutionPth(objval, inNvars, vars, vals);
}

ScipParaDiffSubproblem*
ScipParaCommPth::createScipParaDiffSubproblem(
         SCIP *scip,
         ScipParaSolver *scipParaSolver,
         int nNewBranchVars,
         SCIP_VAR **newBranchVars,
         SCIP_Real *newBranchBounds,
         SCIP_BOUNDTYPE *newBoundTypes
         )
{
    return new ScipParaDiffSubproblemPth(
         scip,
         scipParaSolver,
         nNewBranchVars,
         newBranchVars,
         newBranchBounds,
         newBoundTypes
         );
}

ScipParaInitialStat*
ScipParaCommPth::createScipParaInitialStat(
         SCIP *scip
         )
{
    return new ScipParaInitialStatPth(
         scip
         );
}

ScipParaInitialStat*
ScipParaCommPth::createScipParaInitialStat(
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
    return new ScipParaInitialStatPth(
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
ScipParaCommPth::createScipParaRacingRampUpParamSet(
         int inTerminationCriteria,
         int inNNodesLeft,
         double inTimeLimit,
         int inScipRacingParamSeed,
         int inPermuteProbSeed,
         int inGenerateBranchOrderSeed,
         ScipDiffParamSet *inScipDiffParamSet
         )
{
    return new ScipParaRacingRampUpParamSetPth(
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
ScipParaCommPth::createScipDiffParamSet()
{
    return new ScipDiffParamSetPth();
}

ScipDiffParamSet *
ScipParaCommPth::createScipDiffParamSet(
        SCIP *scip
        )
{
    return new ScipDiffParamSetPth(scip);
}
