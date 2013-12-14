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

/**@file    scipParaSolver.cpp
 * @brief   ParaSolver extension for SCIP: Parallelized solver implementation for SCIP.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cfloat>
#include <typeinfo>
#include "ug/paraComm.h"
#include "ug/paraNode.h"
#include "ug/paraInstance.h"
#include "ug/paraSolver.h"
#include "ug/paraSolution.h"
#include "ug/paraInitialStat.h"
#include "objscip/objscip.h"
#include "scipParaObjMessageHdlr.h"
#include "scipParaObjCommPointHdlr.h"
#include "scipParaObjBranchRule.h"
#include "scipParaInitialStat.h"
#include "scipParaRacingRampUpParamSet.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

using namespace ParaSCIP;

/*
 * Callback methods of conflict handler
 */
#define CONFLICTHDLR_NAME      "conflictCollector"
#define CONFLICTHDLR_DESC      "conflict handler to collect conflicts"
#define CONFLICTHDLR_PRIORITY  +100000000
static
SCIP_DECL_CONFLICTEXEC(conflictExecCollector)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real lhs;
   int i;

   assert(conflicthdlr != NULL);
   assert(strcmp(SCIPconflicthdlrGetName(conflicthdlr), CONFLICTHDLR_NAME) == 0);
   assert(bdchginfos != NULL || nbdchginfos == 0);
   assert(result != NULL);

   /* don't process already resolved conflicts */
   if( resolved )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* create array of variables and coefficients: sum_{i \in P} x_i - sum_{i \in N} x_i >= 1 - |N| */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbdchginfos) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nbdchginfos) );
   lhs = 1.0;
   for( i = 0; i < nbdchginfos; ++i )
   {
      assert(bdchginfos != NULL);

      vars[i] = SCIPbdchginfoGetVar(bdchginfos[i]);

      /* we can only treat binary variables */
      /**@todo extend linear conflict constraints to some non-binary cases */
      if( !SCIPvarIsBinary(vars[i]) )
         break;

      /* check whether the variable is fixed to zero (P) or one (N) in the conflict set */
      if( SCIPbdchginfoGetNewbound(bdchginfos[i]) < 0.5 )
         vals[i] = 1.0;
      else
      {
         vals[i] = -1.0;
         lhs -= 1.0;
      }
   }

   if( i == nbdchginfos )
   {
      std::list<LocalNodeInfoPtr> *conflictConsList = reinterpret_cast<std::list<LocalNodeInfoPtr> *>(SCIPconflicthdlrGetData(conflicthdlr));
      LocalNodeInfo *localNodeInfo = new LocalNodeInfo;
      localNodeInfo->linearRhs = SCIPinfinity(scip);
      localNodeInfo->nLinearCoefs = nbdchginfos;
      localNodeInfo->idxLinearCoefsVars = new int[nbdchginfos];
      localNodeInfo->linearCoefs = new double[nbdchginfos];
      for( i = 0; i < nbdchginfos; ++i )
      {
         SCIP_VAR *transformVar = vars[i];
         SCIP_Real scalar = vals[i];
         SCIP_Real constant = 0.0;
         if( SCIPvarGetOrigvarSum(&transformVar, &scalar, &constant ) ==  SCIP_INVALIDDATA )
            break;
         assert(transformVar != NULL);
         lhs -= constant;
         localNodeInfo->idxLinearCoefsVars[i] = SCIPvarGetIndex(transformVar);
         localNodeInfo->linearCoefs[i] = scalar;
      }
      if( i == nbdchginfos )
      {
         localNodeInfo->linearLhs = lhs;
         conflictConsList->push_back(localNodeInfo);
      }
      else
      {
         delete [] localNodeInfo->idxLinearCoefsVars;
         delete [] localNodeInfo->linearCoefs;
         delete localNodeInfo;
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

void
ScipParaSolver::setWinnerRacingParams(
      UG::ParaRacingRampUpParamSet *inRacingParams   /**< inRacingParams == 0 is winner solver */
      )
{


   SCIP_CALL_ABORT( SCIPresetParams(scip) );
   /* don't catch control+c */
   SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "misc/catchctrlc", FALSE) );

   if( paraParams->getBoolParamValue(UG::SetAllDefaultsAfterRacing) )
   {
      SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
      SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
      SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
      if( inRacingParams )
      {
         ScipParaRacingRampUpParamSet *scipRacingParams = dynamic_cast< ScipParaRacingRampUpParamSet * >(inRacingParams);
         SCIP_CALL_ABORT( SCIPsetIntParam(scip, "misc/permutationseed", scipRacingParams->getPermuteProbSeed()) );
      }
   }
   else
   {
      ScipParaRacingRampUpParamSet *scipRacingParams = dynamic_cast< ScipParaRacingRampUpParamSet * >(inRacingParams);
      setRacingParams(scipRacingParams, true);
   }
}

void
ScipParaSolver::setRacingParams(
      UG::ParaRacingRampUpParamSet *inRacingParams,
      bool winnerParam
      )
{
   ScipParaRacingRampUpParamSet *scipRacingParams = dynamic_cast< ScipParaRacingRampUpParamSet * >(inRacingParams);

   if( paraParams->getBoolParamValue(UG::ProvingRun) )
   {
      if( !winnerParam && scipRacingParams->getScipDiffParamSet() )
      {
         scipRacingParams->getScipDiffParamSet()->setParametersInScip(scip);
      }
   }
   else
   {
      int nHeuristics;
      if( paraParams->getBoolParamValue(UG::SetAllDefaultsAfterRacing ))
      {
         nHeuristics = scipRacingParams->getScipRacingParamSeed() % 2;
      }
      else
      {
         nHeuristics = scipRacingParams->getScipRacingParamSeed() % 4;
      }
      int nPresolving = (scipRacingParams->getScipRacingParamSeed()/4) % 4;
      int nSeparating = (scipRacingParams->getScipRacingParamSeed()/(4*4)) % 4;

      switch( nHeuristics )
      {
         case 0:
         {
            SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
            break;
         }
         case 1:
         {
            SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );
            break;
         }
         case 2:
         {
            SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_FAST, TRUE) );
            break;
         }
         case 3:
         {
            SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
            break;
         }
         default:
            THROW_LOGICAL_ERROR1("invalid nHeuristics");
      }

      switch( nPresolving )
      {
         case 0:
         {
            SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
            break;
         }
         case 1:
         {
            SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );
            break;
         }
         case 2:
         {
            SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );
            break;
         }
         case 3:
         {
            SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
            break;
         }
         default:
            THROW_LOGICAL_ERROR1("invalid nPresolving");
      }

      switch( nSeparating )
      {
         case 0:
         {
            SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
            break;
         }
         case 1:
         {
            if( paraParams->getBoolParamValue(UG::NoAggressiveSeparatorInRacing) )
            {
               SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
               SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );
            }
            else
            {
               SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );
            }
            break;
         }
         case 2:
         {
            SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_FAST, TRUE) );
            break;
         }
         case 3:
         {
            SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
            break;
         }
         default:
            THROW_LOGICAL_ERROR1("invalid nSeparating");
      }
   }

   assert(SCIPgetStage(scip) <= SCIP_STAGE_TRANSFORMED);
   SCIP_CALL_ABORT( SCIPsetIntParam(scip, "misc/permutationseed", scipRacingParams->getPermuteProbSeed()) );

   if( !winnerParam && scipRacingParams->getPermuteProbSeed() >= 64 )  // after all parameters tested, random branchig variable selection
   {
      SCIP_CALL_ABORT( SCIPsetIntParam(scip, "branching/random/maxdepth", 2) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scip, "branching/random/priority", 100000) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scip, "branching/random/seed", scipRacingParams->getGenerateBranchOrderSeed()) );
   }
#ifdef _DEBUG_DET
   writeSubproblem();
#endif
}

void
ScipParaSolver::createSubproblem(
      )
{
   UG::ParaNode *paraNode = getCurrentNode();
   assert(paraNode);
   double dualBoundValue = paraNode->getDualBoundValue();
   ScipParaDiffSubproblem *scipParaDiffSubproblem = dynamic_cast< ScipParaDiffSubproblem* >(paraNode->getDiffSubproblem());
   SCIP_VAR **orgVars = SCIPgetOrigVars(scip);  // variables are indexed by index
   int nOrg = SCIPgetNOrigVars(scip);       // the number of original variables
   if( scipParaDiffSubproblem )
   {
      for(int v = 0; v <  scipParaDiffSubproblem->getNBoundChanges(); v++)
      {
         if( scipParaDiffSubproblem->getIndex(v) < nOrg )
         {
            if( scipParaDiffSubproblem->getBoundType(v) == SCIP_BOUNDTYPE_LOWER )
            {
               SCIP_CALL_ABORT(
                     SCIPchgVarLbGlobal(
                           scip,
                           orgVars[scipParaDiffSubproblem->getIndex(v)],
                           scipParaDiffSubproblem->getBranchBound(v) )
                     );
               assert(SCIPisEQ(scip,SCIPvarGetLbGlobal(orgVars[scipParaDiffSubproblem->getIndex(v)]),scipParaDiffSubproblem->getBranchBound(v)));
            }
            else if (scipParaDiffSubproblem->getBoundType(v) == SCIP_BOUNDTYPE_UPPER)
            {
               SCIP_CALL_ABORT(SCIPchgVarUbGlobal(
                     scip,
                     orgVars[scipParaDiffSubproblem->getIndex(v)],
                     scipParaDiffSubproblem->getBranchBound(v) )
               );
               assert(SCIPisEQ(scip,SCIPvarGetUbGlobal(orgVars[scipParaDiffSubproblem->getIndex(v)]),scipParaDiffSubproblem->getBranchBound(v)));
            }
            else
            {
               THROW_LOGICAL_ERROR2("Invalid bound type: type = ", static_cast<int>(scipParaDiffSubproblem->getBoundType(v))) ;
            }
         }
         else
         {
            std::cout << "fixing branching variable index = " << scipParaDiffSubproblem->getIndex(v) << " is omitted!" << std::endl;
         }
      }

      if( scipParaDiffSubproblem->getNLinearConss() > 0 )
      {
         addedConss = new SCIP_CONS*[scipParaDiffSubproblem->getNLinearConss()];
      }

      for(int c = 0; c < scipParaDiffSubproblem->getNLinearConss() ; c++ )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int nVars = scipParaDiffSubproblem->getNLinearCoefs(c);

         /* create array of variables and coefficients */
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vars, nVars) );
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vals, nVars) );

         for( int v = 0; v < nVars; ++v )
         {
            vars[v] = orgVars[scipParaDiffSubproblem->getIdxLinearCoefsVars(c,v)];
            vals[v] = scipParaDiffSubproblem->getLinearCoefs(c,v);
         }

         SCIP_CONS* cons;
         char consname[SCIP_MAXSTRLEN];

         /* create a constraint */
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cli%d", c);
         SCIP_CALL_ABORT( SCIPcreateConsLinear(scip, &cons, consname, nVars, vars, vals,
               scipParaDiffSubproblem->getLinearLhs(c), scipParaDiffSubproblem->getLinearRhs(c),
               TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
         /** only a constraint whose "enforce is TRUE can be written in transformed problem */

         /* add constraint to SCIP */
         SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
         SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
         addedConss[c] = cons;
         /* free temporary memory */
         SCIPfreeBufferArray(scip, &vals);
         SCIPfreeBufferArray(scip, &vars);
      }
   }
   int addingConsParam = getParaParamSet()->getIntParamValue(UG::AddDualBoundCons);
   addedDualCons = 0;
   if( addingConsParam != 0 )
   {
      if( ( addingConsParam == 1 && !SCIPisGT(scip, paraNode->getDualBoundValue(), paraNode->getInitialDualBoundValue()) )
            || addingConsParam == 2 || addingConsParam == 3 )
      {
         SCIP_CONS* cons;
         int nvars = SCIPgetNVars(scip);
         SCIP_VAR **vars = SCIPgetVars(scip);
         SCIP_Real* vals = new SCIP_Real[nvars];
         for(int v = 0; v < nvars; ++v )
         {
            vals[v] = SCIPvarGetObj(vars[v]);
         }
         SCIP_CALL_ABORT( SCIPcreateConsLinear(scip, &cons, "objective", nvars, vars, vals, dualBoundValue, SCIPinfinity(scip),
               TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
         assert( SCIPisEQ( scip, SCIPgetTransObjscale(scip), 1.0 ) );
         assert( SCIPisZero( scip, SCIPgetTransObjoffset(scip) ) );

         /** try to automatically convert a linear constraint into a more specific and more specialized constraint */

         /* add constraint to SCIP */
         SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
         SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
         addedDualCons = cons;
      }
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM)
   {
      SCIP_CALL_ABORT( SCIPtransformProb(scip));
   }

   if( scipParaDiffSubproblem && scipParaDiffSubproblem->getNVarBranchStats() > 0 )
   {
      orgVars = SCIPgetOrigVars(scip);       /* original problem's variables              */
      for( int i = 0; i < scipParaDiffSubproblem->getNVarBranchStats(); i++ )
      {
         SCIP_CALL_ABORT( SCIPinitVarBranchStats(scip, orgVars[scipParaDiffSubproblem->getIdxLBranchStatsVars(i)],
               scipParaDiffSubproblem->getDownpscost(i),
               scipParaDiffSubproblem->getUppscost(i),
               scipParaDiffSubproblem->getDownvsids(i),
               scipParaDiffSubproblem->getUpvsids(i),
               scipParaDiffSubproblem->getDownconflen(i),
               scipParaDiffSubproblem->getUpconflen(i),
               scipParaDiffSubproblem->getDowninfer(i),
               scipParaDiffSubproblem->getUpinfer(i),
               scipParaDiffSubproblem->getDowncutoff(i),
               scipParaDiffSubproblem->getUpcutoff(i)
               )
         );
      }
   }
}

void
ScipParaSolver::freeSubproblem(
      )
{
   if( isRacingStage() && paraParams->getBoolParamValue(UG::RacingStatBranching) )
   {
	  DEF_SCIP_PARA_COMM( scipParaComm, paraComm);
      ScipParaInitialStat *initialStat = scipParaComm->createScipParaInitialStat(scip);
      initialStat->send(paraComm, 0);
      delete initialStat;
   }
   SCIP_CALL_ABORT( SCIPfreeTransform(scip) );

   if( currentNode->getDiffSubproblem() )
   {
      ScipParaDiffSubproblem *scipParaDiffSubproblem = dynamic_cast< ScipParaDiffSubproblem* >(currentNode->getDiffSubproblem());
      if( scipParaDiffSubproblem->getNLinearConss() > 0 )
      {
         for(int c = 0; c < scipParaDiffSubproblem->getNLinearConss() ; c++ )
         {
            SCIP_CALL_ABORT( SCIPdelCons(scip, addedConss[c]) );
         }
      }
   }
   if( addedDualCons )
   {
      SCIP_CALL_ABORT( SCIPdelCons(scip, addedDualCons) );
   }

   SCIP_VAR **orgVars = SCIPgetOrigVars(scip);  // variables are indexed by index
   int n = SCIPgetNOrigVars(scip);       // the number of original variables
   assert( n == nOrgVars );
   for( int v = 0; v < n; v++ )
   {
      SCIP_CALL_ABORT( SCIPchgVarLbGlobal( scip,orgVars[v], orgVarLbs[v] ) );
      SCIP_CALL_ABORT( SCIPchgVarUbGlobal( scip,orgVars[v], orgVarUbs[v] ) );
   }

}

void
ScipParaSolver::solve(
      )
{
   /** set instance specific parameters */
   // if( currentNode->isRootNode() && !(paraParams->getBoolParamValue(UseRootNodeCuts)) )
   // Probably, in order to avoid root node settings twice. Once root nodes is solved in LC
   // when UseRootNodeCuts is specified. However, for racing ramp-up, this is too bad.
   //  So, I changed the specification. The root node parameter settings is applied twice now
   SCIP_CALL_ABORT( SCIPresetParams(scip) ); 
   if( currentNode->isRootNode() )
   {
      scipDiffParamSetRoot->setParametersInScip(scip);
      if( paraParams->getBoolParamValue(UG::NoSolverPresolvingAtRoot) ||
            ( paraParams->getBoolParamValue(UG::NoSolverPresolvingAtRootDefaultSet) &&
                  isRacingStage() &&
                  paraComm->getRank() == 1 )    // rank 1 should be all default
            )
      {
         SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
      }
   }
   else
   {
      if( paraParams->getBoolParamValue(UG::NoSolverPresolvingAtRoot) )
      {
         SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
      }
      scipDiffParamSet->setParametersInScip(scip);
   }

   /** set cutoff value */
   SCIP_CALL_ABORT( SCIPsetObjlimit(scip, globalBestIncumbentValue) );
   /** solve */
   if( paraParams->getBoolParamDefaultValue(UG::TransferConflicts) )
   {
      assert(conflictConsList->size() == 0);
   }
   SCIP_RETCODE ret = SCIPsolve(scip);
   if( ret != SCIP_OKAY )
   {
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
      SCIPprintError(ret, NULL);
#else
      SCIPprintError(ret);
#endif
      THROW_LOGICAL_ERROR1("SCIP terminate with NO SCIP_OKAY");
   }

   // Notification message has to complete
   if( notificationProcessed )
   {
      waitNotificationIdMessage();
   }
   // Then, solver status should be checked
   SCIP_STATUS status = SCIPgetStatus(scip);
   if( status == SCIP_STATUS_OPTIMAL )   // when sub-MIP is solved at root node, the solution may not be saved
   {
      SCIP_SOL *sol = SCIPgetBestSol(scip);
      int nVars = SCIPgetNOrigVars(scip);
      SCIP_VAR **vars = SCIPgetOrigVars(scip);
      SCIP_Real *vals = new SCIP_Real[nVars];
      SCIP_CALL_ABORT( SCIPgetSolVals(scip, sol, nVars, vars, vals) );
      DEF_SCIP_PARA_COMM( scipParaComm, paraComm);
      saveIfImprovedSolutionWasFound(
               scipParaComm->createScipParaSolution(
                     SCIPgetSolOrigObj(scip, sol),
                     nVars,
                     vars,
                     vals
                     )
      );
      delete [] vals;
      sendLocalSolution();
   }
   else
   {
      if( status == SCIP_STATUS_MEMLIMIT  )
      {
         std::cout << "Warning: SCIP was interrupted because the memory limit was reached" << std::endl;
      }
   }
   if( conflictConsList && conflictConsList->size() > 0 )
   {
      int nConfilcts = conflictConsList->size();
      for(int i = 0; i < nConfilcts; i++ )
      {
         assert(!conflictConsList->empty());
         LocalNodeInfo *info= conflictConsList->front();
         conflictConsList->pop_front();
         if( info->linearCoefs ) delete[] info->linearCoefs;
         if( info->idxLinearCoefsVars ) delete[] info->idxLinearCoefsVars;
         delete info;
      }
   }
}

long long
ScipParaSolver::getNNodesSolved(
      )
{
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVED  )
   {
      return SCIPgetNTotalNodes(scip);
   }
   else
   {
      return 0;
   }
}

int
ScipParaSolver::getNNodesLeft(
      )
{
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVED  )
   {
      return SCIPgetNNodesLeft(scip);
   }
   else
   {
      if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVING && SCIPgetStage(scip) <= SCIP_STAGE_INITSOLVE )
      {
         return 1;
      }
      else
      {
         return 0;
      }
   }
}

double
ScipParaSolver::getDualBoundValue(
      )
{
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_INITSOLVE )
   {
      return currentNode->getDualBoundValue();
   }
   else
   {
      return SCIPgetDualbound(scip);
   }
}


ScipParaSolver::ScipParaSolver(
      int argc,
      char **argv,
      UG::ParaComm     *comm,
      UG::ParaParamSet *inParaParamSet,
      UG::ParaInstance *inParaInstance,
      UG::ParaDeterministicTimer *inDetTimer
      ) : ParaSolver(argc, argv, comm, inParaParamSet, inParaInstance, inDetTimer),
      messagehdlr(0), logfile(0), originalParamSet(0), conflictConsList(0), userPlugins(0),
      nOrgVars(0), orgVarLbs(0), orgVarUbs(0), addedConss(0), addedDualCons(0)
{
   char* logname = NULL;
   SCIP_Bool quiet = false;

   ScipParaInstance *scipParaInstance = dynamic_cast< ScipParaInstance *>(paraInstance);

   /* Initialize the SCIP environment */
   /*********
    * Setup *
    *********/
   /* initialize SCIP */
   SCIP_CALL_ABORT( SCIPcreate(&scip) );
   /* include default SCIP plugins */
   SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scip) );
   /** user include plugins */
   includeUserPlugins(scip);
   /* include communication point handler */
   SCIP_CALL_ABORT( SCIPincludeObjEventhdlr(scip, new ScipParaObjCommPointHdlr(paraComm, this), TRUE) );
   /* include branch rule plugins */
   SCIP_CALL_ABORT( SCIPincludeObjBranchrule(scip, new ScipParaObjBranchRule(this), TRUE) );

   if( inParaParamSet->getBoolParamValue(UG::TransferConflicts) )
   {
      conflictConsList = new std::list<LocalNodeInfoPtr>;
      SCIP_CONFLICTHDLRDATA *conflictHdrData = reinterpret_cast< SCIP_CONFLICTHDLRDATA * >(conflictConsList);
      /* create conflict handler to collects conflicts */
#if SCIP_VERSION == 211 && SCIP_SUBVERSION == 0
      SCIP_CALL_ABORT( SCIPincludeConflicthdlr(scip, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
            NULL, NULL, NULL, NULL, NULL, NULL, conflictExecCollector, conflictHdrData) );
#else
      SCIP_CALL_ABORT( SCIPincludeConflicthdlrBasic(scip, NULL, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
            conflictExecCollector, conflictHdrData) );
#endif
   }

   /********************
    * Parse parameters *
    ********************/
   for( int i = 3; i < argc; ++i )   /** the first argument is runtime parameter file for ParaSCIP */
   {
      if( strcmp(argv[i], "-l") == 0 )
      {
         i++;
         if( i < argc )
            logname = argv[i];
         else
         {
            THROW_LOGICAL_ERROR1("missing log filename after parameter '-l'");
         }
      }
      else if( strcmp(argv[i], "-q") == 0 )
         quiet = true;
      // other arguments are omitted in Solver
   }

   /***********************************
    * create log file message handler *
    ***********************************/
   if( paraParams->getBoolParamValue(UG::Quiet) )
   {
      // SCIP_CALL_ABORT( SCIPsetMessagehdlr(NULL) );
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
      SCIP_CALL_ABORT( SCIPcreateObjMessagehdlr(&messagehdlr, new ScipParaObjMessageHdlr(paraComm, NULL, TRUE, FALSE), TRUE) );
      SCIP_CALL_ABORT( SCIPsetMessagehdlr(messagehdlr) );
#else
      SCIPsetMessagehdlrQuiet(scip, TRUE);
#endif
   }
   else
   {
      if( logname != NULL || quiet  )
      {
         if( logname != NULL )
         {
            std::ostringstream os;
            os << logname << comm->getRank();
            logfile = fopen(os.str().c_str(), "a"); // append to log file */
            if( logfile == NULL )
            {
               THROW_LOGICAL_ERROR3("cannot open log file <", logname, "> for writing");
            }
         }
         SCIP_CALL_ABORT( SCIPcreateObjMessagehdlr(&messagehdlr, new ScipParaObjMessageHdlr(paraComm, logfile, quiet, FALSE), TRUE) );
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
         SCIP_CALL_ABORT( SCIPsetMessagehdlr(messagehdlr) );
#else
         SCIP_CALL_ABORT( SCIPsetMessagehdlr(scip, messagehdlr) );
#endif
      }
   }

   DEF_SCIP_PARA_COMM( scipParaComm, paraComm );
   scipDiffParamSetRoot = scipParaComm->createScipDiffParamSet();
   scipDiffParamSetRoot->bcast(comm, 0);    /** this bcast is sent as SolverInitializationMessage */
   scipDiffParamSet = scipParaComm->createScipDiffParamSet();
   scipDiffParamSet->bcast(comm, 0);    /** this bcast is sent as SolverInitializationMessage */
   int tempIsWarmStarted;
   comm->bcast(&tempIsWarmStarted, 1, UG::ParaINT, 0);
   warmStarted = (tempIsWarmStarted == 1);
   comm->bcast(&globalBestIncumbentValue, 1, UG::ParaDOUBLE, 0);
   if( paraParams->getBoolParamValue( UG::DistributeBestPrimalSolution ) )
   {
      int solutionExists = 0;
      comm->bcast(&solutionExists, 1, UG::ParaINT, 0);
      if( solutionExists )
      {
         globalBestIncumbentSolution = paraComm->createParaSolution();
         globalBestIncumbentSolution->bcast(comm, 0);
      }
   }

   /** set parameters for SCIP: this values are reseted before solving */
   scipDiffParamSet->setParametersInScip(scip);
   SCIP_Real epsilon;
   SCIP_CALL_ABORT( SCIPgetRealParam(scip, "numerics/epsilon", &epsilon));
   eps = epsilon;

   char *settingsNameLC = 0;
   for( int i = 3; i < argc; ++i )   /** the first argument is runtime parameter file for ParaSCIP */
                                     /** the second argument is problem file name */
   {
      if( strcmp(argv[i], "-sl") == 0 )
      {
         i++;
         if( i < argc )
         {
            settingsNameLC = argv[i];
            break;
         }
         else
         {
            std::cerr << "missing settings filename after parameter '-sl'" << std::endl;
            exit(1);
         }
      }
   }

   /** create problem */
   scipParaInstance->createProblem(scip, 
                                   paraParams->getIntParamValue(UG::InstanceTransferMethod),
                                   paraParams->getBoolParamValue(UG::NoPreprocessingInLC),
                                   settingsNameLC
                                  );
   saveOrgProblemBounds();
   if( paraParams->getBoolParamValue(UG::CheckEffectOfRootNodePreprocesses) )
   {
      /* initialize SCIP to check root solvability */
      SCIP_CALL_ABORT( SCIPcreate(&scipToCheckEffectOfRootNodeProcesses) );
      /* include default SCIP plugins */
      SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scipToCheckEffectOfRootNodeProcesses) );
      /* include scipParaConshdlr plugins */
      scipParaInstance->createProblem(scipToCheckEffectOfRootNodeProcesses, 
                                      paraParams->getIntParamValue(UG::InstanceTransferMethod),
                                      paraParams->getBoolParamValue(UG::NoPreprocessingInLC),
                                      settingsNameLC
                                     );
      scipDiffParamSet->setParametersInScip(scipToCheckEffectOfRootNodeProcesses);
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "presolving/maxrestarts", 0) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "presolving/maxrounds", 0) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scipToCheckEffectOfRootNodeProcesses, "constraints/linear/presolpairwise", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scipToCheckEffectOfRootNodeProcesses, "constraints/and/presolpairwise", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scipToCheckEffectOfRootNodeProcesses, "constraints/logicor/presolpairwise", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scipToCheckEffectOfRootNodeProcesses, "constraints/setppc/presolpairwise", FALSE) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "propagating/probing/maxprerounds", 0) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "heuristics/feaspump/freq", -1) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "heuristics/rens/freq", -1) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "separating/maxcutsroot", 100) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "separating/maxroundsroot", 5) );
   }
   delete paraInstance;
   paraInstance = 0;
}

ScipParaSolver::ScipParaSolver(
      int argc,
      char **argv,
      UG::ParaComm     *comm,
      UG::ParaParamSet *inParaParamSet,
      UG::ParaInstance *inParaInstance,
      UG::ParaDeterministicTimer *inDetTimer,
      bool thread
      ) : ParaSolver(argc, argv, comm, inParaParamSet, inParaInstance, inDetTimer),
      messagehdlr(0), logfile(0), originalParamSet(0), conflictConsList(0), userPlugins(0),
      nOrgVars(0), orgVarLbs(0), orgVarUbs(0), addedConss(0), addedDualCons(0)
{
   assert(thread);  // This is a constructor for threads parallel version

   char* logname = NULL;
   SCIP_Bool quiet = false;

   ScipParaInstance *scipParaInstance = dynamic_cast<ScipParaInstance *>(paraInstance);

   /* Initialize the SCIP environment */
   /*********
    * Setup *
    *********/
   if( paraParams->getIntParamValue(UG::InstanceTransferMethod) == 0 )
   {
      /* copy SCIP environment */
      scip = scipParaInstance->getScip();
      SCIP_CALL_ABORT( SCIPresetParams(scip) );  // if LC parameter settings are applied,
                                                 // it is necessary to reset them
   }
   else
   {
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
      SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scip) );
   }

   /* don't catch control+c */
   SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "misc/catchctrlc", FALSE) );

   /* include communication point handler */
   SCIP_CALL_ABORT( SCIPincludeObjEventhdlr(scip, new ScipParaObjCommPointHdlr(paraComm, this), TRUE) );
   /* include branch rule plugins */
   SCIP_CALL_ABORT( SCIPincludeObjBranchrule(scip, new ScipParaObjBranchRule(this), TRUE) );

   if( inParaParamSet->getBoolParamValue(UG::TransferConflicts) )
   {
      conflictConsList = new std::list<LocalNodeInfoPtr>;
      SCIP_CONFLICTHDLRDATA *conflictHdrData = reinterpret_cast< SCIP_CONFLICTHDLRDATA * >(conflictConsList);
      /* create conflict handler to collects conflicts */
#if SCIP_VERSION == 211 && SCIP_SUBVERSION == 0
      SCIP_CALL_ABORT( SCIPincludeConflicthdlr(scip, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
            NULL, NULL, NULL, NULL, NULL, NULL, conflictExecCollector, conflictHdrData) );
#else
      SCIP_CALL_ABORT( SCIPincludeConflicthdlrBasic(scip, NULL, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
            conflictExecCollector, conflictHdrData) );
#endif
   }

   /********************
    * Parse parameters *
    ********************/
   for( int i = 3; i < argc; ++i )   /** the first argument is runtime parameter file for ParaSCIP */
   {
      if( strcmp(argv[i], "-l") == 0 )
      {
         i++;
         if( i < argc )
            logname = argv[i];
         else
         {
            THROW_LOGICAL_ERROR1("missing log filename after parameter '-l'");
         }
      }
      else if( strcmp(argv[i], "-q") == 0 )
         quiet = true;
      // other arguments are omitted in Solver
   }

   /***********************************
    * create log file message handler *
    ***********************************/
   if( paraParams->getBoolParamValue(UG::Quiet) )
   {
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
      SCIP_CALL_ABORT( SCIPcreateObjMessagehdlr(&messagehdlr, new ScipParaObjMessageHdlr(paraComm, NULL, TRUE, FALSE), TRUE) );
      SCIP_CALL_ABORT( SCIPsetMessagehdlr(messagehdlr) );
#else
      SCIPsetMessagehdlrQuiet(scip, TRUE );
#endif
   }
   else
   {
      if( logname != NULL || quiet )
      {
         paraComm->lockApp();    // if solver runs as thread, this lock is necessary
         if( logname != NULL )
         {
            std::ostringstream os;
            os << logname << comm->getRank();
            logfile = fopen(os.str().c_str(), "a"); // append to log file */
            if( logfile == NULL )
            {
               THROW_LOGICAL_ERROR3("cannot open log file <", logname, "> for writing");
            }
         }
         SCIP_CALL_ABORT( SCIPcreateObjMessagehdlr(&messagehdlr, new ScipParaObjMessageHdlr(paraComm, logfile, quiet, TRUE), TRUE) );
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
         SCIP_CALL_ABORT( SCIPsetMessagehdlr(messagehdlr) );
#else
         SCIP_CALL_ABORT( SCIPsetMessagehdlr(scip, messagehdlr) );
#endif
         paraComm->unlockApp();
      }
   }

   DEF_SCIP_PARA_COMM( scipParaComm, paraComm );
   scipDiffParamSetRoot = scipParaComm->createScipDiffParamSet();
   scipDiffParamSetRoot->bcast(comm, 0);    /** this bcast is sent as SolverInitializationMessage */
   scipDiffParamSet = scipParaComm->createScipDiffParamSet();
   scipDiffParamSet->bcast(comm, 0);    /** this bcast is sent as SolverInitializationMessage */
   int tempIsWarmStarted;
   comm->bcast(&tempIsWarmStarted, 1, UG::ParaINT, 0);
   warmStarted = (tempIsWarmStarted == 1);
   comm->bcast(&globalBestIncumbentValue, 1, UG::ParaDOUBLE, 0);

   if( paraParams->getBoolParamValue(UG::NoUpperBoundTransferInRacing) )
   {
      int solutionExists = 0;
      paraComm->bcast(&solutionExists, 1, UG::ParaINT, 0);
   }
   else
   {
      if( paraParams->getBoolParamValue( UG::DistributeBestPrimalSolution ) )
      {
         int solutionExists = 0;
         comm->bcast(&solutionExists, 1, UG::ParaINT, 0);
         if( solutionExists )
         {
            globalBestIncumbentSolution = paraComm->createParaSolution();
            globalBestIncumbentSolution->bcast(comm, 0);
         }
      }
   }

   /** set parameters for SCIP: this values are reseted before solving */
   scipDiffParamSet->setParametersInScip(scip);
   SCIP_Real epsilon;
   SCIP_CALL_ABORT( SCIPgetRealParam(scip, "numerics/epsilon", &epsilon));
   eps = epsilon;

   char *settingsNameLC = 0;
   for( int i = 3; i < argc; ++i )   /** the first argument is runtime parameter file for ParaSCIP */
                                     /** the second argument is problem file name */
   {
      if( strcmp(argv[i], "-sl") == 0 )
      {
         i++;
         if( i < argc )
         {
            settingsNameLC = argv[i];
            break;
         }
         else
         {
            std::cerr << "missing settings filename after parameter '-sl'" << std::endl;
            exit(1);
         }
      }
   }

   /** create problem */
   scipParaInstance->createProblem(scip, 
                                   paraParams->getIntParamValue(UG::InstanceTransferMethod),
                                   paraParams->getBoolParamValue(UG::NoPreprocessingInLC),
                                   settingsNameLC
                                  );

   saveOrgProblemBounds();
   if( paraParams->getBoolParamValue(UG::CheckEffectOfRootNodePreprocesses) )
   {
      /* initialize SCIP to check root solvability */
      SCIP_CALL_ABORT( SCIPcreate(&scipToCheckEffectOfRootNodeProcesses) );
      /* include default SCIP plugins */
      SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scipToCheckEffectOfRootNodeProcesses) );
      /* include scipParaConshdlr plugins */
      scipParaInstance->createProblem(scipToCheckEffectOfRootNodeProcesses, 
                                      paraParams->getIntParamValue(UG::InstanceTransferMethod),
                                      paraParams->getBoolParamValue(UG::NoPreprocessingInLC),
                                      settingsNameLC
                                     );
      scipDiffParamSet->setParametersInScip(scipToCheckEffectOfRootNodeProcesses);
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "presolving/maxrestarts", 0) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "presolving/maxrounds", 0) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scipToCheckEffectOfRootNodeProcesses, "constraints/linear/presolpairwise", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scipToCheckEffectOfRootNodeProcesses, "constraints/and/presolpairwise", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scipToCheckEffectOfRootNodeProcesses, "constraints/logicor/presolpairwise", FALSE) );
      SCIP_CALL_ABORT( SCIPsetBoolParam(scipToCheckEffectOfRootNodeProcesses, "constraints/setppc/presolpairwise", FALSE) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "propagating/probing/maxprerounds", 0) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "heuristics/feaspump/freq", -1) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "heuristics/rens/freq", -1) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "separating/maxcutsroot", 100) );
      SCIP_CALL_ABORT( SCIPsetIntParam(scipToCheckEffectOfRootNodeProcesses, "separating/maxroundsroot", 5) );
   }
   delete paraInstance;
   paraInstance = 0;

}

ScipParaSolver::~ScipParaSolver(
      )
{
   /** delete scip diff param sets */
   if( scipDiffParamSetRoot ) delete scipDiffParamSetRoot;
   if( scipDiffParamSet ) delete scipDiffParamSet;
   if( userPlugins ) delete userPlugins;

   /** reset message handler */
   // message handler is mangaed within scip. It is freed at SCIPfree
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
   if( messagehdlr )
   {
      SCIP_CALL_ABORT( SCIPsetDefaultMessagehdlr() );
      SCIP_CALL_ABORT( SCIPfreeObjMessagehdlr(&messagehdlr) );
   }
#endif

   /* free SCIP environment */
   if( paraParams->getBoolParamValue(UG::CheckEffectOfRootNodePreprocesses) )
   {
      SCIP_CALL_ABORT( SCIPfree(&scipToCheckEffectOfRootNodeProcesses) );
   }
   SCIP_CALL_ABORT( SCIPfree(&scip) );

   /** close log file */
   if( logfile ) fclose(logfile);

   if( conflictConsList && conflictConsList->size() > 0 )
   {
      int nConfilcts = conflictConsList->size();
      for(int i = 0; i < nConfilcts; i++ )
      {
         assert(!conflictConsList->empty());
         LocalNodeInfo *info= conflictConsList->front();
         conflictConsList->pop_front();
         if( info->linearCoefs ) delete[] info->linearCoefs;
         if( info->idxLinearCoefsVars ) delete[] info->idxLinearCoefsVars;
         delete info;
      }
   }

   if( orgVarLbs ) delete [] orgVarLbs;
   if( orgVarUbs ) delete [] orgVarUbs;
}

void
ScipParaSolver::writeCurrentNodeProblem(
      const std::string& filename
      )
{
   FILE *file = fopen(filename.c_str(),"a");
   if( !file )
   {
      std::cout << "file : " << filename << "cannot open." << std::endl;
      abort();
   }
   SCIP_CALL_ABORT( SCIPprintTransProblem(scip, file, "lp", FALSE) );
}

void
ScipParaSolver::solveToCheckEffectOfRootNodePreprocesses(
      )
{
   SCIP *backupScip = scip;
   scip = scipToCheckEffectOfRootNodeProcesses;
   createSubproblem();
   /** set cutoff value */
   SCIP_CALL_ABORT( SCIPsetObjlimit(scip, globalBestIncumbentValue) );
   /** solve */
   SCIP_CALL_ABORT( SCIPsolve(scip) );
   nSolvedWithNoPreprocesses = SCIPgetNNodes(scip);
   freeSubproblem();
   scip = backupScip;
}

void
ScipParaSolver::tryNewSolution(
      UG::ParaSolution *sol
      )
{

   if( SCIPgetStage(scip) <= SCIP_STAGE_TRANSFORMING || SCIPgetStage(scip) >= SCIP_STAGE_SOLVED  )
   {
      THROW_LOGICAL_ERROR1("invalid tyrNewSolution");
   }

   ScipParaSolution *tempSol = dynamic_cast< ScipParaSolution * >(sol);
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   if( SCIPcreateOrigSol(scip, &newsol, 0) != SCIP_OKAY )  // SCIP_CALL_ABORT ????
   {
      return;
   }

   SCIP_VAR** vars = SCIPgetOrigVars(scip);
   int i;
   for(i = 0; i < tempSol->getNVars(); i++ )
   {
      /* skip variables which are definitely not in our scip */
      if( tempSol->indexAmongSolvers(i) >= SCIPgetNOrigVars(scip) ) continue;
      if( tempSol->indexAmongSolvers(i) > ( tempSol->getNVars() - 1 ) ) break;
      SCIP_CALL_ABORT( SCIPsetSolVal(scip, newsol, vars[tempSol->indexAmongSolvers(i)], tempSol->getValues()[i]) );
   }
   if( i != tempSol->getNVars() )
   {
      /** the given solution should be generated in original space,
       * therefore the solution values cannot use for ParaSCIP
       */
      SCIP_CALL_ABORT( SCIPfreeSol(scip, &newsol) );
      // delete tempSol;  // this case, DO NOT DELETE tempSol.
      return;
   }

#if SCIP_VERSION == 211 && SCIP_SUBVERSION == 0
   if( SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMED || SCIPgetStage(scip) == SCIP_STAGE_PRESOLVED  )
#else
   if( SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMED ||
		   SCIPgetStage(scip) == SCIP_STAGE_PRESOLVED ||
		   SCIPgetStage(scip) == SCIP_STAGE_INITPRESOLVE )
#endif
   {
      SCIP_Bool success;
      SCIP_CALL_ABORT( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, &success) );
   }
   else
   {
      SCIP_CALL_ABORT( SCIPheurPassSolTrySol(scip, SCIPfindHeur(scip, "trysol"), newsol) );
      SCIP_CALL_ABORT( SCIPfreeSol(scip, &newsol) );
   }
}

void
ScipParaSolver::setLightWeightRootNodeProcess(
      )
{
   lightWeightRootNodeComputation = true;
   if( !originalParamSet )
   {
	  DEF_SCIP_PARA_COMM( scipParaComm, paraComm );
      originalParamSet = scipParaComm->createScipDiffParamSet(scip);;
   }
   SCIP_CALL_ABORT( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_FAST, TRUE) );
   SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );
   SCIP_CALL_ABORT( SCIPsetSeparating(scip, SCIP_PARAMSETTING_FAST, TRUE) );
}

void
ScipParaSolver::setOriginalRootNodeProcess(
      )
{
   assert(originalParamSet);
   originalParamSet->setParametersInScip(scip);
   lightWeightRootNodeComputation = false;
}

void
ScipParaSolver::writeSubproblem(
      )
{
   std::string subcipprefix("SolverLp");
   std::string subcipfilename;
   std::ostringstream oss;
   oss << subcipprefix << paraComm->getRank();
   subcipfilename = oss.str();
   subcipfilename += ".lp";
   SCIP_CALL_ABORT( SCIPwriteOrigProblem(scip, subcipfilename.c_str(), "lp", FALSE) );
   subcipfilename += "t";
   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL_ABORT( SCIPwriteTransProblem(scip, subcipfilename.c_str(), "lp", FALSE) );
   }
   char name[SCIP_MAXSTRLEN];
   (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "SolverLp%d.set", paraComm->getRank());
   SCIP_CALL_ABORT( SCIPwriteParams(scip, name, TRUE, FALSE) );
}

void
ScipParaSolver::saveOrgProblemBounds(
      )
{
   assert(paraInstance);
   ScipParaInstance *scipParaInstance = dynamic_cast<ScipParaInstance*>(paraInstance);
   nOrgVars = paraInstance->getNVars();
   orgVarLbs = new SCIP_Real[nOrgVars];
   orgVarUbs = new SCIP_Real[nOrgVars];
   for( int v = 0; v < nOrgVars; v++ )
   {
      orgVarLbs[v] = scipParaInstance->getVarLb(v);
      orgVarUbs[v] = scipParaInstance->getVarUb(v);
   }
}

