/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_orig.c
 * @brief  branching rule for original problem in GCG
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_orig.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "cons_origbranch.h"
#include "branch_relpsprob.h"
#include "scip/cons_linear.h"
#include "type_branchgcg.h"
#include "pub_gcgvar.h"

#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_restartdfs.h"
#include "scip/branch_allfullstrong.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"


#define BRANCHRULE_NAME          "orig"
#define BRANCHRULE_DESC          "branching for the original program in generic column generation"
#define BRANCHRULE_PRIORITY      100
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_ENFORCEBYCONS FALSE
#define DEFAULT_MOSTFRAC      FALSE
#define DEFAULT_USEPSEUDO     TRUE
#define DEFAULT_USEPSSTRONG   FALSE

/** branching data for branching decisions */
struct GCG_BranchData
{
   SCIP_VAR*             origvar;            /**< original variable on which the branching is done */
   SCIP_BOUNDTYPE        boundtype;          /**< type of the new bound of original variable */
   SCIP_Real             newbound;           /**< new lower/upper bound of the original variable */
   SCIP_Real             oldbound;           /**< old lower/upper bound of the pricing variable */
   SCIP_Real             oldvalue;           /**< old value of the original variable */
   SCIP_Real             olddualbound;       /**< dual bound before the branching was performed */
   SCIP_CONS*            cons;               /**< constraint that enforces the branching restriction in the original
                                              *   problem, or NULL if this is done by variable bounds */
};


/** branches on a given variable */
static
SCIP_RETCODE branchVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer of the orig branching rule */
   SCIP_VAR*             branchvar,          /**< variable to branch on */
   SCIP_Real             solval              /**< value of the variable in the current solution */
   )
{
   /* data for b&b child creation */
   SCIP_NODE* childup;
   SCIP_NODE* childdown;
   SCIP_CONS* origbranchup;
   SCIP_CONS* origbranchdown;
   GCG_BRANCHDATA* branchupdata;
   GCG_BRANCHDATA* branchdowndata;
   char upname[SCIP_MAXSTRLEN];
   char downname[SCIP_MAXSTRLEN];

   /* parameter data */
   SCIP_Bool enforcebycons;

   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(branchvar != NULL);


   /* get values of parameters */
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/enforcebycons", &enforcebycons) );

   SCIPdebugMessage("Branching on var %s with value %g in current solution\n", SCIPvarGetName(branchvar), solval);

   /* create the b&b-tree child-nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childup, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdown, 0.0, SCIPgetLocalTransEstimate(scip)) );

   /** @todo use block memory here */
   /* create the branch data for the childs and assign the values */
   SCIP_CALL( SCIPallocMemory(scip, &(branchupdata)) );
   SCIP_CALL( SCIPallocMemory(scip, &(branchdowndata)) );

   branchupdata->origvar = branchvar;
   branchupdata->oldvalue = solval;
   branchupdata->olddualbound = SCIPgetLocalLowerbound(GCGrelaxGetMasterprob(scip));
   branchupdata->boundtype = SCIP_BOUNDTYPE_LOWER;
   branchupdata->newbound = SCIPceil(scip, solval);
   branchupdata->oldbound = SCIPvarGetLbLocal(branchvar);

   branchdowndata->origvar = branchvar;
   branchdowndata->oldvalue = solval;
   branchdowndata->olddualbound = SCIPgetLocalLowerbound(GCGrelaxGetMasterprob(scip));
   branchdowndata->boundtype = SCIP_BOUNDTYPE_UPPER;
   branchdowndata->newbound = SCIPfloor(scip, solval);
   branchdowndata->oldbound = SCIPvarGetUbLocal(branchvar);


   (void) SCIPsnprintf(upname, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchupdata->origvar),
      ">=", branchupdata->newbound);
   (void) SCIPsnprintf(downname, SCIP_MAXSTRLEN, "%s %s %f", SCIPvarGetName(branchdowndata->origvar),
      "<=", branchdowndata->newbound);

   /* enforce branching decision by a constraint rather than by bound changes */
   if( enforcebycons )
   {
      /* enforce new bounds by linear constraints */
      SCIP_CONS* consup;
      SCIP_CONS* consdown;

      /* create corresponding constraints */
      SCIP_CALL( SCIPcreateConsLinear(scip, &consup, upname, 0, NULL, NULL,
            SCIPceil(scip, solval), SCIPinfinity(scip),
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcreateConsLinear(scip, &consdown, downname, 0, NULL, NULL,
            -1.0 * SCIPinfinity(scip), SCIPfloor(scip, solval),
            TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPaddCoefLinear(scip, consup, branchvar, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, consdown, branchvar, 1.0) );

      /* add constraints to nodes */
      SCIP_CALL( SCIPaddConsNode(scip, childup, consup, NULL) );
      SCIP_CALL( SCIPaddConsNode(scip, childdown, consdown, NULL) );

      branchupdata->cons = consup;
      branchdowndata->cons = consdown;
   }
   else
   {
      /* enforce new bounds by setting variable bounds */
      SCIP_CALL( SCIPchgVarUbNode(scip, childdown, branchvar, solval) );
      SCIP_CALL( SCIPchgVarLbNode(scip, childup, branchvar, solval) );

      branchupdata->cons = NULL;
      branchdowndata->cons = NULL;
   }

   /* create the origbranch constraints */
   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchup, upname, childup,
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchupdata) );
   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchdown, downname, childdown,
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchdowndata) );

   /* add constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childup, origbranchup, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdown, origbranchdown, NULL) );

   /* store bound change of variables that were directly transferred to the master problem */
   if( !enforcebycons && GCGvarGetBlock(branchvar) == -1 )
   {
      SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, origbranchup, branchdowndata->origvar,
            branchupdata->boundtype, branchupdata->newbound) );
      SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, origbranchdown, branchdowndata->origvar,
            branchdowndata->boundtype, branchdowndata->newbound) );
   }

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchup) );
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchdown) );

   return SCIP_OKAY;
}


/** branching method for relaxation solutions */
static
SCIP_RETCODE branchExtern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< pointer of the orig branching rule */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching call */
   )
{
   int i;

   /* parameter data */
   SCIP_Bool mostfrac;
   SCIP_Bool usepseudocosts;
   SCIP_Bool usepsstrong;

   /* branching candidates */
   SCIP_VAR** branchcands;
   SCIP_Real* branchcandsscore;
   SCIP_Real* branchcandssol;
   int nbranchcands;
   int npriobranchcands;

   /* values for choosing the variable to branch on */
   SCIP_VAR* branchvar;
   SCIP_Real solval;
   SCIP_Real maxfrac;
   SCIP_Real frac;
   SCIP_Real maxpsscore;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPisRelaxSolValid(scip));

   SCIPdebugMessage("Execrel method of orig branching\n");

   *result = SCIP_DIDNOTRUN;

   /* get values of parameters */
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/mostfrac", &mostfrac) );
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/usepseudocosts", &usepseudocosts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "branching/orig/usepsstrong", &usepsstrong) );

   /* get the branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &branchcands, &branchcandssol, &branchcandsscore, &nbranchcands,
         &npriobranchcands, NULL, NULL, NULL) );

   branchvar = NULL;
   solval = 0.0;

   maxfrac = 0.0;
   maxpsscore = -1.0;

   if( usepsstrong )
   {
      SCIP_CALL( SCIPgetRelpsprobBranchVar(scip, branchcands, branchcandssol, branchcandsscore, npriobranchcands,
            npriobranchcands, result, &branchvar) );
      assert(branchvar != NULL || *result == SCIP_CUTOFF);
      assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF);

      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;

      solval = SCIPgetRelaxSolVal(scip, branchvar);
   }

   /* branch on an integer variable belonging to a unique block with fractional value */
   if( branchvar == NULL )
      for( i = 0; i < npriobranchcands; i++ )
      {
         assert(GCGvarIsOriginal(branchcands[i]));

         /* variable belongs to no block */
         if( GCGvarGetBlock(branchcands[i]) == -1 )
            continue;

         /* block is not unique (non-linking variables) */
         if( !GCGvarIsLinking(branchcands[i]) && GCGrelaxGetNIdenticalBlocks(scip, GCGvarGetBlock(branchcands[i])) != 1 )
            continue;

         /* check that blocks of linking variable are unique */
         if( GCGvarIsLinking(branchcands[i]) )
         {
            int nvarblocks;
            int* varblocks;
            SCIP_Bool unique;
            int j;

            nvarblocks = GCGlinkingVarGetNBlocks(branchcands[i]);
            SCIP_CALL( SCIPallocBufferArray(scip, &varblocks, nvarblocks) );
            SCIP_CALL( GCGlinkingVarGetBlocks(branchcands[i], nvarblocks, varblocks) );

            unique = TRUE;
            for( j = 0; j < nvarblocks; ++j )
               if( GCGrelaxGetNIdenticalBlocks(scip, varblocks[j]) != 1 )
                  unique = FALSE;

            SCIPfreeBufferArray(scip, &varblocks);

            if( !unique )
               continue;
         }

         /* use pseudocost variable selection rule */
         if( usepseudocosts )
         {
            /* select the variable, if its pseudocost are higher than the ones of the currently saved variable */
            if( SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]) > maxpsscore )
            {
               branchvar = branchcands[i];
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               maxpsscore = SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]);
            }
         }
         /* use most fractional variable selection rule */
         else
         {
            /* compute the fractionality */
            frac = MIN( branchcandsscore[i], 1.0 - branchcandsscore[i] );
            assert(frac > 0);

            /* fractionality is higher than that of the current highest fractionality */
            if( frac >= maxfrac )
            {
               SCIPdebugMessage("Var %s has fractional value in current solution: %f\n", SCIPvarGetName(branchcands[i]), branchcandssol[i]);
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               branchvar = branchcands[i];
               /* if we do not look for the most fractional variable, but for the first fractional variable,
                * we can stop here since we found a variable to branch on */
               if( !mostfrac )
                  break;
            }
         }
      }

   /* we did not find a variable to branch on so far, so we look for an integer variable that belongs to no block
    * but was directly transferred to the master problem and which has fractional value in the current solution */
   if( branchvar == NULL )
   {
      for( i = 0; i < npriobranchcands; i++ )
      {
         assert(GCGvarIsOriginal(branchcands[i]));

         /* continue if variable belongs to a block */
         if( GCGvarGetBlock(branchcands[i]) != -1 )
            continue;

         /* use pseudocost variable selection rule */
         if( usepseudocosts )
         {
            if( SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]) > maxpsscore )
            {
               branchvar = branchcands[i];
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               maxpsscore = SCIPgetVarPseudocostScore(scip, branchcands[i], branchcandssol[i]);
            }
         }
         /* use most fractional variable selection rule */
         else
         {
            /* compute fractionality */
            frac = MIN( branchcandsscore[i], 1.0 - branchcandsscore[i] );
            assert(frac > 0);

            if( frac >= maxfrac )
            {
               SCIPdebugMessage("Var %s has fractional value in current solution: %f\n",
                  SCIPvarGetName(branchcands[i]), branchcandssol[i]);
               solval = SCIPgetRelaxSolVal(scip, branchcands[i]);
               branchvar = branchcands[i];
               /* if we do not look for the most fractional variable, but for the first fractional variable,
                * we stop here since we found a variable to branch on */
               if( !mostfrac )
                  break;
            }
         }
      }
   }

   if( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(branchvar != NULL);

   SCIP_CALL( branchVar(scip, branchrule, branchvar, solval) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/* copy default SCIP branching rules to allow solving restrictions of the original problem as a subSCIP without
 * Dantzig-Wolfe decomposition
 */
static
SCIP_RETCODE GCGincludeOriginalCopyPlugins(
   SCIP* scip
   )
{
   SCIP_CALL( SCIPincludeNodeselEstimate(scip) );
   SCIP_CALL( SCIPincludeNodeselHybridestim(scip) );
   SCIP_CALL( SCIPincludeNodeselRestartdfs(scip) );
   SCIP_CALL( SCIPincludeBranchruleAllfullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleFullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleInference(scip) );
   SCIP_CALL( SCIPincludeBranchruleMostinf(scip) );
   SCIP_CALL( SCIPincludeBranchruleLeastinf(scip) );
   SCIP_CALL( SCIPincludeBranchrulePscost(scip) );
   SCIP_CALL( SCIPincludeBranchruleRandom(scip) );
   SCIP_CALL( SCIPincludeBranchruleRelpscost(scip) );
   return SCIP_OKAY;
}

/*
 * Callback methods for enforcing branching constraints
 */

#define branchDeactiveMasterOrig NULL
#define branchPropMasterOrig NULL

/** callback activation method */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterOrig)
{
   SCIP* origscip;
   SCIP_CONS* mastercons;

   assert(scip != NULL);
   assert(branchdata != NULL);

   /* branching restrictions are enforced by variable bounds, this is done automatically, so we can abort here */
   if( branchdata->cons == NULL )
      return SCIP_OKAY;

   assert(branchdata->origvar != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   SCIPdebugMessage("branchActiveMasterOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=" ), branchdata->newbound);

   /* transform constraint to the master variable space */
   SCIP_CALL( GCGrelaxTransOrigToMasterCons(origscip, branchdata->cons, &mastercons) );
   assert(mastercons != NULL);

   /* add constraint to the master problem */
   SCIP_CALL( SCIPaddConsNode(scip, SCIPgetCurrentNode(scip), mastercons, NULL) );

   /* constraint was added locally to the node where it is needed, so we do not need to care about this
    * at the next activation of the node and can set the constraint pointer to NULL */
   SCIP_CALL( SCIPreleaseCons(scip, &branchdata->cons) );
   branchdata->cons = NULL;

   return SCIP_OKAY;
}

/** callback solved method */
static
GCG_DECL_BRANCHMASTERSOLVED(branchMasterSolvedOrig)
{
   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->origvar != NULL);

   SCIPdebugMessage("branchMasterSolvedOrig: %s %s %f\n", SCIPvarGetName(branchdata->origvar),
      ( branchdata->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=" ), branchdata->newbound);

   if( !SCIPisInfinity(scip, newlowerbound) && SCIPgetStage(GCGrelaxGetMasterprob(scip)) == SCIP_STAGE_SOLVING
      && SCIPisRelaxSolValid(GCGrelaxGetMasterprob(scip)) )
   {
      SCIP_CALL( SCIPupdateVarPseudocost(scip, branchdata->origvar,
            SCIPgetRelaxSolVal(scip, branchdata->origvar) - branchdata->oldvalue,
            newlowerbound - branchdata->olddualbound, 1.0) );
   }

   return SCIP_OKAY;
}

/** callback deletion method for branching data */
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteOrig)
{
   assert(scip != NULL);
   assert(branchdata != NULL);

   SCIPdebugMessage("branchDataDeleteOrig: %s %s %f\n", SCIPvarGetName((*branchdata)->origvar),
      ( (*branchdata)->boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=" ), (*branchdata)->newbound);

   /* release constraint */
   if( (*branchdata)->cons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*branchdata)->cons) );
   }

   SCIPfreeMemory(scip, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpOrig)
{  /*lint --e{715}*/
   SCIPdebugMessage("Execlp method of orig branching\n");

   if( SCIPgetNExternBranchCands(scip) > 0 )
   {
      assert(SCIPisRelaxSolValid(scip));
      SCIP_CALL( branchExtern(scip, branchrule, result) );
   }

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextOrig)
{  /*lint --e{715}*/
   SCIP_CALL( branchExtern(scip, branchrule, result) );

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitOrig)
{
   assert(branchrule != NULL);

   SCIP_CALL( GCGrelaxIncludeBranchrule(scip, branchrule, branchActiveMasterOrig,
         branchDeactiveMasterOrig, branchPropMasterOrig, branchMasterSolvedOrig, branchDataDeleteOrig) );

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsOrig)
{  /*lint --e{715}*/
   int i;

   /* branching candidates */
   SCIP_VAR** branchcands;
   int nbranchcands;
   int npriobranchcands;

   /* values for choosing the variable to branch on */
   SCIP_VAR* branchvar;
   SCIP_Real solval;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of orig branching\n");

   *result = SCIP_DIDNOTRUN;
   if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) > SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* get the branching candidates */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &branchcands, &nbranchcands, &npriobranchcands) );

   branchvar = NULL;
   solval = 0.0;

   /* branch on an integer variable belonging to a unique block with fractional value */
   for( i = 0; i < npriobranchcands; i++ )
   {
      assert(GCGvarIsOriginal(branchcands[i]));

      /* variable belongs to no block or the block is not unique */
      if( GCGvarGetBlock(branchcands[i]) <= -1 || GCGrelaxGetNIdenticalBlocks(scip, GCGvarGetBlock(branchcands[i])) != 1 )
         continue;

      branchvar = branchcands[i];
      assert(SCIPvarGetUbLocal(branchvar) - SCIPvarGetLbLocal(branchvar) > 0.8);
      solval = SCIPvarGetLbLocal(branchvar) + 0.5;
   }

   /* we did not find a variable to branch on so far, so we look for an unfixed linking variable or an integer variable
    * that belongs to no block but was directly transferred to the master problem
    */
   if( branchvar == NULL )
   {
      for( i = 0; i < npriobranchcands; i++ )
      {
         assert(GCGvarIsOriginal(branchcands[i]));

         /* continue if variable belongs to a block */
         if( GCGvarGetBlock(branchcands[i]) > -1 )
            continue;

         /* check that blocks of linking variable are unique */
         if( GCGvarIsLinking(branchcands[i]) )
         {
            int nvarblocks;
            int* varblocks;
            SCIP_Bool unique;
            int j;

            nvarblocks = GCGlinkingVarGetNBlocks(branchcands[i]);
            SCIP_CALL( SCIPallocBufferArray(scip, &varblocks, nvarblocks) );
            SCIP_CALL( GCGlinkingVarGetBlocks(branchcands[i], nvarblocks, varblocks) );

            unique = TRUE;
            for( j = 0; j < nvarblocks; ++j )
               if( GCGrelaxGetNIdenticalBlocks(scip, varblocks[j]) != 1 )
                  unique = FALSE;

            SCIPfreeBufferArray(scip, &varblocks);

            if( !unique )
               continue;
         }

         branchvar = branchcands[i];
         assert(SCIPvarGetUbLocal(branchvar) - SCIPvarGetLbLocal(branchvar) > 0.8);
         solval = SCIPvarGetLbLocal(branchvar) + 0.5;
      }
   }

   if( branchvar == NULL )
   {
      SCIPdebugMessage("Original branching rule could not find a variable to branch on!\n");
      return SCIP_OKAY;
   }

   assert(branchvar != NULL);

   SCIP_CALL( branchVar(scip, branchrule, branchvar, solval) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

static
SCIP_DECL_BRANCHCOPY(branchCopyOrig)
{
   assert(scip != NULL);
   assert(branchrule != NULL);

   SCIPdebugMessage("orig copy called.\n");
   SCIP_CALL( GCGincludeOriginalCopyPlugins(scip) );

   return SCIP_OKAY;
}

/* define not used callback as NULL*/
#define branchFreeOrig NULL
#define branchExitOrig NULL
#define branchInitsolOrig NULL
#define branchExitsolOrig NULL


/*
 * branching specific interface methods
 */

/** creates the branching on original variable branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleOrig(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchCopyOrig,
         branchFreeOrig, branchInitOrig, branchExitOrig, branchInitsolOrig, branchExitsolOrig,
         branchExeclpOrig, branchExecextOrig, branchExecpsOrig, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/orig/enforcebycons",
         "should bounds on variables be enforced by constraints(TRUE) or by bounds(FALSE)",
         NULL, FALSE, DEFAULT_ENFORCEBYCONS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/orig/mostfrac",
         "should branching be performed on the most fractional variable instead of the first variable?",
         NULL, FALSE, DEFAULT_MOSTFRAC, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/orig/usepseudocosts",
         "should pseudocosts be used to determine the variable on which the branching is performed?",
         NULL, FALSE, DEFAULT_USEPSEUDO, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/orig/usepsstrong",
         "should strong branching with propagation be used to determine the variable on which the branching is performed?",
         NULL, FALSE, DEFAULT_USEPSSTRONG, NULL, NULL) );

   return SCIP_OKAY;
}
