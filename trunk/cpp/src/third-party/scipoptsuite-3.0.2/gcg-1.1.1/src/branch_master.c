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

/**@file   branch_master.c
 * @brief  branching rule for master problem
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_master.h"
#include "cons_origbranch.h"
#include "cons_masterbranch.h"

#include "scip/nodesel_bfs.h"
#include "scip/nodesel_dfs.h"
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


#define BRANCHRULE_NAME          "master"
#define BRANCHRULE_DESC          "branching for generic column generation master"
#define BRANCHRULE_PRIORITY      1000000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

/** includes all plugins in the master copy.
 *
 * This method makes SCIPcopy work for heuristics in the master
 */
static
SCIP_RETCODE GCGincludeMasterCopyPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludeNodeselBfs(scip) );
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );
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
 * Callback methods
 */

/** copy method for master branching rule */
static
SCIP_DECL_BRANCHCOPY(branchCopyMaster)
{
   assert(branchrule != NULL);
   assert(scip != NULL);
   SCIPdebugMessage("pricer copy called.\n");
   SCIP_CALL( GCGincludeMasterCopyPlugins(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#define branchFreeMaster NULL


/** initialization method of branching rule (called after problem was transformed) */
#define branchInitMaster NULL


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitMaster NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolMaster NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolMaster NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMaster)
{  /*lint --e{715}*/
   SCIP_NODE* child1;
   SCIP_NODE* child2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;

   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of master branching\n");

   /* create two child-nodes of the current node in the b&b-tree and add the masterbranch constraints */
   SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &child2, 0.0, SCIPgetLocalTransEstimate(scip)) );

   SCIP_CALL( GCGcreateConsMasterbranch(scip, &cons1, child1, GCGconsMasterbranchGetActiveCons(scip)) );
   SCIP_CALL( GCGcreateConsMasterbranch(scip, &cons2, child2, GCGconsMasterbranchGetActiveCons(scip)) );

   SCIP_CALL( SCIPaddConsNode(scip, child1, cons1, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, child2, cons2, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons2) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** branching execution method relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextMaster)
{  /*lint --e{715}*/
   SCIPdebugMessage("Execext method of master branching\n");

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsMaster)
{  /*lint --e{715}*/
   SCIP_NODE* child1;
   SCIP_NODE* child2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;

   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of master branching\n");

   /* create two child-nodes of the current node in the b&b-tree and add the masterbranch constraints */
   SCIP_CALL( SCIPcreateChild(scip, &child1, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &child2, 0.0, SCIPgetLocalTransEstimate(scip)) );

   SCIP_CALL( GCGcreateConsMasterbranch(scip, &cons1, child1, GCGconsMasterbranchGetActiveCons(scip)) );
   SCIP_CALL( GCGcreateConsMasterbranch(scip, &cons2, child2, GCGconsMasterbranchGetActiveCons(scip)) );

   SCIP_CALL( SCIPaddConsNode(scip, child1, cons1, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, child2, cons2, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons2) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/*
 * branching specific interface methods
 */

/** creates the master branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMaster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create inference branching rule data */
   branchruledata = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchCopyMaster, branchFreeMaster, branchInitMaster, branchExitMaster, branchInitsolMaster,
         branchExitsolMaster, branchExeclpMaster, branchExecextMaster, branchExecpsMaster,
         branchruledata) );

   return SCIP_OKAY;
}
