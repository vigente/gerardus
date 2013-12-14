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

/**@file   heur_xprins.c
 * @brief  Extreme Point RINS
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "heur_xprins.h"
#include "pub_gcgvar.h"
#include "relax_gcg.h"
#include "gcgplugins.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"


#define HEUR_NAME             "xprins"
#define HEUR_DESC             "Extreme Point RINS"
#define HEUR_DISPCHAR         'Y'
#define HEUR_PRIORITY         -1100600
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE

#define DEFAULT_EQUALITYRATE  0.5           /**< minimum percentage of coincidence of relaxation and extreme pts     */
#define DEFAULT_MAXNODES      1000LL        /**< maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01          /**< factor by which xprins should at least improve the incumbent        */
#define DEFAULT_MINNODES      200LL         /**< minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINFIXINGRATE 0.5           /**< minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_NODESOFS      200LL         /**< number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1           /**< subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_NUSEDPTS      -1            /**< number of extreme pts per block that will be taken into account     */
#define DEFAULT_NWAITINGNODES 200LL         /**< number of nodes without incumbent change heuristic should wait      */
#define DEFAULT_RANDOMIZATION FALSE         /**< should the choice which sols to take be randomized?                 */
#define DEFAULT_DONTWAITATROOT FALSE        /**< should the nwaitingnodes parameter be ignored at the root node?     */
#define DEFAULT_USELPROWS     FALSE         /**< should subproblem be created out of the rows in the LP rows,
                                             * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE          /**< if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                             * of the original scip be copied to constraints of the subscip
                                             */




/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Real             equalityrate;       /**< minimum percentage of coincidence of relaxation and extreme pts   */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem               */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem               */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes        */
   SCIP_Longint          usednodes;          /**< nodes already used by xprins in earlier calls                     */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem     */
   int                   nusedpts;           /**< number of extreme pts per block that will be taken into account   */

   SCIP_Longint          nwaitingnodes;      /**< number of nodes without incumbent change heuristic should wait    */
   unsigned int          nfailures;          /**< number of failures since last successful call                     */
   SCIP_Longint          nextnodenumber;     /**< number of BnB nodes at which crossover should be called next      */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed     */
   SCIP_Real             minimprove;         /**< factor by which xprins should at least improve the incumbent      */
   SCIP_Bool             randomization;      /**< should the choice which sols to take be randomized?               */
   SCIP_Bool             dontwaitatroot;     /**< should the nwaitingnodes parameter be ignored at the root node?   */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?      */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?
                                              */
   unsigned int          randseed;           /**< seed value for random number generator                            */
};




/*
 * Local methods
 */


/** for each block, select extreme points (represented by mastervars) to be compared to the relaxation solution */
static
SCIP_RETCODE selectExtremePoints(
   SCIP*                 scip,               /**< original SCIP data structure                                    */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data                                           */
   int*                  selection,          /**< indices of selected extreme points                              */
   SCIP_Bool*            success             /**< pointer to store whether the process was successful             */
   )
{
   SCIP* masterprob;
   int nblocks;

   SCIP_VAR** mastervars;
   int nmastervars;

   int nusedpts;
   int block;
   int nidentblocks;
   SCIP_Real value;
   SCIP_Real* selvalue;

   int i;
   int j;
   int k;

   assert(scip != NULL);

   /* get master problem, number of blocks and extreme points per block */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   /* get number of blocks */
   nblocks = GCGrelaxGetNPricingprobs(scip);

   /* get variables of the master problem */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   /* get number of extreme points per block */
   nusedpts = heurdata->nusedpts;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &selvalue, nblocks * nusedpts) );

   *success = FALSE;

   /* loop over all given master variables */
   for( i = 0; i < nmastervars; i++ )
   {
      SCIP_VAR* mastervar;

      mastervar = mastervars[i];
      assert(GCGvarIsMaster(mastervar));

      /* get block information and solution value */
      block = GCGvarGetBlock(mastervar);
      value = SCIPgetSolVal(masterprob, NULL, mastervar);

      /** @todo handle infinite master solution values */
      assert(!SCIPisInfinity(scip, value));

      /* ignore irrelevant extreme points */
      if( SCIPisZero(scip, value) )
         continue;

      /* ignore rays
       * @todo do it smarter */
      if( GCGmasterVarIsRay(mastervar) )
         continue;

      /* variables belonging to no block are not treated here */
      if( block == -1 )
         continue;

      /* ignore "empty" master variables, i.e. variables representing the zero vector */
//      if( norigvars == 0 )
//         continue;

      /* get number of blocks that are identical to this block */
      assert(block >= 0);
      nidentblocks = GCGrelaxGetNIdenticalBlocks(scip, block);

      value = value / nidentblocks;

      /* check if the extreme point is good enough to be inserted in the selection */
      for( j = block * nusedpts; j < (block + 1) * nusedpts; ++j )
      {
         /* if the extreme point is better than a point in the selection
          * or there are < nusedpts, insert it */
         if( selection[j] == -1 || SCIPisGT(scip, value, selvalue[j]) )
         {
//            SCIPdebugMessage("insert new point: block %d, mastervar %d, value %g, pos %d\n",
//                  block+1, i, value, j % nusedpts);
            for( k = (block + 1) * nusedpts - 1; k > j; --k )
            {
//               SCIPdebugMessage("  shift point %d from pos %d to pos %d\n", selection[k-1], (k-1) % nusedpts, k % nusedpts);
               selection[k] = selection[k-1];
               selvalue[k] = selvalue[k-1];
            }
            selection[j] = i;
            selvalue[j] = value;
            break;
         }
      }
   }

   *success = TRUE;

   SCIPfreeBufferArray(scip, &selvalue);

   return SCIP_OKAY;
}


/** select extreme points (represented by mastervars) randomly */
static
SCIP_RETCODE selectExtremePointsRandomized(
   SCIP*                 scip,               /**< original SCIP data structure                                    */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data                                           */
   int*                  selection,          /**< indices of selected extreme points                              */
   SCIP_Bool*            success             /**< pointer to store whether the process was successful             */
   )
{
   SCIP* masterprob;
   int nblocks;

   SCIP_VAR** mastervars;
   int nmastervars;

   int nusedpts;         /* number of extreme points per block to be chosen        */
   int* npts;            /* for each block, the number of available extreme points */
   int* blockpts;        /* all points of a block which to be considered           */
   SCIP_Real* ptvals;    /* solution values of extreme points in master problem    */
   int lastpt;           /* the worst extreme point possible to choose             */

   int i;
   int j;
   int k;

   assert(scip != NULL);

   /* get master problem, number of blocks and extreme points per block */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   /* get number of blocks */
   nblocks = GCGrelaxGetNPricingprobs(scip);

   /* get variables of the master problem */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   /* get number of extreme points per block */
   nusedpts = heurdata->nusedpts;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &npts, nblocks) );

   *success = TRUE;

   /* check whether we have enough points per block to perform a randomization */
   for( i = 0; i < nblocks; ++i )
      npts[i] = 0;
   for( i = 0; i < nmastervars; ++i )
   {
      SCIP_VAR* mastervar;
      SCIP_Real solval;
      int block;

      mastervar = mastervars[i];
      solval = SCIPgetSolVal(masterprob, NULL, mastervar);
      block = GCGvarGetBlock(mastervar);

      if( block >= 0 && !SCIPisFeasZero(scip, solval) )
         ++npts[block];
   }
   for( i = 0; i < nblocks; ++i )
      if( GCGrelaxIsPricingprobRelevant(scip, i) && npts[i] <= nusedpts )
               *success = FALSE;

   /* do not randomize if there are not enough points available */
   if( !*success )
   {
      SCIPdebugMessage(" -> not enough extreme points available for randomization.\n");

      /* free memory */
      SCIPfreeBufferArray(scip, &npts);

      return SCIP_OKAY;
   }

   *success = FALSE;

   /* perform randomization: for each block, select a set of extreme points to be considered */
   for( i = 0; i < nblocks; ++i )
   {
      int blockrep;

      SCIP_CALL( SCIPallocBufferArray(scip, &blockpts, npts[i]) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ptvals, npts[i]) );

      /* get representative of this block */
      blockrep = GCGrelaxGetBlockRepresentative(scip, i);
      assert(blockrep >= 0 && blockrep <= i);

      /* get all relevant extreme points for this block */
      k = 0;
      for( j = 0; j < nmastervars; ++j )
      {
         SCIP_VAR* mastervar;
         SCIP_Real solval;
         int block;

         mastervar = mastervars[j];
         solval = SCIPgetSolVal(masterprob, NULL, mastervar);
         block = GCGvarGetBlock(mastervar);

         if( block == blockrep && !SCIPisFeasZero(scip, solval) )
         {
            assert(k < npts[blockrep]);
            blockpts[k] = j;
            ++k;
         }
      }
      assert(k == npts[blockrep]);

      /* sort the extreme points */
      SCIPsortRealInt(ptvals, blockpts, npts[blockrep]);
      lastpt = npts[blockrep];

      /* perform a random selection for this block */
      for( k = 0; k < nusedpts; ++k )
      {
         int idx;
         int selidx;

         idx = SCIPgetRandomInt(nusedpts-k-1, lastpt-1, &heurdata->randseed);
         selidx = i * nusedpts + k;
         selection[selidx] = blockpts[idx];
         lastpt = idx;
      }

      SCIPfreeBufferArray(scip, &ptvals);
      SCIPfreeBufferArray(scip, &blockpts);
   }

   *success = TRUE;

   /* free memory */
   SCIPfreeBufferArray(scip, &npts);

   return SCIP_OKAY;
}


/** initialize the subSCIP instance: copy SCIP to subSCIP, set the parameters */
static
SCIP_RETCODE initializeSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure                                  */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                        */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                               */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data                                         */
   SCIP_Longint          nstallnodes,        /**< node limit for subproblem                                     */
   SCIP_Real             timelimit,          /**< time limit for subproblem                                     */
   SCIP_Real             memorylimit,        /**< memory limit for subproblem                                   */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_VAR** vars;
   int nvars;

   SCIP_Real cutoff;                         /* objective cutoff for the subproblem                 */
   SCIP_Real upperbound;

   int i;

   char probname[SCIP_MAXSTRLEN];
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to subSCIP variables      */

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );

   /* copy the SCIP instance to the subSCIP */

   /* copy all plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* get name of the original problem and add the string "_extremeptsub" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_extremeptsub", SCIPgetProbName(scip));

   /* create the subproblem */
   SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* copy all variables */
   SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );

   /* if the lp rows are not used, also copy the constraints */
   if( !heurdata->uselprows )
   {
      SCIP_Bool valid;
      valid = FALSE;

      SCIP_CALL( SCIPcopyConss(scip, subscip, varmapfw, NULL, TRUE, FALSE, &valid) );
      if( heurdata->copycuts )
      {
         /** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, TRUE, NULL) );
      }
      SCIPdebugMessage("Copying the SCIP constraints was %s complete.\n", valid ? "" : "not ");
   }

   /* get the subproblem variables */
   for( i = 0; i < nvars; i++ )
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* setup parameters of subSCIP */
   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nstallnodes) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(scip, "estimate") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(scip, "inference") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useinflp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useboundlp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) );

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
      {
         cutoff = (1-heurdata->minimprove)*SCIPgetUpperbound(scip) + heurdata->minimprove*SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound ( scip ) >= 0 )
            cutoff = ( 1 - heurdata->minimprove ) * SCIPgetUpperbound ( scip );
         else
            cutoff = ( 1 + heurdata->minimprove ) * SCIPgetUpperbound ( scip );
      }
      cutoff = MIN(upperbound, cutoff );
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );
   }

   *success = TRUE;

   return SCIP_OKAY;
}


/** fix variables; for each variable, we evaluate the percentage of extreme points in which it has the same value
 *  as in the relaxation solution and fix it if the percentage exceeds a certain value */
static SCIP_RETCODE fixVariables(
   SCIP*                 scip,               /**< original SCIP data structure                                  */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                        */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                               */
   int*                  selection,          /**< pool of solutions crossover will use                          */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data                                         */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP* masterprob;                         /* master problem                         */
   SCIP_VAR** mastervars;                    /* master variables                       */
   int nmastervars;                          /* number of master variables             */

   SCIP_VAR** vars;                          /* original scip variables                */
   SCIP_Real fixingrate;                     /* percentage of variables that are fixed */

   int nblocks;                              /* number of blocks                                   */
   int nusedpts;                             /* number of extreme points per block                 */
   int nvars;                                /* number of original variables                       */
   int nbinvars;                             /* number of binary variables in the original problem */
   int nintvars;                             /* number of general integer variables                */

   int* neqpts;                              /* for each original variable, count the number of
						                              points where it has the same value as the relaxation solution */
   int* npts;                                /* for each block, the number of extreme points                  */
   int fixingcounter;                        /* count how many original variables are fixed                   */
   int zerocounter;                          /* count how many variables are fixed to zero                    */

   int i;
   int j;
   int k;
   int l;
   int selidx;

   /* get master problem and its variables */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nblocks = GCGrelaxGetNPricingprobs(scip);
   nusedpts = heurdata->nusedpts;
   fixingcounter = 0;
   zerocounter = 0;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &neqpts, nbinvars + nintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &npts, nbinvars + nintvars) );

   /* count the number of considered extreme points per block */
   if( nusedpts <= 0 )
   {
      assert(nusedpts == 0 || nusedpts == -1);

      for( i = 0; i < nblocks; ++i )
         npts[i] = 0;
      for( j = 0; j < nmastervars; ++j )
      {
         SCIP_VAR* mastervar;
         int block;

         mastervar = mastervars[j];
         assert(mastervar != NULL);
         block = GCGvarGetBlock(mastervar);
         if( block >= 0 )
            if( nusedpts == -1 || !SCIPisZero(scip, SCIPgetSolVal(masterprob, NULL, mastervar)) )
               ++npts[block];
      }
   }
   else
   {
      assert(selection != NULL);

      for( i = 0; i < nblocks; ++i )
      {
         for( j = 0; j < nusedpts; ++j )
         {
            selidx = i * nusedpts + j;

            if( selection[selidx] == -1 )
               break;
         }
         npts[i] = j;
      }
   }

   /* initialize counters; if the relaxation solution value is different to zero, we count upwards,
    * otherwise downwards */
   for( i = 0; i < nbinvars + nintvars; ++i )
   {
      SCIP_VAR* var;
      int block;
      SCIP_Real solval;

      /* get variable, block and relaxation value */
      var = vars[i];
      assert(var != NULL);
      block = GCGvarGetBlock(var);
      assert(block >= -2);
      solval = SCIPgetRelaxSolVal(scip, var);

      if( !SCIPisZero(scip, solval) || block == -1 )
         neqpts[i] = 0;
      else if( block == -2 )
      {
         SCIP_VAR** linkingpricingvars;

         assert(GCGvarIsLinking(var));
         linkingpricingvars = GCGlinkingVarGetPricingVars(var);

         neqpts[i] = 0;
         for( j = 0; j < nblocks; ++j )
            if( linkingpricingvars[j] != NULL )
               neqpts[i] += npts[j];
      }
      else
      {
         assert(block >= 0);
         neqpts[i] = npts[block];
      }
   }

   /* compare the relaxation solution to the extreme points */
   if( nusedpts <= 0 )
   {
      assert(nusedpts == 0 || nusedpts == -1);

      for( i = 0; i < nmastervars; ++i )
      {
         SCIP_VAR* mastervar;
         int block;
         SCIP_VAR** origvars;
         SCIP_Real* origvals;
         int norigvars;

         /* get master variable and block */
         mastervar = mastervars[i];
         assert(GCGvarIsMaster(mastervar));
         block = GCGvarGetBlock(mastervar);

         /* ignore copied original variables */
         if( block == -1 )
            continue;
         assert(block >= 0);

         /* ignore master variable if it is zero and only the nonzeroes should be considered */
         if( nusedpts == 0 && SCIPisZero(scip, SCIPgetSolVal(masterprob, NULL, mastervar)) )
            continue;

         /* get extreme point */
         origvars = GCGmasterVarGetOrigvars(mastervar);
         origvals = GCGmasterVarGetOrigvals(mastervar);
         norigvars = GCGmasterVarGetNOrigvars(mastervar);

         /* compare the solution values;
          * count downwards if relaxation solution is zero, count upwards otherwise */
         for( k = 0; k < norigvars; ++k )
         {
            SCIP_VAR* pricingvar;
            SCIP_VAR** pricingorigvars;
            int npricingorigvars;

            if( SCIPvarGetType(origvars[k]) > SCIP_VARTYPE_INTEGER )
               continue;

            /* get the corresponding pricing variable */
            if( GCGvarIsLinking(origvars[k]) )
            {
               SCIP_VAR** linkingpricingvars;

               linkingpricingvars = GCGlinkingVarGetPricingVars(origvars[k]);
               pricingvar = linkingpricingvars[block];
            }
            else
               pricingvar = GCGoriginalVarGetPricingVar(origvars[k]);
          
            assert(pricingvar != NULL);
            assert(GCGvarIsPricing(pricingvar));

            /* get all origvars represented by the current origvar */
            pricingorigvars = GCGpricingVarGetOrigvars(pricingvar);
            npricingorigvars = GCGpricingVarGetNOrigvars(pricingvar);
            assert(pricingorigvars != NULL);
            assert(npricingorigvars >= 0);

            for( l = 0; l < npricingorigvars; ++l )
            {
               int idx;
               SCIP_Real solval;

               idx = SCIPvarGetProbindex(pricingorigvars[l]);
               assert(idx < nbinvars + nintvars);
               solval = SCIPgetRelaxSolVal(scip, pricingorigvars[l]);

               if( SCIPisZero(scip, solval) )
               {
                  if( !SCIPisZero(scip, origvals[k]) )
                     --neqpts[idx];
               }
               else
               {
                  if( SCIPisEQ(scip, solval, origvals[k]) )
                     ++neqpts[idx];
               }
            }
         }
      }
   }
   else
   {
      assert(selection != NULL);

      for( i = 0; i < nblocks; ++i )
      {
         /* ignore blocks represented by others */
         if( !GCGrelaxIsPricingprobRelevant(scip, i) )
            continue;

         /* compare the relaxation solution to the selected extreme points */
         for( j = 0; j < nusedpts; ++j )
         {
            selidx = i * nusedpts + j;
            if( selection[selidx] != -1 )
            {
               SCIP_VAR* mastervar;
               SCIP_VAR** origvars;
               SCIP_Real* origvals;
               int norigvars;

               /* get master variable */
               mastervar = mastervars[selection[selidx]];
               assert(mastervar != NULL);
               assert(GCGvarIsMaster(mastervar));
               assert(GCGvarGetBlock(mastervar) == i);

               /* get extreme point */
               origvars = GCGmasterVarGetOrigvars(mastervar);
               origvals = GCGmasterVarGetOrigvals(mastervar);
               norigvars = GCGmasterVarGetNOrigvars(mastervar);

               /* compare the solution values;
                * count downwards if relaxation solution is zero, count upwards otherwise */
               for( k = 0; k < norigvars; ++k )
               {
                  SCIP_VAR* pricingvar;
                  SCIP_VAR** pricingorigvars;
                  int npricingorigvars;

                  if( SCIPvarGetType(origvars[k]) > SCIP_VARTYPE_INTEGER )
                     continue;

                  /* get the corresponding pricing variable */
                  if( GCGvarIsLinking(origvars[k]) )
                  {
                     SCIP_VAR** linkingpricingvars;

                     linkingpricingvars = GCGlinkingVarGetPricingVars(origvars[k]);
                     pricingvar = linkingpricingvars[i];
                  }
                  else
                     pricingvar = GCGoriginalVarGetPricingVar(origvars[k]);
          
                  assert(pricingvar != NULL);
                  assert(GCGvarIsPricing(pricingvar));

                  /* get all origvars represented by the current origvar */
                  pricingorigvars = GCGpricingVarGetOrigvars(pricingvar);
                  npricingorigvars = GCGpricingVarGetNOrigvars(pricingvar);
                  assert(pricingorigvars != NULL);
                  assert(npricingorigvars >= 0);

                  for( l = 0; l < npricingorigvars; ++l )
                  {
                     int idx;
                     SCIP_Real solval;

                     idx = SCIPvarGetProbindex(pricingorigvars[l]);
                     assert(idx < nbinvars + nintvars);
                     solval = SCIPgetRelaxSolVal(scip, pricingorigvars[l]);

                     if( SCIPisZero(scip, solval) )
                     {
                        if( !SCIPisZero(scip, origvals[k]) )
                           --neqpts[idx];
                     }
                     else
                     {
                        if( SCIPisEQ(scip, solval, origvals[k]) )
                           ++neqpts[idx];
                     }
                  }
               }
            }
         }
      }
   }

   /* try to fix the binary and general integer variables */
   for( i = 0; i < nbinvars + nintvars; ++i )
   {
      SCIP_VAR* var;
      int block;                             /* current block we are working in                    */
      SCIP_Real solval;

      var = vars[i];
      assert(GCGvarIsOriginal(var));
      block = GCGvarGetBlock(var);
      solval = SCIPgetRelaxSolVal(scip, var);

      /* if the variable is represented by another one, it is not treated here */
      if( block >= 0 && !GCGrelaxIsPricingprobRelevant(scip, block) )
         continue;

      /* we still need to treat variables belonging to no block (as they did not appear in any extreme point) */
      /* if the variable belongs to no block, fix it in a RENS-like fashion */
      if( block == -1 )
      {
         if( SCIPisFeasIntegral(scip, solval) )
         {
            /* fix variable to current relaxation solution if it is integral,
             * use exact integral value, if the variable is only integral within numerical tolerances
             */
            solval = SCIPfloor(scip, solval + 0.5);
            SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], solval) );
            SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], solval) );

            fixingcounter++;
            if( SCIPisZero(scip, solval) )
               zerocounter++;
         }
      }
      else
      {
         SCIP_Real quoteqpts;

         assert(block == -2 || block >= 0);

         /* evaluate percentage of extreme points having the same variable value as the relaxation solution */
         if( block >= 0 )
         {
            assert(neqpts[i] <= npts[block]);
            quoteqpts = (SCIP_Real) neqpts[i] / (SCIP_Real) MAX(npts[block],1);

            SCIPdebugMessage("Variable %s: %d/%d (%.2f percent) extreme points identical to relaxation solution.\n",
                        SCIPvarGetName(var), neqpts[i], npts[block], quoteqpts * 100);
         }
         else
         {
            SCIP_VAR** linkingpricingvars;
            int ntotalpts;

            assert(GCGvarIsLinking(var));
            linkingpricingvars = GCGlinkingVarGetPricingVars(var);

            ntotalpts = 0;
            for( j = 0; j < nblocks; ++j )
               if( linkingpricingvars[j] != NULL )
                  ntotalpts += npts[j];

            assert(neqpts[i] <= ntotalpts);
            quoteqpts = (SCIP_Real) neqpts[i] / (SCIP_Real) MAX(ntotalpts,1);

            SCIPdebugMessage("Variable %s: %d/%d (%.2f percent) extreme points identical to relaxation solution.\n",
                        SCIPvarGetName(var), neqpts[i], ntotalpts, quoteqpts * 100);
         }

         /* the variable can be fixed if the relaxation value is shared by enough extreme points */
         if( quoteqpts >= heurdata->equalityrate )
         {
            SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], solval) );
            SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], solval) );

            fixingcounter++;
            if( SCIPisZero(scip, solval) )
               zerocounter++;
         }
      }
   }

   fixingrate = (SCIP_Real) fixingcounter / (SCIP_Real) (MAX(nbinvars + nintvars, 1));

   SCIPdebugMessage("subSCIP: %i out of %i (%.2f percent) variables have been fixed.\n", fixingcounter, nbinvars + nintvars, fixingrate * 100.0);
   SCIPdebugMessage("subSCIP: %i out of %i (%.2f percent) fixed variables are zero.\n", zerocounter, fixingcounter,
         (SCIP_Real)zerocounter / MAX((SCIP_Real)fixingcounter,1.0) * 100.0);

   /* if all variables were fixed or amount of fixed variables is insufficient, skip residual part of
    * subproblem creation ans abort immediately */
   if( fixingcounter == nbinvars + nintvars || fixingrate < heurdata->minfixingrate )
   {
      *success = FALSE;
      SCIPdebugMessage("Fixing of variables was not successful - fixing rate %f percent.\n", fixingrate * 100.0);
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &neqpts);
   SCIPfreeBufferArray(scip, &npts);

   return SCIP_OKAY;
}


/** creates the rows of the subproblem by copying LP rows of the SCIP instance;
 *  only used if the uselprows parameter is TRUE */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< original SCIP data structure                                  */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                        */
   SCIP_VAR**            subvars             /**< the variables of the subproblem                               */
   )
{
   SCIP_ROW** rows;                          /* original scip rows                       */
   SCIP_CONS* cons;                          /* new constraint                           */
   SCIP_VAR** consvars;                      /* new constraint's variables               */
   SCIP_COL** cols;                          /* original row's columns                   */

   SCIP_Real constant;                       /* constant added to the row                */
   SCIP_Real lhs;                            /* left hand side of the row                */
   SCIP_Real rhs;                            /* left right side of the row               */
   SCIP_Real* vals;                          /* variables' coefficient values of the row */

   int nrows;
   int nnonz;
   int i;
   int j;

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* copy all rows to linear constraints */
   for( i = 0; i < nrows; i++ )
   {
      /* ignore rows that are only locally valid */
      if( SCIProwIsLocal(rows[i]) )
         continue;

      /* get the row's data */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIProwGetRhs(rows[i]) - constant;
      vals = SCIProwGetVals(rows[i]);
      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);

      assert(lhs <= rhs);

      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ )
         consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );

      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< crossover heuristic structure                       */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   int*                  solindex,           /**< index of the solution                               */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );
   *solindex = SCIPsolGetIndex(newsol);

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   if( *success )
   {
      SCIPdebugMessage("Extreme Point RINS: new solution added.\n");
   }

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** updates heurdata after a run of crossover */
static
void updateFailureStatistic(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data                               */
   )
{
   /* increase number of failures, calculate next node at which crossover should be called and update actual solutions */
   heurdata->nfailures++;
   heurdata->nextnodenumber = (heurdata->nfailures <= 25
      ? SCIPgetNNodes(scip) + 100*(2LL << heurdata->nfailures) /*lint !e703*/
      : SCIP_LONGINT_MAX);
}



/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyXprins NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeXprins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitXprins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->randseed = 0;
   heurdata->nfailures = 0;
   heurdata->nextnodenumber = 0;

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitXprins NULL

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolXprins NULL

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolXprins NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecXprins)
{  /*lint --e{715}*/

   SCIP* masterprob;
   SCIP_HEURDATA* heurdata;

   SCIP* subscip;
   SCIP_VAR** subvars;

   SCIP_Real memorylimit;                    /* memory limit for the subproblem                     */
   SCIP_Real timelimit;                      /* time limit for the subproblem                       */
   SCIP_Longint nstallnodes;                 /* node limit for the subproblem                       */

   int* selection;
   int nblocks;

   SCIP_Bool success;

   int i;

#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   nblocks = GCGrelaxGetNPricingprobs(scip);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DELAYED;

   /* do not execute the heuristic on invalid relaxation solutions
    * (which is the case if the node has been cut off)
    */
   if( !SCIPisRelaxSolValid(scip) )
   {
      SCIPdebugMessage("skipping Extreme Point RINS: invalid relaxation solution\n");
      return SCIP_OKAY;
   }

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetStage(masterprob) > SCIP_STAGE_SOLVING || SCIPgetLPSolstat(masterprob) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("skipping Extreme Point RINS: master LP not solved to optimality.\n");
      return SCIP_OKAY;
   }

   assert(SCIPhasCurrentNodeLP(masterprob));

   /* if heuristic should be delayed, wait until certain number of nodes is reached */
   if( SCIPgetNNodes(scip) < heurdata->nextnodenumber )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* only continue with some fractional variables */
   if( SCIPgetNExternBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward Crossover if it succeeded often */
   nstallnodes = (SCIP_Longint)
                              (nstallnodes * (1.0 + 2.0*(SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur)+1.0)));

   /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMessage("skipping Extreme Point RINS: nstallnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( timelimit < 10.0 || memorylimit <= 0.0 )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPdebugMessage("Executing Extreme Point RINS heuristic ...\n");

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, SCIPgetNVars(scip)) );
   selection = NULL;

   if( heurdata->nusedpts > 0 )
   {
      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &selection, nblocks * heurdata->nusedpts) );

      /* initialize empty selection */
      for( i = 0; i < nblocks * heurdata->nusedpts; ++i )
         selection[i] = -1;

      /* for each block, select extreme points (represented by master variables) to perform RINS */
      success = FALSE;
      if( heurdata->randomization )
      {
         SCIPdebugMessage("selecting extreme points randomly...\n");
         SCIP_CALL( selectExtremePointsRandomized(scip, heurdata, selection, &success) );
      }
      if( !heurdata->randomization || !success )
      {
         SCIPdebugMessage("selecting extreme points...\n");
         SCIP_CALL( selectExtremePoints(scip, heurdata, selection, &success) );
      }

      /* do not execute heuristic if no new selection of extreme points was found */
      if( !success )
      {
         SCIPdebugMessage("no proper selection could be created - aborting heuristic.\n");

         updateFailureStatistic(scip, heurdata);

         /* free memory */
         SCIPfreeBufferArray(scip, &selection);
         SCIPfreeBufferArray(scip, &subvars);

         return SCIP_OKAY;
      }
   }

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( initializeSubproblem(scip, subscip, subvars, heurdata, nstallnodes, timelimit, memorylimit, &success) );

   /* fix variables the variables of the subproblem */
   SCIP_CALL( fixVariables(scip, subscip, subvars, selection, heurdata, &success) );

   /* if creation of subscip was aborted (e.g. due to number of fixings), free subscip and abort */
   if( !success )
   {
      /* if creation was aborted due to number of fixings, free the already created subproblem */
      if( SCIPgetStage(subscip) != SCIP_STAGE_INIT )
      {
         int nbinvars;
         int nintvars;
         SCIP_CALL( SCIPgetVarsData(scip, NULL, NULL, &nbinvars, &nintvars, NULL, NULL) );
//            for( i = 0; i < nbinvars + nintvars; i++ )
//            {
//               SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
//            }
         SCIP_CALL( SCIPfreeTransform(subscip) );
      }

      /* free memory */
      if( heurdata->nusedpts > 0 )
         SCIPfreeBufferArray(scip, &selection);
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );

      /* this run will be counted as a failure since no new solution tuple could be generated or the neighborhood of the
       * solution was not fruitful in the sense that it was too big
       */
      updateFailureStatistic(scip, heurdata);

      return SCIP_OKAY;
   }

   /* if enough variables could be fixed, create rows of the subproblem */
   if( heurdata->uselprows )
   {
      SCIP_CALL( createRows(scip, subscip, subvars) );
   }

   *result = SCIP_DIDNOTFIND;

   /* solve the subproblem */
   SCIPdebugMessage("subSCIP: Solving... (node limit = %lld, time limit = %.2g)\n", nstallnodes, timelimit);

   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is catched and a warning is printed, only in debug mode, SCIP will stop.
    */
#ifdef NDEBUG
   retstat = SCIPsolve(subscip);
   if( retstat != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while solving subMIP in Extreme Point RINS heuristic; subSCIP terminated with code <%d>\n",
            retstat);
   }
#else
   SCIP_CALL( SCIPsolve(subscip) );
#endif

   heurdata->usednodes += SCIPgetNNodes(subscip);

   /* check, whether a solution was found */
   success = FALSE;
   if( SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;
      int solindex;                             /* index of the solution created by crossover          */

      SCIPdebugMessage("  -> found %i feasible solution(s).\n", SCIPgetNSols(subscip));

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      solindex = -1;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &solindex, &success) );
      }

      if( success )
         *result = SCIP_FOUNDSOL;
      else
         updateFailureStatistic(scip, heurdata);
   }
   else
   {
      /* if no new solution was found, run was a failure */
      updateFailureStatistic(scip, heurdata);
      SCIPdebugMessage(" -> no subMIP solution found - subSCIP status is %d\n", SCIPgetStatus(subscip));
   }

   /* free subproblem */
   SCIP_CALL( SCIPfreeTransform(subscip) );
   //      for( i = 0; i < nvars; i++ )
   //      {
   //         SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
   //      }
   SCIP_CALL( SCIPfree(&subscip) );

   /* free memory */
   if( heurdata->nusedpts > 0 )
      SCIPfreeBufferArray(scip, &selection);
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the xprins primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurXprins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create extreme point based fixing primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyXprins, heurFreeXprins, heurInitXprins, heurExitXprins,
         heurInitsolXprins, heurExitsolXprins, heurExecXprins,
         heurdata) );

   /* add extreme point based fixing primal heuristic parameters */

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/equalityrate",
	 "minimum percentage of coincidence of relaxation and extreme pts",
         &heurdata->equalityrate, FALSE, DEFAULT_EQUALITYRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/nusedpts",
         "number of extreme pts per block that will be taken into account",
         &heurdata->nusedpts, FALSE, DEFAULT_NUSEDPTS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
         "minimum percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which crossover should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/randomization",
         "should the choice which sols to take be randomized?",
         &heurdata->randomization, TRUE, DEFAULT_RANDOMIZATION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/dontwaitatroot",
         "should the nwaitingnodes parameter be ignored at the root node?",
         &heurdata->dontwaitatroot, TRUE, DEFAULT_DONTWAITATROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
