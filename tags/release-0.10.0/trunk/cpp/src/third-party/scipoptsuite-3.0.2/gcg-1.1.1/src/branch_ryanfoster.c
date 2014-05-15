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

/**@file   branch_ryanfoster.c
 * @brief  branching rule for original problem in GCG implementing the Ryan and Foster branching scheme
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_ryanfoster.h"
#include "relax_gcg.h"
#include "cons_origbranch.h"
#include "pricer_gcg.h"
#include "scip/cons_varbound.h"
#include "type_branchgcg.h"
#include "pub_gcgvar.h"

#define BRANCHRULE_NAME          "ryanfoster"
#define BRANCHRULE_DESC          "ryan and foster branching in generic column generation"
#define BRANCHRULE_PRIORITY      10
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0


/** branching data for branching decisions */
struct GCG_BranchData
{
   SCIP_VAR*             var1;               /**< first original variable on which the branching is done */
   SCIP_VAR*             var2;               /**< second original variable on which the branching is done */
   SCIP_Bool             same;               /**< should each master var contain either both or none of the vars? */
   int                   blocknr;            /**< number of the block in which branching was performed */
   SCIP_CONS*            pricecons;          /**< constraint enforcing the branching restriction in the pricing problem */
};


/*
 * Callback methods for enforcing branching constraints
 */

/** callback activation method */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterRyanfoster)
{
   SCIP* origscip;
   SCIP* pricingscip;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->var1 != NULL);
   assert(branchdata->var2 != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   pricingscip = GCGrelaxGetPricingprob(origscip, branchdata->blocknr);
   assert(pricingscip != NULL);

   SCIPdebugMessage("branchActiveMasterRyanfoster: %s(%s, %s)\n", ( branchdata->same ? "same" : "differ" ),
      SCIPvarGetName(branchdata->var1), SCIPvarGetName(branchdata->var2));

   assert(GCGvarIsOriginal(branchdata->var1));
   /** @todo it is not clear if linking variables interfere with ryan foster branching */
   assert(GCGvarGetBlock(branchdata->var1) == branchdata->blocknr);

   assert(GCGvarIsOriginal(branchdata->var2));
   assert(GCGvarGetBlock(branchdata->var2) == branchdata->blocknr);

   /* create corresponding constraint in the pricing problem, if not yet created */
   if( branchdata->pricecons == NULL )
   {
      if( branchdata->same )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "same(%s, %s)", SCIPvarGetName(branchdata->var1),
            SCIPvarGetName(branchdata->var2));

         SCIP_CALL( SCIPcreateConsVarbound(pricingscip,
               &(branchdata->pricecons), name, GCGoriginalVarGetPricingVar(branchdata->var1),
               GCGoriginalVarGetPricingVar(branchdata->var2), -1.0, 0.0, 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }
      else
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "differ(%s, %s)", SCIPvarGetName(branchdata->var1),
            SCIPvarGetName(branchdata->var2));

         SCIP_CALL( SCIPcreateConsVarbound(pricingscip,
               &(branchdata->pricecons), name, GCGoriginalVarGetPricingVar(branchdata->var1),
               GCGoriginalVarGetPricingVar(branchdata->var2), 1.0, -SCIPinfinity(scip), 1.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }
   }
   /* add constraint to the pricing problem that enforces the branching decision */
   SCIP_CALL( SCIPaddCons(pricingscip, branchdata->pricecons) );

   return SCIP_OKAY;
}

/** callback deactivation method */
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterRyanfoster)
{
   SCIP* origscip;
   SCIP* pricingscip;

   assert(scip != NULL);
    assert(branchdata != NULL);
   assert(branchdata->var1 != NULL);
   assert(branchdata->var2 != NULL);
   assert(branchdata->pricecons != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   pricingscip = GCGrelaxGetPricingprob(origscip, branchdata->blocknr);
   assert(pricingscip != NULL);

   SCIPdebugMessage("branchDeactiveMasterRyanfoster: %s(%s, %s)\n", ( branchdata->same ? "same" : "differ" ),
      SCIPvarGetName(branchdata->var1), SCIPvarGetName(branchdata->var2));

   /* remove constraint from the pricing problem that enforces the branching decision */
   assert(branchdata->pricecons != NULL);
   SCIP_CALL( SCIPdelCons(pricingscip, branchdata->pricecons) );

   return SCIP_OKAY;
}

/** callback propagation method */
static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterRyanfoster)
{
   SCIP_VAR** vars;
   SCIP_Real val1;
   SCIP_Real val2;
   int nvars;
   int propcount;
   int i;
   int j;

   assert(scip != NULL);
   assert(branchdata != NULL);
   assert(branchdata->var1 != NULL);
   assert(branchdata->var2 != NULL);
   assert(branchdata->pricecons != NULL);

   assert(GCGpricerGetOrigprob(scip) != NULL);

   SCIPdebugMessage("branchPropMasterRyanfoster: %s(%s, %s)\n", ( branchdata->same ? "same" : "differ" ),
      SCIPvarGetName(branchdata->var1), SCIPvarGetName(branchdata->var2));

   *result = SCIP_DIDNOTFIND;

   propcount = 0;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* iterate over all master variables */
   for( i = 0; i < nvars; i++ )
   {
      int norigvars;
      SCIP_Real* origvals;
      SCIP_VAR** origvars;

      origvars = GCGmasterVarGetOrigvars(vars[i]);
      origvals = GCGmasterVarGetOrigvals(vars[i]);
      norigvars = GCGmasterVarGetNOrigvars(vars[i]);

      /* only look at variables not fixed to 0 */
      if( !SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[i])) )
      {
         assert(GCGvarIsMaster(vars[i]));

         /* if variable belongs to a different block than the branching restriction, we do not have to look at it */
         if( branchdata->blocknr != GCGvarGetBlock(vars[i]) )
            continue;

         /* save the values of the original variables for the current master variable */
         val1 = 0.0;
         val2 = 0.0;
         for( j = 0; j < norigvars; j++ )
         {
            if( origvars[j] == branchdata->var1 )
            {
               val1 = origvals[j];
               continue;
            }
            if( origvars[j] == branchdata->var2 )
            {
               val2 = origvals[j];
            }
         }

         /* if branching enforces that both original vars are either both contained or none of them is contained
          * and the current master variable has different values for both of them, fix the variable to 0 */
         if( branchdata->same && !SCIPisEQ(scip, val1, val2) )
         {
            SCIP_CALL( SCIPchgVarUb(scip, vars[i], 0.0) );
            propcount++;
         }
         /* if branching enforces that both original vars must be in different mastervars, fix all
          * master variables to 0 that contain both */
         if( !branchdata->same && SCIPisEQ(scip, val1, 1.0) && SCIPisEQ(scip, val2, 1.0) )
         {
            SCIP_CALL( SCIPchgVarUb(scip, vars[i], 0.0) );
            propcount++;
         }
      }
   }

   SCIPdebugMessage("Finished propagation of branching decision constraint: %s(%s, %s), %d vars fixed.\n",
      ( branchdata->same ? "same" : "differ" ), SCIPvarGetName(branchdata->var1), SCIPvarGetName(branchdata->var2), propcount);

   if( propcount > 0 )
   {
      *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}

/** callback deletion method for branching data*/
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteRyanfoster)
{
   assert(scip != NULL);
   assert(branchdata != NULL);

   SCIPdebugMessage("branchDataDeleteRyanfoster: %s(%s, %s)\n", ( (*branchdata)->same ? "same" : "differ" ),
      SCIPvarGetName((*branchdata)->var1), SCIPvarGetName((*branchdata)->var2));

   /* release constraint that enforces the branching decision */
   if( (*branchdata)->pricecons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(GCGrelaxGetPricingprob(scip, (*branchdata)->blocknr),
            &(*branchdata)->pricecons) );
   }

   SCIPfreeMemory(scip, branchdata);
   *branchdata = NULL;

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** for the two given original variables, create two Ryan&Foster branching nodes, one for same, one for differ */
static
SCIP_RETCODE createChildNodesRyanfoster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_VAR*             ovar1,              /**< first original variable */
   SCIP_VAR*             ovar2,              /**< second original variable */
   int                   blocknr             /**< number of the pricing block */
   )
{
   SCIP_NODE* childsame;
   SCIP_NODE* childdiffer;
   SCIP_CONS* origbranchsame;
   SCIP_CONS* origbranchdiffer;
   SCIP_VAR* pricingvar1;
   SCIP_VAR* pricingvar2;
   GCG_BRANCHDATA* branchsamedata;
   GCG_BRANCHDATA* branchdifferdata;
   char samename[SCIP_MAXSTRLEN];
   char differname[SCIP_MAXSTRLEN];

   SCIP_VAR** origvars1;
   SCIP_VAR** origvars2;
   int norigvars1;
   int v;

   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(ovar1 != NULL);
   assert(ovar2 != NULL);
   assert(GCGvarIsOriginal(ovar1));
   assert(GCGvarIsOriginal(ovar2));

   SCIPdebugMessage("Ryanfoster branching rule: branch on original variables %s and %s!\n",
      SCIPvarGetName(ovar1), SCIPvarGetName(ovar2));

   /* create the b&b-tree child-nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdiffer, 0.0, SCIPgetLocalTransEstimate(scip)) );

   /* allocate branchdata for same child and store information */
   SCIP_CALL( SCIPallocMemory(scip, &branchsamedata) );
   branchsamedata->var1 = ovar1;
   branchsamedata->var2 = ovar2;
   branchsamedata->same = TRUE;
   branchsamedata->blocknr = blocknr;
   branchsamedata->pricecons = NULL;

   /* allocate branchdata for differ child and store information */
   SCIP_CALL( SCIPallocMemory(scip, &branchdifferdata) );
   branchdifferdata->var1 = ovar1;
   branchdifferdata->var2 = ovar2;
   branchdifferdata->same = FALSE;
   branchdifferdata->blocknr = blocknr;
   branchdifferdata->pricecons = NULL;

   /* define names for origbranch constraints */
   (void) SCIPsnprintf(samename, SCIP_MAXSTRLEN, "same(%s, %s)", SCIPvarGetName(branchsamedata->var1),
      SCIPvarGetName(branchsamedata->var2));
   (void) SCIPsnprintf(differname, SCIP_MAXSTRLEN, "differ(%s, %s)", SCIPvarGetName(branchsamedata->var1),
      SCIPvarGetName(branchsamedata->var2));

   /* create origbranch constraints */
   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchsame, samename, childsame,
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchsamedata) );
   SCIP_CALL( GCGcreateConsOrigbranch(scip, &origbranchdiffer, differname, childdiffer,
         GCGconsOrigbranchGetActiveCons(scip), branchrule, branchdifferdata) );

   /* add and release constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childsame, origbranchsame, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdiffer, origbranchdiffer, NULL) );
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchsame) );
   SCIP_CALL( SCIPreleaseCons(scip, &origbranchdiffer) );

   pricingvar1 = GCGoriginalVarGetPricingVar(branchdifferdata->var1);
   pricingvar2 = GCGoriginalVarGetPricingVar(branchdifferdata->var2);
   assert(GCGvarIsPricing(pricingvar1));
   assert(GCGvarIsPricing(pricingvar2));
   assert(GCGvarGetBlock(pricingvar1) == GCGvarGetBlock(pricingvar2));
   assert(GCGpricingVarGetNOrigvars(pricingvar1) == GCGpricingVarGetNOrigvars(pricingvar2));

   norigvars1 = GCGpricingVarGetNOrigvars(pricingvar1);
   assert(norigvars1 == GCGpricingVarGetNOrigvars(pricingvar2));

   origvars1 = GCGpricingVarGetOrigvars(pricingvar1);
   origvars2 = GCGpricingVarGetOrigvars(pricingvar2);

   /* add branching decision as varbound constraints to original problem */
   for( v = 0; v < norigvars1; v++ )
   {
      SCIP_CONS* origcons;

      assert(GCGvarGetBlock(origvars1[v]) == GCGvarGetBlock(origvars2[v]));

      /* create constraint for same-child */
      SCIP_CALL( SCIPcreateConsVarbound(scip, &origcons, samename, origvars1[v], origvars2[v],
            -1.0, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      /* add cons locally to the problem and release it */
      SCIP_CALL( SCIPaddConsNode(scip, childsame, origcons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &origcons) );

      /* create constraint for differ-child */
      SCIP_CALL( SCIPcreateConsVarbound(scip, &origcons, differname, origvars1[v], origvars2[v],
            1.0, -SCIPinfinity(scip), 1.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      /* add cons locally to the problem and release it */
      SCIP_CALL( SCIPaddConsNode(scip, childdiffer, origcons, NULL) );
      SCIP_CALL( SCIPreleaseCons(scip, &origcons) );
   }

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpRyanfoster)
{  /*lint --e{715}*/
   SCIPdebugMessage("Execlp method of ryanfoster branching\n");

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextRyanfoster)
{  /*lint --e{715}*/
   SCIP* masterscip;

   SCIP_Bool feasible;
   SCIP_Bool contained;

   SCIP_VAR** branchcands;
   int nbranchcands;

   int v1;
   int v2;
   int o1;
   int o2;
   int j;

   SCIP_VAR* mvar1;
   SCIP_VAR* mvar2;
   SCIP_VAR* ovar1;
   SCIP_VAR* ovar2;

   SCIP_VAR** origvars1;
   SCIP_VAR** origvars2;
   int norigvars1;
   int norigvars2;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execrel method of ryanfoster branching\n");

   *result = SCIP_DIDNOTRUN;

   /* do not perform Ryan & Foster branching if we have neither a set partitioning nor a set covering structure */
   if( !GCGrelaxIsMasterSetCovering(scip) && !GCGrelaxIsMasterSetPartitioning(scip) )
   {
      SCIPdebugMessage("Not executing Ryan&Foster branching, master is neither set covering nor set partitioning\n");
      return SCIP_OKAY;
   }

   /* check whether the current original solution is integral */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), TRUE, TRUE, TRUE, TRUE, &feasible) );
#else
   SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), FALSE, TRUE, TRUE, TRUE, &feasible) );
#endif
   if( feasible )
   {
      SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
         SCIPgetSolOrigObj(scip, GCGrelaxGetCurrentOrigSol(scip)));

      *result = SCIP_CUTOFF;

      return SCIP_OKAY;
   }

   /* the current original solution is not integral, now we have to branch;
    * first, get the master problem and all variables of the master problem
    */
   masterscip = GCGrelaxGetMasterprob(scip);
   SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL) );

   /* now search for two (fractional) columns mvar1, mvar2 in the master and 2 original variables ovar1, ovar2
    * s.t. mvar1 contains both ovar1 and ovar2 and mvar2 contains ovar1, but not ovar2
    */
   ovar1 = NULL;
   ovar2 = NULL;
   mvar1 = NULL;
   mvar2 = NULL;
   feasible = FALSE;

   /* select first fractional column (mvar1) */
   for( v1 = 0; v1 < nbranchcands && !feasible; v1++ )
   {
      mvar1 = branchcands[v1];
      assert(GCGvarIsMaster(mvar1));

      origvars1 = GCGmasterVarGetOrigvars(mvar1);
      norigvars1 = GCGmasterVarGetNOrigvars(mvar1);

      /* select first original variable ovar1, that should be contained in both master variables */
      for( o1 = 0; o1 < norigvars1 && !feasible; o1++ )
      {
         ovar1 = origvars1[o1];
         /* if we deal with a trivial variable, skip it */
         if( SCIPisZero(scip, GCGmasterVarGetOrigvals(mvar1)[o1]) )
            continue;

         /* mvar1 contains ovar1, look for mvar2 which constains ovar1, too */
         for( v2 = v1+1; v2 < nbranchcands && !feasible; v2++ )
         {
            mvar2 = branchcands[v2];
            assert(GCGvarIsMaster(mvar2));

            origvars2 = GCGmasterVarGetOrigvars(mvar2);
            norigvars2 = GCGmasterVarGetNOrigvars(mvar2);

            /* check whether ovar1 is contained in mvar2, too */
            contained = FALSE;
            for( j = 0; j < norigvars2; j++ )
            {
               /* if we deal with a trivial variable, skip it */
               if( SCIPisZero(scip, GCGmasterVarGetOrigvals(mvar2)[j]) )
                  continue;

               if( origvars2[j] == ovar1 )
               {
                  contained = TRUE;
                  break;
               }
            }

            /* mvar2 does not contain ovar1, so look for another mvar2 */
            if( !contained )
               continue;

            /* mvar2 also contains ovar1, now look for ovar2 contained in mvar1, but not in mvar2 */
            for( o2 = 0; o2 < norigvars1; o2++ )
            {
               /* if we deal with a trivial variable, skip it */
               if( !SCIPisZero(scip, GCGmasterVarGetOrigvals(mvar1)[o2]) )
                  continue;

               ovar2 = origvars1[o2];
               if( ovar2 == ovar1 )
                  continue;

               /* check whether ovar2 is contained in mvar2, too */
               contained = FALSE;
               for( j = 0; j < norigvars2; j++ )
               {
                  /* if we deal with a trivial variable, skip it */
                  if( !SCIPisZero(scip, GCGmasterVarGetOrigvals(mvar2)[j]) )
                     continue;

                  if( origvars2[j] == ovar2 )
                  {
                     contained = TRUE;
                     break;
                  }
               }

               /* ovar2 should be contained in mvar1 but not in mvar2, so look for another ovar2,
                * if the current one is contained in mvar2
                */
               if( contained )
                  continue;

               /* if we arrive here, ovar2 is contained in mvar1 but not in mvar2, so everything is fine */
               feasible = TRUE;
               break;
            }

            /* we did not find an ovar2 contained in mvar1, but not in mvar2,
             * now look for one contained in mvar2, but not in mvar1
             */
            if( !feasible )
            {
               for( o2 = 0; o2 < norigvars2; o2++ )
               {
                  /* if we deal with a trivial variable, skip it */
                  if( SCIPisZero(scip, GCGmasterVarGetOrigvals(mvar2)[o2]) )
                     continue;

                  ovar2 = origvars2[o2];

                  if( ovar2 == ovar1 )
                     continue;

                  contained = FALSE;
                  for( j = 0; j < norigvars1; j++ )
                  {
                     /* if we deal with a trivial variable, skip it */
                     if( SCIPisZero(scip, GCGmasterVarGetOrigvals(mvar1)[j]) )
                        continue;
                     if( origvars1[j] == ovar2 )
                     {
                        contained = TRUE;
                        break;
                     }
                  }

                  /* ovar2 should be contained in mvar2 but not in mvar1, so look for another ovar2,
                   * if the current one is contained in mvar1
                   */
                  if( contained )
                     continue;

                  /* if we arrive here, ovar2 is contained in mvar2 but not in mvar1, so everything is fine */
                  feasible = TRUE;
                  break;
               }
            }
         }
      }
   }

   if( !feasible )
   {
      SCIPdebugMessage("Ryanfoster branching rule could not find variables to branch on!\n");
      return SCIP_OKAY;
   }

   /* create the two child nodes in the branch-and-bound tree */
   SCIP_CALL( createChildNodesRyanfoster(scip, branchrule, ovar1, ovar2, GCGvarGetBlock(mvar1)) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsRyanfoster)
{  /*lint --e{715}*/
   SCIP_CONS** origbranchconss;
   GCG_BRANCHDATA* branchdata;
   SCIP_VAR** branchcands;
   SCIP_VAR* ovar1;
   SCIP_VAR* ovar2;
   SCIP_Bool feasible;
   int norigbranchconss;
   int nbranchcands;
   int o1;
   int o2;
   int c;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of ryanfoster branching\n");

   *result = SCIP_DIDNOTRUN;

   /* do not perform Ryan & Foster branching if we have neither a set partitioning nor a set covering structure */
   if( !GCGrelaxIsMasterSetCovering(scip) || !GCGrelaxIsMasterSetPartitioning(scip) )
   {
      SCIPdebugMessage("Not executing Ryanfoster branching, master is neither set covering nor set partitioning\n");
      return SCIP_OKAY;
   }

   /* get unfixed variables and stack of active origbranchconss */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &branchcands, NULL, &nbranchcands) );
   GCGconsOrigbranchGetStack(scip, &origbranchconss, &norigbranchconss);

   ovar1 = NULL;
   ovar2 = NULL;
   feasible = FALSE;

   /* select first original variable ovar1 */
   for( o1 = 0; o1 < nbranchcands && !feasible; ++o1 )
   {
      ovar1 = branchcands[o1];

      /* select second original variable o2 */
      for( o2 = o1 + 1; o2 < nbranchcands; ++o2 )
      {
         ovar2 = branchcands[o2];

         assert(ovar2 != ovar1);

         /* check whether we already branched on this combination of variables */
         for( c = 0; c < norigbranchconss; ++c )
         {
            if( GCGconsOrigbranchGetBranchrule(origbranchconss[c]) == branchrule )
               continue;

            branchdata = GCGconsOrigbranchGetBranchdata(origbranchconss[c]);

            if( (branchdata->var1 == ovar1 && branchdata->var2 == ovar2 )
               || (branchdata->var1 == ovar2 && branchdata->var2 == ovar1) )
            {
               break;
            }
         }

         /* we did not break, so there is no active origbranch constraint with both variables */
         if( c == norigbranchconss )
         {
            feasible = TRUE;
            break;
         }
      }
   }

   if( !feasible )
   {
      SCIPdebugMessage("Ryanfoster branching rule could not find variables to branch on!\n");
      return SCIP_OKAY;
   }

   /* create the two child nodes in the branch-and-bound tree */
   SCIP_CALL( createChildNodesRyanfoster(scip, branchrule, ovar1, ovar2, GCGvarGetBlock(ovar1)) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitRyanfoster)
{
   assert(branchrule != NULL);

   SCIP_CALL( GCGrelaxIncludeBranchrule(scip, branchrule, branchActiveMasterRyanfoster,
         branchDeactiveMasterRyanfoster, branchPropMasterRyanfoster, NULL, branchDataDeleteRyanfoster) );

   return SCIP_OKAY;
}



/* define not used callback as NULL*/
#define branchCopyRyanfoster NULL
#define branchFreeRyanfoster NULL
#define branchExitRyanfoster NULL
#define branchInitsolRyanfoster NULL
#define branchExitsolRyanfoster NULL


/*
 * branching specific interface methods
 */

/** creates the Ryan-Foster LP braching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleRyanfoster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create branching rule data */
   branchruledata = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchCopyRyanfoster,
         branchFreeRyanfoster, branchInitRyanfoster, branchExitRyanfoster, branchInitsolRyanfoster,
         branchExitsolRyanfoster, branchExeclpRyanfoster, branchExecextRyanfoster, branchExecpsRyanfoster,
         branchruledata) );

   return SCIP_OKAY;
}
