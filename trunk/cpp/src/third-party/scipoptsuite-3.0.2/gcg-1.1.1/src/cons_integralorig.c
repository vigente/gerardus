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

/**@file   cons_integralorig.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for enforcing integrality of the transferred master solution in the original problem
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "cons_integralorig.h"
#include "pricer_gcg.h"
#include "cons_masterbranch.h"
#include "pub_gcgvar.h"

#define CONSHDLR_NAME          "integralorig"
#define CONSHDLR_DESC          "integrality constraint"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY      1000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY     1000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/*
 * Callback methods
 */

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpIntegralOrig)
{  /*lint --e{715}*/
   SCIP* origprob;
   SCIP_VAR** origvars;
   int norigvars;
   SCIP_Real solval;
   SCIP_Bool discretization;
   int v;
   int i;

   SCIP_NODE* child1;
   SCIP_NODE* child2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);
   assert(result != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   SCIPdebugMessage("LP solution enforcing method of integralorig constraint\n");

   *result = SCIP_FEASIBLE;

   /* if we use the discretization approach, we do not have to check for integrality of the solution in the
    * original variable space, we obtain it by enforcing integrality of the master solution*/
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( discretization )
   {
      return SCIP_OKAY;
   }

   origvars = SCIPgetOrigVars(origprob);
   norigvars = SCIPgetNOrigVars(origprob);

   /* check for each integral original variable whether it has a fractional value */
   for( v = 0; v < norigvars; v++ )
   {
      SCIP_Real* mastervals;
      SCIP_VAR** mastervars;
      int nmastervars;

      if( SCIPvarGetType(origvars[v]) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      solval = 0;
      assert(GCGvarIsOriginal(origvars[v]));

      mastervals = GCGoriginalVarGetMastervals(origvars[v]);
      mastervars = GCGoriginalVarGetMastervars(origvars[v]);
      nmastervars = GCGoriginalVarGetNMastervars(origvars[v]);

      for( i = 0; i < nmastervars; i++ )
      {
         solval += mastervals[i] * SCIPgetSolVal(scip, NULL, mastervars[i]);
      }
      /* create two children if a variable with fractional value is found */
      if( !SCIPisFeasIntegral(scip, solval) )
      {
         /* create the b&b-tree child-nodes of the current node */
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
   }

   return SCIP_OKAY;

}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsIntegralOrig)
{  /*lint --e{715}*/
   SCIP* origprob;
   SCIP_Bool discretization;

   SCIP_NODE* child1;
   SCIP_NODE* child2;
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);
   assert(result != NULL);

   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   *result = SCIP_FEASIBLE;

   /* if we use the discretization approach, we do not have to check for integrality of the solution in the
    * original variable space, we obtain it by enforcing integrality of the master solution*/
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( discretization )
   {
      return SCIP_OKAY;
   }

   assert(SCIPgetNPseudoBranchCands(origprob) > 0);

   /* create the b&b-tree child-nodes of the current node */
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


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckIntegralOrig)
{  /*lint --e{715}*/
   SCIP* origprob;
   SCIP_VAR** origvars;
   int norigvars;
   SCIP_Real solval;
   SCIP_Bool discretization;
   int v;
   int i;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   SCIPdebugMessage("Check method of integralorig constraint\n");

   *result = SCIP_FEASIBLE;

   /* if we use the discretization approach, we do not have to check for integrality of the solution in the
    * original variable space, we obtain it by enforcing integrality of the master solution*/
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( discretization )
   {
      return SCIP_OKAY;
   }

   origvars = SCIPgetOrigVars(origprob);
   norigvars = SCIPgetNOrigVars(origprob);

   /* check for each integral original variable whether it has a fractional value */
   for( v = 0; v < norigvars && *result == SCIP_FEASIBLE; v++ )
   {
      SCIP_Real* mastervals;
      SCIP_VAR** mastervars;
      int nmastervars;

      if( SCIPvarGetType(origvars[v]) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      solval = 0;
      assert(GCGvarIsOriginal(origvars[v]));

      mastervals = GCGoriginalVarGetMastervals(origvars[v]);
      mastervars = GCGoriginalVarGetMastervars(origvars[v]);
      nmastervars = GCGoriginalVarGetNMastervars(origvars[v]);

      for( i = 0; i < nmastervars; i++ )
      {
         solval += mastervals[i] * SCIPgetSolVal(scip, sol, mastervars[i]);
      }
      if( !SCIPisFeasIntegral(scip, solval) )
      {
         *result = SCIP_INFEASIBLE;

         if( printreason )
         {
            SCIPinfoMessage(scip, NULL, "violation: integrality condition of variable <%s> = %.15g\n",
               SCIPvarGetName(origvars[v]), solval);
         }
      }
   }

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockIntegralOrig)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/* define not used callback as NULL */
#define conshdlrCopyIntegralOrig NULL
#define consFreeIntegralOrig NULL
#define consInitIntegralOrig NULL
#define consExitIntegralOrig NULL
#define consInitpreIntegralOrig NULL
#define consExitpreIntegralOrig NULL
#define consInitsolIntegralOrig NULL
#define consExitsolIntegralOrig NULL
#define consDeleteIntegralOrig NULL
#define consTransIntegralOrig NULL
#define consInitlpIntegralOrig NULL
#define consSepalpIntegralOrig NULL
#define consSepasolIntegralOrig NULL
#define consPropIntegralOrig NULL
#define consPresolIntegralOrig NULL
#define consRespropIntegralOrig NULL
#define consActiveIntegralOrig NULL
#define consDeactiveIntegralOrig NULL
#define consEnableIntegralOrig NULL
#define consDisableIntegralOrig NULL
#define consDelvarsIntegralOrig NULL
#define consPrintIntegralOrig NULL
#define consCopyIntegralOrig NULL
#define consParseIntegralOrig NULL
#define consGetVarsIntegralOrig NULL
#define consGetNVarsIntegralOrig NULL


/*
 * constraint specific interface methods
 */

/** creates the handler for integrality constraint and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrIntegralOrig(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create integral constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         SCIP_PROPTIMING_ALWAYS,
         conshdlrCopyIntegralOrig, consFreeIntegralOrig, consInitIntegralOrig, consExitIntegralOrig,
         consInitpreIntegralOrig, consExitpreIntegralOrig, consInitsolIntegralOrig, consExitsolIntegralOrig,
         consDeleteIntegralOrig, consTransIntegralOrig, consInitlpIntegralOrig,
         consSepalpIntegralOrig, consSepasolIntegralOrig, consEnfolpIntegralOrig, consEnfopsIntegralOrig, consCheckIntegralOrig,
         consPropIntegralOrig, consPresolIntegralOrig, consRespropIntegralOrig, consLockIntegralOrig,
         consActiveIntegralOrig, consDeactiveIntegralOrig,
         consEnableIntegralOrig, consDisableIntegralOrig,
         consDelvarsIntegralOrig, consPrintIntegralOrig, consCopyIntegralOrig, consParseIntegralOrig,
         consGetVarsIntegralOrig, consGetNVarsIntegralOrig,
         conshdlrdata) );

   return SCIP_OKAY;
}
