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

/**@file   cons_masterbranch.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 * @author Martin Bergner
 */

#include <assert.h>
#include <string.h>

#include "cons_masterbranch.h"
#include "scip/cons_linear.h"
#include "cons_origbranch.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "pub_gcgvar.h"

/*#define CHECKPROPAGATEDVARS*/

/* constraint handler properties */
#define CONSHDLR_NAME          "masterbranch"
#define CONSHDLR_DESC          "store branching decision at nodes of the tree constraint handler"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              * propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */


#define EVENTHDLR_NAME         "origvarbound"
#define EVENTHDLR_DESC         "event handler for origvarbound event"


/** constraint data for masterbranch constraints */
struct SCIP_ConsData
{
   int                   propagatedvars;     /**< number of Vars that existed, the last time, the related node was propagated,
                                                  used to determine whether the constraint should be repropagated */
   SCIP_Bool             needprop;           /**< should the constraint be propagated? */
   SCIP_Bool             created;
   SCIP_NODE*            node;               /**< the node at which the cons is sticking */
   SCIP_CONS*            parentcons;         /**< the masterbranch constraint of the parent node */
   SCIP_CONS*            child1cons;         /**< the masterbranch constraint of the first child node */
   SCIP_CONS*            child2cons;         /**< the masterbranch constraint of the second child node */
   SCIP_CONS*            probingtmpcons;     /**< pointer to save the second child if the child2cons pointer is overwritten in probing mode */
   SCIP_CONS*            origcons;           /**< the corresponding origbranch cons in the original program */

   GCG_BRANCHDATA*       branchdata;         /**< branching data stored by the branching rule at the corresponding origcons constraint
                                              *   containing information about the branching restrictions */
   SCIP_BRANCHRULE*      branchrule;         /**< branching rule that created the corresponding node in the original problem and imposed
                                              *   branching restrictions */

   SCIP_VAR**            boundchgvars;       /**< variables of bound changes stored at the current node */
   SCIP_Real*            newbounds;          /**< new bounds for the bound changes stored at the current node */
   SCIP_Real*            oldbounds;          /**< old bounds for the bound changes stored at the current node */
   SCIP_BOUNDTYPE*       boundtypes;         /**< types of the bound changes stored at the current node */

   int*                  nboundchangestreated; /**< number of bound changes of the nodes on the way from the current node to
                                                *   the root node that are treated so far */
   int                   nboundchanges;      /**< number of bound changes */
   int                   nbranchingchanges;  /**< number of bound changes due to branching (<= nboundchanges) */
   int                   nactivated;         /**< number of times the constraint was activated so far */
   char*                 name;               /**< name of the constraint */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_CONS**           stack;              /**< stack for storing active constraints */
   int                   nstack;             /**< number of elements on the stack */
   int                   maxstacksize;       /**< maximum size of the stack */
   SCIP_VAR**            pendingvars;        /**< variables corresponding to pending bound changes (global bound changes) */
   SCIP_BOUNDTYPE*       pendingbndtypes;    /**< types of the pending bound changes (global bound changes) */
   SCIP_Real*            pendingnewbnds;     /**< new bounds corresponding to pending bound changes (global bound changes) */
   SCIP_Real*            pendingoldbnds;     /**< old bounds corresponding to pending bound changes (global bound changes) */
   int                   npendingbnds;       /**< number of pending bound changes (global bound changes) */
   SCIP_Bool             pendingbndsactivated; /**< were pending bound changes already activated? */
   int                   maxpendingbnds;     /**< size of the array corresponding to pending bound changes */
   SCIP_Bool             enforceproper;      /**< should proper variables be enforced? */
};


/*
 * Local methods
 */

/** the function initializes the condata data structure */
static
SCIP_RETCODE GCGconsMasterbranchCreateConsData(
   SCIP*                 scip,               /**< SCIP data structure*/
   SCIP*                 origscip,           /**< Original problem SCIP data structure*/
   SCIP_CONSDATA*        consdata,           /**< consdata to initialize*/
   SCIP_CONSHDLRDATA*    conshdlrData,       /**< constraint handler data */
   SCIP_CONS*            cons                /**< constraint for which the consdata is created */
)
{
#ifdef SCIP_DEBUG
   SCIP_CONS* origcons_parent;
   SCIP_CONS* parent_origcons;
#endif

   SCIP_CONS* origcons;
   SCIP_CONSDATA* parentdata;
   int i;
   int j;

   SCIP_DOMCHG* domchg;
   SCIP_BOUNDCHG* boundchg;

   SCIP_CONSDATA* stackconsdata;
   assert(scip != NULL);
   assert(origscip != NULL);
   assert(consdata != NULL);
   assert(conshdlrData != NULL);
   assert(cons != NULL);

   origcons = GCGconsOrigbranchGetActiveCons(origscip);
   assert(origcons != NULL);

   consdata->origcons = origcons;
   consdata->branchrule = GCGconsOrigbranchGetBranchrule(origcons);
   consdata->branchdata = GCGconsOrigbranchGetBranchdata(origcons);
   GCGconsOrigbranchSetMastercons(origcons, cons);

   SCIP_ALLOC( BMSduplicateBlockMemoryArray(SCIPblkmem(scip), &consdata->name, SCIPconsGetName(consdata->origcons),
         strlen(SCIPconsGetName(consdata->origcons))+1) );

#ifdef SCIP_DEBUG
      if( consdata->parentcons != NULL )
         parent_origcons = SCIPconsGetData(consdata->parentcons)->origcons;
      else
         parent_origcons = NULL;

      if( consdata->origcons != NULL )
         origcons_parent = GCGconsOrigbranchGetParentcons(consdata->origcons);
      else
         origcons_parent = NULL;

      SCIPdebugMessage("cons: %s, origcons: %s, parent: %s => %s\n", SCIPconsGetName(cons), consdata->origcons == NULL? "NULL" : SCIPconsGetName( consdata->origcons ),
         parent_origcons == NULL? "NULL" :  SCIPconsGetName(parent_origcons), origcons_parent == NULL? "NULL" : SCIPconsGetName(origcons_parent) );
#endif

   assert(SCIPgetCurrentNode(scip) == consdata->node || consdata->node == SCIPgetRootNode(scip));
   assert((SCIPgetNNodesLeft(scip)+SCIPgetNNodes(scip) == 1) == (consdata->node == SCIPgetRootNode(scip)));
   assert(SCIPnodeGetDepth(GCGconsOrigbranchGetNode(consdata->origcons)) == SCIPnodeGetDepth(consdata->node));
   assert(consdata->parentcons != NULL || SCIPnodeGetDepth(consdata->node) == 0);

   assert(consdata->parentcons == NULL ||
      SCIPconsGetData(consdata->parentcons)->origcons == GCGconsOrigbranchGetParentcons(consdata->origcons));

   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->nboundchangestreated, conshdlrData->nstack+1) );

   domchg = SCIPnodeGetDomchg(GCGconsOrigbranchGetNode(origcons));
   consdata->nboundchanges = SCIPdomchgGetNBoundchgs(domchg);
   consdata->nboundchangestreated[conshdlrData->nstack] = consdata->nboundchanges;

   if( consdata->nboundchanges > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->boundchgvars, consdata->nboundchanges) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->boundtypes, consdata->nboundchanges) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->newbounds, consdata->nboundchanges) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->oldbounds, consdata->nboundchanges) );
   }

   consdata->nbranchingchanges = 0;

   for( i = 0; i < consdata->nboundchanges; i++ )
   {
      boundchg = SCIPdomchgGetBoundchg(domchg, i);

      consdata->boundchgvars[i] = SCIPboundchgGetVar(boundchg);
      consdata->newbounds[i] = SCIPboundchgGetNewbound(boundchg);
      consdata->boundtypes[i] = SCIPboundchgGetBoundtype(boundchg);

      if( SCIPboundchgGetBoundchgtype(boundchg) == SCIP_BOUNDCHGTYPE_BRANCHING )
      {
         consdata->nbranchingchanges++;
         assert(consdata->nbranchingchanges == i+1);
      }
   }

   consdata->created = TRUE;
   consdata->needprop = TRUE;

   assert((consdata->parentcons == NULL) == (conshdlrData->nstack == 0));
   if( consdata->parentcons != NULL )
   {
      parentdata = SCIPconsGetData(consdata->parentcons);

      assert(consdata->parentcons == conshdlrData->stack[conshdlrData->nstack-1]);
      assert(SCIPconsGetData(conshdlrData->stack[0])->parentcons == NULL);

      /* check whether bound changes were added in nodes on the path to the current node after first activation */
      for( i = 1; i < conshdlrData->nstack; i++ )
      {
         stackconsdata = SCIPconsGetData(conshdlrData->stack[i]);
         domchg = SCIPnodeGetDomchg(GCGconsOrigbranchGetNode(stackconsdata->origcons));

         assert(SCIPdomchgGetNBoundchgs(domchg) >= parentdata->nboundchangestreated[i]);

         if( SCIPdomchgGetNBoundchgs(domchg) != parentdata->nboundchangestreated[i] )
         {
            int diff;

            diff = SCIPdomchgGetNBoundchgs(domchg) - parentdata->nboundchangestreated[i];

            SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->boundchgvars, consdata->nboundchanges + diff) );
            SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->boundtypes, consdata->nboundchanges + diff) );
            SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->newbounds, consdata->nboundchanges + diff) );
            SCIP_CALL( SCIPreallocMemoryArray(scip, &consdata->oldbounds, consdata->nboundchanges + diff) );

            /* add bound changes to the boundchanges array */
            for( j = 0; j < SCIPdomchgGetNBoundchgs(domchg); j++ )
            {
               int bndchgindex;
               boundchg = SCIPdomchgGetBoundchg(domchg, j);

               if( j < stackconsdata->nboundchangestreated[i] )
               {
                  assert(stackconsdata->boundchgvars[j] == SCIPboundchgGetVar(boundchg)
                     && SCIPisEQ(scip, stackconsdata->newbounds[j], SCIPboundchgGetNewbound(boundchg))
                     && stackconsdata->boundtypes[j] == SCIPboundchgGetBoundtype(boundchg));
                  continue;
               }
               if( j < parentdata->nboundchangestreated[i] )
                  continue;

               bndchgindex = consdata->nboundchanges + j - parentdata->nboundchangestreated[i];

               consdata->boundchgvars[bndchgindex] = SCIPboundchgGetVar(boundchg);
               consdata->newbounds[bndchgindex] = SCIPboundchgGetNewbound(boundchg);
               consdata->boundtypes[bndchgindex] = SCIPboundchgGetBoundtype(boundchg);

            }
            consdata->nboundchanges += diff;
         }
         consdata->nboundchangestreated[i] = SCIPdomchgGetNBoundchgs(domchg);
      }
   }

   return SCIP_OKAY;
}

/** local method performs the bound change on pricing variables */
static
SCIP_RETCODE GCGresetPricingVarBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origscip,           /**< original problem SCIP data structure */
   SCIP_VAR*             pricingvar,         /**< variable on which the bound change should be performed */
   SCIP_CONSDATA*        consdata,           /**< constraint data structure */
   int                   i,                  /**< index of the information in the constraint data structure */
   int                   blocknr             /**< number of the pricing problem */
   )
{
   assert( scip != NULL );
   assert(origscip != NULL);
   assert(pricingvar != NULL);
   assert(consdata != NULL);
   assert(consdata->created);
   assert(consdata->boundchgvars != NULL);
   assert(consdata->boundtypes != NULL);
   assert(consdata->newbounds != NULL);
   assert(consdata->oldbounds != NULL);
   assert(blocknr >= 0 && blocknr < GCGrelaxGetNPricingprobs(origscip));

   /* lower bound was changed */
   if( consdata->boundtypes[i] == SCIP_BOUNDTYPE_LOWER )
   {
      assert( SCIPisEQ(scip, SCIPvarGetLbLocal(pricingvar), consdata->newbounds[i])
           || SCIPisLE(scip, SCIPvarGetLbLocal(pricingvar), SCIPvarGetLbGlobal(consdata->boundchgvars[i])));

      if( SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->boundchgvars[i]), consdata->newbounds[i]) )
         return SCIP_OKAY;

      if( SCIPisGT(scip, SCIPvarGetLbGlobal(consdata->boundchgvars[i]), consdata->oldbounds[i]) )
      {
         SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, blocknr), pricingvar, SCIPvarGetLbGlobal(consdata->boundchgvars[i])) );
         SCIPdebugMessage("relaxed lower bound of pricing var %s from %g to global bound %g (%s)\n",
            SCIPvarGetName(pricingvar), consdata->newbounds[i], SCIPvarGetLbGlobal(consdata->boundchgvars[i]), consdata->name);
      }
      else
      {
         SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, blocknr), pricingvar, consdata->oldbounds[i]) );
         SCIPdebugMessage("relaxed lower bound of pricing var %s from %g to %g (%s)\n",
            SCIPvarGetName(pricingvar), consdata->newbounds[i], consdata->oldbounds[i], consdata->name);
      }
   }
   /* upper bound was changed */
   else
   {
      assert( SCIPisEQ(scip, SCIPvarGetUbLocal(pricingvar), consdata->newbounds[i])
           || SCIPisGE(scip, SCIPvarGetUbLocal(pricingvar), SCIPvarGetUbGlobal(consdata->boundchgvars[i])));

      if( SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->boundchgvars[i]), consdata->newbounds[i]) )
         return SCIP_OKAY;

      if( SCIPisLT(scip, SCIPvarGetUbGlobal(consdata->boundchgvars[i]), consdata->oldbounds[i]) )
      {
         SCIP_CALL( SCIPchgVarUb(GCGrelaxGetPricingprob(origscip, blocknr), pricingvar, SCIPvarGetUbGlobal(consdata->boundchgvars[i])) );
         SCIPdebugMessage("relaxed upper bound of pricing var %s from %g to global bound %g (%s)\n",
            SCIPvarGetName(pricingvar), consdata->newbounds[i], SCIPvarGetUbGlobal(consdata->boundchgvars[i]), consdata->name);
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(GCGrelaxGetPricingprob(origscip, blocknr), pricingvar, consdata->oldbounds[i]) );
         SCIPdebugMessage("relaxed upper bound of pricing var %s from %g to %g (%s)\n",
            SCIPvarGetName(pricingvar), consdata->newbounds[i], consdata->oldbounds[i], consdata->name);
      }
   }
   return SCIP_OKAY;
}

/** local method performs the bound change on pricing variables */
static
SCIP_RETCODE GCGtightenPricingVarBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origscip,           /**< original problem SCIP data structure */
   SCIP_VAR*             pricingvar,         /**< variable on which the bound change should be performed */
   SCIP_CONSDATA*        consdata,           /**< constraint data structure */
   int                   i,                  /**< index of the information in the constraint data structure */
   int                   blocknr             /**< number of the pricing problem */
   )
{
   assert( scip != NULL );
   assert(origscip != NULL);
   assert(pricingvar != NULL);
   assert(consdata != NULL);
   assert(consdata->created);
   assert(consdata->boundchgvars != NULL);
   assert(consdata->boundtypes != NULL);
   assert(consdata->newbounds != NULL);
   assert(consdata->oldbounds != NULL);
   assert(blocknr >= 0 && blocknr < GCGrelaxGetNPricingprobs(origscip));

   /* lower bound was changed */
   if( consdata->boundtypes[i] == SCIP_BOUNDTYPE_LOWER )
   {
      consdata->oldbounds[i] = SCIPvarGetLbLocal(pricingvar);

      SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, blocknr), pricingvar, consdata->newbounds[i]) );
      SCIPdebugMessage("tightened lower bound of var %s from %g to %g\n",
         SCIPvarGetName(pricingvar), consdata->oldbounds[i], consdata->newbounds[i]);
   }
   /* upper bound was changed */
   else
   {
      consdata->oldbounds[i] = SCIPvarGetUbLocal(pricingvar);

      SCIP_CALL( SCIPchgVarUb(GCGrelaxGetPricingprob(origscip, blocknr), pricingvar, consdata->newbounds[i]) );
      SCIPdebugMessage("tightened upper bound of var %s from %g to %g\n",
         SCIPvarGetName(pricingvar), consdata->oldbounds[i], consdata->newbounds[i]);
   }
   return SCIP_OKAY;
}

/* local method to add a global bound change to the pending bound changes array */
static
SCIP_RETCODE GCGconsMasterbranchAddPendingBndChg(
   SCIP*                 scip,
   SCIP_VAR*             var,
   SCIP_BOUNDTYPE        boundtype,
   SCIP_Real             oldbound,
   SCIP_Real             newbound
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert((conshdlrData->npendingbnds > 0) || conshdlrData->pendingbndsactivated);

   /* realloc memory if needed */
   if( conshdlrData->npendingbnds >= conshdlrData->maxpendingbnds )
   {
      conshdlrData->maxpendingbnds = conshdlrData->npendingbnds+5;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->pendingvars), conshdlrData->maxpendingbnds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->pendingbndtypes), conshdlrData->maxpendingbnds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->pendingoldbnds), conshdlrData->maxpendingbnds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->pendingnewbnds), conshdlrData->maxpendingbnds) );
   }

   /* store pending bound change */
   conshdlrData->pendingvars[conshdlrData->npendingbnds] = var;
   conshdlrData->pendingbndtypes[conshdlrData->npendingbnds] = boundtype;
   conshdlrData->pendingoldbnds[conshdlrData->npendingbnds] = oldbound;
   conshdlrData->pendingnewbnds[conshdlrData->npendingbnds] = newbound;
   conshdlrData->npendingbnds++;
   conshdlrData->pendingbndsactivated = FALSE;

   return SCIP_OKAY;
}

/*
 * Callback methods
 */


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   SCIPdebugMessage("freeing masterbranch constraint handler\n");

   /* free constraint handler storage */
   assert(conshdlrData->stack == NULL);
   SCIPfreeMemory(scip, &conshdlrData);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitMasterbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   SCIPdebugMessage("consInitMasterbranch()\n");

   /* prepare stack */
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrData->stack, conshdlrData->maxstacksize) );
   conshdlrData->nstack = 0;

   /* prepare pending bound changes */
   conshdlrData->npendingbnds = 0;
   conshdlrData->maxpendingbnds = 5;
   conshdlrData->pendingbndsactivated = TRUE;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrData->pendingvars), conshdlrData->maxpendingbnds) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrData->pendingbndtypes), conshdlrData->maxpendingbnds) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrData->pendingoldbnds), conshdlrData->maxpendingbnds) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(conshdlrData->pendingnewbnds), conshdlrData->maxpendingbnds) );

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolMasterbranch)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   SCIPdebugMessage("consInitsolMasterbranch()\n");

   /* create masterbranch constraint for the root node */

   SCIP_CALL( GCGcreateConsMasterbranch(scip, &cons, NULL, NULL) );
   GCGconsOrigbranchSetMastercons(GCGconsOrigbranchGetActiveCons(GCGpricerGetOrigprob(scip)), cons);

   conshdlrData->nstack = 1;
   conshdlrData->stack[0] = cons;

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitMasterbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->nstack == 1 || SCIPgetNNodes(scip) == 0);

   SCIPdebugMessage("exiting masterbranch constraint handler\n");

   /* free stack */
   SCIPfreeMemoryArray(scip, &(conshdlrData->stack));
   SCIPfreeMemoryArray(scip, &(conshdlrData->pendingvars));
   SCIPfreeMemoryArray(scip, &(conshdlrData->pendingbndtypes));
   SCIPfreeMemoryArray(scip, &(conshdlrData->pendingoldbnds));
   SCIPfreeMemoryArray(scip, &(conshdlrData->pendingnewbnds));

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteMasterbranch)
{
   SCIP_CONSDATA* consdata2;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   SCIPdebugMessage("Deleting masterbranch constraint: <%s>.\n", (*consdata)->name);

   /* set the mastercons pointer of the corresponding origcons to NULL */
   if( (*consdata)->origcons != NULL )
      GCGconsOrigbranchSetMastercons((*consdata)->origcons, NULL);

   /* set the pointer in the parent node to NULL */
   if( (*consdata)->parentcons != NULL )
   {
      consdata2 = SCIPconsGetData((*consdata)->parentcons);
      if( consdata2->child1cons == cons )
      {
         consdata2->child1cons = NULL;
      }
      else
      {
         if( consdata2->child2cons == cons )
         {
            consdata2->child2cons = NULL;

            if( SCIPinProbing(scip) )
            {
               consdata2->child2cons = consdata2->probingtmpcons;
               consdata2->probingtmpcons = NULL;
            }
         }
         else
         {
            assert(SCIPinProbing(scip));
            assert(consdata2->probingtmpcons == cons);
            assert(SCIPisLE(scip, SCIPgetCutoffbound(scip), SCIPgetNodeLowerbound(scip, (*consdata)->node)));

            consdata2->probingtmpcons = NULL;
         }
      }
   }
   /* the node should not have children anymore */
   assert((*consdata)->child1cons == NULL);
   assert((*consdata)->child2cons == NULL);

   /* delete branchdata, if the corresponding origcons was already deleted, otherwise, it will be deleted by the
    * corresponding origbranch constraint */
   if( (*consdata)->origcons == NULL && (*consdata)->branchdata != NULL )
   {
      SCIP_CALL( GCGrelaxBranchDataDelete(GCGpricerGetOrigprob(scip), (*consdata)->branchrule, &(*consdata)->branchdata) );
   }

   /* delete array with bound changes */
   if( (*consdata)->nboundchanges > 0 )
   {
      SCIPfreeMemoryArray(scip, &(*consdata)->oldbounds);
      SCIPfreeMemoryArray(scip, &(*consdata)->newbounds);
      SCIPfreeMemoryArray(scip, &(*consdata)->boundtypes);
      SCIPfreeMemoryArray(scip, &(*consdata)->boundchgvars);
   }

   if( (*consdata)->nboundchangestreated != NULL )
   {
      SCIPfreeMemoryArray(scip, &(*consdata)->nboundchangestreated);
   }

   /* free constraint data */
   if( (*consdata)->name != NULL )
   {
      BMSfreeBlockMemoryArray(SCIPblkmem(scip), &(*consdata)->name, strlen((*consdata)->name)+1);
   }
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveMasterbranch)
{
   SCIP* origscip;
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA* consdata;
   int i;


   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert((consdata->node == NULL) == (consdata->parentcons == NULL) );

   if( consdata->node == NULL )
      consdata->node = SCIPgetRootNode(scip);

   assert(consdata->node != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   consdata->nactivated++;

   SCIPdebugMessage("Activating ");
   /* if the node is activated the first time, we first have to setup the constraint data */
   if( !consdata->created )
   {
      SCIPdebugPrintf("for the first time\n");
      SCIP_CALL( GCGconsMasterbranchCreateConsData(scip, origscip, consdata, conshdlrData, cons) );
      assert(consdata->created);
   }
   else
      SCIPdebugPrintf("\n");

   /* the node has to be repropagated if new variables were created after the node was left the last time
    * or if new bound changes on directly transferred variables were found */
   assert(GCGpricerGetNPricedvars(scip) >= consdata->propagatedvars);
   if( GCGpricerGetNPricedvars(scip) > consdata->propagatedvars || GCGconsOrigbranchGetNPropBoundChgs(origscip, consdata->origcons) > 0 )
   {
      consdata->needprop = TRUE;
      SCIP_CALL( SCIPrepropagateNode(scip, consdata->node) );
   }

   if( consdata->nboundchanges - consdata->nboundchangestreated[conshdlrData->nstack] > 0 )
      SCIPdebugMessage("added %d boundchanges from previous nodes!\n", consdata->nboundchanges - consdata->nboundchangestreated[conshdlrData->nstack]);

   /* put constraint on the stack */
   if( conshdlrData->nstack >= conshdlrData->maxstacksize )
   {
      conshdlrData->maxstacksize = 2*(conshdlrData->maxstacksize);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->stack), conshdlrData->maxstacksize) );
      SCIPdebugMessage("reallocating Memory for stack! %d --> %d\n", conshdlrData->maxstacksize/2, conshdlrData->maxstacksize);
   }

   conshdlrData->stack[conshdlrData->nstack] = cons;
   (conshdlrData->nstack)++;

   SCIPdebugMessage("Activating masterbranch constraint: <%s> [stack size: %d], needprop = %u.\n",
      consdata->name, conshdlrData->nstack, consdata->needprop);

   /* apply global bound changes in the original problem to the master problem */
   if( !conshdlrData->pendingbndsactivated )
   {
      assert(conshdlrData->npendingbnds > 0);
      for( i = 0; i < conshdlrData->npendingbnds; i++ )
      {
         /* This should not have an effect on linking variables */
         assert(GCGvarIsMaster(conshdlrData->pendingvars[i]) || GCGvarIsPricing(conshdlrData->pendingvars[i]));

         if( GCGvarIsMaster(conshdlrData->pendingvars[i]) )
         {
            if( conshdlrData->pendingbndtypes[i] == SCIP_BOUNDTYPE_LOWER )
            {
               if( !SCIPisEQ(scip, SCIPvarGetLbGlobal(conshdlrData->pendingvars[i]), conshdlrData->pendingnewbnds[i]) )
               {
                  SCIPdebugMessage("Global lower bound of var <%s> set to %g\n", SCIPvarGetName(conshdlrData->pendingvars[i]),
                     conshdlrData->pendingnewbnds[i]);
                  SCIP_CALL( SCIPchgVarLbGlobal(scip, conshdlrData->pendingvars[i], conshdlrData->pendingnewbnds[i]) );
               }
            }
            else
            {
               if( !SCIPisEQ(scip, SCIPvarGetUbGlobal(conshdlrData->pendingvars[i]), conshdlrData->pendingnewbnds[i]) )
               {
                  SCIPdebugMessage("Global upper bound of var <%s> set to %g\n", SCIPvarGetName(conshdlrData->pendingvars[i]),
                     conshdlrData->pendingnewbnds[i]);
                  SCIP_CALL( SCIPchgVarUbGlobal(scip, conshdlrData->pendingvars[i], conshdlrData->pendingnewbnds[i]) );
               }
            }
         }
         else
         {
            /* this is a global boundchange on a variable that belongs to a block,
             * we have to adjust the bound of the corresponding variable in the pricing problem */
            if( conshdlrData->pendingbndtypes[i] == SCIP_BOUNDTYPE_LOWER )
            {
               SCIP_CALL( SCIPchgVarLb(GCGrelaxGetPricingprob(origscip, GCGvarGetBlock(conshdlrData->pendingvars[i]) ),
                     conshdlrData->pendingvars[i], conshdlrData->pendingnewbnds[i]) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarUbGlobal(GCGrelaxGetPricingprob(origscip, GCGvarGetBlock(conshdlrData->pendingvars[i]) ),
                     conshdlrData->pendingvars[i], conshdlrData->pendingnewbnds[i]) );
            }
         }
      }
      conshdlrData->pendingbndsactivated = TRUE;
   }


   /* apply local bound changes in the original problem to the pricing problems */
   for( i = 0; i < consdata->nboundchanges; i++ )
   {
      assert(GCGvarIsOriginal(consdata->boundchgvars[i]));
      assert(GCGvarGetBlock(consdata->boundchgvars[i]) < GCGrelaxGetNPricingprobs(origscip));

      /* if variable belongs to no block, skip it here because the bound changes are treated in the propagation */
      if( GCGvarGetBlock(consdata->boundchgvars[i]) == -1 )
         continue;

      else if( GCGvarGetBlock(consdata->boundchgvars[i]) >= 0 )
      {
         SCIPdebugMessage("adjusting bound of pricing var %s\n", SCIPvarGetName(consdata->boundchgvars[i]));
         /* set corresponding bound in the pricing problem */
         SCIP_CALL(GCGtightenPricingVarBound(scip, origscip,
               GCGoriginalVarGetPricingVar(consdata->boundchgvars[i]), consdata, i, GCGvarGetBlock(consdata->boundchgvars[i])));
      }
      else if( GCGvarGetBlock(consdata->boundchgvars[i]) == -2 )
      {
         int j;
         int npricingprobs;
         SCIP_VAR** pricingvars;

         npricingprobs = GCGrelaxGetNPricingprobs(origscip);
         pricingvars = GCGlinkingVarGetPricingVars(consdata->boundchgvars[i]);
         SCIPdebugMessage("adjusting bound of linking pricing var %s\n", SCIPvarGetName(consdata->boundchgvars[i]));
         /* set corresponding bound in the pricing problem */
         for( j = 0; j < npricingprobs; ++j )
         {
            if( pricingvars[j] == NULL )
               continue;

            SCIP_CALL( GCGtightenPricingVarBound(scip, origscip, pricingvars[j], consdata, i, j) );
         }
      }
      else
      {
         SCIPerrorMessage("blocknr = %d is not valid! This is a serious error!", GCGvarGetBlock(consdata->boundchgvars[i]));
         SCIPABORT();
      }
   }


   /* call branching specific activation method */
   if( consdata->branchrule != NULL )
   {
      SCIP_CALL( GCGrelaxBranchActiveMaster(origscip, consdata->branchrule, consdata->branchdata) );
   }

   return SCIP_OKAY;
}

/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveMasterbranch)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA* consdata;
   SCIP* origscip;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL || conshdlrData->nstack == 1);
   assert(conshdlrData->nstack > 0);
   assert(conshdlrData->nstack == 1 || cons == conshdlrData->stack[conshdlrData->nstack-1]);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->created);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   if( !conshdlrData->pendingbndsactivated )
   {
      SCIPdebugMessage("We need repropagation\n");
      consdata->needprop = TRUE;
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      consdata->propagatedvars = GCGpricerGetNPricedvars(scip);

   /* remove constraint from the stack */
   (conshdlrData->nstack)--;


   SCIPdebugMessage("Deactivating masterbranch constraint: <%s> [stack size: %d].\n",
      consdata->name, conshdlrData->nstack);

   /* undo local bound changes in the original problem to the pricing problems */
   for( i = consdata->nboundchanges - 1; i >= 0; i-- )
   {
      int blocknr;
      blocknr = GCGvarGetBlock(consdata->boundchgvars[i]);
      assert(GCGvarIsOriginal(consdata->boundchgvars[i]));
      assert(blocknr < GCGrelaxGetNPricingprobs(origscip));

      /* if variable belongs to no block, local bound in master was set, is reset automatically */
      if( blocknr == -1 )
         continue;

      else if( blocknr >= 0 )
      {
         assert(GCGrelaxGetPricingprob(origscip, GCGvarGetBlock(consdata->boundchgvars[i])) != NULL);

         /* reset corresponding bound in the pricing problem */

         SCIP_CALL(GCGresetPricingVarBound(scip, origscip,
               GCGoriginalVarGetPricingVar(consdata->boundchgvars[i]), consdata, i, blocknr));
      }
      else if( blocknr == -2 )
      {
         int j;
         SCIP_VAR** pricingvars;
         int npricingprobs;

         /* if the variable is linking, we have to perform the same step as above for every existing block*/
         assert(GCGvarIsLinking(consdata->boundchgvars[i]));
         pricingvars = GCGlinkingVarGetPricingVars(consdata->boundchgvars[i]);
         npricingprobs = GCGrelaxGetNPricingprobs(origscip);

         /* reset corresponding bound in the pricing problem */
         /* lower bound was changed */
         for( j = 0; j < npricingprobs; ++j )
         {
            assert(GCGrelaxGetPricingprob(origscip, j) != NULL);
            if( pricingvars[j] == NULL )
               continue;

            assert(GCGrelaxGetPricingprob(origscip, j) != NULL);

            /* reset corresponding bound in the pricing problem */
            SCIP_CALL( GCGresetPricingVarBound(scip, origscip, pricingvars[j], consdata, i, j) );
         }
      }
      else
      {
         SCIPerrorMessage("blocknr = %d is not valid! This is a serious error!", blocknr);
         SCIPABORT();
      }
   }

   /* call branching specific deactivation method */
   if( consdata->branchrule != NULL )
   {
      SCIP_CALL( GCGrelaxBranchDeactiveMaster(GCGpricerGetOrigprob(scip), consdata->branchrule, consdata->branchdata) );
   }

   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropMasterbranch)
{  /*lint --e{715}*/
   SCIP* origscip;
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Real val;

   int propcount;
   int i;
   int j;
   int k;

   SCIP_VAR** vars;
   int nvars;
   int nboundchanges;

   SCIP_VAR** propvars;                      /**< original variable for which the propagation found domain reductions */
   SCIP_BOUNDTYPE* propboundtypes;           /**< type of the domain new bound found by propagation */
   SCIP_Real* propbounds;                    /**< new lower/upper bound of the propagated original variable */
   int npropbounds;
   SCIP_VAR* mastervar;


   assert(conshdlr != NULL);
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   *result = SCIP_DIDNOTRUN;

   /* the constraint data of the cons related to the current node */
   cons = conshdlrData->stack[conshdlrData->nstack-1];
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !consdata->needprop && GCGconsOrigbranchGetNPropBoundChgs(origscip, consdata->origcons) == 0 )
   {
      SCIPdebugMessage("No propagation of masterbranch constraint needed: <%s>, stack size = %d.\n",
         consdata->name, conshdlrData->nstack);

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   vars = GCGpricerGetPricedvars(scip);
   nvars = GCGpricerGetNPricedvars(scip);

   SCIPdebugMessage("Starting propagation of masterbranch constraint: <%s>, stack size = %d, newvars = %d, npendingbnds = %d, npropbounds = %d.\n",
      consdata->name, conshdlrData->nstack, nvars - consdata->propagatedvars, conshdlrData->npendingbnds, GCGconsOrigbranchGetNPropBoundChgs(origscip, consdata->origcons));

   *result = SCIP_DIDNOTFIND;

   propcount = 0;

   /* propagate all bound changes or only the branching bound changes, depending on the setting for the enforcement of proper variables */
   nboundchanges = (conshdlrData->enforceproper ? consdata->nboundchanges : consdata->nbranchingchanges);

   assert((conshdlrData->npendingbnds > 0) || conshdlrData->pendingbndsactivated);

   /* iterate over all master variables and apply global bound changes */
   if( conshdlrData->npendingbnds > 0 && conshdlrData->pendingbndsactivated )
   {
      for( i = 0; i < nvars; i++ )
      {
         SCIP_Bool ismastervariablerelevant;
         SCIP_VAR** origvars;
         SCIP_Real* origvals;
         int norigvars;
         int blocknr;
         blocknr = GCGvarGetBlock(vars[i]);

         assert(GCGvarIsMaster(vars[i]));
         norigvars = GCGmasterVarGetNOrigvars(vars[i]);
         origvars = GCGmasterVarGetOrigvars(vars[i]);
         origvals = GCGmasterVarGetOrigvals(vars[i]);

         assert(blocknr < GCGrelaxGetNPricingprobs(origscip));
         assert(norigvars >= 0);
         assert(origvars != NULL || norigvars == 0);

         /* only look at master variables not globally fixed to zero that belong to a block */
         ismastervariablerelevant = !SCIPisFeasZero(scip, SCIPvarGetUbGlobal(vars[i]));
         ismastervariablerelevant = ismastervariablerelevant && (norigvars > 0);
         ismastervariablerelevant = ismastervariablerelevant && (blocknr >= 0 || GCGvarIsLinking(origvars[0])); /*lint !e613*/
         if( !ismastervariablerelevant )
            continue;

         /* iterate over global bound changes that were not yet checked for the master variables */
         for( k = 0; k < conshdlrData->npendingbnds; k++ )
         {
            SCIP_Bool ismastervarrelevant;
            SCIP_VAR** pricingvars;
            int bndchgblocknr;
            SCIP_VAR** bndchgorigvars;

            assert(!GCGvarIsOriginal(conshdlrData->pendingvars[k]));

            bndchgblocknr = GCGvarGetBlock(conshdlrData->pendingvars[k]);
            if( GCGvarIsMaster(conshdlrData->pendingvars[k]) )
               bndchgorigvars = GCGmasterVarGetOrigvars(conshdlrData->pendingvars[k]);
            else if( GCGvarIsPricing(conshdlrData->pendingvars[k]) )
               bndchgorigvars = GCGpricingVarGetOrigvars(conshdlrData->pendingvars[k]);
            else
            {
               SCIPerrorMessage("Variable %s is not pricing nor master.\n", SCIPvarGetName(conshdlrData->pendingvars[k]));
               assert(GCGvarIsMaster(conshdlrData->pendingvars[k]) || GCGvarIsPricing(conshdlrData->pendingvars[k]));
               bndchgorigvars = NULL;
            }
            assert(bndchgblocknr < GCGrelaxGetNPricingprobs(origscip));
            assert(bndchgorigvars != NULL);

            /* LINK: mb: this might work  */
            /* the boundchange was performed on a variable in another block, continue */

            /* if we are not dealing with a linking variable, skip the master variable if it is useless */
            ismastervarrelevant = (bndchgblocknr == blocknr);

            /* if we are dealing with a linking master variable but it has nothing to do with the
             * boundchangevar's block, skip it, too */
            if( origvars != NULL )
            {
               if( GCGvarIsLinking(origvars[0]) )
               {
                  pricingvars = GCGlinkingVarGetPricingVars(origvars[0]);
                  ismastervarrelevant = ismastervarrelevant || (pricingvars[bndchgblocknr] != NULL);
               }
            }

            if( !ismastervarrelevant )
               continue;

            assert(bndchgorigvars[0] != NULL);
            /* val is the value of the branching variable in the current mastervar,
             * we set it to 0.0, since variables with 0 coefficient are not stored in the origvars array,
             * if we do not find the branching variable in this array, it has value 0.0 */
            val = 0.0;

            for( j = 0; j < norigvars; j++ )
            {
               /* Make sure that the original variable and the master variable belong to the same block
                * or that, in case of linking variables, the linking variable is in that block */
               assert(GCGvarGetBlock(origvars[j]) == blocknr || (GCGisLinkingVarInBlock(origvars[j], blocknr))); /*lint !e613*/

               /* check whether the original variable contained in the master variable equals the variable
                * on which the current branching was performed */
               if( origvars[j] == bndchgorigvars[0] ) /*lint !e613*/
               {
                  val = origvals[j];
                  break;
               }
            }

            /* if the variable contains a part of the branching variable that violates the bound,
             * fix the master variable to 0 */

            /* branching imposes new lower bound */
            if( conshdlrData->pendingbndtypes[k] == SCIP_BOUNDTYPE_LOWER &&
               SCIPisFeasLT(scip, val, conshdlrData->pendingnewbnds[k]) )
            {
               SCIP_CALL( SCIPchgVarUbGlobal(scip, vars[i], 0.0) );
               propcount++;
               break;
            }
            /* branching imposes new upper bound */
            if( conshdlrData->pendingbndtypes[k] == SCIP_BOUNDTYPE_UPPER &&
               SCIPisFeasGT(scip, val, conshdlrData->pendingnewbnds[k]) )
            {
               SCIP_CALL( SCIPchgVarUbGlobal(scip, vars[i], 0.0) );
               propcount++;
               break;
            }
         }
      }
      conshdlrData->pendingbndsactivated = TRUE;
      conshdlrData->npendingbnds = 0;

      SCIPdebugMessage("Finished handling of pending global bound changes: %d changed bounds\n", propcount);
   }

   /* iterate over all master variables created after the current node was left the last time */
   for( i = consdata->propagatedvars; i < nvars; i++ )
   {
      SCIP_VAR** origvars;
      int norigvars;
      SCIP_Real* origvals;
      int blocknr;

      assert(GCGvarIsMaster(vars[i]));
      blocknr = GCGvarGetBlock(vars[i]);
      assert(blocknr >= -1 && blocknr < GCGrelaxGetNPricingprobs(origscip));

      origvals = GCGmasterVarGetOrigvals(vars[i]);
      norigvars = GCGmasterVarGetNOrigvars(vars[i]);
      origvars = GCGmasterVarGetOrigvars(vars[i]);
      /** @todo check if this really works with linking variables */

      /* only look at variables not already fixed to 0 or that belong to no block */
      if( (SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[i]))) && blocknr >= 0 )
         continue;

      /* the variable was copied from original to master */
      /** @todo This code might never be executed as the vars array only contains variables generated
        * during pricing. We might want to check that with an assert */
      if( blocknr == -1 )
      {
         /* iterate over bound changes performed at the current node's equivalent in the original tree */
         for( k = 0; k < nboundchanges; k++ )
         {
            assert(SCIPisFeasEQ(scip, origvals[0], 1.0));
            if( origvars[0] == consdata->boundchgvars[k] )
            {
               /* branching imposes new lower bound */
               if( consdata->boundtypes[k] == SCIP_BOUNDTYPE_LOWER
                  && SCIPisGT(scip, consdata->newbounds[k], SCIPvarGetLbLocal(vars[i])) )
               {
                  SCIP_CALL( SCIPchgVarLb(scip, vars[i], consdata->newbounds[k]) );
                  propcount++;
               }
               /* branching imposes new upper bound */
               if( consdata->boundtypes[k] == SCIP_BOUNDTYPE_UPPER
                  && SCIPisLT(scip, consdata->newbounds[k], SCIPvarGetUbLocal(vars[i])) )
               {
                  SCIP_CALL( SCIPchgVarUb(scip, vars[i], consdata->newbounds[k]) );
                  propcount++;
               }
            }
         }
      }
      else
      {
         /* iterate over bound changes performed at the current node's equivalent in the original tree */
         for( k = 0; k < nboundchanges; k++ )
         {
#ifdef SCIP_DEBUG
            SCIP_Bool contained = FALSE;
            SCIP_Bool handled = FALSE;
#endif
            int bndchgblocknr;

            /* get the block the original variable is in */
            bndchgblocknr = GCGvarGetBlock(consdata->boundchgvars[k]);
            assert(GCGvarIsOriginal(consdata->boundchgvars[k]));
            assert(bndchgblocknr < GCGrelaxGetNPricingprobs(origscip));

            /* ignore master variables that contain no original variables */
            /** @todo move this statement one for loop higher? */
            if( origvars == NULL || origvars[0] == NULL )
               continue;

            /* the boundchange was performed on a variable in another block, continue */
            if( (!GCGvarIsLinking(consdata->boundchgvars[k]) && bndchgblocknr != blocknr) ||
               (GCGvarIsLinking(consdata->boundchgvars[k]) && !GCGisLinkingVarInBlock(consdata->boundchgvars[k], blocknr)) )
               continue;

            assert(bndchgblocknr != -1);

            /* val is the value of the branching variable in the current mastervar,
             * we set it to 0.0, since variables with 0 coefficient are not stored in the origvars array,
             * if we do not find the branching variable in this array, it has value 0.0 */
            val = 0.0;

            /* iterate over all original variables contained in the current master variable */
            for( j = 0; j < norigvars; j++ )
            {
               assert(GCGvarGetBlock(origvars[j]) == blocknr || GCGisLinkingVarInBlock(origvars[j], blocknr));

               /* check whether the original variable contained in the master variable equals the variable
                * on which the current branching was performed */
               if( origvars[j] == consdata->boundchgvars[k] )
               {
#ifdef SCIP_DEBUG
                  contained = TRUE;
#endif
                  val = origvals[j];
                  break;
               }
            }

            /* if the variable contains a part of the branching variable that violates the bound,
             * fix the master variable to 0 */

            /* branching imposes new lower bound */
            if( consdata->boundtypes[k] == SCIP_BOUNDTYPE_LOWER && SCIPisFeasLT(scip, val, consdata->newbounds[k]) )
            {
               SCIPdebugMessage("Changing lower bound of var %s\n", SCIPvarGetName(vars[i]));
               SCIP_CALL( SCIPchgVarUb(scip, vars[i], 0.0) );
               propcount++;
#ifdef SCIP_DEBUG
               handled = TRUE;
#endif
               break;
            }
            /* branching imposes new upper bound */
            if( consdata->boundtypes[k] == SCIP_BOUNDTYPE_UPPER && SCIPisFeasGT(scip, val, consdata->newbounds[k]) )
            {
               SCIPdebugMessage("Changing upper bound of var %s\n", SCIPvarGetName(vars[i]));
               SCIP_CALL( SCIPchgVarUb(scip, vars[i], 0.0) );
               propcount++;


#ifdef SCIP_DEBUG
               handled = TRUE;
#endif
               break;
            }
#ifdef SCIP_DEBUG
            if( contained || !handled )
            {
               SCIPdebugMessage("orig var %s is contained in %s but not handled val = %f \n", SCIPvarGetName(consdata->boundchgvars[k]), SCIPvarGetName(vars[i]), val);

            }
            assert(j == norigvars || contained);
#endif
         }
      }
   }
   SCIPdebugMessage("Finished propagation of newly created variables: %d changed bounds\n", propcount);

   /* get local bound changes on variables directly transferred to the master problem and apply them */
   SCIP_CALL( GCGconsOrigbranchGetPropBoundChgs(origscip, consdata->origcons, &propvars, &propboundtypes,
         &propbounds, &npropbounds) );
   for( i = 0; i < npropbounds; i++ )
   {

      assert(GCGvarIsOriginal(propvars[i]));
      assert(GCGvarGetBlock(propvars[i]) < 0); /** @todo this might lead to an error with linking variables*/
      assert(GCGoriginalVarGetNMastervars(propvars[i]) >= 1);
      mastervar = GCGoriginalVarGetMastervars(propvars[i])[0];

      if( propboundtypes[i] == SCIP_BOUNDTYPE_LOWER )
      {
         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(mastervar), propbounds[i]) )
         {
            SCIP_CALL( SCIPchgVarLb(scip, mastervar, propbounds[i]) );
            propcount++;
            SCIPdebugMessage("changed lb of var %s locally to %g\n", SCIPvarGetName(propvars[i]), propbounds[i]);
         }
      }
      else
      {
         if( !SCIPisEQ(scip, SCIPvarGetUbLocal(mastervar), propbounds[i]) )
         {
            SCIP_CALL( SCIPchgVarUb(scip, mastervar, propbounds[i]) );
            propcount++;
            SCIPdebugMessage("changed ub of var %s locally to %g\n", SCIPvarGetName(propvars[i]), propbounds[i]);
         }
      }
   }

   SCIPdebugMessage("Finished propagation of %d stored propagated bounds: %d vars fixed.\n", npropbounds, propcount);

   /* call branching rule specific propagation method */
   if( consdata->branchrule != NULL )
   {
      /** @todo count number of propagations */
      SCIP_CALL( GCGrelaxBranchPropMaster(GCGpricerGetOrigprob(scip), consdata->branchrule, consdata->branchdata, result) );
   }

   if( *result != SCIP_CUTOFF )
      if( propcount > 0 )
         *result = SCIP_REDUCEDDOM;

   consdata->needprop = FALSE;
   consdata->propagatedvars = GCGpricerGetNPricedvars(scip);

   return SCIP_OKAY;
}


/** enforcement method for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpMasterbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** enforcement method for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsMasterbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** check method for solutions */
static
SCIP_DECL_CONSCHECK(consCheckMasterbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** variable lock method */
static
SCIP_DECL_CONSLOCK(consLockMasterbranch)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/* define not used callbacks as NULL */
#define conshdlrCopyMasterbranch NULL
#define consPresolMasterbranch NULL
#define consRespropMasterbranch NULL
#define consExitsolMasterbranch NULL
#define consInitpreMasterbranch NULL
#define consExitpreMasterbranch NULL
#define consTransMasterbranch NULL
#define consInitlpMasterbranch NULL
#define consSepalpMasterbranch NULL
#define consSepasolMasterbranch NULL
#define consEnableMasterbranch NULL
#define consDisableMasterbranch NULL
#define consPrintMasterbranch NULL
#define consDelvarsMasterbranch NULL
#define consCopyMasterbranch NULL
#define consParseMasterbranch NULL
#define consGetVarsMasterbranch NULL
#define consGetNVarsMasterbranch NULL


/*
 * Callback methods of event handler
 */

#define eventCopyOrigvarbound NULL
#define eventFreeOrigvarbound NULL
#define eventExitOrigvarbound NULL

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitOrigvarbound)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}



/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolOrigvarbound)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int i;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, vars[i], SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_BOUNDCHANGED,
            eventhdlr, NULL, NULL) );
   }

   return SCIP_OKAY;
}


#define eventExitsolOrigvarbound NULL
#define eventDeleteOrigvarbound NULL

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecOrigvarbound)
{  /*lint --e{715}*/
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR* var;
   SCIP_Real oldbound;
   SCIP_Real newbound;
   int i;
   int blocknr;
#ifdef SCIP_DEBUG
   SCIP_Bool handled = FALSE;
#endif
   SCIP_VAR** mastervars;
#ifndef NDEBUG
   int nmastervars;
   SCIP_Real* mastervals;
#endif

   eventtype = SCIPeventGetType(event);
   var = SCIPeventGetVar(event);
   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);

   SCIPdebugMessage("eventexec: eventtype = 0x%x, var = %s, oldbound = %f, newbound = %f\n", eventtype, SCIPvarGetName(var), oldbound, newbound);

   if( !GCGrelaxIsInitialized(scip) )
   {
      assert(SCIPvarGetData(var) == NULL);
      SCIPdebugMessage("Ignoring since in presolving / propagating.\n");
      return SCIP_OKAY;
   }

   assert(GCGvarIsOriginal(var));
   blocknr = GCGvarGetBlock(var);

   mastervars = GCGoriginalVarGetMastervars(var);
#ifndef NDEBUG
   nmastervars = GCGoriginalVarGetNMastervars(var);
   mastervals = GCGoriginalVarGetMastervals(var);
#endif

   /* deal with variables present in the pricing */
   if( blocknr >= 0 && GCGrelaxIsPricingprobRelevant(scip, blocknr) )
   {
      SCIPdebugMessage("Pricing var!\n");
      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip),
               GCGoriginalVarGetPricingVar(var), SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
      }
      if( (eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip),
               GCGoriginalVarGetPricingVar(var), SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
      }
   }
   /* deal with variables appearing in the master only */
   if( blocknr == -1 && SCIPgetStage(GCGrelaxGetMasterprob(scip)) >= SCIP_STAGE_SOLVING )
   {
      assert(nmastervars == 1);
      assert(mastervals[0] == 1);
      assert(mastervars[0] != NULL);

      SCIPdebugMessage("Master var!\n");
      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip),
               mastervars[0], SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
      }
      if( (eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip),
               mastervars[0], SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
      }
      if( (eventtype & SCIP_EVENTTYPE_LBTIGHTENED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, GCGconsOrigbranchGetActiveCons(scip), var,
               SCIP_BOUNDTYPE_LOWER, newbound) );
      }
      if( (eventtype & SCIP_EVENTTYPE_UBTIGHTENED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, GCGconsOrigbranchGetActiveCons(scip), var,
               SCIP_BOUNDTYPE_UPPER, newbound) );

         /** @todo do we also have to iterate over the pricing problems or is this handled elsewhere? */
      }
   }
   /* deal with linking variables */
   if( blocknr == -2 )
   {

      int npricingprobs;
      SCIP_VAR** pricingvars;

      SCIPdebugMessage("Linking var!\n");
      pricingvars = GCGlinkingVarGetPricingVars(var);
      npricingprobs = GCGrelaxGetNPricingprobs(scip);

      assert(nmastervars >= 1);
      assert(mastervals[0] == 1);
      assert(mastervars[0] != NULL);

      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
      {
         if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) >= SCIP_STAGE_SOLVING )
         {
#ifdef SCIP_DEBUG
            handled = TRUE;
#endif
            /* add the bound change in the master */
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(mastervars[0]), oldbound));
            SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip),
                  mastervars[0], SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
         }

         /* add the bound change to the pricing problems */
         for( i = 0; i < npricingprobs; ++i )
         {
            if( pricingvars[i] == NULL )
               continue;
#ifdef SCIP_DEBUG
            handled = TRUE;
#endif
            SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip),
                  pricingvars[i], SCIP_BOUNDTYPE_LOWER, oldbound, newbound) );
         }
      }
      if( (eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0 )
      {
         if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) >= SCIP_STAGE_SOLVING )
         {
#ifdef SCIP_DEBUG
            handled = TRUE;
#endif
            /* add the bound change in the master */
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(mastervars[0]), oldbound));
            SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip),
                  mastervars[0], SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
         }

         /* add the bound change to the pricing problems */
         for( i = 0; i < npricingprobs; ++i )
         {
            if( pricingvars[i] == NULL )
               continue;
#ifdef SCIP_DEBUG
            handled = TRUE;
#endif
            SCIP_CALL( GCGconsMasterbranchAddPendingBndChg(GCGrelaxGetMasterprob(scip),
                  pricingvars[i], SCIP_BOUNDTYPE_UPPER, oldbound, newbound) );
         }

      }

      /* store tightened bounds as prop bound changes */
      if( (eventtype & SCIP_EVENTTYPE_LBTIGHTENED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, GCGconsOrigbranchGetActiveCons(scip), var,
               SCIP_BOUNDTYPE_LOWER, newbound) );
      }
      if( (eventtype & SCIP_EVENTTYPE_UBTIGHTENED) != 0 )
      {
#ifdef SCIP_DEBUG
         handled = TRUE;
#endif
         SCIP_CALL( GCGconsOrigbranchAddPropBoundChg(scip, GCGconsOrigbranchGetActiveCons(scip), var,
               SCIP_BOUNDTYPE_UPPER, newbound) );
      }
   }
#ifdef SCIP_DEBUG
   if( !handled )
   {
      SCIPdebugMessage("Effectively ignoring this change\n");
   }
#endif
   return SCIP_OKAY;
}


/*
 * interface methods
 */


/** creates the handler for masterbranch constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrMasterbranch(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   SCIPdebugMessage("Including masterbranch constraint handler.\n");

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrData) );
   conshdlrData->stack = NULL;
   conshdlrData->nstack = 0;
   conshdlrData->maxstacksize = 25;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         SCIP_PROPTIMING_ALWAYS,
         conshdlrCopyMasterbranch, consFreeMasterbranch, consInitMasterbranch, consExitMasterbranch,
         consInitpreMasterbranch, consExitpreMasterbranch, consInitsolMasterbranch, consExitsolMasterbranch,
         consDeleteMasterbranch, consTransMasterbranch, consInitlpMasterbranch,
         consSepalpMasterbranch, consSepasolMasterbranch, consEnfolpMasterbranch, consEnfopsMasterbranch, consCheckMasterbranch,
         consPropMasterbranch, consPresolMasterbranch, consRespropMasterbranch, consLockMasterbranch,
         consActiveMasterbranch, consDeactiveMasterbranch,
         consEnableMasterbranch, consDisableMasterbranch,
         consDelvarsMasterbranch, consPrintMasterbranch, consCopyMasterbranch, consParseMasterbranch,
         consGetVarsMasterbranch, consGetNVarsMasterbranch,
         conshdlrData) );

   /* create event handler data */
   eventhdlrdata = NULL;

   /* include event handler into original SCIP */
   SCIP_CALL( SCIPincludeEventhdlr(GCGpricerGetOrigprob(scip), EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventCopyOrigvarbound, eventFreeOrigvarbound, eventInitOrigvarbound, eventExitOrigvarbound,
         eventInitsolOrigvarbound, eventExitsolOrigvarbound, eventDeleteOrigvarbound, eventExecOrigvarbound,
         eventhdlrdata) );

   SCIP_CALL( SCIPaddBoolParam(GCGpricerGetOrigprob(scip), "relaxing/gcg/enforceproper",
         "should propagated bound changes in the original be enforced in the master (only proper vars)?",
         &conshdlrData->enforceproper, FALSE, TRUE, NULL, NULL) );



   return SCIP_OKAY;
}


/** creates and captures a masterbranch constraint */
SCIP_RETCODE GCGcreateConsMasterbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   SCIP_NODE*            node,               /**< node at which the constraint should be created */
   SCIP_CONS*            parentcons          /**< parent constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(node != NULL || parentcons == NULL);
   if( node != NULL )
      assert((parentcons == NULL) == (SCIPnodeGetDepth(node) == 0));
   else
      assert(parentcons == NULL);

   /* find the masterbranch constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->propagatedvars = 0;
   consdata->needprop = TRUE;

   consdata->node = node;
   consdata->parentcons = parentcons;
   consdata->child1cons = NULL;
   consdata->child2cons = NULL;
   consdata->probingtmpcons = NULL;
   consdata->created = FALSE;
   consdata->origcons = NULL;
   consdata->name = NULL;

   consdata->branchrule = NULL;
   consdata->branchdata = NULL;

   consdata->boundchgvars = NULL;
   consdata->boundtypes = NULL;
   consdata->newbounds = NULL;
   consdata->oldbounds = NULL;
   consdata->nboundchangestreated = NULL;
   consdata->nboundchanges = 0;
   consdata->nactivated = 0;


   SCIPdebugMessage("Creating masterbranch constraint with parent %p.\n", parentcons);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, "masterbranch", conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

   if( parentcons != NULL )
   {
      SCIP_CONSDATA* parentdata;

      parentdata = SCIPconsGetData(parentcons);
      assert(parentdata != NULL);

      if( parentdata->child1cons == NULL )
      {
         parentdata->child1cons = *cons;
      }
      else
      {
         assert(parentdata->child2cons == NULL || SCIPinProbing(scip));

         /* store the second child in case we are in probing and have to overwrite it */
         if( SCIPinProbing(scip) )
         {
            assert(parentdata->probingtmpcons == NULL);
            parentdata->probingtmpcons = parentdata->child2cons;
         }

         parentdata->child2cons = *cons;
      }
   }

   return SCIP_OKAY;
}




/* ----------------------------------- external methods -------------------------- */

/** returns the masterbranch constraint of the current node */
SCIP_CONS* GCGconsMasterbranchGetActiveCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return NULL;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);

   assert(conshdlrData->stack[conshdlrData->nstack-1] != NULL);
   return conshdlrData->stack[conshdlrData->nstack-1];
}


/** returns the stack and the number of elements on it */
void GCGconsMasterbranchGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   *stack = conshdlrData->stack;
   *nstackelements = conshdlrData->nstack;
}

/** returns the number of elements on the stack */
int GCGconsMasterbranchGetNStackelements(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return -1;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   return conshdlrData->nstack;
}

/** returns the branching data for a given masterbranch constraint */
GCG_BRANCHDATA* GCGconsMasterbranchGetBranchdata(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchdata;
}

/** returns the node in the B&B tree at which the given masterbranch constraint is sticking */
SCIP_NODE* GCGconsMasterbranchGetNode(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->node;
}

/** returns the masterbranch constraint of the B&B father of the node at which the
    given masterbranch constraint is sticking */
SCIP_CONS* GCGconsMasterbranchGetParentcons(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->parentcons;
}

/** returns the masterbranch constraint of the first child of the node at which the
    given masterbranch constraint is sticking */
SCIP_CONS* GCGconsMasterbranchGetChild1cons(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->child1cons;
}

/** returns the masterbranch constraint of the second child of the node at which the
    given masterbranch constraint is sticking */
SCIP_CONS* GCGconsMasterbranchGetChild2cons(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->child2cons;
}

/** returns the origbranch constraint of the node in the original program corresponding to the node
    at which the given masterbranch constraint is sticking */
SCIP_CONS* GCGconsMasterbranchGetOrigcons(
   SCIP_CONS*            cons                /**< constraint pointer */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->origcons;
}

/** sets the origbranch constraint of the node in the master program corresponding to the node
    at which the given masterbranchbranch constraint is sticking */
void GCGconsMasterbranchSetOrigcons(
   SCIP_CONS*            cons,               /**< constraint pointer */
   SCIP_CONS*            origcons            /**< original origbranch constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->origcons == NULL || origcons == NULL);

   consdata->origcons = origcons;
}


/** checks the consistency of the masterbranch constraints in the problem */
void GCGconsMasterbranchCheckConsistency(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   int nconss;
#ifndef NDEBUG
   int i;
   SCIP_CONS** conss;
   SCIP_CONSDATA* consdata;
#endif

   if( scip == NULL )
      return;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      SCIPABORT();
#ifdef NDEBUG
      return;
#endif
   }
   nconss = SCIPconshdlrGetNConss(conshdlr);
#ifndef NDEBUG
   conss = SCIPconshdlrGetConss(conshdlr);

   for( i = 0; i < nconss; i++ )
   {
#ifdef SCIP_DEBUG
      SCIP_CONS* parent_origcons;
      SCIP_CONS* origcons_parent;
#endif

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      assert((consdata->parentcons == NULL) == ((consdata->node == NULL) || (SCIPnodeGetDepth(consdata->node) == 0)));

      assert(consdata->origcons == NULL || consdata->created);
      assert(consdata->parentcons == NULL || SCIPconsGetData(consdata->parentcons)->child1cons == conss[i]
         || SCIPconsGetData(consdata->parentcons)->child2cons == conss[i]
         || ( SCIPinProbing(scip) && SCIPconsGetData(consdata->parentcons)->probingtmpcons == conss[i]));
      assert(consdata->child1cons == NULL || SCIPconsGetData(consdata->child1cons)->parentcons == conss[i]);
      assert(consdata->child2cons == NULL || SCIPconsGetData(consdata->child2cons)->parentcons == conss[i]);
      assert(consdata->probingtmpcons == NULL || SCIPinProbing(scip));
      assert(consdata->probingtmpcons == NULL || SCIPconsGetData(consdata->probingtmpcons)->parentcons == conss[i]);
      assert(consdata->origcons == NULL ||
         GCGconsOrigbranchGetMastercons(consdata->origcons) == conss[i]);

#ifdef SCIP_DEBUG
      if( consdata->parentcons != NULL )
         parent_origcons = SCIPconsGetData(consdata->parentcons)->origcons;
      else
         parent_origcons = NULL;

      if( consdata->origcons != NULL )
         origcons_parent = GCGconsOrigbranchGetParentcons(consdata->origcons);
      else
         origcons_parent = NULL;

      SCIPdebugMessage("cons: %s (node %p), origcons: %s, parent %s: %s => %s\n",
         SCIPconsGetName(conss[i]),
         consdata->node,
         consdata->origcons == NULL? "NULL" : SCIPconsGetName(consdata->origcons),
         consdata->parentcons == NULL? "NULL" : SCIPconsGetName(consdata->parentcons),
         parent_origcons == NULL? "NULL" :  SCIPconsGetName(parent_origcons),
         origcons_parent == NULL? "NULL" : SCIPconsGetName(origcons_parent) );
#endif

   }
#endif

   SCIPdebugMessage("checked consistency of %d masterbranch constraints, all ok!\n", nconss);
}

/** adds initial constraint to root node */
SCIP_RETCODE SCIPconsMasterbranchAddRootCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS* cons;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("masterbranch constraint handler not found\n");
      return SCIP_ERROR;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   cons = conshdlrData->stack[0];
   conshdlrData->stack[0] = NULL;
   assert(conshdlrData->nstack == 1);
   conshdlrData->nstack = 0;

   SCIP_CALL( SCIPaddConsNode(scip, SCIPgetRootNode(scip), cons, SCIPgetRootNode(scip)) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}
