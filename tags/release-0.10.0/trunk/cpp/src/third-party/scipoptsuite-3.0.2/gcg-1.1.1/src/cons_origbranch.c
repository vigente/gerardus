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

/**@file   cons_origbranch.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 */

/* #define CHECKCONSISTENCY */

#include <assert.h>
#include <string.h>
#include "cons_origbranch.h"
#include "scip/cons_linear.h"
#include "relax_gcg.h"
#include "cons_masterbranch.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "origbranch"
#define CONSHDLR_DESC          "store branching decision at nodes of the tree constraint handler"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              * propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/** constraint data for branch orig constraints */
struct SCIP_ConsData
{
   SCIP_NODE*            node;               /**< the node at which the cons is sticking */
   SCIP_CONS*            parentcons;         /**< the origbranch constraint of the parent node */
   SCIP_CONS*            child1cons;         /**< the origbranch constraint of the first child node */
   SCIP_CONS*            child2cons;         /**< the origbranch constraint of the second child node */
   SCIP_CONS*            probingtmpcons;     /**< pointer to save the second child if the child2cons pointer is overwritten in probing mode */
   SCIP_CONS*            mastercons;         /**< the masterbranch constraint of the corresponding node
                                              *   in the master program */
   GCG_BRANCHDATA*       branchdata;         /**< branching data stored by the branching rule containing information
                                              *   about the branching restrictions */
   SCIP_BRANCHRULE*      branchrule;         /**< branching rule that created the corresponding node and imposed
                                              *   branching restrictions */
   SCIP_VAR**            propvars;           /**< original variable for which the propagation found domain reductions */
   SCIP_BOUNDTYPE*       propboundtypes;     /**< type of the new bound found by propagation */
   SCIP_Real*            propbounds;         /**< new lower/upper bound of the propagated original variable */
   int                   npropbounds;        /**< number of propagation bounds stored */
   int                   maxpropbounds;      /**< size of propvars, propboundtypes, and propbounds arrays */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_CONS**           stack;              /**< stack for storing active constraints */
   int                   nstack;             /**< number of elements on the stack */
   int                   maxstacksize;       /**< maximum size of the stack */
   SCIP_CONS*            rootcons;           /**< constraint in the root node */
};

/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   SCIPdebugMessage("freeing branch orig constraint handler\n");

   /* free constraint handler storage */
   assert(conshdlrData->stack == NULL);
   if( conshdlrData->rootcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conshdlrData->rootcons) );
   }

   SCIPfreeMemory(scip, &conshdlrData);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   /* prepare stack */
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrData->stack, conshdlrData->maxstacksize) );
   assert( conshdlrData->nstack >= 0 );

   /* check consistency */
   if( conshdlrData->rootcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conshdlrData->rootcons) );
      conshdlrData->rootcons = NULL;
      --(conshdlrData->nstack);
   }
   GCGconsOrigbranchCheckConsistency(scip);

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->nstack <= 1);
   SCIPdebugMessage("exiting solution process branch orig constraint handler\n");

   /* free stack */
   SCIPfreeMemoryArray(scip, &conshdlrData->stack);
   conshdlrData->stack = NULL;

   return SCIP_OKAY;
}

/** exit method of constraint handler (called before problem is free transformed) */
static
SCIP_DECL_CONSEXIT(consExitOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(conshdlr != NULL);
   assert(scip != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   SCIPdebugMessage("exiting transformed branch orig constraint handler\n");

   if( conshdlrdata->rootcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conshdlrdata->rootcons) );
      conshdlrdata->rootcons = NULL;
   }
   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSDATA* parentdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   SCIPdebugMessage("Deleting branch orig constraint: <%s>.\n", SCIPconsGetName(cons));

   /* set the origcons pointer of the corresponding mastercons to NULL */
   if( (*consdata)->mastercons != NULL )
      GCGconsMasterbranchSetOrigcons((*consdata)->mastercons, NULL);

   /* set the pointer in the parent constraint to NULL */
   if( (*consdata)->parentcons != NULL )
   {
      parentdata = SCIPconsGetData((*consdata)->parentcons);
      if( parentdata->child1cons == cons )
         parentdata->child1cons = NULL;
      else if( parentdata->probingtmpcons == cons )
      {
         assert(SCIPinProbing(scip));
         parentdata->probingtmpcons = NULL;
      }
      else
      {
         assert(parentdata->child2cons == cons);
         parentdata->child2cons = NULL;
         if( SCIPinProbing(scip) )
         {
            parentdata->child2cons = parentdata->probingtmpcons;
            parentdata->probingtmpcons = NULL;
         }
      }
   }
   /* no child nodes may exist */
   assert((*consdata)->child1cons == NULL);
   assert((*consdata)->child2cons == NULL);

   /* delete branchdata, if no mastercons is linked, which would still need the branchdata
    * otherwise, the mastercons deletes the branchdata when it is deleted itself */
   if( (*consdata)->mastercons == NULL && (*consdata)->branchdata != NULL )
   {
      SCIP_CALL( GCGrelaxBranchDataDelete(scip, (*consdata)->branchrule, &(*consdata)->branchdata) );
   }

   /* free propagation domain changes arrays */
   if( (*consdata)->maxpropbounds > 0 )
   {
      SCIPfreeMemoryArray(scip, &((*consdata)->propvars));
      SCIPfreeMemoryArray(scip, &((*consdata)->propboundtypes));
      SCIPfreeMemoryArray(scip, &((*consdata)->propbounds));
   }

   /* free constraint data */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   assert(SCIPconsGetData(cons) != NULL);

   if( SCIPconsGetData(cons)->node == NULL )
      SCIPconsGetData(cons)->node = SCIPgetRootNode(scip);

   SCIPdebugMessage("Activating branch orig constraint: <%s>[stack size: %d].\n", SCIPconsGetName(cons),
      conshdlrData->nstack+1);

   /* put constraint on the stack */
   if( conshdlrData->nstack >= conshdlrData->maxstacksize )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrData->stack), 2*(conshdlrData->maxstacksize)) );
      conshdlrData->maxstacksize = 2*(conshdlrData->maxstacksize);
      SCIPdebugMessage("reallocating Memory for stack! %d --> %d\n", conshdlrData->maxstacksize/2, conshdlrData->maxstacksize);
   }

   /* put constraint on the stack */
   assert(conshdlrData->stack != NULL);
   conshdlrData->stack[conshdlrData->nstack] = cons;
   ++(conshdlrData->nstack);

   return SCIP_OKAY;
}



/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveOrigbranch)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL || conshdlrData->nstack <= 1);
   assert(conshdlrData->nstack <= 1 || cons == conshdlrData->stack[conshdlrData->nstack-1]);

   assert(SCIPconsGetData(cons) != NULL);

   SCIPdebugMessage("Deactivating branch orig constraint: <%s> [stack size: %d].\n",
      SCIPconsGetName(cons), conshdlrData->nstack-1);

   /* remove constraint from the stack */
   if( conshdlrData->nstack > 0 )
      --(conshdlrData->nstack);

   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropOrigbranch)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** lp solution enforcement method */
static
SCIP_DECL_CONSENFOLP(consEnfolpOrigbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** pseudo solution enforcement method */
static
SCIP_DECL_CONSENFOPS(consEnfopsOrigbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** solution check method */
static
SCIP_DECL_CONSCHECK(consCheckOrigbranch)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** variable lock method */
static
SCIP_DECL_CONSLOCK(consLockOrigbranch)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/* define not used callbacks as NULL */
#define consPresolOrigbranch NULL
#define consRespropOrigbranch NULL
#define consInitOrigbranch NULL
#define consInitpreOrigbranch NULL
#define consExitpreOrigbranch NULL
#define consTransOrigbranch NULL
#define consInitlpOrigbranch NULL
#define consSepalpOrigbranch NULL
#define consSepasolOrigbranch NULL
#define consEnableOrigbranch NULL
#define consDisableOrigbranch NULL
#define consDelvarsOrigbranch NULL
#define consPrintOrigbranch NULL
#define consCopyOrigbranch NULL
#define consParseOrigbranch NULL
#define consGetVarsOrigbranch NULL
#define consGetNVarsOrigbranch NULL
/*
 * interface methods
 */


/** creates the handler for origbranch constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrOrigbranch(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;

   SCIPdebugMessage("Including branch orig constraint handler.\n");

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrData) );
   conshdlrData->stack = NULL;
   conshdlrData->nstack = 0;
   conshdlrData->maxstacksize = 25;
   conshdlrData->rootcons = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         SCIP_PROPTIMING_ALWAYS,
         consCopyOrigbranch, consFreeOrigbranch, consInitOrigbranch, consExitOrigbranch,
         consInitpreOrigbranch, consExitpreOrigbranch, consInitsolOrigbranch, consExitsolOrigbranch,
         consDeleteOrigbranch, consTransOrigbranch, consInitlpOrigbranch,
         consSepalpOrigbranch, consSepasolOrigbranch, consEnfolpOrigbranch, consEnfopsOrigbranch, consCheckOrigbranch,
         consPropOrigbranch, consPresolOrigbranch, consRespropOrigbranch, consLockOrigbranch,
         consActiveOrigbranch, consDeactiveOrigbranch,
         consEnableOrigbranch, consDisableOrigbranch,
         consDelvarsOrigbranch, consPrintOrigbranch, consCopyOrigbranch, consParseOrigbranch,
         consGetVarsOrigbranch, consGetNVarsOrigbranch,
         conshdlrData) );

   return SCIP_OKAY;
}


/** creates and captures a origbranch constraint */
SCIP_RETCODE GCGcreateConsOrigbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_NODE*            node,               /**< the node to which this origbranch constraint belongs */
   SCIP_CONS*            parentcons,         /**< origbranch constraint associated with the father node */
   SCIP_BRANCHRULE*      branchrule,         /**< the branching rule that created the b&b node the constraint belongs to */
   GCG_BRANCHDATA*       branchdata          /**< branching data storing information about the branching restrictions at the
                                              *   corresponding node */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert((parentcons == NULL) == (node == NULL));

   /* find the origbranch constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   /* initialize the fields in the constraint data */
   consdata->parentcons = parentcons;
   consdata->node = node;
   consdata->child1cons = NULL;
   consdata->child2cons = NULL;
   consdata->probingtmpcons = NULL;
   consdata->mastercons = NULL;
   consdata->branchrule = branchrule;
   consdata->branchdata = branchdata;
   consdata->npropbounds = 0;
   consdata->maxpropbounds = 0;
   consdata->propvars = NULL;
   consdata->propboundtypes = NULL;
   consdata->propbounds = NULL;

   SCIPdebugMessage("Creating branch orig constraint: <%s>.\n", name);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

   /* store the pointer to the new constraint in the origbranch constraint of the parent node */
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

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
SCIP_CONS* GCGconsOrigbranchGetActiveCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return NULL;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);

   return conshdlrData->stack[conshdlrData->nstack-1];
}


/** returns the stack and the number of elements on it */
void GCGconsOrigbranchGetStack(
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
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   *stack = conshdlrData->stack;
   *nstackelements = conshdlrData->nstack;

}

/** returns the branching data for a given origbranch constraint */
GCG_BRANCHDATA* GCGconsOrigbranchGetBranchdata(
   SCIP_CONS*            cons                /**< origbranch constraint for which the branching data is requested */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchdata;
}

/** returns the branchrule for a given origbranch constraint */
SCIP_BRANCHRULE* GCGconsOrigbranchGetBranchrule(
   SCIP_CONS*            cons                /**< origbranch constraint for which the branchrule is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->branchrule;
}

/** returns the node in the B&B tree at which the given origbranch constraint is sticking */
SCIP_NODE* GCGconsOrigbranchGetNode(
   SCIP_CONS*            cons                /**< origbranch constraint for which the corresponding node is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->node;
}

/** returns the origbranch constraint of the B&B father of the node at which the
    given origbranch constraint is sticking */
SCIP_CONS* GCGconsOrigbranchGetParentcons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the origbranch constraint of
                                              *   the father node is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->parentcons;
}

/** returns the origbranch constraint of the first child of the node at which the
    given origbranch constraint is sticking */
SCIP_CONS* GCGconsOrigbranchGetChild1cons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the origbranch constraint of
                                              *   the first child node is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->child1cons;
}

/** returns the origbranch constraint of the second child of the node at which the
    given origbranch constraint is sticking */
SCIP_CONS* GCGconsOrigbranchGetChild2cons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the origbranch constraint of
                                              *   the second child node is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->child2cons;
}

/** sets the masterbranch constraint of the node in the master program corresponding to the node
    at which the given origbranchbranch constraint is sticking */
void GCGconsOrigbranchSetMastercons(
   SCIP_CONS*            cons,               /**< origbranch constraint for which the masterbranch constraint should be set */
   SCIP_CONS*            mastercons          /**< masterbranch constraint corresponding to the given origbranch constraint */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   //   if( consdata->mastercons != mastercons )
   //   assert(consdata->mastercons == NULL || mastercons == NULL);

   consdata->mastercons = mastercons;
}

/** returns the masterbranch constraint of the node in the master program corresponding to the node
    at which the given origbranchbranch constraint is sticking */
SCIP_CONS* GCGconsOrigbranchGetMastercons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the corresponding masterbranch
                                              *   constraint is requested */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->mastercons != NULL);

   return consdata->mastercons;
}


/** checks the consistency of the origbranch constraints in the problem */
void GCGconsOrigbranchCheckConsistency(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
#ifdef CHECKCONSISTENCY

   SCIP_CONSHDLR*     conshdlr;

#ifndef NDEBUG
   SCIP_CONS** conss;
   int nconss;
   int i;
   SCIP_CONSDATA* consdata;
#endif

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return;
   }
#ifndef NDEBUG
   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   for( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      assert(consdata->node != NULL);
      assert((consdata->parentcons == NULL) == (SCIPnodeGetDepth(consdata->node) == 0));
      assert(consdata->parentcons == NULL || SCIPconsGetData(consdata->parentcons)->child1cons == conss[i]
         || SCIPconsGetData(consdata->parentcons)->child2cons == conss[i]
         || ( SCIPinProbing(scip) && SCIPconsGetData(consdata->parentcons)->probingtmpcons == conss[i]));
      assert(consdata->child1cons == NULL || SCIPconsGetData(consdata->child1cons)->parentcons == conss[i]);
      assert(consdata->child2cons == NULL || SCIPconsGetData(consdata->child2cons)->parentcons == conss[i]);
      assert(consdata->probingtmpcons == NULL || SCIPinProbing(scip));
      assert(consdata->probingtmpcons == NULL || SCIPconsGetData(consdata->probingtmpcons)->parentcons == conss[i]);
      assert(consdata->mastercons == NULL ||
         GCGconsMasterbranchGetOrigcons(consdata->mastercons) == conss[i]);
   }
#endif
#endif
}

/** adds a bound change on an original variable found by propagation in the original problem
 *  to the given origbranch constraint so that it will be transferred to the master problem */
SCIP_RETCODE GCGconsOrigbranchAddPropBoundChg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< origbranch constraint to which the bound change is added */
   SCIP_VAR*             var,                /**< variable on which the bound change was performed */
   SCIP_BOUNDTYPE        boundtype,          /**< bound type of the bound change */
   SCIP_Real             newbound            /**< new bound of the variable after the bound change */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* realloc the arrays, if needed */
   if( consdata->npropbounds >= consdata->maxpropbounds )
   {
      consdata->maxpropbounds = consdata->npropbounds+5;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(consdata->propvars), consdata->maxpropbounds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(consdata->propboundtypes), consdata->maxpropbounds) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(consdata->propbounds), consdata->maxpropbounds) );
   }

   SCIPdebugMessage("Bound change stored at branch orig constraint: <%s>.\n", SCIPconsGetName(cons));

   /* store the new bound change */
   consdata->propvars[consdata->npropbounds] = var;
   consdata->propboundtypes[consdata->npropbounds] = boundtype;
   consdata->propbounds[consdata->npropbounds] = newbound;
   consdata->npropbounds++;

   /* mark the corresponding master node to be repropagated */
   if( consdata->mastercons != NULL )
   {
      SCIP_CALL( SCIPrepropagateNode(GCGrelaxGetMasterprob(scip), GCGconsMasterbranchGetNode(consdata->mastercons)) );
   }

   return SCIP_OKAY;
}

/** returns the array of bound changes on original variables found by propagation in the original problem
 *  at the node corresponding to the given origbranch constraint */
SCIP_RETCODE GCGconsOrigbranchGetPropBoundChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< origbranch constraint for which the bound changes are requested */
   SCIP_VAR***           vars,               /**< pointer to store array of variables corresponding to the bound changes */
   SCIP_BOUNDTYPE**      boundtypes,         /**< pointer to store array of the types of the bound changes */
   SCIP_Real**           newbounds,          /**< pointer to store array of the new bounds */
   int*                  npropbounds         /**< pointer to store the number of bound changes stored at the constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *vars = consdata->propvars;
   *boundtypes = consdata->propboundtypes;
   *newbounds = consdata->propbounds;
   *npropbounds = consdata->npropbounds;

   consdata->npropbounds = 0;

   return SCIP_OKAY;
}

/** returns the number of bound changes on original variables found by propagation in the original problem
 *  at the node corresponding to the given origbranch constraint */
int GCGconsOrigbranchGetNPropBoundChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< origbranch constraint for which the bound changes are requested */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->npropbounds;
}

/** adds initial constraint to root node */
SCIP_RETCODE SCIPconsOrigbranchAddRootCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("origbranch constraint handler not found\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   if( conshdlrdata->rootcons == NULL )
   {
      SCIP_CALL( GCGcreateConsOrigbranch(scip, &cons, "root-origbranch", NULL, NULL, NULL, NULL) );

      SCIP_CALL( SCIPaddConsNode(scip, SCIPgetRootNode(scip), cons, SCIPgetRootNode(scip)) );
      conshdlrdata->rootcons = cons;
   }

   /* check consistency */
   GCGconsOrigbranchCheckConsistency(scip);

   return SCIP_OKAY;
}
