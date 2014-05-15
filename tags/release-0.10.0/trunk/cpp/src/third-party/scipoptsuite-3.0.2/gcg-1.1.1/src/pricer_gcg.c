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

/**@file   pricer_gcg.c
 * @brief  pricer for generic column generation
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @author Alexander Gross
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"

#include "pricer_gcg.h"
#include "sepa_master.h"

#include "relax_gcg.h"
#include "struct_solver.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"
#include "cons_masterbranch.h"

#define PRICER_NAME            "gcg"
#define PRICER_DESC            "pricer for gcg"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

#define DEFAULT_MAXVARSROUNDFARKAS       10         /**< maximal number of variables per farkas pricing round */
#define DEFAULT_MAXVARSROUNDREDCOSTROOT  100        /**< maximal number of variables per reduced cost pricing round at root node */
#define DEFAULT_MAXVARSROUNDREDCOST      100        /**< maximal number of variables per reduced cost pricing round */
#define DEFAULT_MAXSUCCESSFULMIPSREDCOST INT_MAX    /**< maximal number of successful MIP solves */
#define DEFAULT_MAXROUNDSREDCOST         INT_MAX    /**< maximal number of reduced cost pricing rounds */
#define DEFAULT_MAXSOLSPROB              INT_MAX    /**< maximal number of solution per pricing problem*/
#define DEFAULT_USEHEURPRICING           FALSE      /**< should heuristic pricing be used */
#define DEFAULT_ABORTPRICINGINT          TRUE       /**< should the pricing be aborted when integral */
#define DEFAULT_ABORTPRICINGGAP          0.00       /**< gap at which the pricing is aborted */
#define DEFAULT_SUCCESSFULMIPSREL        1.0        /**< factor of successful mips to be solved */
#define DEFAULT_MIPSRELREDCOSTROOT       1.0        /**< factor of reduced cost pricing MIPs to be solved at root node */
#define DEFAULT_MIPSRELREDCOST           1.0        /**< factor of reduced cost pricing MIPs to be solver*/
#define DEFAULT_MIPSRELFARKAS            1.0        /**< factor of farkas pricing MIPs to be solved */
#define DEFAULT_DISPINFOS                FALSE      /**< should additional information be displayed */
#define DEFAULT_SORTING                  2          /**< default sorting method for pricing mips
                                                     *    0 :   order of pricing problems
                                                     *    1 :   according to dual solution of convexity constraint
                                                     *    2 :   according to reliability from previous round)
                                                     */

#define EVENTHDLR_NAME         "probdatavardeleted"
#define EVENTHDLR_DESC         "event handler for variable deleted event"

/** small macro to simplify printing pricer information */
#define GCGpricerPrintInfo(scip,pricerdata, ...) do { \
   if( pricerdata->dispinfos ) { \
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,__VA_ARGS__);\
   } else {\
      SCIPdebugMessage(__VA_ARGS__); \
   }\
   }while( FALSE )

#define PRICER_STAT_ARRAYLEN_TIME 1024                /**< length of the array for Time histogram representation */
#define PRICER_STAT_BUCKETSIZE_TIME 10                /**< size of the buckets for Time histogram representation */
#define PRICER_STAT_ARRAYLEN_VARS 1024                /**< length of the array for foundVars histogram representation */
#define PRICER_STAT_BUCKETSIZE_VARS 1                 /**< size of the buckets for foundVars histogram representation */

/*
 * Data structures
 */


/** variable pricer data */
struct SCIP_PricerData
{
   int                   npricingprobs;      /**< number of pricing problems */
   SCIP**                pricingprobs;       /**< pointers to the pricing problems */
   SCIP_Real*            dualsolconv;        /**< array of dual solutions for the convexity constraints */
   SCIP*                 origprob;           /**< the original program */
   SCIP_Real*            solvals;            /**< solution values of variables in the pricing problems */
   int*                  npointsprob;        /**< number of variables representing points created by the pricing probs */
   int*                  nraysprob;          /**< number of variables representing rays created by the pricing probs */
   SCIP_Longint          currnodenr;         /**< current node number in the masterproblem*/
   SCIP_HASHMAP*         mapcons2idx;        /**< hashmap mapping constraints to their index in the conss array */
   SCIP_Real*            score;              /**< score of the pricing problem problems */
   int*                  permu;              /**< current permutation of the pricing problems */
   int                   npricingprobsnotnull; /**< number of non-Null pricing problems*/

   SCIP_VAR**            pricedvars;         /**< array of all priced variables */
   int                   npricedvars;        /**< number of priced variables */
   int                   maxpricedvars;      /**< maximal number of priced variables */

   /** variables used for statistics */
   SCIP_CLOCK*           redcostclock;       /**< time for reduced cost pricing */
   SCIP_CLOCK*           farkasclock;        /**< time for farkas pricing */
   SCIP_CLOCK*           freeclock;          /**< time for freeing pricing problems */
   SCIP_CLOCK*           transformclock;     /**< time for transforming pricing problems */
   int                   solvedsubmipsoptimal; /**< number of optimal pricing runs */
   int                   solvedsubmipsheur;  /**< number of heuristical pricing runs*/
   int                   calls;              /**< number of total pricing calls */
   int                   farkascalls;        /**< number of farkas pricing calls */
   int                   redcostcalls;       /**< number of reduced cost pricing calls */

   /* solver data */
   GCG_SOLVER**          solvers;            /**< pricing solvers array */
   int                   nsolvers;           /**< number of pricing solvers */

   /* event handler */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler */

   /** parameter values */
   SCIP_VARTYPE          vartype;            /**< vartype of created master variables */
   int                   maxvarsroundfarkas; /**< maximal number of variables per farkas round */
   int                   maxvarsroundredcost; /**< maximal number of variables per reduced cost round */
   int                   maxvarsroundredcostroot;  /**< maximal number of variables per round in the root node */
   int                   maxsuccessfulmipsredcost; /**< maximal number of sucessful pricing runs for redcost pricing */
   int                   maxroundsredcost;   /**< maximal number of reduced cost rounds */
   int                   maxsolsprob;        /**< maximal number of solutions per pricing problem */
   int                   nroundsredcost;     /**< number of reduced cost rounds */
   int                   sorting;            /**< how should pricing problems be sorted */
   SCIP_Bool             useheurpricing;     /**< should heuristic pricing be used */
   SCIP_Bool             abortpricingint;    /**< should the pricing be aborted on integral solutions */
   SCIP_Bool             dispinfos;          /**< should pricing information be displayed*/
   SCIP_Real             successfulmipsrel;  /**< Factor of successful MIPs solved until pricing be aborted */
   SCIP_Real             mipsrelredcost;     /**< Factor of successful reduced cost MIPs solved until pricing aborted */
   SCIP_Real             mipsrelredcostroot; /**< Factor of successful reduced cost MIPs solved until pricing aborted at root node */
   SCIP_Real             mipsrelfarkas;      /**< Factor of successful farkas MIPs solved until pricing aborted */
   SCIP_Real             abortpricinggap;    /**< Gap at which pricing should be aborted */


   /** statistics */
   int                   oldvars;            /**< Vars of last pricing iteration */
   int*                  farkascallsdist;    /**< Calls of each farkas pricing problem */
   int*                  farkasfoundvars;    /**< Found vars of each farkas pricing problem */
   double*               farkasnodetimedist; /**< Time spend in each farkas pricing problem */

   int*                  redcostcallsdist;   /**< Calls of each redcost pricing problem */
   int*                  redcostfoundvars;   /**< Found vars of each redcost pricing problem */
   double*               redcostnodetimedist; /**< Time spend in each redcost pricing problem */

   int*                  nodetimehist;       /**< Histogram of nodetime distribution */
   int*                  foundvarshist;      /**< Histogram of foundvars distribution */

   double                rootnodedegeneracy; /**< degeneracy of the root node */
   double*               nodedegeneracy;     /**< degeneracy of the remaining nodes */
   double                avgnodedegeneracy;  /**< average degeneray of all nodes */
   int                   nnodes;             /**< number of nodes handled so far */
   int                   maxnnodes;          /**< maximal number of nodes to handle */
   SCIP_NODE*            lastnode;           /**< last handled node */
};


/*
 * Callback methods of event handler
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
#define eventFreeVardeleted NULL

/** initialization method of event handler (called after problem was transformed) */
#define eventInitVardeleted NULL

/** deinitialization method of event handler (called before transformed problem is freed) */
#define eventExitVardeleted NULL

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#define eventInitsolVardeleted NULL

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#define eventExitsolVardeleted NULL

/** frees specific event data */
#define eventDeleteVardeleted NULL

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecVardeleted)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_VAR** origvars;
   int i;

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARDELETED);
   var = SCIPeventGetVar(event);
   assert(var != NULL);

   SCIPdebugMessage("remove master variable %s from pricerdata and corresponding original variables\n", SCIPvarGetName(var));

   assert(GCGvarIsMaster(var));
   origvars = GCGmasterVarGetOrigvars(var);
   assert(origvars != NULL);

   /* remove master variable from corresponding pricing original variables */
   for( i = 0; i < GCGmasterVarGetNOrigvars(var); ++i )
   {
      SCIP_CALL( GCGoriginalVarRemoveMasterVar(scip, origvars[i], var) );
   }

   /* remove variable from array of stored priced variables */
   for( i = 0; i < pricerdata->npricedvars; ++i )
   {
      if( pricerdata->pricedvars[i] == var )
      {
         /* drop vardeleted event on variable */
         SCIP_CALL( SCIPdropVarEvent(scip, pricerdata->pricedvars[i], SCIP_EVENTTYPE_VARDELETED,
               pricerdata->eventhdlr, NULL, -1) );

         SCIP_CALL( SCIPreleaseVar(scip, &(pricerdata->pricedvars[i])) );
         (pricerdata->npricedvars)--;
         pricerdata->pricedvars[i] = pricerdata->pricedvars[pricerdata->npricedvars];

         break;
      }
   }
   assert(i <= pricerdata->npricedvars);
#ifndef NDEBUG
   for( ; i < pricerdata->npricedvars; ++i )
   {
      assert(pricerdata->pricedvars[i] != var);
   }
#endif

   return SCIP_OKAY;
}


/*
 * Local methods
 */

/** returns TRUE or FALSE, depending whether we are in the root node or not */
static SCIP_Bool isRootNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   return (SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip));
}


/** ensures size of pricedvars array */
static
SCIP_RETCODE ensureSizePricedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< Pricerdata data structure */
   int                   size                /**< needed size */
   )
{
   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricerdata->pricedvars != NULL);

   if( pricerdata->maxpricedvars < size )
   {
      int oldsize;

      oldsize = pricerdata->maxpricedvars;
      pricerdata->maxpricedvars = SCIPcalcMemGrowSize(scip, size);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(pricerdata->pricedvars), oldsize, pricerdata->maxpricedvars) );
   }
   assert(pricerdata->maxpricedvars >= size);

   return SCIP_OKAY;
}


/** ensures size of solvers array */
static
SCIP_RETCODE ensureSizeSolvers(
   SCIP*                 scip,               /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata          /**< Pricerdata data structure  */
   )
{
   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));

   if( pricerdata->nsolvers == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->solvers), 1) ); /*lint !e506*/
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(pricerdata->solvers), pricerdata->nsolvers+1) );
   }

   return SCIP_OKAY;
}
#ifdef ENABLESTATISTICS
/** ensures size of nodes array */
static
SCIP_RETCODE ensureSizeAvgnodedegeneracy(
   SCIP*                 scip,               /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata          /**< Pricerdata data structure  */
   )
{
   int memgrowsize;
   assert(scip != NULL);
   assert(pricerdata != NULL);

   if( pricerdata->maxnnodes > pricerdata->nnodes )
      return SCIP_OKAY;

   memgrowsize = SCIPcalcMemGrowSize(scip, pricerdata->nnodes+1);
   SCIP_CALL( SCIPreallocMemoryArray(scip, &(pricerdata->nodedegeneracy), memgrowsize) );

   pricerdata->maxnnodes = memgrowsize;

   assert(pricerdata->maxnnodes > pricerdata->nnodes);
   return SCIP_OKAY;
}
#endif

/** gets the NodeTimeDistribution in the form of a histogram */
static
void GCGpricerGetNodeTimeHistogram(
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   SCIP_Real             time                /**< time the pricingproblem needed */
   )
{
   int i;

   /* 1000* because mapping milliseconds on the index i */
   i = 1000*time/PRICER_STAT_BUCKETSIZE_TIME; /*lint !e524 */

   if( i >= PRICER_STAT_ARRAYLEN_TIME )
   {
      i = PRICER_STAT_ARRAYLEN_TIME-1;
   }
   pricerdata->nodetimehist[i]++;

}


/** gets the FoundVarsDistribution in form of a histogram */
static
void GCGpricerGetFoundVarsHistogram(
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   int                   foundvars           /**< foundVars in pricingproblem */
   )
{
   int i;
   i = foundvars/PRICER_STAT_BUCKETSIZE_VARS;

   if( i >= PRICER_STAT_ARRAYLEN_VARS )
   {
      i = PRICER_STAT_ARRAYLEN_VARS-1;
   }
   pricerdata->foundvarshist[i]++;

}


/** gets the statistics of the pricingprobs like calls, foundvars and time */
static
void GCGpricerCollectStatistic(
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   GCG_PRICETYPE         type,               /**< type of pricing: optimal or heuristic */
   int                   probindex,          /**< index of the pricingproblem */
   SCIP_Real             time                /**< time the pricingproblem needed */
   )
{
   int foundvars;

   foundvars = pricerdata->npricedvars - pricerdata->oldvars;

   if( type == GCG_PRICETYPE_FARKAS )
   {

      pricerdata->farkascallsdist[probindex]++; /*Calls*/
      pricerdata->farkasfoundvars[probindex] += foundvars;
      pricerdata->farkasnodetimedist[probindex] += time;   /*Time*/

   }
   else if( type == GCG_PRICETYPE_REDCOST )
   {

      pricerdata->redcostcallsdist[probindex]++;
      pricerdata->redcostfoundvars[probindex] += foundvars;
      pricerdata->redcostnodetimedist[probindex] += time;

   }
   pricerdata->oldvars = pricerdata->npricedvars;

   GCGpricerGetNodeTimeHistogram(pricerdata, time);
   GCGpricerGetFoundVarsHistogram(pricerdata, foundvars);

}


/** frees all solvers */
static
SCIP_RETCODE solversFree(
   SCIP*                 scip,               /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata          /**< Pricerdata data structure  */
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverfree != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverfree(scip, pricerdata->solvers[i]) );

         BMSfreeMemoryArray(&pricerdata->solvers[i]->name);
         BMSfreeMemoryArray(&pricerdata->solvers[i]->description);

         SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->solvers[i]->optfarkasclock)) );
         SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->solvers[i]->optredcostclock)) );
         SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->solvers[i]->heurfarkasclock)) );
         SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->solvers[i]->heurredcostclock)) );

         SCIPfreeMemory(scip, &(pricerdata->solvers[i]));
      }
   }

   return SCIP_OKAY;
}

/** calls the init method on all solvers */
static
SCIP_RETCODE solversInit(
   SCIP*                 scip,               /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata          /**< Pricerdata data structure  */
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverinit != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverinit(scip, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

/** calls the exit method on all solvers */
static
SCIP_RETCODE solversExit(
   SCIP*                 scip,               /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata          /**< Pricerdata data structure  */
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverexit != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverexit(scip, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

/** calls the initsol method on all solvers */
static
SCIP_RETCODE solversInitsol(
   SCIP*                 scip,               /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata          /**< Pricerdata data structure  */
   )
{
   int i;

   if( pricerdata->npricingprobs == 0)
      return SCIP_OKAY;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverinitsol != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverinitsol(scip, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

/** calls the exitsol method of all solvers */
static
SCIP_RETCODE solversExitsol(
   SCIP*                 scip,               /**< SCIP data structure        */
   SCIP_PRICERDATA*      pricerdata          /**< Pricerdata data structure  */
   )
{
   int i;

   assert(scip != NULL);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      if( pricerdata->solvers[i]->solverexitsol != NULL )
      {
         SCIP_CALL( pricerdata->solvers[i]->solverexitsol(scip, pricerdata->solvers[i]) );
      }
   }

   return SCIP_OKAY;
}

#ifdef ENABLESTATISTICS
/** returns the gegeneracy of the masterproblem */
static
SCIP_RETCODE computeCurrentDegeneracy(
   SCIP*                 scip,               /**< SCIP data structure */
   double*               degeneracy          /**< pointer to store degeneracy */
   )
{
   int ncols;
   int nrows;
   int i;
   int count;
   int countz;
   int colindex;
   double currentVal;
   int* indizes;
   SCIP_COL** cols;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(degeneracy != NULL);

   *degeneracy = 0.0;
   ncols = SCIPgetNLPCols(scip);
   nrows = SCIPgetNLPRows(scip);
   cols = SCIPgetLPCols(scip);

   SCIP_CALL( SCIPallocMemoryArray(scip, &indizes, ncols+nrows) );

   for( i = 0; i < ncols+nrows; i++ )
   {
      indizes[i] = 0;
   }

   /* gives indices of Columns in Basis and indices of vars in Basis */
   SCIP_CALL( SCIPgetLPBasisInd(scip, indizes) );

   countz = 0;
   count = 0;

   for( i = 0; i < nrows; i++ )
   {
      colindex = indizes[i];
      /* is column if >0 it is column in basis, <0 is for row */
      if( colindex > 0 )
      {
         var = SCIPcolGetVar(cols[colindex]);

         currentVal = SCIPgetSolVal(scip, NULL, var);

         if( SCIPisEQ(scip, currentVal, 0.0) )
            countz++;

         count++;
      }
   }

   /* Degeneracy in % */
   if( count > 0 )
      *degeneracy = ((double)countz / count)*100;

   assert(*degeneracy <= 100);

   SCIPfreeMemoryArray(scip, &indizes);

   return SCIP_OKAY;
}
#endif

/** solves a specific pricing problem */
static
SCIP_RETCODE solvePricingProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   int                   prob,               /**< index of pricing problem */
   GCG_PRICETYPE         pricetype,          /**< type of pricing: optimal or heuristic */
   SCIP_Real*            lowerbound,         /**< dual bound returned by pricing problem */
   SCIP_VAR****          solvars,            /**< solution variables for found solutions */
   SCIP_Real***          solvals,            /**< values for solution variables for each solution */
   int**                 nsolvars,           /**< number of non-zero variables for each solution */
   SCIP_Bool**           solisray,           /**< array to indicate whether solution is a ray */
   int*                  nsols,              /**< number of solutions */
   SCIP_STATUS*          status              /**< solution status of the pricing problem */
   )
{
   int i;
   SCIP_Real timelimit;

   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricerdata->pricingprobs[prob] != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   *status = SCIP_STATUS_UNKNOWN;

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      SCIP_CLOCK* clock;
      int* calls;

      if( pricetype == GCG_PRICETYPE_FARKAS )
      {
         clock = pricerdata->solvers[i]->optfarkasclock;
         calls = &(pricerdata->solvers[i]->optfarkascalls);
      }
      else
      {
         clock = pricerdata->solvers[i]->optredcostclock;
         calls = &(pricerdata->solvers[i]->optredcostcalls);
      }

      /* get time limit */
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      if( !SCIPisInfinity(scip, timelimit) && timelimit - SCIPgetSolvingTime(scip) < 0 )
      {
         *nsols = 0;
         *status = SCIP_STATUS_TIMELIMIT;
      }

      if( pricerdata->solvers[i]->solversolve != NULL )
      {
         SCIP_CALL( SCIPstartClock(scip, clock) );

         SCIP_CALL( pricerdata->solvers[i]->solversolve(scip, pricerdata->solvers[i],
               pricerdata->pricingprobs[prob], prob, lowerbound,
               solvars, solvals, nsolvars, solisray, nsols, status) );

         SCIP_CALL( SCIPstopClock(scip, clock) );

         if( pricetype == GCG_PRICETYPE_FARKAS && *status != SCIP_STATUS_UNKNOWN )

         if( *status != SCIP_STATUS_UNKNOWN )
            (*calls)++;

         if( *status == SCIP_STATUS_OPTIMAL || *status == SCIP_STATUS_UNBOUNDED )
         {
            GCGpricerCollectStatistic(pricerdata, pricetype, prob,
                          SCIPgetSolvingTime(pricerdata->pricingprobs[prob]));
            break;
         }

      }
   }

   return SCIP_OKAY;
}

/** solves the specific pricing problem heuristically */
static
SCIP_RETCODE solvePricingProblemHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   int                   prob,               /**< index of pricing problem */
   GCG_PRICETYPE         pricetype,          /**< type of pricing: optimal or heuristic */
   SCIP_Real*            lowerbound,         /**< dual bound returned by pricing problem */
   SCIP_VAR****          solvars,            /**< solution variables for found solutions */
   SCIP_Real***          solvals,            /**< values for solution variables for each solution */
   int**                 nsolvars,           /**< number of non-zero variables for each solution */
   SCIP_Bool**           solisray,           /**< array to indicate whether solution is a ray */
   int*                  nsols,              /**< number of solutions */
   SCIP_STATUS*          status              /**< solution status of the pricing problem */
   )
{
   int i;
   SCIP_Real timelimit;

   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricerdata->pricingprobs[prob] != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);

   *status = SCIP_STATUS_UNKNOWN;

   for( i = 0; i < pricerdata->nsolvers; i++ )
   {
      SCIP_CLOCK* clock;
      int* calls;

      if( pricetype == GCG_PRICETYPE_FARKAS )
      {
         clock = pricerdata->solvers[i]->heurfarkasclock;
         calls = &(pricerdata->solvers[i]->heurfarkascalls);
      }
      else
      {
         clock = pricerdata->solvers[i]->heurredcostclock;
         calls = &(pricerdata->solvers[i]->heurredcostcalls);
      }

      /* get time limit */
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      if( !SCIPisInfinity(scip, timelimit) && timelimit - SCIPgetSolvingTime(scip) < 1 )
      {
         *nsols = 0;
         *status = SCIP_STATUS_TIMELIMIT;
      }

      if( pricerdata->solvers[i]->solversolve != NULL )
      {
         SCIP_CALL( SCIPstartClock(scip, clock) );

         SCIP_CALL( pricerdata->solvers[i]->solversolveheur(scip, pricerdata->solvers[i],
               pricerdata->pricingprobs[prob], prob, lowerbound,
               solvars, solvals, nsolvars, solisray, nsols, status) );

         SCIP_CALL( SCIPstopClock(scip, clock) );

         if( *status != SCIP_STATUS_UNKNOWN )
            (*calls)++;

         if( *status == SCIP_STATUS_OPTIMAL || *status == SCIP_STATUS_UNBOUNDED )
            break;
      }
   }

   return SCIP_OKAY;
}

/** computes the pricing problem objectives
 *  @todo this method could use more parameters as it is private
 */
static
SCIP_RETCODE setPricingObjs(
   SCIP*                 scip,               /**< SCIP data structure            */
   GCG_PRICETYPE         pricetype           /**< Farkas or Reduced cost pricing */
   )
{
   SCIP* origprob;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS** origconss;
   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_VAR** probvars;
   int nprobvars;

   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;
   SCIP_COL** cols;
   SCIP_Real* consvals;
   SCIP_Real dualsol;

   SCIP_VAR** consvars;
   int nconsvars;

   int i;
   int j;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   /* get the constraints of the master problem and the corresponding constraints in the original problem */
   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);
   origconss = GCGrelaxGetLinearOrigMasterConss(origprob);

   /* set objective value of all variables in the pricing problems to 0 (for farkas pricing) /
    * to the original objective of the variable (for redcost pricing)
    */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if( pricerdata->pricingprobs[i] == NULL )
         continue;
      probvars = SCIPgetVars(pricerdata->pricingprobs[i]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[i]);

      for( j = 0; j < nprobvars; j++ )
      {
         if( pricetype == GCG_PRICETYPE_FARKAS )
         {
            SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], 0.0) );
         }
         else
         {
            SCIP_VAR* origvar;
            origvar = GCGpricingVarGetOrigvars(probvars[j])[0];

            assert(GCGvarGetBlock(probvars[j]) == i);

            if( GCGvarIsLinking(origvar) )
            {
               SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], 0.0) );
            }
            else
            {
               assert( GCGvarGetBlock(origvar) == i);
               SCIP_CALL( SCIPchgVarObj(pricerdata->pricingprobs[i], probvars[j], SCIPvarGetObj(origvar)) );
            }
         }
      }
   }

   /* compute reduced cost for linking variable constraints and update objectives in the pricing problems */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if( pricerdata->pricingprobs[i] == NULL )
         continue;
      probvars = SCIPgetVars(pricerdata->pricingprobs[i]);
      nprobvars = SCIPgetNVars(pricerdata->pricingprobs[i]);

      for( j = 0; j < nprobvars; j++ )
      {
         SCIP_VAR* origvar;
         SCIP_CONS** linkconss;
#ifndef NDEBUG
         SCIP_VAR** pricingvars;
#endif
         assert(GCGvarIsPricing(probvars[j]));
         assert(GCGvarGetBlock(probvars[j]) == i);

         origvar = GCGpricingVarGetOrigvars(probvars[j])[0];
         if( !GCGvarIsLinking(origvar) )
            continue;

#ifndef NDEBUG
         pricingvars = GCGlinkingVarGetPricingVars(origvar);
#endif
         linkconss = GCGlinkingVarGetLinkingConss(origvar);
         assert(pricingvars[i] == probvars[j]);
         assert(linkconss[i] != NULL);

         /* redcost pricing */
         if( pricetype == GCG_PRICETYPE_REDCOST )
            dualsol = SCIPgetDualsolLinear(scip, linkconss[i]);
         /* farkas pricing */
         else
         {
            assert(pricetype == GCG_PRICETYPE_FARKAS);
            dualsol = SCIPgetDualfarkasLinear(scip, linkconss[i]);
         }

         /* add dual solution value to the pricing variable:
          * lambda variables get coef -1 in linking constraints --> add dualsol
          */
         SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[i], probvars[j], dualsol) );
      }
   }

   /* compute reduced cost and update objectives in the pricing problems */
   for( i = 0; i < nmasterconss; i++ )
   {
      /* redcost pricing */
      if( pricetype == GCG_PRICETYPE_REDCOST )
         dualsol = SCIPgetDualsolLinear(scip, masterconss[i]);
      /* farkas pricing */
      else
      {
         assert(pricetype == GCG_PRICETYPE_FARKAS);
         dualsol = SCIPgetDualfarkasLinear(scip, masterconss[i]);
      }
      if( !SCIPisZero(scip, dualsol) )
      {
#ifdef PRINTDUALSOLS
         SCIPdebugMessage("mastercons <%s> dualsol: %g\n", SCIPconsGetName(masterconss[i]), dualsol);
#endif

         /* for all variables in the constraint, modify the objective of the corresponding variable in a pricing problem */
         consvars = SCIPgetVarsLinear(origprob, origconss[i]);
         consvals = SCIPgetValsLinear(origprob, origconss[i]);
         nconsvars = SCIPgetNVarsLinear(origprob, origconss[i]);
         for( j = 0; j < nconsvars; j++ )
         {
            int blocknr;
            blocknr = GCGvarGetBlock(consvars[j]);
            assert(GCGvarIsOriginal(consvars[j]));
            /* nothing to be done if variable belongs to redundant block or variable was directly transferred to the master
             * or variable is linking variable (which means, the directly transferred copy is part of the master cons)
             */
            if( blocknr >= 0 && pricerdata->pricingprobs[blocknr] != NULL )
            {
               assert(GCGoriginalVarGetPricingVar(consvars[j]) != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[blocknr],
                     GCGoriginalVarGetPricingVar(consvars[j]), -1.0 * dualsol * consvals[j]) );
            }
         }
      }
   }

   /* get the cuts of the master problem and the corresponding cuts in the original problem */
   mastercuts = GCGsepaGetMastercuts(scip);
   nmastercuts = GCGsepaGetNMastercuts(scip);
   origcuts = GCGsepaGetOrigcuts(scip);

   assert(mastercuts != NULL);
   assert(origcuts != NULL);
   assert(GCGsepaGetNOrigcuts(scip) == nmastercuts);

   /* compute reduced cost and update objectives in the pricing problems */
   for( i = 0; i < nmastercuts; i++ )
   {
      /* farkas pricing */
      if( pricetype == GCG_PRICETYPE_REDCOST )
         dualsol = SCIProwGetDualsol(mastercuts[i]);
      /* redcost pricing */
      else
      {
         assert(pricetype == GCG_PRICETYPE_FARKAS);
         dualsol = SCIProwGetDualfarkas(mastercuts[i]);
      }
      if( !SCIPisZero(scip, dualsol) )
      {
         /* get columns and vals of the cut */
         nconsvars = SCIProwGetNNonz(origcuts[i]);
         cols = SCIProwGetCols(origcuts[i]);
         consvals = SCIProwGetVals(origcuts[i]);

         /* get the variables corresponding to the columns in the cut */
         SCIP_CALL( SCIPallocMemoryArray(scip, &consvars, nconsvars) );
         for( j = 0; j < nconsvars; j++ )
            consvars[j] = SCIPcolGetVar(cols[j]);

         /* for all variables in the cut, modify the objective of the corresponding variable in a pricing problem */
         for( j = 0; j < nconsvars; j++ )
         {
            int blocknr;
            blocknr = GCGvarGetBlock(consvars[j]);
            assert(GCGvarIsOriginal(consvars[j]));
            /* nothing to be done if variable belongs to redundant block or
             * variable was directly transferred to the master
             * or variable is linking variable (which means, the directly transferred copy is part of the master cut) */
            if( blocknr >= 0 && pricerdata->pricingprobs[blocknr] != NULL )
            {
               assert(GCGoriginalVarGetPricingVar(consvars[j]) != NULL);
               /* modify the objective of the corresponding variable in the pricing problem */
               SCIP_CALL( SCIPaddVarObj(pricerdata->pricingprobs[blocknr],
                     GCGoriginalVarGetPricingVar(consvars[j]), -1.0 * dualsol * consvals[j]) );
            }
         }
         SCIPfreeMemoryArray(scip, &consvars);
      }
   }

   /* get dual solutions / farkas values of the convexity constraints */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      assert( GCGrelaxIsPricingprobRelevant(origprob, i) == (GCGrelaxGetConvCons(origprob, i) != NULL) );
      if( !GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         pricerdata->dualsolconv[i] = -1.0 * SCIPinfinity(scip);
         continue;
      }
      if( pricetype == GCG_PRICETYPE_REDCOST )
         pricerdata->dualsolconv[i] = SCIPgetDualsolLinear(scip, GCGrelaxGetConvCons(origprob, i));
      else
      {
         assert(pricetype == GCG_PRICETYPE_FARKAS);
         pricerdata->dualsolconv[i] = SCIPgetDualfarkasLinear(scip, GCGrelaxGetConvCons(origprob, i));
      }
#ifdef PRINTDUALSOLS
      if( GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         SCIPdebugMessage("convcons <%s> dualsol: %g\n", SCIPconsGetName(GCGrelaxGetConvCons(origprob, i)), pricerdata->dualsolconv[i]);
      }
#endif
   }

   return SCIP_OKAY;
}

/** add master variable to all constraints */
static
SCIP_RETCODE addVariableToMasterconstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< Pricerdata data structure */
   SCIP_VAR*             newvar,             /**< The new variable to add */
   int                   prob,               /**< number of the pricing problem the solution belongs to */
   SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars            /**< number of variables in array solvars */
   )
{
   int i;
   int c;
   int idx;

   SCIP_CONS** masterconss;
   int nmasterconss;
   SCIP_Real* mastercoefs;
   SCIP_CONS* linkcons;

   nmasterconss = GCGrelaxGetNMasterConss(pricerdata->origprob);
   masterconss = GCGrelaxGetMasterConss(pricerdata->origprob);

   SCIP_CALL( SCIPallocBufferArray(scip, &mastercoefs, nmasterconss) );
   BMSclearMemoryArray(mastercoefs, nmasterconss);

   /* compute coef of the variable in the master constraints */
   for( i = 0; i < nsolvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         SCIP_CONS** linkconss;
         SCIP_VAR** origvars;
         SCIP_Real* coefs;
         int ncoefs;

         assert(GCGvarIsPricing(solvars[i]));
         origvars = GCGpricingVarGetOrigvars(solvars[i]);
         assert(GCGvarIsOriginal(origvars[0]));

         coefs = GCGoriginalVarGetCoefs(origvars[0]);
         ncoefs = GCGoriginalVarGetNCoefs(origvars[0]);
         assert(!SCIPisInfinity(scip, solvals[i]));

         /* original variable is a linking variable, just add it to the linkcons */
         if( GCGvarIsLinking(origvars[0]) )
         {
#ifndef NDEBUG
            SCIP_VAR** pricingvars;
            pricingvars = GCGlinkingVarGetPricingVars(origvars[0]);
#endif
            linkconss = GCGlinkingVarGetLinkingConss(origvars[0]);

            assert(pricingvars[prob] == solvars[i]);
            assert(linkconss[prob] != NULL);
            SCIP_CALL( SCIPaddCoefLinear(scip, linkconss[prob], newvar, -solvals[i]) );
            continue;
         }

         /* for each coef, add coef * solval to the coef of the new variable for the corresponding constraint */
         for( c = 0; c < ncoefs; c++ )
         {
            linkconss = GCGoriginalVarGetMasterconss(origvars[0]);
            assert(!SCIPisZero(scip, coefs[c]));
            SCIP_CALL( SCIPgetTransformedCons(scip, linkconss[c], &linkcons) );

            idx = (int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, linkcons); /*lint !e507*/
            assert(0 <= idx && idx < nmasterconss);
            assert(masterconss[idx] == linkcons);
            mastercoefs[idx] += coefs[c] * solvals[i];
         }

      }
   }

   /* add the variable to the master constraints */
   for( i = 0; i < nmasterconss; i++ )
   {
      if( !SCIPisZero(scip, mastercoefs[i]) )
      {
         assert(!SCIPisInfinity(scip, mastercoefs[i]) && !SCIPisInfinity(scip, -mastercoefs[i]));
         SCIP_CALL( SCIPaddCoefLinear(scip, masterconss[i], newvar, mastercoefs[i]) );
      }
   }

   SCIPfreeBufferArray(scip, &mastercoefs);
   return SCIP_OKAY;
}



/** add variable with computed coefficients to the master cuts */
static
SCIP_RETCODE addVariableToMastercuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             newvar,             /**< The new variable to add */
   int                   prob,               /**< number of the pricing problem the solution belongs to */
   SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars            /**< number of variables in array solvars */
   )
{
   SCIP_ROW** mastercuts;
   int nmastercuts;
   SCIP_ROW** origcuts;

   SCIP_COL** cols;
   SCIP_Real conscoef;
   SCIP_VAR* var;
   SCIP_Real* consvals;

   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(newvar != NULL);
   assert(solvars != NULL);
   assert(solvals != NULL);

   /* get the cuts of the master problem and the corresponding cuts in the original problem */
   mastercuts = GCGsepaGetMastercuts(scip);
   nmastercuts = GCGsepaGetNMastercuts(scip);
   origcuts = GCGsepaGetOrigcuts(scip);

   assert(mastercuts != NULL);
   assert(origcuts != NULL);
   assert(GCGsepaGetNOrigcuts(scip) == nmastercuts);

   /* compute coef of the variable in the cuts and add it to the cuts */
   for( i = 0; i < nmastercuts; i++ )
   {
      if( !SCIProwIsInLP(mastercuts[i]) )
         continue;

      /* get columns of the cut and their coefficients */
      cols = SCIProwGetCols(origcuts[i]);
      consvals = SCIProwGetVals(origcuts[i]);

      conscoef = 0;

      for( j = 0; j < SCIProwGetNNonz(origcuts[i]); j++ )
      {
         int blocknr;
         var = SCIPcolGetVar(cols[j]);
         blocknr = GCGvarGetBlock(var);
         assert(GCGvarIsOriginal(var));

         /* if the belongs to the same block and is no linking variable, update the coef */
         if( blocknr == prob )
            for( k = 0; k < nsolvars; k++ )
               if( solvars[k] == GCGoriginalVarGetPricingVar(var) )
               {
                  conscoef += ( consvals[j] * solvals[k] );
                  break;
               }
      }

      if( !SCIPisZero(scip, conscoef) )
         SCIP_CALL( SCIPaddVarToRow(scip , mastercuts[i], newvar, conscoef) );
   }

   return SCIP_OKAY;
}

/** adds new variable to the end of the priced variables array */
static
SCIP_RETCODE addVariableToPricedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< pricer data structure */
   SCIP_VAR*             newvar              /**< variable to add */
   )
{
   SCIP_CALL( ensureSizePricedvars(scip, pricerdata, pricerdata->npricedvars + 1) );
   pricerdata->pricedvars[pricerdata->npricedvars] = newvar;
   pricerdata->npricedvars++;

   return SCIP_OKAY;
}

/** creates a new master variable corresponding to the given solution and problem */
static
SCIP_RETCODE createNewMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            solvars,            /**< array of variables with non-zero value in the solution of the pricing problem */
   SCIP_Real*            solvals,            /**< array of values in the solution of the pricing problem for variables in array solvars*/
   int                   nsolvars,           /**< number of variables in array solvars */
   SCIP_Bool             solisray,           /**< is the solution a ray? */
   int                   prob,               /**< number of the pricing problem the solution belongs to */
   SCIP_Bool             force,              /**< should the given variable be added also if it has non-negative reduced cost? */
   SCIP_Bool*            added,              /**< pointer to store whether the variable was successfully added */
   SCIP_VAR**            addedvar            /**< pointer to store the created variable */
   )
{
   SCIP* origprob;
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   char varname[SCIP_MAXSTRLEN];

   SCIP_Real objcoeff;
   SCIP_VAR* newvar;

   SCIP_VARDATA* vardata;
   long long int nodenumber;

   SCIP_Real objvalue;
   SCIP_Real redcost;
   SCIP_Real gap;
   SCIP_Real origgap;
   int i;

   assert(scip != NULL);
   assert(solvars != NULL);
   assert(solvals != NULL);
   assert(nsolvars >= 0);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   if( addedvar != NULL )
      *addedvar = NULL;

   objvalue = 0.0;
   redcost = 0.0;

   if( !force )
   {
      /* compute the objective function value of the solution */
      for( i = 0; i < nsolvars; i++ )
         objvalue += solvals[i] * SCIPvarGetObj(solvars[i]);

      /* compute reduced cost of variable (i.e. subtract dual solution of convexity constraint, if solution corresponds to a point) */
      redcost = ( solisray ? objvalue : objvalue - pricerdata->dualsolconv[prob]);

      if( !SCIPisSumNegative(scip, redcost) )
      {
         SCIPdebugMessage("var with redcost %g (objvalue = %g, dualsol =%g) was not added\n", redcost, objvalue, pricerdata->dualsolconv[prob]);
         *added = FALSE;

         return SCIP_OKAY;
      }
      SCIPdebugMessage("found var with redcost %g (objvalue = %g, dualsol =%g)\n", redcost, objvalue, pricerdata->dualsolconv[prob]);
   }
   else
   {
      SCIPdebugMessage("force var (objvalue = %g, dualsol =%g)\n",  objvalue, pricerdata->dualsolconv[prob]);
   }

   *added = TRUE;

   /* compute objective coefficient of the variable */
   objcoeff = 0;
   for( i = 0; i < nsolvars; i++ )
   {
      if( !SCIPisZero(scip, solvals[i]) )
      {
         SCIP_VAR* origvar;

         assert(GCGvarIsPricing(solvars[i]));
         origvar = GCGpricingVarGetOrigvars(solvars[i])[0];

         /* original variable is linking variable --> directly transferred master variable got the full obj,
          * priced-in variables get no objective value for this origvar */
         if( GCGvarIsLinking(origvar) )
            continue;

         /* add quota of original variable's objcoef to the master variable's coef */
         objcoeff += solvals[i] * SCIPvarGetObj(origvar);
      }

   }



   if( SCIPisInfinity(scip, objcoeff) )
   {
      SCIPwarningMessage(scip, "variable with infinite objective value found in pricing, change objective to SCIPinfinity()/2\n");
      objcoeff = SCIPinfinity(scip) / 2;
   }

   if( solisray )
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "r_%d_%d", prob, pricerdata->nraysprob[prob]);
      pricerdata->nraysprob[prob]++;
   }
   else
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "p_%d_%d", prob, pricerdata->npointsprob[prob]);
      pricerdata->npointsprob[prob]++;
   }

   SCIP_CALL( GCGcreateMasterVar(scip, pricerdata->pricingprobs[prob], &newvar, varname, objcoeff,
         pricerdata->vartype, solisray, prob, nsolvars, solvals, solvars));

   SCIPvarMarkDeletable(newvar);

   SCIP_CALL( SCIPcatchVarEvent(scip, newvar, SCIP_EVENTTYPE_VARDELETED,
         pricerdata->eventhdlr, NULL, NULL) );


   /* add variable */
   if( !force )
   {
      SCIP_CALL( SCIPaddPricedVar(scip, newvar, pricerdata->dualsolconv[prob] - objvalue) );
   }
   else
   {
      SCIP_CALL( SCIPaddVar(scip, newvar) );
   }

   SCIP_CALL( addVariableToPricedvars(scip, pricerdata, newvar) );
   SCIP_CALL( addVariableToMasterconstraints(scip, pricerdata, newvar, prob, solvars, solvals, nsolvars) );
   SCIP_CALL( addVariableToMastercuts(scip, newvar, prob, solvars, solvals, nsolvars) );

   /* add variable to convexity constraint */
   if( !solisray )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, GCGrelaxGetConvCons(origprob, prob), newvar, 1.0) );
   }

   if( addedvar != NULL )
   {

      *addedvar = newvar;
   }
   nodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(origprob));
   vardata = SCIPvarGetData(newvar);
   GCGsetCreationNode(origprob, vardata, nodenumber);
   GCGsetCreationTime(origprob, vardata, SCIPgetSolvingTime(scip));
   GCGsetIteration(origprob, vardata, SCIPgetNLPIterations(scip));

   origgap = SCIPgetGap(origprob);
   gap = SCIPgetGap(scip);
   GCGsetGap(origprob, vardata, MIN(origgap, gap));
   GCGsetRedcost(origprob, vardata, redcost);


   return SCIP_OKAY;
}


/**
 * check whether pricing can be aborted:
 * if objective value is always integral and the current node's current
 * lowerbound rounded up equals the current lp objective value rounded
 * up we don't need to continue pricing since the best possible feasible
 * solution must have at least this value
 */
static
SCIP_Bool canPricingBeAborted(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata          /**< pricerdata data structure */
   )
{
   SCIP_Bool canabort;
   assert(scip != NULL);
   assert(pricerdata != NULL);
   canabort = FALSE;
   if( pricerdata->abortpricingint && SCIPisObjIntegral(scip )
      && SCIPisEQ(scip, SCIPceil(scip, SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip))), SCIPceil(scip, SCIPgetLPObjval(scip))) /* && SCIPgetNNodes(scip) > 1 ??????*/)
   {
      GCGpricerPrintInfo(scip, pricerdata,
            "pricing aborted due to integral objective: node LB = %g, LP obj = %g\n",
            SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip)), SCIPgetLPObjval(scip));

      canabort = TRUE;
   }
   if( pricerdata->abortpricinggap > 0 )
   {
      SCIP_Real gap;
      gap = (SCIPgetLPObjval(scip) - SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip)))/SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip));
      gap = ABS(gap);

      if( gap < pricerdata->abortpricinggap )
      {
         GCGpricerPrintInfo(scip, pricerdata,
               "pricing aborted due to small gap: node LB = %g, LP obj = %g, gap = %g\n",
               SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip)), SCIPgetLPObjval(scip), gap);
         canabort = TRUE;
      }
   }

   return canabort;
}

/** sorts pricing problems according to their score */
static
void sortPricingProblemsByScore(
   SCIP_PRICERDATA*      pricerdata          /**< pricerdata data structure */
)
{
   int i;
   assert(pricerdata != NULL);
    /** @todo sort w.r.t. other measures? Don't sort in Farkas pricing? Randomized? */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      pricerdata->permu[i] = i;
      switch( pricerdata->sorting )
      {
      case 1:
         pricerdata->score[i] = pricerdata->dualsolconv[i];
         break;
      case 2:
         pricerdata->score[i] = -(0.2 * pricerdata->npointsprob[i] + pricerdata->nraysprob[i]);
         break;
      default:
         pricerdata->score[i] = 0.0;
      }
   }

   if( pricerdata->sorting > 0 )
      SCIPsortDownRealInt(pricerdata->score, pricerdata->permu, pricerdata->npricingprobs);
}

/** indicates whether heuristic pricing can be aborted */
static
SCIP_Bool abortHeuristicPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   GCG_PRICETYPE         pricetype,          /**< type of pricing */
   int                   nfoundvars,         /**< number of found variables */
   int                   solvedmips,         /**< number of solved mips */
   int                   successfulmips      /**< number of successfully solved mips so far */
   )
{

   assert( pricetype == GCG_PRICETYPE_FARKAS || pricetype == GCG_PRICETYPE_REDCOST );

   if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      return !((nfoundvars < pricerdata->maxvarsroundredcost)
         && successfulmips < pricerdata->maxsuccessfulmipsredcost
         && successfulmips < pricerdata->successfulmipsrel * pricerdata->npricingprobsnotnull
         && (nfoundvars == 0 ||
            solvedmips < pricerdata->mipsrelredcost * pricerdata->npricingprobsnotnull ));
   }
   else
   {
      assert(pricetype == GCG_PRICETYPE_FARKAS);
      return !(nfoundvars < pricerdata->maxvarsroundfarkas
         && (nfoundvars == 0 || solvedmips < pricerdata->mipsrelfarkas * pricerdata->npricingprobsnotnull));
   }
}

/** set subproblem timelimit */
static
SCIP_RETCODE subproblemSetTimelimit(
   SCIP*                 scip,               /**< SCIP data structure*/
   SCIP*                 pricingscip,        /**< SCIP of the pricingproblem */
   int                   prob,               /**< number of the pricing problem */
   SCIP_Real*            timelimit           /**< pointer to store timelimit */
   )
{
   /* set time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", timelimit) );
   if( !SCIPisInfinity(scip, *timelimit) )
   {
      if( *timelimit - SCIPgetSolvingTime(scip) > 0 )
      {
         SCIP_CALL( SCIPsetRealParam(pricingscip, "limits/time", *timelimit - SCIPgetSolvingTime(scip)) );
         SCIPdebugMessage("Tilim for pricing %d is %f\n", prob, *timelimit- SCIPgetSolvingTime(scip));
      }
      else
      {
         SCIPdebugMessage("Tilim for pricing %d is < 0\n", prob);
      }
   }
   return SCIP_OKAY;
}

/** free pricing problems */
static
SCIP_RETCODE freePricingProblems(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata          /**< pricerdata data structure */
   )
{
   int j;
   assert(pricerdata != NULL);
   assert(pricerdata->pricingprobs != NULL);

   for( j = 0; j < pricerdata->npricingprobs; j++ )
      if( pricerdata->pricingprobs[j] != NULL
         && SCIPgetStage(pricerdata->pricingprobs[j]) > SCIP_STAGE_PROBLEM)
         {
            SCIP_CALL( SCIPstartClock(scip, pricerdata->freeclock) );
            SCIP_CALL( SCIPfreeTransform(pricerdata->pricingprobs[j]) );
            SCIP_CALL( SCIPstopClock(scip, pricerdata->freeclock) );
         }

   return SCIP_OKAY;
}

/** performs heuristic pricing */
static
SCIP_RETCODE performHeuristicPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   GCG_PRICETYPE         pricetype,          /**< type of the pricing */
   int*                  nfoundvars,         /**< pointer to store the number of found vars */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   int i;
   int j;
   int prob;
   int solvedmips;
   int successfulmips;
   int nfoundvarsprob;
   SCIP_Real pricinglowerbound;

   SCIP_VAR*** solvars;
   SCIP_Real** solvals;
   int* nsolvars;
   int nsols;
   SCIP_Bool* solisray;
   SCIP_STATUS status;
   SCIP_Real timelimit;

   SCIPdebugMessage("heuristical pricing\n");
   solvedmips = 0;
   successfulmips = 0;


   /* solve the pricing MIPs heuristically and check whether solutions
    * corresponding to variables with negative reduced costs where found
    */
   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      if( abortHeuristicPricing(scip, pricerdata, pricetype, *nfoundvars, solvedmips, successfulmips) )
         break;

      prob = pricerdata->permu[i];

      if( pricerdata->pricingprobs[prob] == NULL )
         continue;

      SCIP_CALL( subproblemSetTimelimit(scip, pricerdata->pricingprobs[prob], prob, &timelimit) );

      if( timelimit - SCIPgetSolvingTime(scip) <= 0 )
      {
         if( result != NULL )
            *result = SCIP_DIDNOTRUN;

         return SCIP_OKAY;
      }

      /* set objective limit, such that only solutions with negative reduced costs are accepted */
      SCIP_CALL( SCIPsetObjlimit(pricerdata->pricingprobs[prob], pricerdata->dualsolconv[prob]) );

      pricerdata->solvedsubmipsheur++;
      solvedmips++;

      SCIP_CALL( solvePricingProblemHeur(scip, pricerdata, prob, pricetype, &pricinglowerbound, &solvars, &solvals,
            &nsolvars, &solisray, &nsols, &status) );

      nfoundvarsprob = 0;

      for( j = 0; j < nsols && nfoundvarsprob <= pricerdata->maxsolsprob &&
              (pricetype == GCG_PRICETYPE_REDCOST || *nfoundvars < pricerdata->maxvarsroundfarkas)
              && (pricetype == GCG_PRICETYPE_FARKAS || *nfoundvars < pricerdata->maxvarsroundredcost); j++ )
      {
         SCIP_Bool added;
         /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
         SCIP_CALL( createNewMasterVar(scip, solvars[j], solvals[j], nsolvars[j], solisray[j], prob,
               FALSE, &added, NULL) );

         if( added )
         {
            ++(*nfoundvars);
            nfoundvarsprob++;

            if( nfoundvarsprob == 1 )
               successfulmips++;
         }
      }
   }

   /* free the pricingproblems if they exist and need to be freed */
   SCIP_CALL( freePricingProblems(scip, pricerdata) );

   return SCIP_OKAY;
}

static
/** returns TRUE if optimal pricing can be aborted, FALSE otherwise*/
SCIP_Bool abortOptimalPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   GCG_PRICETYPE         pricetype,          /**< type of pricing*/
   int                   nfoundvars,         /**< number of variables found so far */
   int                   solvedmips,         /**< number of MIPS solved so far */
   int                   successfulmips      /**< number of sucessful mips solved so far */
   )
{
   SCIP_Bool root;
   root = isRootNode(scip);

   if( pricetype == GCG_PRICETYPE_FARKAS )
      return !(nfoundvars < pricerdata->maxvarsroundfarkas
         && (nfoundvars == 0 || solvedmips < pricerdata->mipsrelfarkas * pricerdata->npricingprobsnotnull));
   else if( pricetype == GCG_PRICETYPE_REDCOST )
      return !(((((nfoundvars < pricerdata->maxvarsroundredcostroot) || !root ) && ((nfoundvars < pricerdata->maxvarsroundredcost) || root)))
            && successfulmips < pricerdata->maxsuccessfulmipsredcost
            && successfulmips < pricerdata->successfulmipsrel * pricerdata->npricingprobsnotnull
            && (nfoundvars == 0 || ( (root || solvedmips < pricerdata->mipsrelredcost * pricerdata->npricingprobsnotnull)
                  && (!root || solvedmips < pricerdata->mipsrelredcostroot * pricerdata->npricingprobsnotnull))));
   return FALSE;

}
static
/** performs optimal pricing */
SCIP_RETCODE performOptimalPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICERDATA*      pricerdata,         /**< pricerdata data structure */
   GCG_PRICETYPE         pricetype,          /**< type of pricing */
   SCIP_RESULT*          result,             /**< result pointer */
   int*                  nfoundvars,         /**< pointer to store number of found variables */
   SCIP_Real*            bestredcost,        /**< pointer to store reduced cost */
   SCIP_Bool*            bestredcostvalid    /**< pointer to store whether the reduced cost returned is valid */
   )
{
   int i;
   int j;
   int prob;
   int solvedmips;
   int successfulmips;
   int nfoundvarsprob;
   SCIP_Real timelimit;
   SCIP* origprob;
   SCIP_Bool added;

   SCIP_Real pricinglowerbound;

   SCIP_VAR*** solvars;
   SCIP_Real** solvals;
   int* nsolvars;
   int nsols;
   SCIP_Bool* solisray;
   SCIP_STATUS status;
   SCIP_Bool root;

   assert(scip != NULL);
   assert(pricerdata != NULL);
   assert(pricetype == GCG_PRICETYPE_FARKAS || result != NULL);

   assert(nfoundvars != NULL);
   assert(bestredcost != NULL);
   assert(bestredcostvalid != NULL);

   SCIPdebugMessage("optimal pricing\n");

   origprob = pricerdata->origprob;
   root = isRootNode(scip);
   solvedmips = 0;
   successfulmips = 0;

   *bestredcost = 0.0;
   *bestredcostvalid = ( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL ? TRUE : FALSE );

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      if( abortOptimalPricing(scip, pricerdata, pricetype, *nfoundvars, solvedmips, successfulmips) )
         break;

      prob = pricerdata->permu[i];

      if( pricerdata->pricingprobs[prob] == NULL )
         continue;

      /* @todo set objective limit, such that only solutions with negative reduced costs are accepted? */
      /* SCIP_CALL( SCIPsetObjlimit(pricerdata->pricingprobs[prob], pricerdata->dualsolconv[prob]) ); */

      SCIP_CALL( subproblemSetTimelimit(scip, pricerdata->pricingprobs[prob], prob, &timelimit) );

      if( timelimit - SCIPgetSolvingTime(scip) <= 0 )
      {
         if( result != NULL )
            *result = SCIP_DIDNOTRUN;

         *bestredcostvalid = FALSE;
         return SCIP_OKAY;
      }

      SCIP_CALL( solvePricingProblem(scip, pricerdata, prob, pricetype, &pricinglowerbound, &solvars, &solvals,
            &nsolvars, &solisray, &nsols, &status) );

      pricerdata->solvedsubmipsoptimal++;
      solvedmips++;
      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL && !SCIPisInfinity(scip, pricinglowerbound) && SCIPgetStatus(pricerdata->pricingprobs[prob]) == SCIP_STATUS_OPTIMAL )
      {
         assert( !SCIPisSumPositive(scip, pricinglowerbound - pricerdata->dualsolconv[prob]) );
      }

      (*bestredcost) += GCGrelaxGetNIdenticalBlocks(origprob, prob) * (pricinglowerbound - pricerdata->dualsolconv[prob]);

      if( status != SCIP_STATUS_OPTIMAL )
      {
         *bestredcostvalid = FALSE;
         if( result != NULL )
            *result = SCIP_DIDNOTRUN;
      }
      nfoundvarsprob = 0;

      for( j = 0; j < nsols && nfoundvarsprob <= pricerdata->maxsolsprob &&
              (pricetype == GCG_PRICETYPE_REDCOST || *nfoundvars < pricerdata->maxvarsroundfarkas)
              && (pricetype == GCG_PRICETYPE_FARKAS || ((*nfoundvars < pricerdata->maxvarsroundredcost || root ) && (*nfoundvars < pricerdata->maxvarsroundredcostroot || !root))); j++ )
      {
         /* create new variable, compute objective function value and add it to the master constraints and cuts it belongs to */
         SCIP_CALL( createNewMasterVar(scip, solvars[j], solvals[j], nsolvars[j], solisray[j], prob,
               FALSE, &added, NULL) );

         if( added )
         {
            ++(*nfoundvars);
            nfoundvarsprob++;

            if( nfoundvarsprob == 1 )
               successfulmips++;
         }
      }
   }

   /** @todo perhaps solve remaining pricing problems, if only few left? */
   /** @todo solve all pricing problems all k iterations? */
   /* this makes sure that if a pricing problem has not been solved, the langrangian bound cannot be calculated */
   for( j = i; j < pricerdata->npricingprobs && bestredcostvalid; j++ )
      if( pricerdata->pricingprobs[pricerdata->permu[j]] != NULL )
         *bestredcostvalid = FALSE;

   /* free the pricingproblems if they exist and need to be freed */
   SCIP_CALL( freePricingProblems(scip, pricerdata) );

   return SCIP_OKAY;
}

/** performs the pricing routine, gets the type of pricing that should be done: farkas or redcost pricing */
static
SCIP_RETCODE performPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer,             /**< the pricer */
   GCG_PRICETYPE         pricetype,          /**< type of the pricing */
   SCIP_RESULT*          result,             /**< result pointer */
   SCIP_Real*            lowerbound          /**< lowerbound pointer */
   )
{
   SCIP_PRICERDATA* pricerdata;
   int nfoundvars;
#ifdef ENABLESTATISTICS
   double degeneracy;
#endif
   SCIP_Real bestredcost;
   SCIP_Bool bestredcostvalid;

   assert(scip != NULL);
   assert(pricer != NULL);

   /* get pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   assert(result != NULL || pricetype == GCG_PRICETYPE_FARKAS);
   assert(lowerbound != NULL || pricetype == GCG_PRICETYPE_FARKAS);

   if( lowerbound != NULL )
      *lowerbound = -SCIPinfinity(scip);

   GCGpricerPrintInfo(scip, pricerdata, "nvars = %d, current LP objval = %g, time = %f, node = %lld\n",
         SCIPgetNVars(scip), SCIPgetLPObjval(scip), SCIPgetSolvingTime(scip), SCIPgetNNodes(scip));

   if( pricetype == GCG_PRICETYPE_REDCOST )
   {
      assert(result != NULL);

      /* terminate early, if applicable */
      if( canPricingBeAborted(scip, pricerdata) )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }

      pricerdata->redcostcalls++;
      *result = SCIP_SUCCESS;
   }
   else if( pricetype == GCG_PRICETYPE_FARKAS )
   {
      pricerdata->farkascalls++;
   }

   pricerdata->calls++;
   nfoundvars = 0;

   /* set objectives of the variables in the pricing sub-MIPs */
   SCIP_CALL( setPricingObjs(scip, pricetype) );

   sortPricingProblemsByScore(pricerdata);

   bestredcost = 0.0;
   bestredcostvalid = FALSE;

   if( pricerdata->useheurpricing )
   {
      SCIP_CALL( performHeuristicPricing(scip, pricerdata, pricetype, &nfoundvars, result) );
   }

   /* if no variables were found so far, solve the pricing MIPs to optimality and check whether
    * solutions corresponding to variables with negative reduced costs where found
    */
   if( nfoundvars == 0 )
   {
      SCIP_CALL( performOptimalPricing(scip, pricerdata, pricetype, result, &nfoundvars, &bestredcost, &bestredcostvalid) );
   }

   if( nfoundvars == 0 && isRootNode(scip) )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );
   }

   if( pricetype == GCG_PRICETYPE_REDCOST && bestredcostvalid )
   {
      assert(lowerbound != NULL);
      GCGpricerPrintInfo(scip, pricerdata, "lower bound = %g, bestredcost = %g\n", SCIPgetLPObjval(scip) + bestredcost, bestredcost);

      *lowerbound = SCIPgetLPObjval(scip) + bestredcost;
   }

   SCIPdebugMessage("%s pricing: found %d new vars\n", (pricetype == GCG_PRICETYPE_REDCOST ? "Redcost" : "Farkas"), nfoundvars);

#ifdef ENABLESTATISTICS
   SCIP_CALL( computeCurrentDegeneracy(scip, &degeneracy) );

   if( pricerdata->lastnode != SCIPgetCurrentNode(scip) )
   {
      pricerdata->lastnode = SCIPgetCurrentNode(scip);

      SCIP_CALL( ensureSizeAvgnodedegeneracy(scip, pricerdata) );
      assert(pricerdata->nnodes < pricerdata->maxnnodes);

      if( pricerdata->nnodes == 1 )
         pricerdata->avgnodedegeneracy = degeneracy;
      else if( pricerdata->nnodes > 2 )
      {
         /* Complicated calculation for numerical stability:
          *     E[\sum_{i=1}^n x_i] = (E[\sum_{i=1}^{n-1} x_i]*(n-1) + x_n)/n
          *     E[\sum_{i=1}^n x_i] = E[\sum_{i=1}^{n-1} x_i]*(n-1)/n + x_n/n
          * <=> E[\sum_{i=1}^n x_i] = E[\sum_{i=1}^{n-1} x_i]-E[\sum_{i=1}^{n-1} x_i]/n + x_n/n
          * <=> E_n = E_{n-1} - E_{n-1}/n + x_n/n
          * <=> E -= E/n - x_n(n
          */
         pricerdata->avgnodedegeneracy -= pricerdata->avgnodedegeneracy/(pricerdata->nnodes-2) - pricerdata->nodedegeneracy[pricerdata->nnodes-1]/(pricerdata->nnodes-1);
      }
   }
   if( pricerdata->lastnode == SCIPgetRootNode(scip) )
      pricerdata->rootnodedegeneracy = degeneracy;

   pricerdata->nodedegeneracy[pricerdata->nnodes] = degeneracy;
   ++(pricerdata->nnodes);
   assert(pricerdata->lastnode == SCIPgetCurrentNode(scip));
#endif
   return SCIP_OKAY;
}

/*
 * Callback methods of variable pricer
 */


/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeGcg)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   /* get pricerdata */
   pricerdata = SCIPpricerGetData(pricer);

   SCIP_CALL( solversFree(scip, pricerdata) );

   SCIPfreeMemoryArray(scip, &pricerdata->solvers);

   /* free memory for pricerdata*/
   if( pricerdata != NULL )
   {
      SCIPfreeMemory(scip, &pricerdata);
   }

   SCIPpricerSetData(pricer, NULL);
   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitGcg)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   /* get pricerdata */
   pricerdata = SCIPpricerGetData(pricer);

   SCIP_CALL( solversInit(scip, pricerdata) );

   return SCIP_OKAY;
}


/** deinitialization method of variable pricer (called before transformed problem is freed) */
static
SCIP_DECL_PRICEREXIT(pricerExitGcg)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   /* get pricerdata */
   pricerdata = SCIPpricerGetData(pricer);

   SCIP_CALL( solversExit(scip, pricerdata) );

   return SCIP_OKAY;
}


/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
static
SCIP_DECL_PRICERINITSOL(pricerInitsolGcg)
{
   SCIP_PRICERDATA* pricerdata;
   int i;
   SCIP* origprob;
   int norigvars;
   SCIP_Bool discretization;
   SCIP_CONS** masterconss;
   int nmasterconss;
   int origverblevel;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   /* at the beginning, the output of the master problem gets the same verbosity level
    * as the output of the original problem */
   SCIP_CALL( SCIPgetIntParam(origprob, "display/verblevel", &origverblevel) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", origverblevel) );

   pricerdata->currnodenr = -1;

   nmasterconss = GCGrelaxGetNMasterConss(origprob);
   masterconss = GCGrelaxGetMasterConss(origprob);

   /* init array containing all pricing problems */
   pricerdata->npricingprobs = GCGrelaxGetNPricingprobs(origprob);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->pricingprobs), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->npointsprob), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nraysprob), pricerdata->npricingprobs) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->farkascallsdist), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->farkasfoundvars), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->farkasnodetimedist), pricerdata->npricingprobs) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcostcallsdist), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcostfoundvars), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->redcostnodetimedist), pricerdata->npricingprobs) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nodetimehist), PRICER_STAT_ARRAYLEN_TIME) ); /*lint !e506*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->foundvarshist), PRICER_STAT_ARRAYLEN_VARS) ); /*lint !e506*/

   BMSclearMemoryArray(pricerdata->nodetimehist, PRICER_STAT_ARRAYLEN_TIME);
   BMSclearMemoryArray(pricerdata->foundvarshist, PRICER_STAT_ARRAYLEN_VARS);

   pricerdata->oldvars=0;

   pricerdata->npricingprobsnotnull = 0;

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {

      pricerdata->farkascallsdist[i] = 0;
      pricerdata->farkasfoundvars[i] = 0;
      pricerdata->farkasnodetimedist[i] = 0;
      pricerdata->redcostcallsdist[i] = 0;
      pricerdata->redcostfoundvars[i] = 0;
      pricerdata->redcostnodetimedist[i]= 0;


      if( GCGrelaxIsPricingprobRelevant(origprob, i) )
      {
         pricerdata->pricingprobs[i] = GCGrelaxGetPricingprob(origprob, i);
         pricerdata->npricingprobsnotnull++;
      }
      else
      {
         pricerdata->pricingprobs[i] = NULL;
      }
      pricerdata->npointsprob[i] = 0;
      pricerdata->nraysprob[i] = 0;
   }

   /* alloc memory for arrays of reduced cost */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->dualsolconv), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->score), pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->permu), pricerdata->npricingprobs) );

   /* alloc memory for solution values of variables in pricing problems */
   norigvars = SCIPgetNOrigVars(origprob);
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->solvals), norigvars) );

   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->redcostclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->farkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->freeclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->transformclock)) );

   pricerdata->solvedsubmipsoptimal = 0;
   pricerdata->solvedsubmipsheur = 0;
   pricerdata->calls = 0;
   pricerdata->redcostcalls = 0;
   pricerdata->farkascalls = 0;

   /* set variable type for master variables */
   SCIP_CALL( SCIPgetBoolParam(origprob, "relaxing/gcg/discretization", &discretization) );
   if( discretization )
   {
      pricerdata->vartype = SCIP_VARTYPE_INTEGER;
   }
   else
   {
      pricerdata->vartype = SCIP_VARTYPE_CONTINUOUS;
   }

   SCIP_CALL( SCIPhashmapCreate(&(pricerdata->mapcons2idx), SCIPblkmem(scip), 10 * nmasterconss +1) );
   for( i = 0; i < nmasterconss; i++ )
   {
      SCIP_CALL( SCIPhashmapInsert(pricerdata->mapcons2idx, masterconss[i], (void*)(size_t)i) );
      assert((int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, masterconss[i]) == i); /*lint !e507*/
   }

   pricerdata->npricedvars = 0;
   pricerdata->maxpricedvars = 50;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pricerdata->pricedvars, pricerdata->maxpricedvars) );

   pricerdata->rootnodedegeneracy = 0.0;
   pricerdata->avgnodedegeneracy = 0.0;
   pricerdata->nnodes = 0;
   pricerdata->maxnnodes = 10;
   pricerdata->lastnode = NULL;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pricerdata->nodedegeneracy), pricerdata->maxnnodes) );

   SCIP_CALL( solversInitsol(scip, pricerdata) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolGcg)
{
   SCIP_PRICERDATA* pricerdata;
   int i;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   SCIPhashmapFree(&(pricerdata->mapcons2idx));

   SCIPfreeMemoryArray(scip, &(pricerdata->pricingprobs));
   SCIPfreeMemoryArray(scip, &(pricerdata->dualsolconv));
   SCIPfreeMemoryArray(scip, &(pricerdata->score));
   SCIPfreeMemoryArray(scip, &(pricerdata->permu));
   SCIPfreeMemoryArray(scip, &(pricerdata->solvals));
   SCIPfreeMemoryArray(scip, &(pricerdata->npointsprob));
   SCIPfreeMemoryArray(scip, &(pricerdata->nraysprob));

   SCIPfreeMemoryArray(scip, &(pricerdata->farkascallsdist));
   SCIPfreeMemoryArray(scip, &(pricerdata->farkasfoundvars));
   SCIPfreeMemoryArray(scip, &(pricerdata->farkasnodetimedist));

   SCIPfreeMemoryArray(scip, &(pricerdata->redcostcallsdist));
   SCIPfreeMemoryArray(scip, &(pricerdata->redcostfoundvars));
   SCIPfreeMemoryArray(scip, &(pricerdata->redcostnodetimedist));

   SCIPfreeMemoryArray(scip, &(pricerdata->nodetimehist));
   SCIPfreeMemoryArray(scip, &(pricerdata->foundvarshist));
   SCIPfreeMemoryArray(scip, &(pricerdata->nodedegeneracy));
   pricerdata->nodetimehist = NULL;
   pricerdata->foundvarshist = NULL;
   pricerdata->nodedegeneracy = NULL;

   for( i = 0; i < pricerdata->npricedvars; i++ )
   {
      SCIP_CALL( SCIPdropVarEvent(scip, pricerdata->pricedvars[i], SCIP_EVENTTYPE_VARDELETED,
            pricerdata->eventhdlr, NULL, -1) );

      SCIP_CALL( SCIPreleaseVar(scip, &pricerdata->pricedvars[i]) );
   }
   SCIPfreeBlockMemoryArray(scip, &pricerdata->pricedvars, pricerdata->maxpricedvars);

   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->redcostclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->farkasclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->freeclock)) );
   SCIP_CALL( SCIPfreeClock(scip, &(pricerdata->transformclock)) );

   SCIP_CALL( solversExitsol(scip, pricerdata) );

   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostGcg)
{
   SCIP_RETCODE retcode;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   assert(pricerdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( pricerdata->redcostcalls == 0 )
   {
      if( pricerdata->farkascalls == 0 )
      {
         SCIP_CALL( SCIPconsMasterbranchAddRootCons(scip) );
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Starting reduced cost pricing...\n");
   }
   /* update number of reduced cost pricing rounds at the current node */
   if( SCIPgetNNodes(scip) == pricerdata->currnodenr )
   {
      pricerdata->nroundsredcost++;
   }
   else
   {
      pricerdata->currnodenr = SCIPgetNNodes(scip);
      pricerdata->nroundsredcost = 0;
   }

   /* if the number of reduced cost pricing rounds at the current node exceeds the limit (and we are not at the root), stop pricing;
    * we always stop pricing, if the maximum number of reduced cost rounds is set to 0
    */
   if( pricerdata->maxroundsredcost == 0 || (pricerdata->nroundsredcost >= pricerdata->maxroundsredcost && pricerdata->currnodenr != 1) )
   {
      SCIPdebugMessage("pricing aborted at node %lld\n", pricerdata->currnodenr);
      return SCIP_OKAY;
   }

   *result = SCIP_SUCCESS;

   /* perform pricing */
   SCIP_CALL( SCIPstartClock(scip, pricerdata->redcostclock) );
   retcode = performPricing(scip, pricer, GCG_PRICETYPE_REDCOST, result, lowerbound);
   SCIP_CALL( SCIPstopClock(scip, pricerdata->redcostclock) );

   return retcode;
}

/** farcas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasGcg)
{
   SCIP_RETCODE retcode;
   SCIP_PRICERDATA* pricerdata;
   SCIP_SOL** origsols;
   int norigsols;
   int i;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   assert(pricerdata != NULL);

   if( pricerdata->redcostcalls == 0 && pricerdata->farkascalls == 0 )
   {
      SCIP_CALL( SCIPconsMasterbranchAddRootCons(scip) );
   }

   /* get solutions from the original problem */
   origsols = SCIPgetSols(pricerdata->origprob);
   norigsols = SCIPgetNSols(pricerdata->origprob);
   assert(norigsols >= 0);

   /* Add already known solutions for the original problem to the master variable space
    * @todo This is just a workaround!
    */
   if( pricerdata->farkascalls == 0 )
   {
      for( i = 0; i < norigsols; ++i )
      {
         assert(origsols[i] != NULL);
         SCIP_CALL( GCGpricerTransOrigSolToMasterVars(scip, origsols[i]) );
      }
   }

   SCIP_CALL( SCIPstartClock(scip, pricerdata->farkasclock) );
   retcode = performPricing(scip, pricer, GCG_PRICETYPE_FARKAS, NULL, NULL);
   SCIP_CALL( SCIPstopClock(scip, pricerdata->farkasclock) );

   return retcode;
}

#define pricerCopyGcg NULL

/*
 * variable pricer specific interface methods
 */

/** creates the GCG variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerGcg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 origprob            /**< SCIP data structure of the original problem */
   )
{
   SCIP_PRICERDATA* pricerdata;

   SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );
   pricerdata->origprob = origprob;

   /* initialize solvers array */
   pricerdata->solvers = NULL;
   pricerdata->nsolvers = 0;
   pricerdata->nodetimehist = NULL;
   pricerdata->foundvarshist = NULL;
   pricerdata->nodedegeneracy = NULL;


   /* include variable pricer */
   SCIP_CALL( SCIPincludePricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerCopyGcg, pricerFreeGcg, pricerInitGcg, pricerExitGcg,
         pricerInitsolGcg, pricerExitsolGcg, pricerRedcostGcg, pricerFarkasGcg,
         pricerdata) );

   /* include event handler into master SCIP */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, eventFreeVardeleted, eventInitVardeleted, eventExitVardeleted,
         eventInitsolVardeleted, eventExitsolVardeleted, eventDeleteVardeleted, eventExecVardeleted,
         NULL) );

   pricerdata->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxsuccessfulmipsredcost",
         "maximal number of pricing mips leading to new variables solved solved in one redcost pricing round",
         &pricerdata->maxsuccessfulmipsredcost, FALSE, DEFAULT_MAXSUCCESSFULMIPSREDCOST, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundredcost",
         "maximal number of variables created in one redcost pricing round",
         &pricerdata->maxvarsroundredcost, FALSE, DEFAULT_MAXVARSROUNDREDCOST, 0, INT_MAX,
         NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundredcostroot",
         "maximal number of variables created in one redcost pricing round at the root node",
         &pricerdata->maxvarsroundredcostroot, FALSE, DEFAULT_MAXVARSROUNDREDCOSTROOT, 0, INT_MAX,
         NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxvarsroundfarkas",
         "maximal number of variables created in one farkas pricing round",
         &pricerdata->maxvarsroundfarkas, FALSE, DEFAULT_MAXVARSROUNDFARKAS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxroundsredcost",
         "maximal number of pricing rounds per node after the root node",
         &pricerdata->maxroundsredcost, FALSE, DEFAULT_MAXROUNDSREDCOST, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/maxsolsprob",
         "maximal number of variables added for each block in a pricinground",
         &pricerdata->maxsolsprob, FALSE, DEFAULT_MAXSOLSPROB, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/useheurpricing",
         "should pricing be performed heuristically before solving the MIPs to optimality?",
         &pricerdata->useheurpricing, TRUE, DEFAULT_USEHEURPRICING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/abortpricingint",
         "should pricing be aborted due to integral objective function?",
         &pricerdata->abortpricingint, TRUE, DEFAULT_ABORTPRICINGINT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/abortpricinggap",
         "should pricing be aborted due to small gap between dual bound and RMP objective?",
         &pricerdata->abortpricinggap, TRUE, DEFAULT_ABORTPRICINGGAP, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/successfulsubmipsrel",
         "part of the submips that are solved and lead to new variables before pricing round is aborted? (1.0 = solve all pricing MIPs)",
         &pricerdata->successfulmipsrel, FALSE, DEFAULT_SUCCESSFULMIPSREL, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/mipsrelredcostroot",
         "part of the submips that are solved before redcost pricing round is aborted at the root node, if variables have been found yed? (1.0 = solve all pricing MIPs)",
         &pricerdata->mipsrelredcostroot, FALSE, DEFAULT_MIPSRELREDCOSTROOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/mipsrelredcost",
         "part of the submips that are solved before redcost pricing round is aborted, if variables have been found yed? (1.0 = solve all pricing MIPs)",
         &pricerdata->mipsrelredcost, FALSE, DEFAULT_MIPSRELREDCOST, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(pricerdata->origprob, "pricing/masterpricer/mipsrelfarkas",
         "part of the submips that are solved before Farkas pricing round is aborted, if variables have been found yed? (1.0 = solve all pricing MIPs)",
         &pricerdata->mipsrelfarkas, FALSE, DEFAULT_MIPSRELFARKAS, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(pricerdata->origprob, "pricing/masterpricer/dispinfos",
         "should additional informations concerning the pricing process be displayed?",
         &pricerdata->dispinfos, FALSE, DEFAULT_DISPINFOS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(pricerdata->origprob, "pricing/masterpricer/sorting",
         "which sorting method should be used to sort the pricing problems (0 = order of pricing problems, 1 = according to dual solution of convexity constraint, 2 = according to reliability from previous round)",
         &pricerdata->sorting, FALSE, DEFAULT_SORTING, 0, 5, NULL, NULL) );

   return SCIP_OKAY;
}

/** returns the pointer to the scip instance representing the original problem */
SCIP* GCGpricerGetOrigprob(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->origprob;
}

/** returns the array of variables that were priced in during the solving process */
SCIP_VAR** GCGpricerGetPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->pricedvars;
}

/** returns the number of variables that were priced in during the solving process */
int GCGpricerGetNPricedvars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->npricedvars;
}


/** adds the given constraint and the given position to the hashmap of the pricer */
SCIP_RETCODE GCGpricerAddMasterconsToHashmap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint that should be added */
   int                   pos                 /**< the position of the constraint in the relaxator's masterconss array */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(pos >= 0);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   SCIP_CALL( SCIPhashmapInsert(pricerdata->mapcons2idx, cons, (void*)(size_t)pos) );
   assert((int)(size_t)SCIPhashmapGetImage(pricerdata->mapcons2idx, cons) == pos); /*lint !e507*/

   SCIPdebugMessage("Added cons %s (%p) to hashmap with index %d\n", SCIPconsGetName(cons), cons, pos);

   return SCIP_OKAY;
}

/** includes a solver into the pricer data */
SCIP_RETCODE GCGpricerIncludeSolver(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of solver */
   const char*           description,        /**< description of solver */
   int                   priority,           /**< priority of solver */
   GCG_DECL_SOLVERSOLVE  ((*solversolve)),   /**< solving method for solver */
   GCG_DECL_SOLVERSOLVEHEUR((*solveheur)),   /**< heuristic solving method for solver */
   GCG_DECL_SOLVERFREE   ((*solverfree)),    /**< free method of solver */
   GCG_DECL_SOLVERINIT   ((*solverinit)),    /**< init method of solver */
   GCG_DECL_SOLVEREXIT   ((*solverexit)),    /**< exit method of solver */
   GCG_DECL_SOLVERINITSOL((*solverinitsol)), /**< initsol method of solver */
   GCG_DECL_SOLVEREXITSOL((*solverexitsol)), /**< exitsol method of solver */
   GCG_SOLVERDATA*       solverdata          /**< solverdata data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   int pos;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   SCIP_CALL( ensureSizeSolvers(scip, pricerdata) );

   /* solvers array is sorted decreasingly wrt. the priority, find right position and shift solvers with smaller priority */
   pos = pricerdata->nsolvers;
   while( pos >= 1 && pricerdata->solvers[pos-1]->priority < priority )
   {
      pricerdata->solvers[pos] = pricerdata->solvers[pos-1];
      pos--;
   }
   SCIP_CALL( SCIPallocMemory(scip, &(pricerdata->solvers[pos])) ); /*lint !e866*/

   SCIP_ALLOC( BMSduplicateMemoryArray(&pricerdata->solvers[pos]->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&pricerdata->solvers[pos]->description, description, strlen(description)+1) );

   pricerdata->solvers[pos]->priority = priority;
   pricerdata->solvers[pos]->solversolve = solversolve;
   pricerdata->solvers[pos]->solversolveheur = solveheur;
   pricerdata->solvers[pos]->solverfree = solverfree;
   pricerdata->solvers[pos]->solverinit = solverinit;
   pricerdata->solvers[pos]->solverexit = solverexit;
   pricerdata->solvers[pos]->solverinitsol = solverinitsol;
   pricerdata->solvers[pos]->solverexitsol = solverexitsol;
   pricerdata->solvers[pos]->solverdata = solverdata;


   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->solvers[pos]->optfarkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->solvers[pos]->optredcostclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->solvers[pos]->heurfarkasclock)) );
   SCIP_CALL( SCIPcreateCPUClock(scip, &(pricerdata->solvers[pos]->heurredcostclock)) );

   pricerdata->solvers[pos]->optfarkascalls = 0;
   pricerdata->solvers[pos]->optredcostcalls = 0;
   pricerdata->solvers[pos]->heurfarkascalls = 0;
   pricerdata->solvers[pos]->heurredcostcalls = 0;

   pricerdata->nsolvers++;

   return SCIP_OKAY;
}

/** returns the solverdata of a solver */
GCG_SOLVERDATA* GCGpricerGetSolverdata(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_SOLVER*           solver              /**< pointer so solver */
   )
{
#ifndef NDEBUG
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
#endif

   assert(scip != NULL);
   assert(solver != NULL);

#ifndef NDEBUG
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);
#endif

   return solver->solverdata;
}

/** sets solver data of specific solver */
void GCGpricerSetSolverdata(
   SCIP*                 scip,               /**< SCIP data structure */
   GCG_SOLVER*           solver,             /**< pointer to solver  */
   GCG_SOLVERDATA*       solverdata          /**< solverdata data structure */
   )
{
#ifndef NDEBUG
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
#endif

   assert(scip != NULL);
   assert(solver != NULL);

#ifndef NDEBUG
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));
   assert(pricerdata->nsolvers > 0);
#endif

   solver->solverdata = solverdata;
}

/** writes out a list of all pricing problem solvers */
void GCGpricerPrintListOfSolvers(
   SCIP*                 scip                /**< SCIP data structure */
   )
{

   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   int nsolvers;
   int i;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   assert((pricerdata->solvers == NULL) == (pricerdata->nsolvers == 0));

   nsolvers = pricerdata->nsolvers;

   SCIPdialogMessage(scip, NULL, " solver               priority description\n --------------       -------- -----------\n");

   for( i = 0; i < nsolvers; ++i )
   {
      SCIPdialogMessage(scip, NULL,  " %-20s", pricerdata->solvers[i]->name);
      SCIPdialogMessage(scip, NULL,  " %8d", pricerdata->solvers[i]->priority);
      SCIPdialogMessage(scip, NULL,  " %s\n", pricerdata->solvers[i]->description);
   }
}

/** prints pricer statistics */
void GCGpricerPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   int i;
   double start;
   double end;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /**@todo add constraint statistics: how many constraints (instead of cuts) have been added? */
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Pricing Solver     : #HeurFarkas  #OptFarkas  #HeurRedcost #OptRedcost Time: HeurFarkas  OptFarkas  HeurRedcost OptRedcost\n");

   for( i = 0; i < pricerdata->nsolvers; ++i )
   {
      GCG_SOLVER* solver;
      solver = pricerdata->solvers[i];
      assert(solver != NULL);

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  %-17.17s:", solver->name);
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " %11d %11d   %11d %11d       %10.2f %10.2f   %10.2f %10.2f \n",
         solver->heurfarkascalls, solver->optfarkascalls,
         solver->heurredcostcalls, solver->optredcostcalls,
         SCIPgetClockTime(scip, solver->heurfarkasclock),
         SCIPgetClockTime(scip, solver->optfarkasclock),
         SCIPgetClockTime(scip, solver->heurredcostclock),
         SCIPgetClockTime(scip, solver->optredcostclock));
   }

   /* print of Pricing Statistics */

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Farkas pricing Statistic:\nno.\t#Calls\t\t#Vars\t\ttime(s)\n");

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d  \t %d \t\t %d \t\t %.2f \n", i, pricerdata->farkascallsdist[i],
         pricerdata->farkasfoundvars[i], pricerdata->farkasnodetimedist[i]);

   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Reduced Cost pricing Statistic:\nno.\t#Calls\t\t#Vars\t\ttime(s)\n");

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d  \t %d \t\t %d \t\t %.2f \n", i, pricerdata->redcostcallsdist[i],
         pricerdata->redcostfoundvars[i], pricerdata->redcostnodetimedist[i]);

   }

   /* print of Histogram Buckets !=0      */

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Histogram Time\n");
   for( i = 0; i < PRICER_STAT_ARRAYLEN_TIME; i++ )
   {
      start = (1.0 * i * PRICER_STAT_BUCKETSIZE_TIME)/1000.0;
      end = start + PRICER_STAT_BUCKETSIZE_TIME/1000.0;

      if( pricerdata->nodetimehist[i] != 0 )
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "From\t%.4f\t-\t%.4f\ts:\t\t%d \n", start, end, pricerdata->nodetimehist[i]);
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Histogram Found Vars\n");

   for( i = 0; i < PRICER_STAT_ARRAYLEN_VARS; i++ )
   {
      start = i * PRICER_STAT_BUCKETSIZE_VARS;
      end = start + PRICER_STAT_BUCKETSIZE_VARS;

      if( pricerdata->foundvarshist[i] != 0 )
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "From\t%.0f\t-\t%.0f\tvars:\t\t%d \n", start, end, pricerdata->foundvarshist[i]);
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Pricing Summary:\n");
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Calls                            : %d\n", pricerdata->calls);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Farkas Pricing Calls             : %d\n", pricerdata->farkascalls);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Farkas Pricing Time              : %f\n", SCIPgetClockTime(scip, pricerdata->farkasclock));
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Reduced Cost Pricing Calls       : %d\n", pricerdata->redcostcalls);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Reduced Cost Pricing Time        : %f\n", SCIPgetClockTime(scip, pricerdata->redcostclock));
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Solved subMIPs Heuristic Pricing : %d\n", pricerdata->solvedsubmipsheur);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Solved subMIPs Optimal Pricing   : %d\n", pricerdata->solvedsubmipsoptimal);
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Time for transformation          : %f\n", SCIPgetClockTime(scip, pricerdata->transformclock));
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Time for freeing subMIPs         : %f\n", SCIPgetClockTime(scip, pricerdata->freeclock));

}


/** transfers a primal solution of the original problem into the master variable space,
 *  i.e. creates one master variable for each block and adds the solution to the master problem  */
SCIP_RETCODE GCGpricerTransOrigSolToMasterVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             origsol             /**< the solution that should be transferred */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_SOL* mastersol;
   SCIP_VAR* newvar;
   SCIP* origprob;
   SCIP_Bool added;
   int prob;
   int i;

   SCIP_VAR** origvars;
   SCIP_Real* origsolvals;
   int norigvars;

   SCIP_VAR*** pricingvars;
   SCIP_Real** pricingvals;
   int* npricingvars;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = GCGpricerGetOrigprob(scip);
   assert(origprob != NULL);

   /* now compute coefficients of the master variables in the master constraint */
   origvars = SCIPgetVars(origprob);
   norigvars = SCIPgetNVars(origprob);

   /* allocate memory for storing variables and solution values from the solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &origsolvals, norigvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pricingvars, pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pricingvals, pricerdata->npricingprobs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &npricingvars, pricerdata->npricingprobs) );

   for( i = 0; i < pricerdata->npricingprobs; i++ )
   {
      npricingvars[i] = 0;
      if( pricerdata->pricingprobs[i] == NULL )
         continue;

      SCIP_CALL( SCIPallocBufferArray(scip, &pricingvars[i], SCIPgetNVars(pricerdata->pricingprobs[i])) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &pricingvals[i], SCIPgetNVars(pricerdata->pricingprobs[i])) ); /*lint !e866*/
   }

   /* get solution values */
   SCIP_CALL( SCIPgetSolVals(scip, origsol, norigvars, origvars, origsolvals) );
   SCIP_CALL( SCIPcreateSol(scip, &mastersol, SCIPgetSolHeur(origprob, origsol)) );

   /* store variables and solutions into arrays */
   for( i = 0; i < norigvars; i++ )
   {
      int blocknr;
      assert(GCGvarIsOriginal(origvars[i]));
      blocknr = GCGvarGetBlock(origvars[i]);
      assert(blocknr < 0 || GCGoriginalVarGetPricingVar(origvars[i]) != NULL);

      if( blocknr >= 0 )
      {
         prob = blocknr;
         if( pricerdata->pricingprobs[prob] == NULL )
            continue;

         if( !SCIPisZero(scip, origsolvals[i]) )
         {
            pricingvars[prob][npricingvars[prob]] = GCGoriginalVarGetPricingVar(origvars[i]);
            pricingvals[prob][npricingvars[prob]] = origsolvals[i];
            npricingvars[prob]++;
         }
      }
      else
      {
         assert((GCGoriginalVarGetNMastervars(origvars[i]) == 1) || (GCGvarIsLinking(origvars[i])));
         assert(GCGoriginalVarGetMastervars(origvars[i])[0] != NULL);
         SCIP_CALL( SCIPsetSolVal(scip, mastersol, GCGoriginalVarGetMastervars(origvars[i])[0], origsolvals[i]) );
      }
   }

   /* create variables in the master problem */
   for( prob = 0; prob < pricerdata->npricingprobs; prob++ )
   {
      if( pricerdata->pricingprobs[prob] == NULL )
      {
         continue;
      }
      SCIP_CALL( createNewMasterVar(scip, pricingvars[prob], pricingvals[prob], npricingvars[prob], FALSE, prob, TRUE, &added, &newvar) );
      assert(added);

      SCIP_CALL( SCIPsetSolVal(scip, mastersol, newvar, 1.0 * GCGrelaxGetNIdenticalBlocks(pricerdata->origprob, prob)) );
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPtrySolFree(scip, &mastersol, TRUE, TRUE, TRUE, TRUE, &added) );
#else
   SCIP_CALL( SCIPtrySolFree(scip, &mastersol, FALSE, TRUE, TRUE, TRUE, &added) );
#endif

   /* free memory for storing variables and solution values from the solution */

   for( i = pricerdata->npricingprobs - 1; i>= 0; i-- )
   {
      if( pricerdata->pricingprobs[i] == NULL )
      {
         continue;
      }

      SCIPfreeBufferArray(scip, &pricingvals[i]);
      SCIPfreeBufferArray(scip, &pricingvars[i]);
   }

   SCIPfreeBufferArray(scip, &npricingvars);
   SCIPfreeBufferArray(scip, &pricingvals);
   SCIPfreeBufferArray(scip, &pricingvars);
   SCIPfreeBufferArray(scip, &origsolvals);

   return SCIP_OKAY;
}


/** create initial master variables */
SCIP_RETCODE GCGpricerCreateInitialMastervars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   int i;
   SCIP* origprob;
   SCIP_VAR** vars;
   int nvars;
   int npricingprobs;
   int v;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   origprob = pricerdata->origprob;
   assert(origprob != NULL);

   npricingprobs = GCGrelaxGetNPricingprobs(origprob);
   assert(npricingprobs >= 0);

   /* for variables in the original problem that do not belong to any block,
    * create the corresponding variable in the master problem
    */
   vars = SCIPgetVars(origprob);
   nvars = SCIPgetNVars(origprob);
   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real* coefs;
      int blocknr;
      int ncoefs;
      SCIP_VAR* var;

      /* var = SCIPvarGetProbvar(vars[v]); */
      var = vars[v];
      blocknr = GCGvarGetBlock(var);
      coefs = GCGoriginalVarGetCoefs(var);
      ncoefs = GCGoriginalVarGetNCoefs(var);

      assert(GCGvarIsOriginal(var));
      if( blocknr < 0 )
      {
         SCIP_CONS** linkconss;
         SCIP_VAR* newvar;

         SCIP_CALL( GCGcreateInitialMasterVar(scip, var, &newvar) );
         SCIP_CALL( SCIPaddVar(scip, newvar) );

         SCIP_CALL( GCGoriginalVarAddMasterVar(scip, var, newvar, 1.0) );

         linkconss = GCGoriginalVarGetMasterconss(var);

         /* add variable in the master to the master constraints it belongs to */
         for( i = 0; i < ncoefs; i++ )
         {
            assert(!SCIPisZero(scip, coefs[i]));
            /*            SCIP_CALL( SCIPgetTransformedCons(scip, linkconss[i], &linkcons) );*/

            SCIP_CALL( SCIPaddCoefLinear(scip, linkconss[i], newvar, coefs[i]) );
         }

         /* we copied a linking variable into the master, add it to the linkcons */
         if( GCGvarIsLinking(var) )
         {
            SCIP_CONS** linkingconss;
            linkingconss = GCGlinkingVarGetLinkingConss(var);

            for( i = 0; i < npricingprobs; i++ )
            {
               if( linkingconss[i] != NULL )
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, linkingconss[i], newvar, 1.0) );
               }
            }
         }

         SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

      }
   }
   return SCIP_OKAY;
}
