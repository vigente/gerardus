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

/**@file   heur_gcgpscostdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the pseudo cost values
 * @author Tobias Achterberg
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_gcgpscostdiving.h"

#include "cons_origbranch.h"
#include "relax_gcg.h"


#define HEUR_NAME             "gcgpscostdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the pseudo cost values"
#define HEUR_DISPCHAR         'p'
#define HEUR_PRIORITY         -1002000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          2
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERPLUNGE
#define HEUR_USESSUBSCIP      FALSE



/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.05 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
//#define DEFAULT_MAXPRICEQUOT       0.05 /**< maximal fraction of pricing rounds compared to node pricing rounds */
//#define DEFAULT_MAXPRICEOFS          10 /**< additional number of allowed pricing rounds (-1: no limit) */
#define DEFAULT_MAXPRICEQUOT       0.00 /**< maximal fraction of pricing rounds compared to node pricing rounds */
#define DEFAULT_MAXPRICEOFS           0 /**< additional number of allowed pricing rounds (-1: no limit) */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */



/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_SOL*             rootsol;            /**< relaxation solution at the root node */
   SCIP_Bool             firstrun;           /**< is the heuristic running for the first time? */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int                   maxlpiterofs;       /**< additional number of allowed LP iterations */
   SCIP_Real             maxpricequot;       /**< maximal fraction of pricing rounds compared to node pricing rounds */
   int                   maxpriceofs;        /**< additional number of allowed pricing rounds (-1: no limit) */
   SCIP_Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Bool             backtrack;          /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Longint          nlpiterations;      /**< LP iterations used in this heuristic */
   int                   npricerounds;       /**< pricing rounds used in this heuristic */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */
};




/*
 * local methods
 */

/** get relaxation solution of root node (in original variables) */
static
SCIP_RETCODE getRootRelaxSol(
      SCIP*        scip,
      SCIP_SOL**   rootsol
      )
{
   SCIP* masterprob;
   SCIP_SOL* masterrootsol;
   SCIP_VAR** mastervars;
   int nmastervars;
   int i;

   masterprob = GCGrelaxGetMasterprob(scip);

   SCIP_CALL( SCIPcreateSol(masterprob, &masterrootsol, NULL) );

   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   for( i = 0; i < nmastervars; i++ )
      SCIP_CALL( SCIPsetSolVal(masterprob, masterrootsol, mastervars[i], SCIPvarGetRootSol(mastervars[i])) );

   SCIP_CALL( GCGrelaxTransformMastersolToOrigsol(scip, masterrootsol, rootsol) );

   SCIP_CALL( SCIPfreeSol(masterprob, &masterrootsol) );

   return SCIP_OKAY;
}

/** for a given variable, calculate its pseudo cost quotient */
static
void calcPscostQuot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             primsol,            /**< primal solution of variable */
   SCIP_Real             rootsolval,         /**< root relaxation solution of variable */
   SCIP_Real             frac,               /**< fractionality of variable */
   int                   rounddir,           /**< -1: round down, +1: round up, 0: select due to pseudo cost values */
   SCIP_Real*            pscostquot,         /**< pointer to store pseudo cost quotient */
   SCIP_Bool*            roundup             /**< pointer to store whether the variable should be rounded up */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   assert(pscostquot != NULL);
   assert(roundup != NULL);
   assert(SCIPisEQ(scip, frac, primsol - SCIPfeasFloor(scip, primsol)));

   /* bound fractions to not prefer variables that are nearly integral */
   frac = MAX(frac, 0.1);
   frac = MIN(frac, 0.9);

   /* get pseudo cost quotient */
   pscostdown = SCIPgetVarPseudocostVal(scip, var, 0.0-frac);
   pscostup = SCIPgetVarPseudocostVal(scip, var, 1.0-frac);
   SCIPdebugMessage("Pseudocosts of variable %s: %g down, %g up\n", SCIPvarGetName(var), pscostdown, pscostup);
   assert(pscostdown >= 0.0 && pscostup >= 0.0);

   /* choose rounding direction */
   if( rounddir == -1 )
      *roundup = FALSE;
   else if( rounddir == +1 )
      *roundup = TRUE;
   else if( primsol < rootsolval - 0.4 )
      *roundup = FALSE;
   else if( primsol > rootsolval + 0.4 )
      *roundup = TRUE;
   else if( frac < 0.3 )
      *roundup = FALSE;
   else if( frac > 0.7 )
      *roundup = TRUE;
   else if( pscostdown < pscostup )
      *roundup = FALSE;
   else
      *roundup = TRUE;

   /* calculate pseudo cost quotient */
   if( *roundup )
      *pscostquot = sqrt(frac) * (1.0+pscostdown) / (1.0+pscostup);
   else
      *pscostquot = sqrt(1.0-frac) * (1.0+pscostup) / (1.0+pscostdown);

   /* prefer decisions on binary variables */
   if( SCIPvarIsBinary(var) )
      (*pscostquot) *= 1000.0;
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#define heurCopyGcgpscostdiving NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeGcgpscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitGcgpscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution and root relaxation solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );
   heurdata->rootsol = NULL;

   /* initialize data */
   heurdata->nlpiterations = 0;
   heurdata->npricerounds = 0;
   heurdata->nsuccess = 0;
   heurdata->firstrun = TRUE;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitGcgpscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution and root relaxation solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );
   if( heurdata->rootsol != NULL )
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->rootsol) );

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolGcgpscostdiving NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolGcgpscostdiving NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecGcgpscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP* masterprob;
   SCIP_HEURDATA* heurdata;
   SCIP_CONS* probingcons;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_NODE* probingnode;
   SCIP_VAR* var;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real lpobj;
   SCIP_Real objval;
   SCIP_Real oldobjval;
   SCIP_Real primsol;
   SCIP_Real rootsolval;
   SCIP_Real frac;
   SCIP_Real bestfrac;
   SCIP_Real pscostquot;
   SCIP_Real bestpscostquot;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   SCIP_Bool bestcandroundup;
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   SCIP_Bool roundup;
   SCIP_Bool lperror;
   SCIP_Bool lpsolved;
   SCIP_Bool cutoff;
   SCIP_Bool feasible;
   SCIP_Bool backtracked;
   SCIP_Longint ncalls;
   SCIP_Longint nsolsfound;
   SCIP_Longint nlpiterations;
   SCIP_Longint maxnlpiterations;
   int maxpricerounds;
   int npricerounds;
   int totalpricerounds;
   int nlpcands;
   int startnlpcands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int bestcand;
   int c;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   /* get master problem */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetStage(masterprob) > SCIP_STAGE_SOLVING ||
         !SCIPhasCurrentNodeLP(masterprob) || SCIPgetLPSolstat(masterprob) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(masterprob) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(masterprob) == SCIPgetNNodes(masterprob) && SCIPgetDepth(masterprob) > 0 )
      return SCIP_OKAY;

   /** @todo for some reason, the heuristic is sometimes called with an invalid relaxation solution;
    *       in that case, don't execute it */
   if( !SCIPisRelaxSolValid(scip) )
   {
      SCIPdebugMessage("not executing GCG pscostdiving: invalid relaxation solution (should not happen!)\n");
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < heurdata->minreldepth*maxdepth || depth > heurdata->maxreldepth*maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNNodeLPIterations(scip) + SCIPgetNNodeLPIterations(masterprob);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = 10*SCIPheurGetNBestSolsFound(heur) + heurdata->nsuccess;
   maxnlpiterations = (SCIP_Longint)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * nlpiterations);
   maxnlpiterations += heurdata->maxlpiterofs;

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + MINLPITER);

   /** @todo limit number of pricing rounds, play with parameters */
   if( heurdata->maxpriceofs > -1 )
   {
      npricerounds = SCIPgetNPriceRounds(masterprob);
      SCIPdebugMessage("GCG pscostdiving - pricing rounds at this node: %d\n", npricerounds);
      maxpricerounds = (int)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxpricequot * npricerounds);
      maxpricerounds += heurdata->maxpriceofs;
   }
   else
      maxpricerounds = -1;

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      if( heurdata->maxdiveubquotnosol > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquotnosol * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquotnosol > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquotnosol * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   else
   {
      if( heurdata->maxdiveubquot > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquot * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquot > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   searchbound = MIN(searchubbound, searchavgbound);
   if( SCIPisObjIntegral(scip) )
      searchbound = SCIPceil(scip, searchbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;


   *result = SCIP_DIDNOTFIND;

   /* if the heuristic is running for the first time, the root relaxation solution needs to be stored */
   if( heurdata->firstrun )
   {
      SCIP_CALL( getRootRelaxSol(scip, &heurdata->rootsol) );
      heurdata->firstrun = FALSE;
   }

   /* start diving */
   SCIP_CALL( SCIPstartProbing(scip) );
   SCIP_CALL( GCGrelaxStartProbing(scip, heur) );

   /* get LP objective value, and fractional variables, that should be integral */
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetRelaxSolObj(scip);
   lpobj = objval;
   SCIP_CALL( SCIPgetExternBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL, NULL, NULL) );

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing GCG pscostdiving heuristic: depth=%d, %d fractionals, dualbound=%g, avgbound=%g, cutoffbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPgetAvgDualbound(scip),
      SCIPretransformObj(scip, SCIPgetCutoffbound(scip)), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   lperror = FALSE;
   cutoff = FALSE;
   divedepth = 0;
   npricerounds = 0;
   totalpricerounds = 0;
//   bestcandmayrounddown = FALSE;
//   bestcandmayroundup = FALSE;
   startnlpcands = nlpcands;
   while( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations && objval < searchbound))
      && !SCIPisStopped(scip) )
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      divedepth++;

      /* choose variable fixing:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, round variable with largest rel. difference of pseudo cost values in corresponding
       *     direction
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - round variable in the objective value direction
       */
      bestcand = -1;
      bestpscostquot = -1.0;
      bestcandmayrounddown = TRUE;
      bestcandmayroundup = TRUE;
      bestcandroundup = FALSE;
      bestfrac = SCIP_INVALID;
      for( c = 0; c < nlpcands; ++c )
      {
         var = lpcands[c];
         mayrounddown = SCIPvarMayRoundDown(var);
         mayroundup = SCIPvarMayRoundUp(var);
         primsol = lpcandssol[c];
         rootsolval = SCIPgetSolVal(scip, heurdata->rootsol, var);
         frac = lpcandsfrac[c];
         if( mayrounddown || mayroundup )
         {
            /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
            if( bestcandmayrounddown || bestcandmayroundup )
            {
               /* choose rounding direction:
                * - if variable may be rounded in both directions, round corresponding to the pseudo cost values
                * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
                *   the current fractional solution
                */
               roundup = FALSE;
               if( mayrounddown && mayroundup )
                  calcPscostQuot(scip, var, primsol, rootsolval, frac, 0, &pscostquot, &roundup);
               else if( mayrounddown )
                  calcPscostQuot(scip, var, primsol, rootsolval, frac, +1, &pscostquot, &roundup);
               else
                  calcPscostQuot(scip, var, primsol, rootsolval, frac, -1, &pscostquot, &roundup);

               /* check, if candidate is new best candidate */
               if( pscostquot > bestpscostquot )
               {
                  bestcand = c;
                  bestpscostquot = pscostquot;
                  bestcandmayrounddown = mayrounddown;
                  bestcandmayroundup = mayroundup;
                  bestcandroundup = roundup;
                  bestfrac = frac;
               }
            }
         }
         else
         {
            /* the candidate may not be rounded: calculate pseudo cost quotient and preferred direction */
            calcPscostQuot(scip, var, primsol, rootsolval, frac, 0, &pscostquot, &roundup);

            /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
            if( bestcandmayrounddown || bestcandmayroundup || pscostquot > bestpscostquot )
            {
               bestcand = c;
               bestpscostquot = pscostquot;
               bestcandmayrounddown = FALSE;
               bestcandmayroundup = FALSE;
               bestcandroundup = roundup;
               bestfrac = frac;
            }
         }
      }
      assert(bestcand != -1);

      /* if all candidates are roundable, try to round the solution */
      if( bestcandmayrounddown || bestcandmayroundup )
      {
         SCIP_Bool success;

         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIPlinkRelaxSol(scip, heurdata->sol) );
         SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

         if( success )
         {
            SCIPdebugMessage("GCG pscostdiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

            /* a rounded solution will only be accepted if its objective value is below the search bound */
            if( SCIPgetSolOrigObj(scip, heurdata->sol) <= searchbound )
            {
               /* try to add solution to SCIP */
               //               SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, &success) );
#ifdef SCIP_DEBUG
               SCIP_CALL( SCIPtrySol(scip, heurdata->sol, TRUE, TRUE, TRUE, TRUE, &success) );
#else
               SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, TRUE, TRUE, TRUE, &success) );
#endif

               /* check, if solution was feasible and good enough */
               if( success )
               {
                  SCIPdebugMessage(" -> solution was feasible and good enough\n");
                  *result = SCIP_FOUNDSOL;
               }
            }
         }
      }

      var = lpcands[bestcand];

      /* if the variable is already fixed, abort diving due to numerical troubles */
      if( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
      {
         SCIPdebugMessage("numerical troubles: selected variable <%s> already fixed to [%g,%g] (solval: %.9f)\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), lpcandssol[bestcand]);
         break;
      }

      backtracked = FALSE;
      do
      {
         /* if the variable is already fixed, numerical troubles may have occured or
          * variable was fixed by propagation while backtracking => Abort diving!
          */
         if( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
         {
            SCIPdebugMessage("Selected variable <%s> already fixed to [%g,%g] (solval: %.9f), diving aborted \n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), lpcandssol[bestcand]);
            cutoff = TRUE;
            break;
         }

         probingnode = SCIPgetCurrentNode(scip);

         /* apply rounding of best candidate */
         if( bestcandroundup == !backtracked )
         {
            /* round variable up */
            SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", pricerounds %d/%d: var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, totalpricerounds, maxpricerounds,
               SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               lpcandssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               SCIPfeasCeil(scip, lpcandssol[bestcand]), SCIPvarGetUbLocal(var));

            SCIP_CALL( GCGcreateConsOrigbranch(scip, &probingcons, "probingcons", probingnode,
                  GCGconsOrigbranchGetActiveCons(scip), NULL, NULL) );
            SCIP_CALL( SCIPaddConsNode(scip, probingnode, probingcons, NULL) );
            SCIP_CALL( SCIPreleaseCons(scip, &probingcons) );
            SCIP_CALL( SCIPchgVarLbProbing(scip, var, SCIPfeasCeil(scip, lpcandssol[bestcand])) );
         }
         else
         {
            /* round variable down */
            SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", pricerounds %d/%d: var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, totalpricerounds, maxpricerounds,
               SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               lpcandssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               SCIPvarGetLbLocal(var), SCIPfeasFloor(scip, lpcandssol[bestcand]));

            SCIP_CALL( GCGcreateConsOrigbranch(scip, &probingcons, "probingcons", probingnode,
                  GCGconsOrigbranchGetActiveCons(scip), NULL, NULL) );
            SCIP_CALL( SCIPaddConsNode(scip, probingnode, probingcons, NULL) );
            SCIP_CALL( SCIPreleaseCons(scip, &probingcons) );
            SCIP_CALL( SCIPchgVarUbProbing(scip, var, SCIPfeasFloor(scip, lpcandssol[bestcand])) );
         }

         /* apply domain propagation */
//         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
         SCIP_CALL( SCIPpropagateProbing(scip, -1, &cutoff, NULL) );
         if( !cutoff )
         {
            /* resolve the diving LP */
            /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
             * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
             */
#ifdef NDEBUG
            SCIP_RETCODE retstat;
            if( maxpricerounds == 0 )
               retstat = GCGrelaxPerformProbing(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &nlpiterations, &lpobj, &lpsolved, &lperror, &cutoff, &feasible);
            else
            {
               retstat = GCGrelaxPerformProbingWithPricing(scip, maxpricerounds == -1 ? -1 : maxpricerounds - totalpricerounds,
                     &nlpiterations, &npricerounds, &lpobj, &lpsolved, &lperror, &cutoff, &feasible);
               npricerounds = 0;
            }
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving LP in GCG coefdiving heuristic; LP solve terminated with code <%d>\n",retstat);
            }
#else
            if( maxpricerounds == 0 )
               SCIP_CALL( GCGrelaxPerformProbing(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &nlpiterations, &lpobj, &lpsolved, &lperror, &cutoff, &feasible) );
            else
            {
               SCIP_CALL( GCGrelaxPerformProbingWithPricing(scip, maxpricerounds == -1 ? -1 : maxpricerounds - totalpricerounds,
                     &nlpiterations, &npricerounds, &lpobj, &lpsolved, &lperror, &cutoff, &feasible) );
               npricerounds = 0;
            }
#endif

            if( lperror || !lpsolved )
               break;

            /* update iteration count */
            heurdata->nlpiterations += nlpiterations;
            heurdata->npricerounds += npricerounds;
            totalpricerounds += npricerounds;

            /* get LP solution status, objective value, and fractional variables, that should be integral */
            lpsolstat = SCIPgetLPSolstat(scip);
//            cutoff = (lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
         }

         /* perform backtracking if a cutoff was detected */
         if( cutoff && !backtracked && heurdata->backtrack )
         {
            SCIPdebugMessage("  *** cutoff detected at level %d - backtracking\n", SCIPgetProbingDepth(scip));
            SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
            SCIP_CALL( SCIPbacktrackProbing(masterprob, SCIPgetProbingDepth(scip)) );
            SCIP_CALL( SCIPnewProbingNode(scip) );
            backtracked = TRUE;
         }
         else
            backtracked = FALSE;
      }
      while( backtracked );

      if( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new objective value */
         oldobjval = objval;
         objval = lpobj;

         /* update pseudo cost values */
         if( SCIPisGT(scip, objval, oldobjval) )
         {
            if( bestcandroundup )
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, var, 1.0-bestfrac,
                     objval - oldobjval, 1.0) );
            }
            else
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, var, 0.0-bestfrac,
                     objval - oldobjval, 1.0) );
            }
         }

         /* get new fractional variables */
         SCIP_CALL( SCIPgetExternBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL, NULL, NULL) );
      }
      SCIPdebugMessage("   -> lpsolstat=%d, objval=%g/%g, nfrac=%d, lpiterations=%"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT"\n",
         lpsolstat, objval, searchbound, nlpcands, heurdata->nlpiterations, maxnlpiterations);
   }

   /* check if a solution has been found */
   /** @todo maybe this is unneccessary since solutions are also added in GCGrelaxUpdateCurrentSol() */
   if( nlpcands == 0 && !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && divedepth > 0 )
   {
      SCIP_Bool success;

      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkRelaxSol(scip, heurdata->sol) );
      SCIPdebugMessage("GCG pscostdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      //      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, &success) );
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, TRUE, TRUE, TRUE, TRUE, &success) );
#else
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, TRUE, TRUE, TRUE, &success) );
#endif

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   SCIP_CALL( SCIPendProbing(scip) );
   SCIP_CALL( GCGrelaxEndProbing(scip) );

   if( *result == SCIP_FOUNDSOL )
      heurdata->nsuccess++;

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") finished GCG pscostdiving heuristic: %d fractionals, dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", pricerounds %d/%d, objval=%g/%g, lpsolstat=%d, cutoff=%u\n",
         SCIPgetNNodes(scip), nlpcands, divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations, totalpricerounds, maxpricerounds,
         SCIPretransformObj(scip, objval), SCIPretransformObj(scip, searchbound), lpsolstat, cutoff);

   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the GCG pscostdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurGcgpscostdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyGcgpscostdiving, heurFreeGcgpscostdiving, heurInitGcgpscostdiving, heurExitGcgpscostdiving,
         heurInitsolGcgpscostdiving, heurExitsolGcgpscostdiving, heurExecGcgpscostdiving,
         heurdata) );

   /* gcgpscostdiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/gcgpscostdiving/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/gcgpscostdiving/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/gcgpscostdiving/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/gcgpscostdiving/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/gcgpscostdiving/maxpricequot",
         "maximal fraction of pricing rounds compared to node pricing rounds",
         &heurdata->maxpricequot, FALSE, DEFAULT_MAXPRICEQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/gcgpscostdiving/maxpriceofs",
         "additional number of allowed pricing rounds (-1: no limit)",
         &heurdata->maxpriceofs, FALSE, DEFAULT_MAXPRICEOFS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/gcgpscostdiving/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveubquot, TRUE, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/gcgpscostdiving/maxdiveavgquot",
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveavgquot, TRUE, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/gcgpscostdiving/maxdiveubquotnosol",
         "maximal UBQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveubquotnosol, TRUE, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/gcgpscostdiving/maxdiveavgquotnosol",
         "maximal AVGQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveavgquotnosol, TRUE, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/gcgpscostdiving/backtrack",
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->backtrack, FALSE, DEFAULT_BACKTRACK, NULL, NULL) );

   return SCIP_OKAY;
}

