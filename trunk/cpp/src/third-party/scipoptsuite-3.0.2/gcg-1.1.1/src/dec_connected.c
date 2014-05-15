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

/**@file   dec_connected.c
 * @ingroup DETECTORS
 * @brief  detector for classical and blockdiagonal problems
 * @author Martin Bergner
 * @todo allow decompositions with only one pricing problem by just removing generalized covering and
 *       partitioning constraints
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "dec_connected.h"
#include "cons_decomp.h"
#include "scip_misc.h"
#include "pub_decomp.h"

/* constraint handler properties */
#define DEC_DETECTORNAME         "connected"    /**< name of detector */
#define DEC_DESC                 "Detector for classical and block diagonal problems" /**< description of detector*/
#define DEC_PRIORITY             0              /**< priority of the constraint handler for separation */
#define DEC_DECCHAR              'C'            /**< display character of detector */

#define DEC_ENABLED              TRUE           /**< should the detection be enabled */
#define DEFAULT_SETPPCINMASTER   TRUE           /**< should the extended structure be detected */
/*
 * Data structures
 */

/** constraint handler data */
struct DEC_DetectorData
{
   SCIP_HASHMAP* constoblock;                /**< hashmap mapping constraints to their associated block */
   SCIP_HASHMAP* vartoblock;                 /**< hashmap mapping variables to their associated block */
   SCIP_Bool blockdiagonal;                  /**< flag to indicate whether the problem is block diagonal */

   SCIP_CLOCK* clock;                        /**< clock to measure detection time */
   int nblocks;                              /**< number of blocks found */

   SCIP_Bool* consismaster;                  /**< boolean array to indicate constraints which should be in the master */
   SCIP_Bool setppcinmaster;                 /**< flag to indicate whether setppc constraints should always be in the master */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/* returns true if the constraint should be a master constraint and false otherwise */
static
SCIP_Bool isConsMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int i;
   int nvars;
   SCIP_Bool relevant = TRUE;
   assert(scip != NULL);
   assert(cons != NULL);

   SCIPdebugMessage("cons %s is ", SCIPconsGetName(cons));

   if( SCIPconsGetType(cons) == setcovering || SCIPconsGetType(cons) == setpartitioning || SCIPconsGetType(cons) == logicor )
   {
      SCIPdebugPrintf("setcov, part or logicor.\n");
      return TRUE;
   }
   nvars = SCIPgetNVarsXXX(scip, cons);
   vars = NULL;
   vals = NULL;
   if( nvars > 0 )
   {
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vals, nvars) );
      SCIP_CALL_ABORT( SCIPgetVarsXXX(scip, cons, vars, nvars) );
      SCIP_CALL_ABORT( SCIPgetValsXXX(scip, cons, vals, nvars) );
   }

   /* check vars and vals for integrality */
   for( i = 0; i < nvars && relevant; ++i )
   {
      assert(vars != NULL);
      assert(vals != NULL);

      if( !SCIPvarIsIntegral(vars[i]) && !SCIPvarIsBinary(vars[i]) )
      {
         SCIPdebugPrintf("(%s is not integral) ", SCIPvarGetName(vars[i]) );
         relevant = FALSE;
      }
      if( !SCIPisEQ(scip, vals[i], 1.0) )
      {
         SCIPdebugPrintf("(coeff for var %s is %.2f != 1.0) ", SCIPvarGetName(vars[i]), vals[i] );
         relevant = FALSE;
      }
   }

   /* free temporary data  */
   SCIPfreeBufferArrayNull(scip, &vals);
   SCIPfreeBufferArrayNull(scip, &vars);

   SCIPdebugPrintf("%s master\n", relevant ? "in" : "not in");
   return relevant;
}


/** looks for connected components in the constraints in detectordata */
static
SCIP_RETCODE findConnectedComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   SCIP_Bool             findextended,       /**< whether the classical structure should be detected */
   SCIP_RESULT*          result              /**< result pointer to indicate success oder failure */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nconss;
   SCIP_VAR** curvars;
   int ncurvars;
   SCIP_CONS* cons;
   SCIP_CONS** conss;

   int i;
   int j;
   int k;
   int tempblock;
   SCIP_Bool masterisempty;
   int* blockrepresentative;
   int nextblock;
   int *vartoblock;
   SCIP_HASHMAP *constoblock;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(result != NULL);

   /* initialize data structures */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   nextblock = 1; /* start at 1 in order to see whether the hashmap has a key*/

   SCIP_CALL( SCIPallocBufferArray(scip, &vartoblock, nvars+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockrepresentative, nconss+1) );
   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss+1) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip), nconss) );

   for( i = 0; i < nvars; ++i )
   {
      vartoblock[i] = -1;
   }

   for( i = 0; i < nconss+1; ++i )
   {
      blockrepresentative[i] = -1;
   }

   blockrepresentative[0] = 0;
   blockrepresentative[1] = 1;
   assert(nconss >= 1);

   masterisempty = findextended;

   /* in a first preprocessing step, indicate which constraints should go in the master */
   if( findextended )
   {
      for( i = 0; i < nconss; ++i )
      {
         detectordata->consismaster[i] = isConsMaster(scip, conss[i]);

         /* mark if there is one constraint in the master */
         masterisempty = masterisempty && !detectordata->consismaster[i];

         /* we look for an extended structure if there is a constraint not in the master! */
         findextended = findextended || !detectordata->consismaster[i];
      }
   }

   if( !findextended )
   {
      for( i = 0; i < nconss; ++i )
         detectordata->consismaster[i] = FALSE;
   }

   /* go through the all constraints */
   for( i = 0; i < nconss; ++i )
   {
      int consblock;

      cons = conss[i];
      assert(cons != NULL);
      if( GCGisConsGCGCons(cons) )
         continue;

      if( detectordata->consismaster[i] )
         continue;

      /* get variables of constraint */
      ncurvars = SCIPgetNVarsXXX(scip, cons);
      curvars = NULL;
      if( ncurvars > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, cons, curvars, ncurvars) );
      }
      assert(ncurvars >= 0);
      assert(ncurvars <= nvars);
      assert(curvars != NULL || ncurvars == 0);

      assert(SCIPhashmapGetImage(constoblock, cons) == NULL);

      /* if there are no variables, put it in the first block, otherwise put it in the next block */
      if( ncurvars == 0 )
         consblock = -1;
      else
         consblock = nextblock;

      /* go through all variables */
      for( j = 0; j < ncurvars; ++j )
      {
         SCIP_VAR* probvar;
         int varindex;
         int varblock;

         assert(curvars != NULL);
         probvar = SCIPvarGetProbvar(curvars[j]);
         assert(probvar != NULL);

         varindex = SCIPvarGetProbindex(probvar);
         assert(varindex >= 0);
         assert(varindex < nvars);

         /** @todo what about deleted variables? */
         /* get block of variable */
         varblock = vartoblock[varindex];

         SCIPdebugMessage("\tVar %s (%d): ", SCIPvarGetName(probvar), varblock);
         /* if variable is assigned to a block, assign constraint to that block */
         if( varblock > -1 && varblock != consblock )
         {
            consblock = MIN(consblock, blockrepresentative[varblock]);
            SCIPdebugPrintf("still in block %d.\n",  varblock);
         }
         else if( varblock == -1 )
         {
            /* if variable is free, assign it to the new block for this constraint */
            varblock = consblock;
            assert(varblock > 0);
            assert(varblock <= nextblock);
            vartoblock[varindex] = varblock;
            SCIPdebugPrintf("new in block %d.\n",  varblock);
         }
         else
         {
            assert((varblock > 0) && (consblock == varblock));
            SCIPdebugPrintf("no change.\n");
         }

         SCIPdebugPrintf("VARINDEX: %d (%d)\n", varindex, vartoblock[varindex]);
      }


      /* if the constraint belongs to a new block, mark it as such */
      if( consblock == nextblock )
      {
         assert(consblock > 0);
         blockrepresentative[consblock] = consblock;
         assert(blockrepresentative[consblock] > 0);
         assert(blockrepresentative[consblock] <= nextblock);
         ++nextblock;
      }

      SCIPdebugMessage("Cons %s will be in block %d (next %d)\n", SCIPconsGetName(cons), consblock, nextblock);
      for( k = 0; k < ncurvars; ++k )
      {
         int curvarindex;
         SCIP_VAR* curprobvar;
         int oldblock;
         assert(curvars != NULL);

         curprobvar = SCIPvarGetProbvar(curvars[k]);
         curvarindex = SCIPvarGetProbindex(curprobvar);
         oldblock = vartoblock[curvarindex];
         assert((oldblock > 0) && (oldblock <= nextblock));
         SCIPdebugMessage("\tVar %s ", SCIPvarGetName(curprobvar));
         if( oldblock != consblock )
         {
            SCIPdebugPrintf("reset from %d to block %d.\n", oldblock, consblock);
            vartoblock[curvarindex] = consblock;
            SCIPdebugPrintf("VARINDEX: %d (%d)\n", curvarindex, consblock);

            if( (blockrepresentative[oldblock] != -1) && (blockrepresentative[oldblock] > blockrepresentative[consblock]) )
            {
               int oldrepr;
               oldrepr = blockrepresentative[oldblock];
               SCIPdebugMessage("\t\tBlock representative from block %d changed from %d to %d.\n", oldblock, blockrepresentative[oldblock], consblock );
               assert(consblock > 0);
               blockrepresentative[oldblock] = consblock;
               if( (oldrepr != consblock) && (oldrepr != oldblock) )
               {
                  blockrepresentative[oldrepr] = consblock;
                  SCIPdebugMessage("\t\tBlock representative from block %d changed from %d to %d.\n", oldrepr, blockrepresentative[oldrepr], consblock );
               }
            }
         }
         else
         {
            SCIPdebugPrintf("will not be changed from %d to %d.\n", oldblock, consblock);
         }
      }

      SCIPfreeBufferArrayNull(scip, &curvars);
      assert(consblock >= 1 || consblock == -1);
      assert(consblock <= nextblock);

      /* store the constraint block */
      if( consblock != -1 )
      {
         SCIPdebugMessage("cons %s in block %d\n", SCIPconsGetName(cons), consblock);
         SCIP_CALL( SCIPhashmapInsert(constoblock, cons, (void*)(size_t)consblock) );
      }
      else
      {
         SCIPdebugMessage("ignoring %s\n", SCIPconsGetName(cons));
      }
   }

   tempblock = 1;

   SCIPdebugPrintf("Blocks: ");

   /* postprocess blockrepresentatives */
   for( i = 1; i < nextblock; ++i )
   {
      /* forward replace the representatives */
      assert(blockrepresentative[i] >= 0);
      assert(blockrepresentative[i] < nextblock);
      if( blockrepresentative[i] != i )
         blockrepresentative[i] = blockrepresentative[blockrepresentative[i]];
      else
      {
         blockrepresentative[i] = tempblock;
         ++tempblock;
      }
      /* It is crucial that this condition holds */
      assert(blockrepresentative[i] <= i);
      SCIPdebugPrintf("%d ", blockrepresentative[i]);
   }
   SCIPdebugPrintf("\n");

   detectordata->nblocks = tempblock-1;

   /* convert temporary data to detectordata */
   for( i = 0; i < nconss; ++i )
   {
      int consblock;

      cons = conss[i];
      if( GCGisConsGCGCons(cons) )
         continue;

      if( detectordata->consismaster[i] )
      {
         SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, cons, (void*) (size_t) (detectordata->nblocks+1)) );
         continue;
      }

      if( !SCIPhashmapExists(constoblock, cons) )
         continue;

      consblock = (int)(size_t) SCIPhashmapGetImage(constoblock, cons); /*lint !e507*/
      assert(consblock > 0);
      consblock = blockrepresentative[consblock];
      assert(consblock < tempblock);
      SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, cons, (void*)(size_t)consblock) );
      SCIPdebugMessage("%d %s\n", consblock, SCIPconsGetName(cons));
   }

   SCIP_CALL( SCIPhashmapCreate(&detectordata->vartoblock, SCIPblkmem(scip), nvars+1) );

   for( i = 0; i < nvars; ++i )
   {
      int varindex;
      int varblock;
      varindex = SCIPvarGetProbindex(SCIPvarGetProbvar(vars[i]));
      assert(varindex >= 0);
      assert(varindex < nvars);

      assert(vartoblock[varindex] < nextblock);
      if( vartoblock[varindex] < 0 )
      {
         SCIP_CALL( SCIPhashmapInsert(detectordata->vartoblock, SCIPvarGetProbvar(vars[i]),
                        (void*)(size_t)(detectordata->nblocks+1)) );
         continue;
      }


      varblock = blockrepresentative[vartoblock[varindex]];
      assert(varblock == -1 || varblock > 0);
      if( varblock > 0 )
      {
         assert(varblock < tempblock);
         SCIPdebugMessage("Var %s in block %d\n", SCIPvarGetName(SCIPvarGetProbvar(vars[i])), varblock-1);
         SCIP_CALL( SCIPhashmapInsert(detectordata->vartoblock, SCIPvarGetProbvar(vars[i]),
               (void*)(size_t)(varblock)) );
      }
   }

   /* free method data */
   SCIPfreeBufferArray(scip, &vartoblock);
   SCIPfreeBufferArray(scip, &blockrepresentative);
   SCIPhashmapFree(&constoblock);

   if( detectordata->nblocks > 1 )
      *result = SCIP_SUCCESS;
   else if( detectordata->nblocks == 1 && findextended && !masterisempty )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/* copy detectordata data to decdecomp */
static
SCIP_RETCODE copyToDecdecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< constraint handler data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp data structure */
   )
{
   SCIP_CONS** conss;
   int nconss;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Bool valid;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomp != NULL);

   assert(DECdecompGetType(decdecomp) == DEC_DECTYPE_UNKNOWN);

   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);


   SCIP_CALL( DECfillOutDecdecompFromHashmaps(scip, decdecomp, detectordata->vartoblock, detectordata->constoblock, detectordata->nblocks, vars, nvars, conss, nconss, &valid, FALSE) );
   assert(valid);

   detectordata->vartoblock = NULL;
   detectordata->constoblock = NULL;

   return SCIP_OKAY;
}

/** destructor of detector to free detector data (called when SCIP is exiting) */
static
DEC_DECL_EXITDETECTOR(exitConnected)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   if( detectordata->clock != NULL )
      SCIP_CALL( SCIPfreeClock(scip, &detectordata->clock) );

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** detection initialization function of detector (called before solving is about to begin) */
static
DEC_DECL_INITDETECTOR(initConnected)
{  /*lint --e{715}*/

   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   detectordata->clock = NULL;
   detectordata->constoblock = NULL;
   detectordata->vartoblock = NULL;
   detectordata->blockdiagonal = FALSE;

   detectordata->nblocks = 0;
   detectordata->consismaster = NULL;

   SCIP_CALL( SCIPcreateClock(scip, &detectordata->clock) );

   return SCIP_OKAY;
}

/** detection function of detector */
static
DEC_DECL_DETECTSTRUCTURE(detectConnected)
{
   int runs;
   int i;
   int nconss;
   SCIP_Bool detectextended;
   *result = SCIP_DIDNOTFIND;
   nconss = SCIPgetNConss(scip);

   runs = detectordata->setppcinmaster ? 2:1;
   detectextended = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &detectordata->consismaster, nconss) );
   *ndecdecomps = 0;

   for( i = 0; i < runs && *result != SCIP_SUCCESS; ++i )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting %s structure:", detectextended ? "set partitioning master":"purely block diagonal" );

      SCIP_CALL( SCIPstartClock(scip, detectordata->clock) );

      SCIP_CALL( findConnectedComponents(scip, detectordata, detectextended, result) );

      SCIP_CALL( SCIPstopClock(scip, detectordata->clock) );

      SCIPdebugMessage("Detection took %fs.\n", SCIPgetClockTime(scip, detectordata->clock));
      if( *result == SCIP_SUCCESS )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found %d blocks.\n", detectordata->nblocks);
         SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, 1) ); /*lint !e506*/
         SCIP_CALL( DECdecompCreate(scip, &((*decdecomps)[0])) );
         SCIP_CALL( copyToDecdecomp(scip, detectordata, (*decdecomps)[0]) );
         detectordata->blockdiagonal = DECdecompGetType((*decdecomps)[0]) == DEC_DECTYPE_DIAGONAL;
         *ndecdecomps = 1;
      }
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " not found.\n");
         SCIPhashmapFree(&detectordata->constoblock);
         SCIPhashmapFree(&detectordata->vartoblock);
      }
      if( detectordata->setppcinmaster == TRUE && *result != SCIP_SUCCESS )
      {
         detectextended = TRUE;
      }
   }

   SCIPfreeBufferArray(scip, &detectordata->consismaster);

   return SCIP_OKAY;
}


/*
 * detector specific interface methods
 */

/** creates the handler for connected constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionConnected(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /* create connected constraint handler data */
   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);

   detectordata->clock = NULL;
   detectordata->constoblock = NULL;
   detectordata->vartoblock = NULL;
   detectordata->blockdiagonal = FALSE;

   detectordata->nblocks = 0;
   detectordata->consismaster = NULL;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, detectordata, detectConnected, initConnected, exitConnected) );

   /* add connected constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/connected/setppcinmaster", "Controls whether SETPPC constraints chould be ignored while detecting and be directly placed in the master", &detectordata->setppcinmaster, FALSE, DEFAULT_SETPPCINMASTER, NULL, NULL) );

   return SCIP_OKAY;
}
