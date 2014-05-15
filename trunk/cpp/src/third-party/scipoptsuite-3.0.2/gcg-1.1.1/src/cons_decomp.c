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

/**@file   cons_decomp.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for structure detection
 * @author Martin Bergner
 *
 * This constraint handler will run all registered structure detectors in
 * increasing priority until the first detector finds a suitable structure.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cons_decomp.h"
#include "dec_connected.h"
#include "relax_gcg.h"
#include "struct_detector.h"
#include "string.h"
#include "scip_misc.h"
#include "scip/clock.h"
#include "pub_gcgvar.h"
#include "pub_decomp.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "decomp"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                          *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA         TRUE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP         TRUE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL       TRUE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define DEFAULT_CREATEBASICDECOMP FALSE /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
/*
 * Data structures
 */

/** score data structure **/
struct SCIP_DecScores
{
   SCIP_Real             borderscore;        /**< score of the border */
   SCIP_Real             densityscore;       /**< score of block densities */
   SCIP_Real             linkingscore;       /**< score related to interlinking blocks */
   SCIP_Real             totalscore;         /**< accumulated score */
};
typedef struct SCIP_DecScores SCIP_DECSCORES;

/** constraint handler data */
struct SCIP_ConshdlrData
{
   DEC_DECOMP**          decdecomps;         /**< array of decomposition structures */
   DEC_DETECTOR**        detectors;          /**< array of structure detectors */
   int*                  priorities;         /**< priorities of the detectors */
   int                   ndetectors;         /**< number of detectors */
   SCIP_CLOCK*           detectorclock;      /**< clock to measure detection time */
   SCIP_Bool             hasrun;             /**< flag to indicate whether we have already detected */
   int                   ndecomps;           /**< number of decomposition structures  */
   SCIP_Bool             createbasicdecomp;  /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** computes the score of the given decomposition */
static
SCIP_RETCODE evaluateDecomposition(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition data structure */
   SCIP_DECSCORES*       score               /**< returns the score of the decomposition */
      )
{
   int matrixarea;
   int borderarea;
   int nvars;
   int nconss;
   int i;
   int j;
   int k;
   /*   int blockarea; */
   SCIP_Real varratio;
   int* nzblocks;
   int nblocks;
   int* nlinkvarsblocks;
   int* nvarsblocks;
   SCIP_Real* blockdensities;
   int* blocksizes;
   SCIP_Real density;

   assert(scip != NULL);
   assert(score != NULL);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);

   nblocks = DECdecompGetNBlocks(decdecomp);

   SCIP_CALL( SCIPallocBufferArray(scip, &nzblocks, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlinkvarsblocks, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockdensities, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocksizes, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nvarsblocks, nblocks) );
   /*
    * 3 Scores
    *
    * - Area percentage (min)
    * - block density (max)
    * - \pi_b {v_b|v_b is linking}/#vb (min)
    */

   /* calculate matrix area */
   matrixarea = nvars*nconss;

   /* calculate slave sizes, nonzeros and linkingvars */
   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CONS** curconss;
      int ncurconss;
      int nvarsblock;
      SCIP_Bool *ishandled;

      SCIP_CALL( SCIPallocBufferArray(scip, &ishandled, nvars) );
      nvarsblock = 0;
      nzblocks[i] = 0;
      nlinkvarsblocks[i] = 0;
      for( j = 0; j < nvars; ++j )
      {
         ishandled[j] = FALSE;
      }
      curconss = DECdecompGetSubscipconss(decdecomp)[i];
      ncurconss = DECdecompGetNSubscipconss(decdecomp)[i];

      for( j = 0; j < ncurconss; ++j )
      {
         SCIP_VAR** curvars;
         SCIP_VAR* var;
         int ncurvars;
         ncurvars = SCIPgetNVarsXXX(scip, curconss[j]);
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, curconss[j], curvars, ncurvars) );

         for( k = 0; k < ncurvars; ++k )
         {
            int block;
            if( !SCIPisVarRelevant(curvars[k]) )
               continue;

            var = SCIPvarGetProbvar(curvars[k]);
            assert(var != NULL);
            if( !SCIPisVarRelevant(var) )
               continue;

            assert(SCIPvarIsActive(var));
            assert(!SCIPvarIsDeleted(var));
            ++(nzblocks[i]);
            if( !SCIPhashmapExists(DECdecompGetVartoblock(decdecomp), var) )
            {
               block = (int)(size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), curvars[k]); /*lint !e507*/
            }
            else
            {
               assert(SCIPhashmapExists(DECdecompGetVartoblock(decdecomp), var));
               block = (int)(size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), var); /*lint !e507*/
            }

            if( block == nblocks+1 && ishandled[SCIPvarGetProbindex(var)] == FALSE )
            {
               ++(nlinkvarsblocks[i]);
            }
            ishandled[SCIPvarGetProbindex(var)] = TRUE;
         }

         SCIPfreeBufferArray(scip, &curvars);
      }

      for( j = 0; j < nvars; ++j )
      {
         if( ishandled[j] )
         {
            ++nvarsblock;
         }
      }

      blocksizes[i] = nvarsblock*ncurconss;
      nvarsblocks[i] = nvarsblock;
      if( blocksizes[i] > 0 )
      {
         blockdensities[i] = 1.0*nzblocks[i]/blocksizes[i];
      }
      else
      {
         blockdensities[i] = 0.0;
      }

      assert(blockdensities[i] >= 0 && blockdensities[i] <= 1.0);
      SCIPfreeBufferArray(scip, &ishandled);
   }

   /* calculate border area */
   borderarea = DECdecompGetNLinkingconss(decdecomp)*nvars+DECdecompGetNLinkingvars(decdecomp)*(nconss-DECdecompGetNLinkingconss(decdecomp));

   /*   blockarea = 0; */
   density = 1E20;
   varratio = 1.0;
   for( i = 0; i < nblocks; ++i )
   {
      /* calculate block area */
      /* blockarea += blocksizes[i]; */


      /* calculate density */
      density = MIN(density, blockdensities[i]);

      /* calculate linking var ratio */
      if( DECdecompGetNLinkingvars(decdecomp) > 0 )
      {
         varratio *= 1.0*nlinkvarsblocks[i]/DECdecompGetNLinkingvars(decdecomp);
      }
      else
      {
         varratio = 0;
      }
   }

   score->linkingscore = (0.5+0.5*varratio);
   score->borderscore = (1.0*(borderarea)/matrixarea);
   score->densityscore = (1-density);

   switch( DECdecompGetType(decdecomp) )
   {
   case DEC_DECTYPE_ARROWHEAD:
      score->totalscore = score->borderscore*score->linkingscore*score->densityscore;
      break;
   case DEC_DECTYPE_BORDERED:
      score->totalscore = score->borderscore*score->linkingscore*score->densityscore;
      break;
   case DEC_DECTYPE_DIAGONAL:
      score->totalscore = 0.0;
      break;
   case DEC_DECTYPE_STAIRCASE:
      SCIPwarningMessage(scip, "Decomposition type is %s, cannot compute score\n", DECgetStrType(DECdecompGetType(decdecomp)));
      score->totalscore = 0.1;
      break;
   case DEC_DECTYPE_UNKNOWN:
      SCIPerrorMessage("Decomposition type is %s, cannot compute score\n", DECgetStrType(DECdecompGetType(decdecomp)));
      assert(FALSE);
      break;
   default:
      SCIPerrorMessage("No rule for this decomposition type, cannot compute score\n");
      assert(FALSE);
      break;
   }

   SCIPfreeBufferArray(scip, &nzblocks);
   SCIPfreeBufferArray(scip, &nlinkvarsblocks);
   SCIPfreeBufferArray(scip, &blockdensities);
   SCIPfreeBufferArray(scip, &blocksizes);
   SCIPfreeBufferArray(scip, &nvarsblocks);

   return SCIP_OKAY;

}


/*
 * Callback methods of constraint handler
 */

#define conshdlrCopyDecomp NULL
#define consInitpreDecomp NULL
#define consExitpreDecomp NULL
#define consDeleteDecomp NULL
#define consTransDecomp NULL
#define consInitlpDecomp NULL
#define consSepalpDecomp NULL
#define consSepasolDecomp NULL
#define consPropDecomp NULL
#define consPresolDecomp NULL
#define consRespropDecomp NULL
#define consActiveDecomp NULL
#define consDeactiveDecomp NULL
#define consEnableDecomp NULL
#define consDisableDecomp NULL
#define consDelvarsDecomp NULL
#define consPrintDecomp NULL
#define consCopyDecomp NULL
#define consParseDecomp NULL
#define consGetVarsDecomp NULL
#define consGetNVarsDecomp NULL

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitDecomp)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->hasrun = FALSE;
   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitDecomp)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(conshdlr != NULL);
   assert(scip != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->ndecomps > 0 )
   {
      int i;
      for( i = 0; i < conshdlrdata->ndecomps; ++i )
      {
         SCIP_CALL( DECdecompFree(scip, &conshdlrdata->decdecomps[i]) );
      }
      SCIPfreeMemoryArray(scip, &conshdlrdata->decdecomps);
      conshdlrdata->decdecomps = NULL;
      conshdlrdata->ndecomps = 0;
   }
   conshdlrdata->hasrun = FALSE;
   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeDecomp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->detectorclock) );

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);
      if( detector->exitDetection != NULL )
      {
         SCIPdebugMessage("Calling exitDetection of %s\n", detector->name);
         SCIP_CALL( (*detector->exitDetection)(scip, detector) );
      }
   }

   if( conshdlrdata->ndecomps > 0 )
   {
      for( i = 0; i < conshdlrdata->ndecomps; ++i )
      {
         SCIP_CALL( DECdecompFree(scip, &conshdlrdata->decdecomps[i]) );
      }
      SCIPfreeMemoryArray(scip, &conshdlrdata->decdecomps);
   }


   SCIPfreeMemoryArray(scip, &conshdlrdata->priorities);
   SCIPfreeMemoryArray(scip, &conshdlrdata->detectors);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->decdecomps);
   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolDecomp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolDecomp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpDecomp)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsDecomp)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckDecomp)
{
   /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}



/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockDecomp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for decomp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create decomp constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   assert(conshdlrdata != NULL);

   conshdlrdata->decdecomps = NULL;
   conshdlrdata->ndecomps = 0;
   conshdlrdata->ndetectors = 0;
   conshdlrdata->priorities = NULL;
   conshdlrdata->detectors = NULL;
   conshdlrdata->hasrun = FALSE;

   SCIP_CALL( SCIPcreateWallClock(scip, &conshdlrdata->detectorclock) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         SCIP_PROPTIMING_AFTERLPNODE, conshdlrCopyDecomp,
         consFreeDecomp, consInitDecomp, consExitDecomp,
         consInitpreDecomp, consExitpreDecomp, consInitsolDecomp, consExitsolDecomp,
         consDeleteDecomp, consTransDecomp, consInitlpDecomp,
         consSepalpDecomp, consSepasolDecomp, consEnfolpDecomp, consEnfopsDecomp, consCheckDecomp,
         consPropDecomp, consPresolDecomp, consRespropDecomp, consLockDecomp,
         consActiveDecomp, consDeactiveDecomp,
         consEnableDecomp, consDisableDecomp,
         consDelvarsDecomp, consPrintDecomp, consCopyDecomp, consParseDecomp,
         consGetVarsDecomp, consGetNVarsDecomp,
         conshdlrdata) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/decomp/createbasicdecomp", "indicates whether to create a decomposition with all constraints in the master if no other specified", &conshdlrdata->createbasicdecomp, FALSE, DEFAULT_CREATEBASICDECOMP, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a decomp constraint */
SCIP_RETCODE SCIPcreateConsDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name                /**< name of constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of decomp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the decomp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("decomp constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         FALSE, FALSE, FALSE, TRUE, FALSE) );

   return SCIP_OKAY;
}

/** sets (and adds) the decomposition structure **/
SCIP_RETCODE SCIPconshdlrDecompAddDecdecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->ndecomps == 0 )
   {
      assert(conshdlrdata->decdecomps == NULL);
      SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->decdecomps, 1) ); /*lint !e506*/
      conshdlrdata->decdecomps[0] = decdecomp;
      conshdlrdata->ndecomps = 1;
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->decdecomps, conshdlrdata->ndecomps+1) );
      conshdlrdata->decdecomps[conshdlrdata->ndecomps] = conshdlrdata->decdecomps[0];
      conshdlrdata->decdecomps[0] = decdecomp;
      conshdlrdata->ndecomps += 1;
   }
   return SCIP_OKAY;
}


/** returns the decomposition structure **/
DEC_DECOMP** SCIPconshdlrDecompGetDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->decdecomps;
}

/** returns the decomposition structure **/
int SCIPconshdlrDecompGetNDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->ndecomps;
}

/** returns the data of the provided detector */
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR*         detector            /**< detector data structure */
   )
{
   assert(detector != NULL);
   return detector->decdata;

}

/** returns the name of the provided detector */
const char* DECdetectorGetName(
   DEC_DETECTOR*         detector            /**< detector data structure */
   )
{
   assert(detector != NULL);
   return detector->name;
}

/** searches for the detector and returns it or returns NULL if detector is not found*/
DEC_DETECTOR* DECfindDetector(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the detector */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);
      if( strcmp(detector->name, name) == 0 )
      {
         return detector;
      }
   }

   return NULL;
}

/** includes the detector */
SCIP_RETCODE DECincludeDetector(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the detector */
   const char            decchar,            /**< display character of the detector */
   const char*           description,        /**< description of the detector */
   int                   priority,           /**< priority of the detector */
   SCIP_Bool             enabled,            /**< whether the detector should be enabled by default */
   DEC_DETECTORDATA*     detectordata,       /**< the associated detector data (or NULL) */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)), /**< the method that will detect the structure (must not be NULL)*/
   DEC_DECL_INITDETECTOR((*initDetector)),   /**< initialization method of detector (or NULL) */
   DEC_DECL_EXITDETECTOR((*exitDetector))    /**< deinitialization method of detector (or NULL) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   DEC_DETECTOR *detector;
   char setstr[SCIP_MAXSTRLEN];
   char descstr[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(name != NULL);
   assert(description != NULL);
   assert(detectStructure != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   assert(detectStructure != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &detector) );
   assert(detector != NULL);

   SCIPdebugMessage("Adding detector %i: %s\n", conshdlrdata->ndetectors+1, name);

#ifndef NDEBUG
   assert(DECfindDetector(scip, name) == NULL);
#endif

   detector->decdata = detectordata;
   detector->name = name;
   detector->description = description;

   detector->detectStructure = detectStructure;

   detector->initDetection = initDetector;
   detector->exitDetection = exitDetector;
   detector->decchar = decchar;

   detector->priority = priority;
   detector->enabled = enabled;

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/enabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabled), FALSE, enabled, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/priority", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "priority of detector <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->priority), FALSE, priority, INT_MIN, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->detectors, conshdlrdata->ndetectors+1) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->priorities, conshdlrdata->ndetectors+1) );

   conshdlrdata->detectors[conshdlrdata->ndetectors] = detector;
   conshdlrdata->ndetectors = conshdlrdata->ndetectors+1;

   return SCIP_OKAY;

}

/** returns the remaning time of scip that the decomposition may use */
SCIP_Real DECgetRemainingTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real timelimit;
   assert(scip != NULL);
   SCIP_CALL_ABORT(SCIPgetRealParam(scip, "limits/time", &timelimit));
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   return timelimit;
}

/** interface method to detect the structure */
SCIP_RETCODE DECdetectStructure(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< Result pointer to indicate whether some structure was found */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_Real* scores;
   int i;
   int j;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "No problem exists, cannot detect structure!\n");

      if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
         conshdlrdata->hasrun = TRUE;

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
      SCIP_CALL( SCIPtransformProb(scip) );

   SCIP_CALL( SCIPresetClock(scip, conshdlrdata->detectorclock) );
   SCIP_CALL( SCIPstartClock(scip, conshdlrdata->detectorclock) );

   if( conshdlrdata->ndecomps == 0 )
   {
      for( i = 0; i < conshdlrdata->ndetectors; ++i )
      {
         DEC_DETECTOR *detector;
         detector = conshdlrdata->detectors[i];
         assert(detector != NULL);
         conshdlrdata->priorities[i] = detector->priority;
      }

      SCIPdebugMessage("Sorting %i detectors\n", conshdlrdata->ndetectors);
      SCIPsortIntPtr(conshdlrdata->priorities, (void**)conshdlrdata->detectors, conshdlrdata->ndetectors);

      SCIPdebugMessage("Trying %d detectors.\n", conshdlrdata->ndetectors);
      for( i = 0; i < conshdlrdata->ndetectors; ++i )
      {
         DEC_DETECTOR* detector;
         DEC_DECOMP** decdecomps;
         int ndecdecomps;

         ndecdecomps = -1;
         detector = conshdlrdata->detectors[i];
         assert(detector != NULL);
         if( !detector->enabled )
            continue;
         if( detector->initDetection != NULL )
         {
            SCIPdebugMessage("Calling initDetection of %s\n", detector->name);
            SCIP_CALL( (*detector->initDetection)(scip, detector) );
         }

         SCIPdebugMessage("Calling detectStructure of %s: ", detector->name);
         SCIP_CALL( (*detector->detectStructure)(scip, detector->decdata, &decdecomps, &ndecdecomps,  result) );
         if( *result == SCIP_SUCCESS )
         {
            assert(ndecdecomps >= 0);
            assert(decdecomps != NULL || ndecdecomps == 0);
            SCIPdebugPrintf("we have %d decompositions!\n", ndecdecomps);
            for( j = 0; j < ndecdecomps; ++j )
            {
               assert(decdecomps != NULL);
               DECdecompSetDetector(decdecomps[j], detector);
            }
            SCIP_CALL( SCIPreallocMemoryArray(scip, &(conshdlrdata->decdecomps), (conshdlrdata->ndecomps+ndecdecomps)) );
            BMScopyMemoryArray(&(conshdlrdata->decdecomps[conshdlrdata->ndecomps]), decdecomps, ndecdecomps); /*lint !e866*/
            SCIPfreeMemoryArray(scip, &decdecomps);
            conshdlrdata->ndecomps += ndecdecomps;
         }
         else
         {
            SCIPdebugPrintf("Failure!\n");
         }
      }
   }
   else
   {
      SCIP_CALL( DECdecompTransform(scip, conshdlrdata->decdecomps[0]) );
   }
   /* evaluate all decompositions and sort them by score */
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, conshdlrdata->ndecomps) );
   for( i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      SCIP_DECSCORES score;
      score.totalscore = 0.0;

      SCIP_CALL( evaluateDecomposition(scip, conshdlrdata->decdecomps[i], &score) );
      scores[i] = score.totalscore;
   }

   SCIPsortRealPtr(scores, (void**)conshdlrdata->decdecomps, conshdlrdata->ndecomps);
   SCIPfreeBufferArray(scip, &scores);

   SCIP_CALL( SCIPstopClock(scip, conshdlrdata->detectorclock) );

   if( conshdlrdata->ndecomps > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Chosen decomposition with %d blocks of type %s.\n",
         DECdecompGetNBlocks(conshdlrdata->decdecomps[0]), DECgetStrType(DECdecompGetType(conshdlrdata->decdecomps[0])));
      GCGsetStructDecdecomp(scip, conshdlrdata->decdecomps[0]);
      *result = SCIP_SUCCESS;
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "No decomposition found!\n");
   }
   SCIPdebugMessage("Detection took %fs\n", SCIPclockGetTime(conshdlrdata->detectorclock));

   /* show that we done our duty */
   conshdlrdata->hasrun = TRUE;

   return SCIP_OKAY;
}

/** write out all known decompositions */
SCIP_RETCODE DECwriteAllDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 extension           /**< extension for decompositions */
   )
{
   int i;
   char name[SCIP_MAXSTRLEN];
   char outname[SCIP_MAXSTRLEN];
   char *pname;
   char decchar;

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   DEC_DECOMP *tmp;
   DEC_DETECTOR *detector;

   assert(scip != NULL);
   assert(extension != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s",  SCIPgetProbName(scip));
   SCIPsplitFilename(name, NULL, &pname, NULL, NULL);

   /** @todo This is a giant hack, but it works quite well */
   tmp = conshdlrdata->decdecomps[0];

   for( i = 0; i < conshdlrdata->ndecomps; ++i )
   {

      conshdlrdata->decdecomps[0] = conshdlrdata->decdecomps[i];

      detector = DECdecompGetDetector(conshdlrdata->decdecomps[i]);
      if( detector == NULL )
         decchar = '?';
      else
         decchar = detector->decchar;

      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%c_%d.%s", pname, decchar, DECdecompGetNBlocks(conshdlrdata->decdecomps[i]), extension);

      SCIP_CALL( SCIPwriteTransProblem(scip, outname, extension, FALSE) );
   }
   conshdlrdata->decdecomps[0] = tmp;

   return SCIP_OKAY;
}

/** returns the best known decomposition, if available and NULL otherwise */
DEC_DECOMP* DECgetBestDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->ndecomps > 0 )
      return conshdlrdata->decdecomps[0];

   else if ( conshdlrdata->createbasicdecomp)
   {
      SCIP_RETCODE retcode;
      DEC_DECOMP* decomp = NULL;
      retcode = DECcreateBasicDecomp(scip, &decomp);
      assert(retcode == SCIP_OKAY);
      assert(decomp != NULL );

      retcode = SCIPconshdlrDecompAddDecdecomp(scip, decomp);
      if( retcode != SCIP_OKAY )
      {
         SCIPerrorMessage("Could not add decomp to cons_decomp!\n");
         return NULL;
      }

      assert(conshdlrdata->ndecomps > 0);
      assert(conshdlrdata->decdecomps[0] != NULL);
      return conshdlrdata->decdecomps[0];
   }

   return NULL;
}

/** writes out a list of all detectors */
void DECprintListOfDetectors(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int ndetectors;
   int i;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ndetectors = conshdlrdata->ndetectors;

   SCIPdialogMessage(scip, NULL, " detector             priority char  description\n --------------       -------- ----  -----------\n");

   for( i = 0; i < ndetectors; ++i )
   {
      SCIPdialogMessage(scip, NULL,  " %-20s", conshdlrdata->detectors[i]->name);
      SCIPdialogMessage(scip, NULL,  " %8d    %c ", conshdlrdata->detectors[i]->priority, conshdlrdata->detectors[i]->decchar);
      SCIPdialogMessage(scip, NULL,  " %s\n", conshdlrdata->detectors[i]->description);
   }
}

/** returns whether the detection has been performed */
SCIP_Bool DEChasDetectionRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->hasrun;
}
