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

/**@file   heur_xpcrossover.c
 * @brief  Extreme Point Crossover
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "heur_xpcrossover.h"
#include "pub_gcgvar.h"
#include "relax_gcg.h"
#include "gcgplugins.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"


#define HEUR_NAME             "xpcrossover"
#define HEUR_DESC             "Extreme Point Crossover"
#define HEUR_DISPCHAR         'X'
#define HEUR_PRIORITY         -1100500
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE

#define DEFAULT_MAXNODES      1000LL        /**< maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01          /**< factor by which crossover should at least improve the incumbent     */
#define DEFAULT_MINNODES      200LL         /**< minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINFIXINGRATE 0.5           /**< minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_NODESOFS      200LL         /**< number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1           /**< subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_NUSEDPTS      2             /**< number of extreme pts per block that will be taken into account     */
#define DEFAULT_NWAITINGNODES 200LL         /**< number of nodes without incumbent change heuristic should wait      */
#define DEFAULT_RANDOMIZATION FALSE         /**< should the choice which sols to take be randomized?                 */
#define DEFAULT_DONTWAITATROOT FALSE        /**< should the nwaitingnodes parameter be ignored at the root node?     */
#define DEFAULT_USELPROWS     FALSE         /**< should subproblem be created out of the rows in the LP rows,
                                             * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE          /**< if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                             * of the original scip be copied to constraints of the subscip
                                             */

#define HASHSIZE_POINTS       11113         /**< size of hash table for extreme point tuples                         */




/*
 * Data structures
 */
typedef struct PointTuple POINTTUPLE;

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem               */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem               */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes        */
   SCIP_Longint          usednodes;          /**< nodes already used by crossover in earlier calls                  */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem     */

   int                   nusedpts;           /**< number of extreme pts per block that will be taken into account   */
   SCIP_Longint          nwaitingnodes;      /**< number of nodes without incumbent change heuristic should wait    */
   unsigned int          nfailures;          /**< number of failures since last successful call                     */
   SCIP_Longint          nextnodenumber;     /**< number of BnB nodes at which crossover should be called next      */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed     */
   SCIP_Real             minimprove;         /**< factor by which crossover should at least improve the incumbent   */
   SCIP_Bool             randomization;      /**< should the choice which sols to take be randomized?               */
   SCIP_Bool             dontwaitatroot;     /**< should the nwaitingnodes parameter be ignored at the root node?   */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?      */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?
                                              */
   unsigned int          randseed;           /**< seed value for random number generator                            */
   SCIP_HASHTABLE*       hashtable;          /**< hashtable used to store the extreme point tuples already used     */
   POINTTUPLE*           lasttuple;          /**< last tuple of extreme points                                      */
};

/** n-tuple of extreme points and their hashkey */
struct PointTuple
{
   int*                  indices;            /**< sorted array of master variable indices                          */
   int                   size;               /**< size of the array (should be heurdata->nusedpts * nblocks)       */
   unsigned int          key;                /**< hashkey of the tuple                                             */
   POINTTUPLE*           prev;               /**< previous extreme point tuple created                             */
};




/*
 * Local methods
 */


/** gets the hash key of an extreme point tuple */
static
SCIP_DECL_HASHGETKEY(hashGetKeyPts)
{  /*lint --e{715}*/
   return elem;
}


/** returns TRUE iff both extreme point tuples are identical */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqPts)
{  /*lint --e{715}*/
   int i;
   int size;

   int* indices1;
   int* indices2;

   indices1 = ((POINTTUPLE*)key1)->indices;
   indices2 = ((POINTTUPLE*)key2)->indices;

   /* there should be two nonempty arrays of the same size */
   assert(indices1 != NULL);
   assert(indices2 != NULL);
   assert(((POINTTUPLE*)key1)->size == ((POINTTUPLE*)key2)->size);

   size  = ((POINTTUPLE*)key1)->size;

   /* compare arrays by components, return TRUE, iff equal */
   for( i = 0; i < size; i++ )
   {
      if( indices1[i] != indices2[i] )
         return FALSE;
   }

   return TRUE;
}


/** returns hashkey of an extreme point tuple */
static
SCIP_DECL_HASHKEYVAL(hashKeyValPts)
{  /*lint --e{715}*/
   return ((POINTTUPLE*)key)->key;
}


/** calculates a hash key for a given tuple of master variable indices */
static
unsigned int calculateHashKey(
   int*                  indices,            /**< indices of master variables                                     */
   int                   size                /**< number of master variables                                      */
   )
{
   int i;
   unsigned int hashkey;

   /* haskey should be (x1+2) * (x2+2) * ... * (xn+2) + (x1+1) + (x2+1) + ... + (xn+1) */
   hashkey = 1;
   for( i = 0; i < size; i++ )
      hashkey *= indices[i] + 2;
   for( i = 0; i < size; i++ )
      hashkey += indices[i] + 1;

   return hashkey;
}


/** creates a new tuple of extreme points (= mastervars) */
static
SCIP_RETCODE createPtTuple(
   SCIP*                 scip,               /**< original SCIP data structure                                    */
   POINTTUPLE**          elem,               /**< tuple of master variables which should be created               */
   int*                  indices,            /**< indices of master variables                                     */
   int                   size,               /**< number of master variables                                      */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data                                           */
   )
{
   /* memory allocation */
   SCIP_CALL( SCIPallocBlockMemory(scip, elem) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*elem)->indices,size) );
   BMScopyMemoryArray((*elem)->indices, indices, size);

   /* data input */
   SCIPsortInt((*elem)->indices, size);
   (*elem)->size = size;
   (*elem)->key = calculateHashKey((*elem)->indices, (*elem)->size);
   (*elem)->prev = heurdata->lasttuple;

   /* update heurdata */
   heurdata->lasttuple = *elem;
   return SCIP_OKAY;
}


/** for each block, select extreme points (represented by mastervars) to be crossed */
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
   SCIP_Real* mastervals;
   int nmastervars;

   int nusedpts;
   int block;
#ifndef NDEBUG
   int nidentblocks;
#endif
   int* blocknrs;
   int* identblock;
   SCIP_Real* blockvalue;
   SCIP_Real value;
   SCIP_Real* selvalue;

   POINTTUPLE* elem;

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

   /* get master LP solution values */
   SCIP_CALL( SCIPallocBufferArray(scip, &mastervals, nmastervars) );
   SCIP_CALL( SCIPgetSolVals(masterprob, NULL, nmastervars, mastervars, mastervals) );

   /* get number of extreme points per block */
   nusedpts = heurdata->nusedpts;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &selvalue, nblocks * nusedpts) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocknrs, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockvalue, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &identblock, nblocks) );

   /* initialize the block values for the pricing problems */
   for( i = 0; i < nblocks; i++ )
   {
      blockvalue[i] = 0.0;
      blocknrs[i] = 0;
      identblock[i] = i;
   }

   *success = FALSE;

   /* loop over all given master variables;
    * this loop treats master variables that have value one or greater
    * (in particular important if blocks are represented by others) */
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
      if( SCIPisFeasZero(scip, value) )
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
#ifndef NDEBUG
      nidentblocks = GCGrelaxGetNIdenticalBlocks(scip, block);
#endif

      while( SCIPisFeasGE(scip, mastervals[i], 1.0) )
      {
         /* insert the extreme point in the selection (should be the only point for this block) */
         j = identblock[block] * nusedpts;
         assert(selection[j] == -1);

         SCIPdebugMessage("insert new point: block %d, mastervar %d, value %g, pos %d\n",
            identblock[block]+1, i, 1.0, j % nusedpts);
         selection[j] = i;
         selvalue[j] = mastervals[i];

         mastervals[i] = mastervals[i] - 1.0;
         blocknrs[block]++;

         /* search the next block to be considered */
         for( j = identblock[block] + 1; j < nblocks; ++j )
            if( GCGrelaxGetBlockRepresentative(scip, j) == block )
            {
               identblock[block] = j;
               break;
            }

#ifndef NDEBUG
         assert(blocknrs[block] >= nidentblocks || j < nblocks);
#endif
      }
   }

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
      if( SCIPisFeasZero(scip, value) )
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
#ifndef NDEBUG
      nidentblocks = GCGrelaxGetNIdenticalBlocks(scip, block);
#endif

      assert(SCIPisFeasGE(scip, mastervals[i], 0.0) && SCIPisFeasLT(scip, mastervals[i], 1.0));

      while( SCIPisFeasPositive(scip, mastervals[i]) )
      {
         value = MIN(mastervals[i], 1.0 - blockvalue[block]);

         /* check if the extreme point is good enough to be inserted in the selection */
         for( j = identblock[block] * nusedpts; j < (identblock[block] + 1) * nusedpts; ++j )
         {
            /* if the extreme point is better than a point in the selection
             * or there are < nusedpts, insert it */
            if( selection[j] == -1 || SCIPisGT(scip, value, selvalue[j]) )
            {
               SCIPdebugMessage("insert new point: block %d, mastervar %d, value %g, pos %d\n",
                     identblock[block]+1, i, value, j % nusedpts);
               for( k = (identblock[block] + 1) * nusedpts - 1; k > j; --k )
               {
                  SCIPdebugMessage("  shift point %d from pos %d to pos %d\n", selection[k-1], (k-1) % nusedpts, k % nusedpts);
                  selection[k] = selection[k-1];
                  selvalue[k] = selvalue[k-1];
               }
               selection[j] = i;
               selvalue[j] = value;
               break;
            }
         }

         mastervals[i] = mastervals[i] - value;
         if( SCIPisFeasZero(scip, mastervals[i]) )
            mastervals[i] = 0.0;
         blockvalue[block] += value;

         /* if the value assigned to the block is equal to 1, this block is full and we consider the next block */
         if( SCIPisFeasGE(scip, blockvalue[block], 1.0) )
         {
            blockvalue[block] = 0.0;
            blocknrs[block]++;

            /* search the next block to be considered */
            for( j = identblock[block] + 1; j < nblocks; ++j )
               if( GCGrelaxGetBlockRepresentative(scip, j) == block )
               {
                  identblock[block] = j;
                  break;
               }

#ifndef NDEBUG
            assert(blocknrs[block] >= nidentblocks || j < nblocks);
#endif
         }
      }
   }

   /* creates an object ready to be inserted into the hashtable */
   SCIP_CALL( createPtTuple(scip, &elem, selection, nusedpts * nblocks, heurdata) );

   /* check whether the set is already in the hashtable, if not, insert it */
   if( !SCIPhashtableExists(heurdata->hashtable, elem) )
   {
      SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
      *success = TRUE;
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &identblock);
   SCIPfreeBufferArray(scip, &blockvalue);
   SCIPfreeBufferArray(scip, &blocknrs);
   SCIPfreeBufferArray(scip, &selvalue);
   SCIPfreeBufferArray(scip, &mastervals);

   return SCIP_OKAY;
}


/** select extreme points (represented by mastervars) to be crossed randomly*/
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
   int iters;            /* iteration counter                                      */

   POINTTUPLE* elem;

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
   iters = 0;

   /* perform at maximum 10 restarts and stop as soon as a new set of solutions is found */
   do
   {
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

      /* creates an object ready to be inserted into the hashtable */
      SCIP_CALL( createPtTuple(scip, &elem, selection, nusedpts * nblocks, heurdata) );

      /* check whether the randomized set is already in the hashtable, if not, insert it */
      if( !SCIPhashtableExists(heurdata->hashtable, elem) )
      {
         SCIP_CALL( SCIPhashtableInsert(heurdata->hashtable, elem) );
         *success = TRUE;
      }
      ++iters;
   }
   while( !*success && iters < 10 );

   /* free memory */
   SCIPfreeBufferArray(scip, &npts);

   return SCIP_OKAY;
}


#ifdef PRINTPOINTS
/** prints selected extreme points to standard output */
static
SCIP_RETCODE printExtremePoints(
   SCIP*                 scip,
   int                   nusedpts,
   int*                  selection
      )
{
   SCIP* masterprob;
   int nblocks;

   SCIP_VAR** mastervars;
   int nmastervars;
   SCIP_VAR* mastervar;
   SCIP_VAR** origvars;
   SCIP_Real* origvals;
   int norigvars;

   int i;
   int j;
   int k;
   int selidx;

   assert(scip != NULL);

   /* get master problem, number of blocks and extreme points per block */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);

   nblocks = GCGrelaxGetNPricingprobs(scip);

   /* get variables of the master problem */
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);
   assert(nmastervars >= 0);

   /* first, print the relaxation solution */
   SCIPdebugPrintf("------------------------------------------------------------\n");
   SCIPdebugPrintf("Current relaxation solution:\n");
   SCIP_CALL( SCIPprintSol(scip, GCGrelaxGetCurrentOrigSol(scip), NULL, FALSE) );
   SCIPdebugPrintf("------------------------------------------------------------\n");

   /* then, print the selected extreme points for each block */
   for( i = 0; i < nblocks; ++i )
   {
      SCIPdebugPrintf("Block %i\n", i+1);
      SCIPdebugPrintf("------------------------------------------------------------\n");

      for( j = 0; j < nusedpts; ++j )
      {
         selidx = i * nusedpts + j;
         if( selection[selidx] != -1 )
         {
            mastervar = mastervars[selection[selidx]];
            assert(GCGvarIsMaster(mastervar));

            origvars = GCGmasterVarGetOrigvars(mastervar);
            origvals = GCGmasterVarGetOrigvals(mastervar);
            norigvars = GCGmasterVarGetNOrigvars(mastervar);

            SCIPdebugPrintf("Extreme point %i, mastervar %s, masterval=%g, index=%d:\n",
                  j+1, SCIPvarGetName(mastervar),
                  SCIPgetSolVal(masterprob, NULL, mastervar), SCIPvarGetProbindex(mastervar));
            for( k = 0; k < norigvars; ++k )
            {
               SCIPdebugPrintf("%-32s % 20.15g \t(obj:%.15g)\n",
                  SCIPvarGetName(origvars[k]), origvals[k], SCIPvarGetObj(origvars[k]));
            }
            SCIPdebugPrintf("------------------------------------------------------------\n");
         }
      }
   }

   return SCIP_OKAY;
}
#endif

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

/** fix those variables which are identical in each extreme point of their blocks */
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

   SCIP_VAR** vars;                          /* original scip variables                */
   SCIP_Real fixingrate;                     /* percentage of variables that are fixed */
   
   int nblocks;                              /* number of blocks                                   */
   int nusedpts;                             /* number of extreme points per block                 */
   int nvars;                                /* number of original variables                       */
   int nbinvars;                             /* number of binary variables in the original problem */
   int nintvars;                             /* number of general integer variables                */

   SCIP_Real* fixvals;                       /* values to which original variables should be fixed    */
   SCIP_Bool* fixable;                       /* for each original variable, remember if it is fixable */
   int* ptcounter;                           /* for each original variable, count in how many extreme points it appears */
   int fixingcounter;                        /* count how many original variables are fixed           */
   int zerocounter;                          /* count how many variables are fixed to zero            */

   int i;
   int j;
   int k;
   int l;
   int idx;

   /* get master problem and its variables */
   masterprob = GCGrelaxGetMasterprob(scip);
   assert(masterprob != NULL);
   SCIP_CALL( SCIPgetVarsData(masterprob, &mastervars, NULL, NULL, NULL, NULL, NULL) );
   assert(mastervars != NULL);

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nblocks = GCGrelaxGetNPricingprobs(scip);
   nusedpts = heurdata->nusedpts;
   fixingcounter = 0;
   zerocounter = 0;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &fixvals, nbinvars + nintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixable, nbinvars + nintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ptcounter, nbinvars + nintvars) );

   /* by default, each original variable can be fixed to zero */
   for( i = 0; i < nbinvars + nintvars; ++i )
   {
      fixable[i] = TRUE;
      fixvals[i] = 0.0;
      ptcounter[i] = 0;
   }

   SCIPdebugMessage("comparing extreme points...\n");

   /* for each block, compare the selected extreme points */
   for( i = 0; i < nblocks; ++i )
   {
      int blockrep;
      SCIP_VAR* mastervar;
      SCIP_VAR** origvars;
      SCIP_Real* origvals;
      int norigvars;
      int selidx;

      /* get the block that represents this block (in case of aggregation) */
      blockrep = GCGrelaxGetBlockRepresentative(scip, i);

      /* at least one extreme point must have been selected */
#ifndef NDEBUG
      selidx = i * nusedpts;
      assert(selection[selidx] != -1);
#endif

      /* compare the selected extreme points, where the first point is the reference point */
      for( j = 0; j < nusedpts; ++j )
      {
         selidx = i * nusedpts + j;
         if( selection[selidx] != -1 )
         {
            /* get master variable */
            mastervar = mastervars[selection[selidx]];
            assert(GCGvarGetBlock(mastervar) == blockrep);

            /* get extreme point */
            origvars = GCGmasterVarGetOrigvars(mastervar);
            origvals = GCGmasterVarGetOrigvals(mastervar);
            norigvars = GCGmasterVarGetNOrigvars(mastervar);

            for( k = 0; k < norigvars; ++k )
            {
               SCIP_VAR* origvar;
               SCIP_VAR* pricingvar;
               SCIP_VAR** pricingorigvars;
               int npricingorigvars;
               SCIP_Bool firstblock;

               if( SCIPvarGetType(origvars[k]) > SCIP_VARTYPE_INTEGER )
                  continue;

               /* get the corresponding pricing variable;
                  check whether this is the first block in which this variable appears;
                  search for the right original variable (in case of aggregation) */
               if( GCGvarIsLinking(origvars[k]) )
               {
                  SCIP_VAR** linkingpricingvars;

                  linkingpricingvars = GCGlinkingVarGetPricingVars(origvars[k]);
#ifndef NDEBUG
                  pricingvar = linkingpricingvars[blockrep];
                  assert(pricingvar != NULL);
                  assert(GCGvarIsPricing(pricingvar));
#endif

                  /* for linking variables, also check whether this is
                     the first block the variable appears in */
                  for( l = 0; l < blockrep; ++l )
                     if( linkingpricingvars[l] != NULL )
                        break;
                  firstblock = (l == blockrep);

                  /** @todo can there be aggregated linking variables? */
                  origvar = origvars[k];
               }
               else
               {
                  pricingvar = GCGoriginalVarGetPricingVar(origvars[k]);
                  assert(pricingvar != NULL);
                  assert(GCGvarIsPricing(pricingvar));

                  firstblock = TRUE;

                  /* search the right original variable (in case of aggregation) */
                  npricingorigvars = GCGpricingVarGetNOrigvars(pricingvar);
                  pricingorigvars = GCGpricingVarGetOrigvars(pricingvar);
                  origvar = NULL;
                  for( l = 0; l < npricingorigvars; ++l )
                     if( GCGvarGetBlock(pricingorigvars[l]) == i )
                     {
                        origvar = pricingorigvars[l];
                        break;
                     }
                  assert(origvar != NULL);
               }

               /* get the variable index */
               idx = SCIPvarGetProbindex(origvar);
               assert(idx < nbinvars + nintvars);

               /* the first extreme point serves as a reference point */
               if( j == 0 && firstblock )
                  fixvals[idx] = origvals[k];
               /* the variable can not be be fixed if its value differs in the extreme points */
               else
                  if( fixable[idx] && !SCIPisEQ(scip, fixvals[idx], origvals[k]) )
                     fixable[idx] = FALSE;

               if( !SCIPisZero(scip, origvals[k]) )
                  ++ptcounter[idx];
            }
         }
      }
   }

   /* try to fix the binary and general integer variables */
   for( i = 0; i < nbinvars + nintvars; ++i )
   {
      SCIP_VAR* var;
      int block;                             /* current block we are working in                    */
      SCIP_Real lb;
      SCIP_Real ub;

      var = vars[i];
      assert(GCGvarIsOriginal(var));
      block = GCGvarGetBlock(var);

      /* we still need to treat variables belonging to no block (as they did not appear in any extreme point) */
      /* if the variable belongs to no block, fix it in a RENS-like fashion */
      if( block == -1 )
      {
         fixvals[i] = SCIPgetRelaxSolVal(scip, var);
         if( SCIPisFeasIntegral(scip, fixvals[i]) )
         {
            /* fix variable to current relaxation solution if it is integral,
             * use exact integral value, if the variable is only integral within numerical tolerances
             */
            fixvals[i] = SCIPfloor(scip, fixvals[i] + 0.5);
         }
         else
            fixable[i] = FALSE;
      }

      /* if the variable is assigned to a block, check whether it was equal in all extreme points */
      else
      {
         int nlinkblocks;

         assert(block >= 0 || block == -2);
         nlinkblocks = GCGvarIsLinking(var) ? GCGlinkingVarGetNBlocks(var) : 1;
         assert(ptcounter[i] <= nusedpts * nlinkblocks);

         /* a variable which has appeared nonzero in some points
          * should have appeared nonzero in all extreme points in order to be fixed */
         if( ptcounter[i] > 0 )
         {
            if( ptcounter[i] < nusedpts * nlinkblocks )
               fixable[i] = FALSE;
         }

         /* if the variable has not appeared in any extreme point, it should be fixed to zero */
#ifdef SCIP_DEBUG
         else
         {
            assert(fixable[i]);
            assert(SCIPisZero(scip, fixvals[i]));
         }
#endif
      }

      /* get variable's global bounds */
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      /* fixing value can be outside transformed global bounds */
      if( fixable[i] && (lb > fixvals[i] || fixvals[i] > ub) )
         fixable[i] = FALSE;

//      SCIPdebugMessage("Trying to fix variable %s: block=%d, fixval=%g, ptcounter=%d, fixable=%d\n",
//                  SCIPvarGetName(var), block+1, fixvals[i], ptcounter[i], fixable[i]);

      /* the variable can be fixed if it has not been marked unfixable and
       *  - it was directly transferred to the master problem and has integer value or
       *  - it appeared zero in all extreme points
       *  - it did not appear zero in some extreme pts and nonzero in other extreme pts
       */
      if( fixable[i] )
      {
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], fixvals[i]) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], fixvals[i]) );

         fixingcounter++;

         if( SCIPisZero(scip, fixvals[i]) )
            zerocounter++;
      }
   }

   fixingrate = (SCIP_Real)fixingcounter / (SCIP_Real)(MAX(nbinvars + nintvars, 1));

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
   SCIPfreeBufferArray(scip, &fixvals);
   SCIPfreeBufferArray(scip, &fixable);
   SCIPfreeBufferArray(scip, &ptcounter);

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
      SCIPdebugMessage("Extreme Point Crossover: new solution added.\n");
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
#define heurCopyXpcrossover NULL

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeXpcrossover)
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
SCIP_DECL_HEURINIT(heurInitXpcrossover)
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
   heurdata->lasttuple = NULL;
   heurdata->nfailures = 0;
   heurdata->nextnodenumber = 0;

   /* initialize hash table */
   SCIP_CALL( SCIPhashtableCreate(&heurdata->hashtable, SCIPblkmem(scip), HASHSIZE_POINTS,
         hashGetKeyPts, hashKeyEqPts, hashKeyValPts, NULL) );
   assert(heurdata->hashtable != NULL );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitXpcrossover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   POINTTUPLE* pointtuple;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   pointtuple = heurdata->lasttuple;

   /* free all pointtuples iteratively */
   while( pointtuple != NULL )
   {
      POINTTUPLE* tmp;
      tmp = pointtuple->prev;
      SCIPfreeBlockMemoryArray(scip, &pointtuple->indices, pointtuple->size);
      SCIPfreeBlockMemory(scip, &pointtuple);
      pointtuple = tmp;
   }

   /* free hash table */
   assert(heurdata->hashtable != NULL );
   SCIPhashtableFree(&heurdata->hashtable);

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolXpcrossover NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolXpcrossover NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecXpcrossover)
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
      SCIPdebugMessage("skipping Extreme Point Crossover: invalid relaxation solution\n");
      return SCIP_OKAY;
   }

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetStage(masterprob) > SCIP_STAGE_SOLVING || SCIPgetLPSolstat(masterprob) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMessage("skipping Extreme Point Crossover: master LP not solved to optimality.\n");
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
      SCIPdebugMessage("skipping Extreme Points Crossover: nstallnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->minnodes);
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

   SCIPdebugMessage("Executing Extreme Point Crossover heuristic ...\n");

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &selection, nblocks * heurdata->nusedpts) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, SCIPgetNVars(scip)) );

   /* initialize empty selection */
   for( i = 0; i < nblocks * heurdata->nusedpts; ++i )
      selection[i] = -1;

   /* for each block, select extreme points (represented by master variables) to perform a crossover */
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

#ifdef PRINTPOINTS
   SCIP_CALL( printExtremePoints(scip, heurdata->nusedpts, selection) );
#endif

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
      SCIP_CALL( SCIPfree(&subscip) );
      SCIPfreeBufferArray(scip, &subvars);
      SCIPfreeBufferArray(scip, &selection);

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
      SCIPwarningMessage(scip, "Error while solving subMIP in Extreme Point Crossover heuristic; subSCIP terminated with code <%d>\n",
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

      SCIPdebugMessage(" -> found %i feasible solution(s).\n", SCIPgetNSols(subscip));

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
   SCIPfreeBufferArray(scip, &selection);
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the Extreme Point Crossover primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurXpcrossover(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create Extreme Point Crossover primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyXpcrossover, heurFreeXpcrossover, heurInitXpcrossover, heurExitXpcrossover,
         heurInitsolXpcrossover, heurExitsolXpcrossover, heurExecXpcrossover,
         heurdata) );

   /* add Extreme Point Crossover primal heuristic parameters */

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
         &heurdata->nusedpts, FALSE, DEFAULT_NUSEDPTS, 2, INT_MAX, NULL, NULL) );

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
