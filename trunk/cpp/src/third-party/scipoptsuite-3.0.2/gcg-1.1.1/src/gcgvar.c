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

/**@file   gcgvar.c
 * @brief  GCG variable access functions
 * @author Martin Bergner
 *
 * @todo capture and release variables stored in other variable's data?
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "pub_gcgvar.h"
#include "struct_vardata.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "scip_misc.h"
#include "scip/cons_linear.h"

#define STARTMAXMASTERVARS 10

/*
 * Vardata methods
 */

/** callback method called when an original GCG variable is deleted */
static
SCIP_DECL_VARDELORIG(GCGvarDelOrig)
{
   /*lint -e715 */
   if( (*vardata)->vartype == GCG_VARTYPE_ORIGINAL )
   {
      if( (*vardata)->blocknr == -2 )
      {
         int nblocks;
         int i;

         nblocks = GCGrelaxGetNPricingprobs(scip);
         assert(nblocks > 0);

         assert((*vardata)->data.origvardata.linkingvardata != NULL);
         if( (*vardata)->data.origvardata.linkingvardata->linkconss != NULL )
         {
            assert((*vardata)->data.origvardata.linkingvardata->pricingvars != NULL);

            for( i = 0; i < nblocks; i++ )
            {
               assert(((*vardata)->data.origvardata.linkingvardata->linkconss[i] == NULL)
                  == ((*vardata)->data.origvardata.linkingvardata->pricingvars[i] == NULL));
            }

            SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.origvardata.linkingvardata->linkconss), nblocks);
            SCIPfreeBlockMemoryArray(scip, &((*vardata)->data.origvardata.linkingvardata->pricingvars), nblocks);
         }
         SCIPfreeMemory(scip, &((*vardata)->data.origvardata.linkingvardata));
         (*vardata)->data.origvardata.linkingvardata = NULL;
      }
      assert((*vardata)->data.origvardata.linkingvardata == NULL);
      assert((*vardata)->data.origvardata.mastervars != NULL);
      assert((*vardata)->data.origvardata.mastervals != NULL);
      SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.mastervars));
      SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.mastervals));
      if( (*vardata)->data.origvardata.ncoefs > 0 )
      {
         assert((*vardata)->data.origvardata.coefs != NULL);
         assert((*vardata)->data.origvardata.masterconss != NULL);
         SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.coefs));
         SCIPfreeMemoryArray(scip, &((*vardata)->data.origvardata.masterconss));
      }
   }
   if( (*vardata)->vartype == GCG_VARTYPE_PRICING )
   {
      assert((*vardata)->data.pricingvardata.norigvars >= 1);
      SCIPfreeMemoryArray(scip, &((*vardata)->data.pricingvardata.origvars));
   }
   assert((*vardata)->vartype != GCG_VARTYPE_MASTER);
   SCIPfreeBlockMemory(scip, vardata);

   return SCIP_OKAY;
}


/** callback method called when a transformed GCG variable is deleted */
static
SCIP_DECL_VARDELTRANS(gcgvardeltrans)
{
   /*lint -e715 */
   assert((*vardata)->vartype == GCG_VARTYPE_MASTER);
   SCIPfreeMemoryArray(scip, &((*vardata)->data.mastervardata.origvals));
   SCIPfreeMemoryArray(scip, &((*vardata)->data.mastervardata.origvars));

   SCIPfreeBlockMemory(scip, vardata);

   return SCIP_OKAY;
}



/** returns TRUE or FALSE whether variable is a pricing variable or not */
SCIP_Bool GCGvarIsPricing(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->vartype == GCG_VARTYPE_PRICING;
}

/** returns TRUE or FALSE whether variable is a master variable or not */
SCIP_Bool GCGvarIsMaster(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->vartype == GCG_VARTYPE_MASTER;
}

/** returns TRUE or FALSE whether variable is a original variable or not */
SCIP_Bool GCGvarIsOriginal(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->vartype == GCG_VARTYPE_ORIGINAL;
}

/** returns TRUE or FALSE whether variable is a linking variable or not */
SCIP_Bool GCGvarIsLinking(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->blocknr == -2;
}

/** returns the pricing var of an original variable */
SCIP_VAR* GCGoriginalVarGetPricingVar(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->data.origvardata.linkingvardata == NULL);
   assert(!GCGvarIsLinking(var));
   if( vardata->data.origvardata.pricingvar != NULL )
      assert(GCGvarIsPricing(vardata->data.origvardata.pricingvar));
   return vardata->data.origvardata.pricingvar;
}

/** returns the pricing var of an original variable */
void GCGoriginalVarSetPricingVar(
   SCIP_VAR*             var,                /**< SCIP variable structure */
   SCIP_VAR*             pricingvar          /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(pricingvar != NULL);
   assert(GCGvarIsOriginal(var));
   assert(GCGvarIsPricing(pricingvar));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.linkingvardata == NULL);
   assert(!GCGvarIsLinking(var));
   vardata->data.origvardata.pricingvar = pricingvar;
}

/** creates the data for all variables of the original program */
SCIP_RETCODE GCGcreateOrigVarsData(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);

   vars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);

   /* loop over the variables in the original problem */
   for( i = 0; i < nvars; i++ )
   {
      assert(vars[i] != NULL);
      SCIP_CALL( GCGorigVarCreateData(scip, vars[i]) );
   }

   return SCIP_OKAY;
}

/** creates the data for a variable of the original program */
SCIP_RETCODE GCGorigVarCreateData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< pointer to variable object */
   )
{
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(var != NULL);
   assert(SCIPvarIsOriginal(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   /* create the vardata and initialize its values */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = GCG_VARTYPE_ORIGINAL;
   vardata->blocknr = -1;
   vardata->data.origvardata.pricingvar = NULL;
   vardata->data.origvardata.coefs = NULL;
   vardata->data.origvardata.masterconss = NULL;
   vardata->data.origvardata.ncoefs = 0;
   vardata->data.origvardata.nmastervars = 0;
   vardata->data.origvardata.maxmastervars = STARTMAXMASTERVARS;
   vardata->data.origvardata.linkingvardata = NULL;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.mastervars),
         vardata->data.origvardata.maxmastervars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.mastervals),
         vardata->data.origvardata.maxmastervars) );

   SCIPvarSetData(var, vardata);
   if( SCIPvarIsOriginal(var) )
   {
      SCIPvarSetDelorigData(var, GCGvarDelOrig);
      if( SCIPvarGetTransVar(var) != NULL )
      {
         SCIPvarSetData(SCIPvarGetProbvar(SCIPvarGetTransVar(var)), vardata);
      }
   }
   else
   {
      assert(SCIPvarIsTransformedOrigvar(var));
      SCIPvarSetDeltransData(var, GCGvarDelOrig);
   }

   return SCIP_OKAY;
}


/** returns the pricing variables of an linking variable */
SCIP_VAR** GCGlinkingVarGetPricingVars(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsLinking(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->pricingvars != NULL);

   return vardata->data.origvardata.linkingvardata->pricingvars;
}

/** sets the pricing var of the corresponding linking variable at the specified position */
void GCGlinkingVarSetPricingVar(
   SCIP_VAR*             origvar,            /**< original variable */
   int                   pricingprobnr,      /**< number of pricing problem */
   SCIP_VAR*             var                 /**< pricing variable */
   )
{
   SCIP_VARDATA* vardata;
   assert(origvar != NULL);
   assert(var != NULL);
   assert(pricingprobnr >= 0);

   assert(GCGvarIsLinking(origvar));
   assert(GCGvarIsPricing(var));

   vardata = SCIPvarGetData(origvar);
   assert(vardata != NULL);
   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->pricingvars != NULL);

   vardata->data.origvardata.linkingvardata->pricingvars[pricingprobnr] = var;
}

/** returns the blocks the linking variable is in */
SCIP_RETCODE GCGlinkingVarGetBlocks(
   SCIP_VAR*             var,                /**< SCIP variable structure */
   int                   nblocks,            /**< size of array blocks */
   int*                  blocks              /**< array to store the blocks of the linking variable */
   )
{
   SCIP_VARDATA* vardata;
   int i;
   int j;

   assert(var != NULL);
   assert(nblocks == 0 || blocks != NULL);

   assert(GCGvarIsLinking(var));
   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->nblocks > 0);

   /* the blocks array must be large enough to hold all block numbers */
   if( nblocks < vardata->data.origvardata.linkingvardata->nblocks )
   {
      SCIPerrorMessage("array too small to store all block numbers!\n");
      return SCIP_INVALIDDATA;
   }
   assert(nblocks >= vardata->data.origvardata.linkingvardata->nblocks);

   /* fill the blocks array */
   j = -1;
   for( i = 0; i < vardata->data.origvardata.linkingvardata->nblocks; ++i )
   {
      /* search the next block the linking variable is contained in */
      do
         ++j;
      while ( vardata->data.origvardata.linkingvardata->pricingvars[j] == NULL );
      blocks[i] = j;
   }

   return SCIP_OKAY;
}

/** returns the number of blocks the linking variable is in */
int GCGlinkingVarGetNBlocks(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);

   assert(GCGvarIsLinking(var));
   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->nblocks > 0);
   return vardata->data.origvardata.linkingvardata->nblocks;
}

/** returns the original var of a pricing variable */
SCIP_VAR* GCGpricingVarGetOriginalVar(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsPricing(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->data.pricingvardata.norigvars >= 0);
   assert(vardata->data.pricingvardata.origvars != NULL);
   assert(vardata->data.pricingvardata.origvars[0] != NULL);
   assert(vardata->blocknr >= 0); /* variable belongs to exactly one block */

   return vardata->data.pricingvardata.origvars[0];
}

/** adds the original var to the pricing variable */
SCIP_RETCODE GCGpricingVarAddOrigVar(
   SCIP*                 scip,               /**< SCIP variable structure */
   SCIP_VAR*             pricingvar,         /**< pricing variable */
   SCIP_VAR*             origvar             /**< original pricing variable */
   )
{
   SCIP_VARDATA* vardata;
   assert(pricingvar != NULL);
   assert(origvar != NULL);
   assert(GCGvarIsPricing(pricingvar));
   assert(GCGvarIsOriginal(origvar));

   vardata = SCIPvarGetData(pricingvar);
   assert(vardata != NULL);
   assert(vardata->data.pricingvardata.norigvars >= 0);
   assert(vardata->data.pricingvardata.origvars != NULL);
   assert(vardata->data.pricingvardata.origvars[0] != NULL);
   assert(vardata->blocknr >= 0); /* variable belongs to exactly one block */
   if( vardata->data.pricingvardata.norigvars >= 1 )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(vardata->data.pricingvardata.origvars),
            vardata->data.pricingvardata.norigvars + 1) );
   }
   vardata->data.pricingvardata.origvars[vardata->data.pricingvardata.norigvars] = origvar;
   vardata->data.pricingvardata.norigvars++;

   return SCIP_OKAY;
}

/** returns the number of master variables the original variable is contained in */
int GCGoriginalVarGetNMastervars(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   assert(vardata->data.origvardata.nmastervars >= 0);

   return vardata->data.origvardata.nmastervars;
}

/** returns the master variables the original variable is contained in */
SCIP_VAR** GCGoriginalVarGetMastervars(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.mastervars != NULL);
   return vardata->data.origvardata.mastervars;
}

/** returns the fraction of master variables the original variable is contained in */
SCIP_Real* GCGoriginalVarGetMastervals(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.mastervals != NULL);
   return vardata->data.origvardata.mastervals;
}

/** returns the coefficients of master constraints the original variable is contained in */
SCIP_Real* GCGoriginalVarGetCoefs(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.coefs != NULL || vardata->data.origvardata.ncoefs == 0 );
   return vardata->data.origvardata.coefs;
}

/** returns the number of coefficients of master constraints the original variable is contained in */
int GCGoriginalVarGetNCoefs(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.coefs != NULL || vardata->data.origvardata.ncoefs == 0 );
   return vardata->data.origvardata.ncoefs;
}

/** sets the number of master variables the original variable is contained in */
void GCGoriginalVarSetNCoefs(
   SCIP_VAR*             var,                /**< SCIP variable structure */
   int                   ncoefs              /**< number of coefficient to set */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(ncoefs >= 0);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.coefs != NULL || vardata->data.origvardata.ncoefs == 0 );
   if( ncoefs == 0 )
      assert(vardata->data.origvardata.coefs == NULL);

   vardata->data.origvardata.ncoefs = ncoefs;
}

/** adds a coefficient of the master variable to the coefs array for the resp. constraint */
SCIP_RETCODE GCGoriginalVarAddCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to add coef */
   SCIP_Real             val,                /**< coefficent to set */
   SCIP_CONS*            cons                /**< constraint the variable is in */
   )
{
   SCIP_VARDATA* vardata;
   assert(scip != NULL);
   assert(var != NULL);
   assert(!SCIPisZero(scip, val));
   assert(cons != NULL);
   assert(GCGvarIsOriginal(var));
   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   if( vardata->data.origvardata.ncoefs == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.coefs), 1) ); /*lint !e506*/
      SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.origvardata.masterconss), 1) ); /*lint !e506*/
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(vardata->data.origvardata.coefs), vardata->data.origvardata.ncoefs+1) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(vardata->data.origvardata.masterconss), vardata->data.origvardata.ncoefs+1) );
   }

   assert(vardata->data.origvardata.coefs != NULL);
   assert(vardata->data.origvardata.masterconss != NULL);

   vardata->data.origvardata.coefs[vardata->data.origvardata.ncoefs] = val;
   vardata->data.origvardata.masterconss[vardata->data.origvardata.ncoefs] = cons;
   vardata->data.origvardata.ncoefs++;

   return SCIP_OKAY;
}


/** returns the fraction of master variables the original variable is contained in */
SCIP_CONS** GCGoriginalVarGetMasterconss(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->data.origvardata.masterconss;
}

/** adds variable to a new block, making a linkingvariable out of it, if necessary */
SCIP_RETCODE GCGoriginalVarAddBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< var that is added to a block */
   int                   newblock            /**< the new block the variable will be in */
   )
{
   SCIP_VARDATA* vardata;
   int nblocks;
   int blocknr;
   assert(scip != NULL);
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   nblocks = GCGrelaxGetNPricingprobs(scip);
   assert(newblock >= 0 && newblock < nblocks);
   blocknr = GCGvarGetBlock(var);
   /* the variable was only in one block so far, so set up the linking variable data */
   if( blocknr > -1 )
   {
      SCIP_CALL( SCIPallocMemory(scip, &vardata->data.origvardata.linkingvardata) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vardata->data.origvardata.linkingvardata->pricingvars, nblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vardata->data.origvardata.linkingvardata->linkconss, nblocks) );
      BMSclearMemoryArray(vardata->data.origvardata.linkingvardata->pricingvars, nblocks);
      BMSclearMemoryArray(vardata->data.origvardata.linkingvardata->linkconss, nblocks);

      /* store old block; store the original variable, it will be exchanged for the correct pricing variable later */
      vardata->data.origvardata.linkingvardata->pricingvars[blocknr] = var;
      vardata->data.origvardata.linkingvardata->nblocks = 1;

      vardata->blocknr = -2;
   }
   assert(GCGvarIsLinking(var));

   /* store new block */
   if( vardata->data.origvardata.linkingvardata->pricingvars[newblock] == NULL )
   {
      assert(vardata->data.origvardata.linkingvardata->linkconss[newblock] == NULL);
      vardata->data.origvardata.linkingvardata->pricingvars[newblock] = var;
      vardata->data.origvardata.linkingvardata->nblocks++;
   }
   assert(vardata->data.origvardata.linkingvardata->nblocks <= nblocks);
   return SCIP_OKAY;
}


/** returns the linking constraints */
SCIP_CONS** GCGlinkingVarGetLinkingConss(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsOriginal(var));
   assert(GCGvarIsLinking(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->linkconss != NULL);
   return vardata->data.origvardata.linkingvardata->linkconss;
}

/** sets the linking constraints */
void GCGlinkingVarSetLinkingCons(
   SCIP_VAR*             var,                /**< variable data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   int                   index               /**< index of pricing problem */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(cons != NULL);
   assert(index >= 0);
   assert(GCGvarIsOriginal(var));
   assert(GCGvarIsLinking(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->linkconss != NULL);
   vardata->data.origvardata.linkingvardata->linkconss[index] = cons;
}

/** returns whether the master variable is a ray */
SCIP_Bool GCGmasterVarIsRay(
   SCIP_VAR*             var                 /**< variable data structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsMaster(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   return vardata->data.mastervardata.isray;
}

/** returns the number of original variables the master variable is contained in */
int GCGmasterVarGetNOrigvars(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsMaster(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.mastervardata.norigvars >= 0);
   assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
   assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);
   assert(vardata->blocknr != -1 || vardata->data.mastervardata.norigvars == 1 );

   return vardata->data.mastervardata.norigvars;
}

/** returns the original variables the master variable is contained in */
SCIP_VAR** GCGmasterVarGetOrigvars(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsMaster(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.mastervardata.origvars != NULL || vardata->data.mastervardata.norigvars == 0);
   assert(vardata->blocknr != -1 || vardata->data.mastervardata.origvars != NULL);
   assert(vardata->blocknr != -1 || vardata->data.mastervardata.origvars[0] != NULL);
   assert(vardata->blocknr != -1 || GCGvarGetBlock(vardata->data.mastervardata.origvars[0]) == -1
      || GCGvarIsLinking(vardata->data.mastervardata.origvars[0]));


   return vardata->data.mastervardata.origvars;
}

/** returns the fraction of original variables the master variable is contained in */
SCIP_Real* GCGmasterVarGetOrigvals(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsMaster(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.mastervardata.origvals != NULL || vardata->data.mastervardata.norigvars == 0);
   return vardata->data.mastervardata.origvals;
}

/** returns the number of original variables the pricing variable is contained in */
int GCGpricingVarGetNOrigvars(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsPricing(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.pricingvardata.norigvars >= 0);
   return vardata->data.pricingvardata.norigvars;
}

/** returns the original variables the pricing variable is contained in */
SCIP_VAR** GCGpricingVarGetOrigvars(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(GCGvarIsPricing(var));

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->data.pricingvardata.origvars != NULL);
   return vardata->data.pricingvardata.origvars;
}

/** returns the block of the variable */
int GCGvarGetBlock(
   SCIP_VAR*             var                 /**< SCIP variable structure */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(vardata->blocknr >= -2);
   return vardata->blocknr;
}

/** sets the block of the variable */
void GCGvarSetBlock(
   SCIP_VAR*             var,                /**< variable to set block for */
   int                   block               /**< block to set */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(block >= -1);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);
   vardata->blocknr = block;
}

/** returns TRUE if the linking variable is in the block, FALSE otherwise */
SCIP_Bool GCGisLinkingVarInBlock(
   SCIP_VAR*             var,                /**< variabel data structure */
   int                   block               /**< pricing problem number */
   )
{
   SCIP_VARDATA* vardata;
   assert(var != NULL);
   assert(block >= 0);

   vardata = SCIPvarGetData(var);
   assert(vardata != NULL);

   assert(GCGvarIsLinking(var));
   assert(GCGvarIsOriginal(var));

   assert(vardata->data.origvardata.linkingvardata != NULL);
   assert(vardata->data.origvardata.linkingvardata->pricingvars != NULL);

   return vardata->data.origvardata.linkingvardata->pricingvars[block] != NULL;

}

/** informs an original variable, that a variable in the master problem was created,
 * that contains a part of the original variable.
 * Saves this information in the original variable's data
 * @todo this method needs a little love
 */
SCIP_RETCODE GCGoriginalVarAddMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< original variable */
   SCIP_VAR*             var,                /**< master variable */
   SCIP_Real             val                 /**< fraction of the original variable */
   )
{
   SCIP_VARDATA* vardata;

   assert(scip != NULL);
   assert(origvar != NULL);
   assert(var != NULL);

   vardata = SCIPvarGetData(origvar);

   assert(vardata != NULL);
   assert(GCGvarIsOriginal(origvar));
   assert(vardata->data.origvardata.mastervars != NULL);
   assert(vardata->data.origvardata.mastervals != NULL);
   assert(vardata->data.origvardata.nmastervars >= 0);
   assert(vardata->data.origvardata.maxmastervars >= vardata->data.origvardata.nmastervars);

   /* realloc mastervars array of the original variable, if needed */
   if( vardata->data.origvardata.maxmastervars == vardata->data.origvardata.nmastervars )
   {
      SCIP_CALL( SCIPreallocMemoryArray(pricerdata->origprob, &(vardata->data.origvardata.mastervars),
            2*vardata->data.origvardata.maxmastervars) );
      SCIP_CALL( SCIPreallocMemoryArray(pricerdata->origprob, &(vardata->data.origvardata.mastervals),
            2*vardata->data.origvardata.maxmastervars) );
      SCIPdebugMessage("mastervars array of var %s resized from %d to %d\n", SCIPvarGetName(origvar),
         vardata->data.origvardata.maxmastervars, 2*vardata->data.origvardata.maxmastervars);
      vardata->data.origvardata.maxmastervars = 2*vardata->data.origvardata.maxmastervars;
   }
   /* add information to the original variable's vardata */
   vardata->data.origvardata.mastervars[vardata->data.origvardata.nmastervars] = var;
   vardata->data.origvardata.mastervals[vardata->data.origvardata.nmastervars] = val;
   vardata->data.origvardata.nmastervars++;

   return SCIP_OKAY;
}

/** informs an original variable, that a variable in the master problem was deleted,
 * that contains a part of the original variable.
 * Update the information in the original variable's data
 * @todo this method needs a little love
 */
SCIP_RETCODE GCGoriginalVarRemoveMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< original variable */
   SCIP_VAR*             var                 /**< master variable */
   )
{
   SCIP_VARDATA* vardata;
   int i;

   assert(scip != NULL);
   assert(origvar != NULL);
   assert(var != NULL);

   vardata = SCIPvarGetData(origvar);

   assert(vardata != NULL);
   assert(GCGvarIsOriginal(origvar));
   assert(vardata->data.origvardata.mastervars != NULL);
   assert(vardata->data.origvardata.mastervals != NULL);
   assert(vardata->data.origvardata.nmastervars > 0);
   assert(vardata->data.origvardata.maxmastervars >= vardata->data.origvardata.nmastervars);

   for( i = 0; i < vardata->data.origvardata.nmastervars; ++i )
   {
      if( vardata->data.origvardata.mastervars[i] == var )
      {
         vardata->data.origvardata.mastervars[i] = vardata->data.origvardata.mastervars[vardata->data.origvardata.nmastervars - 1];
         vardata->data.origvardata.mastervals[i] = vardata->data.origvardata.mastervals[vardata->data.origvardata.nmastervars - 1];
         (vardata->data.origvardata.nmastervars)--;

         break;
      }
   }
   assert(i <= vardata->data.origvardata.nmastervars);
#ifndef NDEBUG
   for( ; i < vardata->data.origvardata.nmastervars; ++i )
   {
      assert(vardata->data.origvardata.mastervars[i] != var);
   }
#endif

   return SCIP_OKAY;
}

/** creates the corresponding pricing variable for the given original variable */
SCIP_RETCODE GCGoriginalVarCreatePricingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             origvar,            /**< original variable */
   SCIP_VAR**            var                 /**< pricing variable */
   )
{
   SCIP_VARDATA* vardata;
   char name[SCIP_MAXSTRLEN];
   int pricingprobnr;
   assert(scip != NULL);
   assert(origvar != NULL);
   assert(var != NULL);
   assert(GCGvarIsOriginal(origvar));
   assert(!GCGvarIsLinking(origvar));
   assert(GCGoriginalVarGetPricingVar(origvar) == NULL);

   /* get the number of the pricing block to which the variable belongs */
   pricingprobnr = GCGvarGetBlock(origvar);

   /* create variable data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = GCG_VARTYPE_PRICING;
   vardata->blocknr = pricingprobnr;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(vardata->data.pricingvardata.origvars), 1) ); /*lint !e506*/
   vardata->data.pricingvardata.origvars[0] = origvar;
   vardata->data.pricingvardata.norigvars = 1;

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pr%d_%s", pricingprobnr, SCIPvarGetName(origvar));
   SCIP_CALL( SCIPcreateVar(scip, var, name, SCIPvarGetLbGlobal(origvar),
         SCIPvarGetUbGlobal(origvar), 0.0, SCIPvarGetType(origvar),
         TRUE, FALSE, GCGvarDelOrig, NULL, NULL, NULL, vardata) );

   return SCIP_OKAY;
}

/** creates the corresponding pricing variable for the given original variable */
SCIP_RETCODE GCGlinkingVarCreatePricingVar(
   SCIP*                 masterscip,         /**< master problem SCIP data structure */
   SCIP*                 pricingscip,        /**< pricing problem SCIP data structure */
   int                   pricingprobnr,      /**< number of the pricing problem */
   SCIP_VAR*             origvar,            /**< original variable */
   SCIP_VAR**            var,                /**< pointer to store new pricing variable */
   SCIP_CONS**           linkcons            /**< constraint linking pricing variables */
   )
{
   SCIP_VARDATA* vardata;
   char name[SCIP_MAXSTRLEN];

   assert(masterscip != NULL);
   assert(pricingscip != NULL);
   assert(pricingprobnr >= 0);
   assert(origvar != NULL);
   assert(GCGvarIsLinking(origvar));
   assert(var != NULL);
   assert(linkcons != NULL);

   /* create variable data */
   SCIP_CALL( SCIPallocBlockMemory(pricingscip, &vardata) );
   vardata->vartype = GCG_VARTYPE_PRICING;
   vardata->blocknr = pricingprobnr;
   SCIP_CALL( SCIPallocMemoryArray(pricingscip, &(vardata->data.pricingvardata.origvars), 1) ); /*lint !e506*/
   vardata->data.pricingvardata.origvars[0] = origvar;
   vardata->data.pricingvardata.norigvars = 1;

   /* create and add variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pr%d_%s", pricingprobnr, SCIPvarGetName(origvar));
   SCIP_CALL( SCIPcreateVar(pricingscip, var, name, SCIPvarGetLbGlobal(origvar),
         SCIPvarGetUbGlobal(origvar), 0.0, SCIPvarGetType(origvar),
         TRUE, FALSE, GCGvarDelOrig, NULL, NULL, NULL, vardata) );

   /* add corresponding linking constraint to the master problem */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "l_%s_%d", SCIPvarGetName(origvar), pricingprobnr);
   SCIP_CALL( SCIPcreateConsLinear(masterscip, linkcons, name, 0, NULL, NULL, 0.0, 0.0,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates the master var and initializes the vardata */
SCIP_RETCODE GCGcreateMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 pricingscip,        /**< pricing problem SCIP data structure */
   SCIP_VAR**            newvar,             /**< pointer to store new master variable */
   char*                 varname,            /**< new variable name */
   SCIP_Real             objcoeff,           /**< new objective coefficient */
   SCIP_VARTYPE          vartype,            /**< new variable type */
   SCIP_Bool             solisray,           /**< indicates whether new variable is a ray */
   int                   prob,               /**< number of pricing problem that created this variable */
   int                   nsolvars,           /**< number of variables in the solution */
   SCIP_Real*            solvals,            /**< values of variables in the solution */
   SCIP_VAR**            solvars             /**< variables with non zero coefficient in the solution */
   )
{
   SCIP_VARDATA* newvardata;
   int i;
   int j;
   SCIP_Bool trivialsol;

   assert(scip != NULL);
   assert(pricingscip != NULL);
   assert(newvar != NULL);
   assert(varname != NULL);
   assert(!SCIPisInfinity(pricingscip, ABS(objcoeff)));
   assert(vartype == SCIP_VARTYPE_INTEGER || vartype == SCIP_VARTYPE_CONTINUOUS);
   assert(prob >= 0);
   assert(nsolvars >= 0);
   assert(solvals != NULL || nsolvars == 0);
   assert(solvars != NULL || nsolvars == 0);

   trivialsol = FALSE;
   /* create data for the new variable in the master problem */
   SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
   newvardata->vartype = GCG_VARTYPE_MASTER;
   newvardata->blocknr = prob;

   /* store whether the variable represents a ray */
   newvardata->data.mastervardata.isray = solisray;

   /* create variable in the master problem */
   SCIP_CALL( SCIPcreateVar(scip, newvar, varname, 0.0, SCIPinfinity(scip), /* GCGrelaxGetNIdenticalBlocks(origprob, prob) */
         objcoeff, vartype, TRUE, TRUE, NULL, NULL, gcgvardeltrans, NULL, newvardata) );

   /* count number of non-zeros */
   newvardata->data.mastervardata.norigvars = 0;

   for( i = 0; i < nsolvars; i++ )
   {
      assert(solvars != NULL);
      assert(solvals != NULL);

      if( !SCIPisZero(scip, solvals[i]) )
      {
         newvardata->data.mastervardata.norigvars++;
      }
   }

   /*
    * if we have not added any original variable to the mastervariable, all coefficients were 0.
    * In that case, we will add all variables in the pricing problem
    */
   if( newvardata->data.mastervardata.norigvars == 0 )
   {
      newvardata->data.mastervardata.norigvars = SCIPgetNOrigVars(pricingscip);
      trivialsol = TRUE;
   }

   if( newvardata->data.mastervardata.norigvars > 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(newvardata->data.mastervardata.origvars), newvardata->data.mastervardata.norigvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(newvardata->data.mastervardata.origvals), newvardata->data.mastervardata.norigvars) );
   }
   else
   {
      newvardata->data.mastervardata.origvars = NULL;
      newvardata->data.mastervardata.origvals = NULL;
   }

   /* number of original variables already saved in mastervardata */
   j = 0;

   /* update variable datas */
   for( i = 0; i < nsolvars && !trivialsol; i++ )
   {
      assert(solvars != NULL);
      assert(solvals != NULL);

      if( !SCIPisZero(scip, solvals[i]) )
      {
         SCIP_VAR* origvar;
         assert(GCGvarIsPricing(solvars[i]));

         origvar = GCGpricingVarGetOrigvars(solvars[i])[0];
         assert(origvar != NULL);

         assert(newvardata->data.mastervardata.origvars != NULL);
         assert(newvardata->data.mastervardata.origvals != NULL);
         assert(GCGvarIsOriginal(origvar));
         /* save in the master problem variable's data the quota of the corresponding original variable */
         newvardata->data.mastervardata.origvars[j] = origvar;
         newvardata->data.mastervardata.origvals[j] = solvals[i];
         /* save the quota in the original variable's data */
         SCIP_CALL( GCGoriginalVarAddMasterVar(scip, origvar, *newvar, solvals[i]) );
         j++;

         /**
          * TODO: TEST
          */

      }
   }
   if( trivialsol )
   {
      SCIP_VAR** pricingvars;
      int npricingvars;

      pricingvars = SCIPgetOrigVars(pricingscip);
      npricingvars = SCIPgetNOrigVars(pricingscip);
      for( j = 0; j < npricingvars; ++j )
      {
         SCIP_VAR* origvar;
         assert(GCGvarIsPricing(pricingvars[j]));

         origvar = GCGpricingVarGetOrigvars(pricingvars[j])[0];
         assert(origvar != NULL);

         assert(newvardata->data.mastervardata.origvars != NULL);
         assert(newvardata->data.mastervardata.origvals != NULL);
         assert(GCGvarIsOriginal(origvar));
         /* save in the master problem variable's data the quota of the corresponding original variable */
         newvardata->data.mastervardata.origvars[j] = origvar;
         newvardata->data.mastervardata.origvals[j] = 0.0;
         /* save the quota in the original variable's data */
         SCIP_CALL( GCGoriginalVarAddMasterVar(scip, origvar, *newvar, 0.0) );
      }
   }
   assert(j == newvardata->data.mastervardata.norigvars);
return SCIP_OKAY;
}

/** creates initial master variables and the vardata */
SCIP_RETCODE GCGcreateInitialMasterVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< original variable */
   SCIP_VAR**            newvar              /**< pointer to store new variable */

   )
{
   SCIP_VARDATA* newvardata;
   int blocknr;

   blocknr = GCGvarGetBlock(var);
   assert( blocknr == -1 || blocknr == -2);

   if( blocknr == -1 )
   {
      SCIPdebugMessage("var %s is in no block - copy it directly to the master\n", SCIPvarGetName(var));
   }
   else
   {
      SCIPdebugMessage("var %s is a linking variable - copy it directly to the master\n", SCIPvarGetName(var));
   }

   /* create vardata */
   SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
   newvardata->vartype = GCG_VARTYPE_MASTER;
   newvardata->blocknr = -1;
   newvardata->data.mastervardata.isray = FALSE;
   newvardata->data.mastervardata.norigvars = 1;

   /* save corresoponding origvar */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(newvardata->data.mastervardata.origvars), 1) ); /*lint !e506*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &(newvardata->data.mastervardata.origvals), 1) ); /*lint !e506*/
   newvardata->data.mastervardata.origvars[0] = var;
   newvardata->data.mastervardata.origvals[0] = 1.0;

   /* create variable in the master problem */
   SCIP_CALL( SCIPcreateVar(scip, newvar, SCIPvarGetName(var),
         SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), SCIPvarGetObj(var), SCIPvarGetType(var),
         TRUE, TRUE, NULL, NULL, gcgvardeltrans, NULL, newvardata) );

   return SCIP_OKAY;
}

/** set creation node of variable */
void GCGsetCreationNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata,            /**< variable data structure */
   long long int         creationNode        /**< node in which the variable is created */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);
   assert(creationNode >= 0);

   vardata->creationnode = creationNode;
}

/** return creation node of variable */
long long int GCGgetCreationNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata             /**< variable data structure */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);

   return vardata->creationnode;
}

/** store creation time */
void GCGsetCreationTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata,            /**< variable data structure */
   SCIP_Real             time                /**< time at which the variable is created */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);
   assert(time >= 0.0);

   vardata->creationtime = time;
}

/** return stored creation time */
SCIP_Real GCGgetCreationTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata             /**< variable data structure */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);

   return vardata->creationtime;
}

/** store iteration */
void GCGsetIteration(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata,            /**< variable data structure */
   SCIP_Longint          iteration           /**< iteration at which the variable is created */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);
   assert(iteration >= 0);

   vardata->iteration = iteration;
}

/** return stored iteration */
SCIP_Longint GCGgetIteration(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata             /**< variable data structure */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);

   return vardata->iteration;

}

/** store gap */
void GCGsetGap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata,            /**< variable data structure */
   SCIP_Real             gap                 /**< present gap when variable is created */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);
   assert(gap >= 0.0);

   vardata->gap = gap;
}

/** return stored gap */
SCIP_Real GCGgetGap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata             /**< variable data structure */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);

   return vardata->gap;
}

/** store reduced cost */
void GCGsetRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata,            /**< variable data structure */
   SCIP_Real             redcost             /**< reduced cost of the variable at creation */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);
   assert(SCIPisLE(scip, redcost, 0.0));
   vardata->redcost = redcost;
}

/** return stored reduced cost */
SCIP_Real GCGgetRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata             /**< variable data structure */
   )
{
   assert(scip != NULL);
   assert(vardata != NULL);

   return vardata->redcost;
}
