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

/**@file   decomp.c
 * @ingroup DECOMP
 * @brief  generic methods for working with different decomposition structures
 * @author Martin Bergner
 *
 * Various methods to work with the decomp structure
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "decomp.h"
#include "pub_decomp.h"
#include "scip/scip.h"
#include "struct_decomp.h"
#include "scip_misc.h"

#include <assert.h>
#include <string.h>

/** fill out subscipvars arrays from the information from vartoblock */
static
SCIP_RETCODE fillOutVarsFromVartoblock(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition structure */
   SCIP_HASHMAP*         vartoblock,         /**< variable to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   SCIP_Bool*            haslinking          /**< returns whether there are linking variables */
   )
{
   SCIP_VAR*** subscipvars;
   int* nsubscipvars;

   SCIP_VAR** linkingvars;
   int nlinkingvars;
   int i;
   SCIP_Bool valid;

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(vartoblock != NULL);
   assert(nblocks >= 0);
   assert(vars != NULL);
   assert(nvars > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );

   nlinkingvars = 0;

   *haslinking = FALSE;

   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars[i], nvars) ); /*lint !e866*/
      nsubscipvars[i] = 0;
   }

   /* handle variables */
   for( i = 0; i < nvars; ++i )
   {
      int block;
      SCIP_VAR* var;

      var = vars[i];
      assert(var != NULL);
      if( !SCIPhashmapExists(vartoblock, var) )
         block = nblocks+1;
      else
      {
         block = (int)(size_t)SCIPhashmapGetImage(vartoblock, var); /*lint !e507*/
      }

      assert(block > 0 && block <= nblocks+1);

      /* if variable belongs to a block */
      if( block <= nblocks )
      {
         SCIPdebugMessage("var %s in block %d.\n", SCIPvarGetName(var), block-1);
         subscipvars[block-1][nsubscipvars[block-1]] = var;
         ++(nsubscipvars[block-1]);
      }
      else /* variable is linking */
      {
         SCIPdebugMessage("var %s is linking.\n", SCIPvarGetName(var));
         assert(block == nblocks+1);
         linkingvars[nlinkingvars] = var;
         ++nlinkingvars;
      }
   }

   if( nlinkingvars > 0 )
   {
      SCIP_CALL( DECdecompSetLinkingvars(scip, decdecomp, linkingvars, nlinkingvars, &valid) );
      assert(valid);
      *haslinking = TRUE;
   }

   for( i = 0; i < nblocks; ++i )
   {
      if( nsubscipvars[i] == 0 )
      {
         SCIPfreeBufferArray(scip, &subscipvars[i]);
         subscipvars[i] = NULL;
      }
   }
   if( nblocks > 0 )
   {
      SCIP_CALL( DECdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars, &valid) );
      assert(valid);
   }
   DECdecompSetVartoblock(decdecomp, vartoblock, &valid);
   assert(valid);
   SCIPfreeBufferArray(scip, &nsubscipvars);

   for( i = 0; i < nblocks; ++i )
   {
     SCIPfreeBufferArrayNull(scip, &subscipvars[i]);
   }

   SCIPfreeBufferArray(scip, &subscipvars);
   SCIPfreeBufferArray(scip, &linkingvars);

   return SCIP_OKAY;
}


/** fill out subscipcons arrays from the information from constoblock */
static
SCIP_RETCODE fillOutConsFromConstoblock(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition structure */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_CONS**           conss,              /**< constraint array */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            haslinking          /**< returns whether there are linking constraints */
   )
{
   SCIP_CONS*** subscipconss;
   int* nsubscipconss;

   SCIP_CONS** linkingconss;
   int nlinkingconss;
   int i;
   SCIP_Bool valid;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(constoblock != NULL);
   assert(nblocks >= 0);
   assert(conss != NULL);
   assert(nconss > 0);

   DECdecompSetConstoblock(decdecomp, constoblock, &valid);
   assert(valid);

   SCIP_CALL( SCIPallocMemoryArray(scip, &linkingconss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );

   *haslinking = FALSE;

   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &subscipconss[i], nconss) ); /*lint !e866*/
      nsubscipconss[i] = 0;
   }

   nlinkingconss = 0;

   /* handle constraints */
   for( i = 0; i < nconss; ++i )
   {
      int block;
      SCIP_CONS* cons;

      cons = conss[i];
      assert(cons != NULL);
      if( !SCIPhashmapExists(decdecomp->constoblock, cons) )
      {
         block = nblocks+1;
         SCIP_CALL( SCIPhashmapInsert(decdecomp->constoblock, cons, (void*) (size_t) block) );
      }
      else
      {
         block = (int)(size_t)SCIPhashmapGetImage(decdecomp->constoblock, cons); /*lint !e507*/
      }

      assert(block > 0 && block <= nblocks+1);

      /* if constraint belongs to a block */
      if( block <= nblocks )
      {
         SCIPdebugMessage("cons %s in block %d.\n", SCIPconsGetName(cons), block-1);
         subscipconss[block-1][nsubscipconss[block-1]] = cons;
         ++(nsubscipconss[block-1]);
      }
      else /* constraint is linking */
      {
         SCIPdebugMessage("cons %s is linking.\n", SCIPconsGetName(cons));
         assert(block == nblocks+1);
         linkingconss[nlinkingconss] = cons;
         ++nlinkingconss;
      }
   }

   if( nlinkingconss > 0 )
   {
      SCIP_CALL( DECdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss, &valid) );
      assert(valid);
      *haslinking = TRUE;
   }
   if( nblocks > 0 )
   {
      SCIP_CALL( DECdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss, &valid) );
      assert(valid);
   }
   SCIPfreeMemoryArray(scip, &linkingconss);
   SCIPfreeBufferArray(scip, &nsubscipconss);

   for( i = 0; i < nblocks; ++i )
   {
     SCIPfreeMemoryArray(scip, &subscipconss[i]);
   }

   SCIPfreeBufferArray(scip, &subscipconss);

   return SCIP_OKAY;
}


const char *DECgetStrType(
   DEC_DECTYPE type
   )
{
   const char * names[] = { "unknown", "arrowhead", "staircase", "diagonal", "bordered" };
   return names[type];
}

/** initializes the decdecomp structure to absolutely nothing */
SCIP_RETCODE DECdecompCreate(
   SCIP*                 scip,               /**< Pointer to the SCIP instance */
   DEC_DECOMP**          decomp              /**< Pointer to the decdecomp instance */
   )
{
   assert(scip != NULL);
   assert(decomp != NULL);
   SCIP_CALL( SCIPallocMemory(scip, decomp) );

   (*decomp)->type = DEC_DECTYPE_UNKNOWN;
   (*decomp)->constoblock = NULL;
   (*decomp)->vartoblock = NULL;
   (*decomp)->subscipvars = NULL;
   (*decomp)->subscipconss = NULL;
   (*decomp)->nsubscipconss = NULL;
   (*decomp)->nsubscipvars = NULL;
   (*decomp)->linkingconss = NULL;
   (*decomp)->nlinkingconss = 0;
   (*decomp)->linkingvars = NULL;
   (*decomp)->nlinkingvars = 0;
   (*decomp)->stairlinkingvars = NULL;
   (*decomp)->nstairlinkingvars = NULL;
   (*decomp)->nblocks = 0;
   (*decomp)->consindex = NULL;
   (*decomp)->varindex = NULL;

   return SCIP_OKAY;
}

/** frees the decdecomp structure */
SCIP_RETCODE DECdecompFree(
   SCIP*                 scip,               /**< pointer to the SCIP instance */
   DEC_DECOMP**          decdecomp           /**< pointer to the decdecomp instance */
   )
{
   DEC_DECOMP* decomp;
   int i;
   int j;

   assert( scip!= NULL );
   assert( decdecomp != NULL);
   decomp = *decdecomp;

   assert(decomp != NULL);

   for( i = 0; i < decomp->nblocks; ++i )
   {
      for( j = 0; j < decomp->nsubscipvars[i]; ++j )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(decomp->subscipvars[i][j])) );
      }
      SCIPfreeMemoryArrayNull(scip, &(decomp->subscipvars[i]));

      for( j = 0; j < decomp->nsubscipconss[i]; ++j )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(decomp->subscipconss[i][j])) );
      }
      SCIPfreeMemoryArray(scip, &decomp->subscipconss[i]);
   }

   for( i = 0; i < decomp->nlinkingvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(decomp->linkingvars[i])) );
   }

   if( decomp->stairlinkingvars != NULL )
      for( i = 0; i < decomp->nblocks-1; ++i )
      {
         for( j = 0; j < decomp->nstairlinkingvars[i]; ++j )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &(decomp->stairlinkingvars[i][j])) );
         }
         SCIPfreeMemoryArray(scip, &decomp->stairlinkingvars[i]);
      }

   /* free hashmaps if they are not NULL */
   if( decomp->constoblock != NULL )
      SCIPhashmapFree(&decomp->constoblock);
   if( decomp->vartoblock != NULL )
      SCIPhashmapFree(&decomp->vartoblock);
   if( decomp->varindex != NULL )
      SCIPhashmapFree(&decomp->varindex);
   if( decomp->consindex != NULL )
      SCIPhashmapFree(&decomp->consindex);

   for( i = 0; i < decomp->nlinkingconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(decomp->linkingconss[i])) );
   }
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->linkingvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->stairlinkingvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->nstairlinkingvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->linkingconss);
   SCIPfreeMemory(scip, decdecomp);

   return SCIP_OKAY;
}

/** sets the type of the decomposition */
void DECdecompSetType(
   DEC_DECOMP*           decdecomp,          /**< pointer to the decdecomp instance */
   DEC_DECTYPE           type,               /**< type of the decomposition */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   assert(decdecomp != NULL);
   switch( type )
   {
   case DEC_DECTYPE_DIAGONAL:
      *valid = decdecomp->nlinkingconss == 0 && decdecomp->linkingconss == NULL;
      *valid = *valid && decdecomp->nlinkingvars == 0 && decdecomp->linkingvars == NULL;
      break;
   case DEC_DECTYPE_ARROWHEAD:
      *valid = TRUE;
      break;
   case DEC_DECTYPE_UNKNOWN:
      *valid = FALSE;
      break;
   case DEC_DECTYPE_BORDERED:
      *valid = decdecomp->nlinkingvars == 0 && decdecomp->linkingvars == NULL;
      break;
   case DEC_DECTYPE_STAIRCASE:
      *valid = decdecomp->nlinkingconss == 0 && decdecomp->linkingconss == NULL;
      break;
   default:
      *valid = FALSE;
      break;
   }

   decdecomp->type = type;
}

/** gets the type of the decomposition */
DEC_DECTYPE DECdecompGetType(
   DEC_DECOMP*           decdecomp           /**< Pointer to the decdecomp instance */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->type;
}


/** sets the presolved flag for decomposition */
void DECdecompSetPresolved(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_Bool             presolved           /**< presolved flag for decomposition */
   )
{
   assert(decdecomp != NULL);

   decdecomp->presolved = presolved;
}

/** gets the presolved flag for decomposition */
SCIP_Bool DECdecompGetPresolved(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->presolved;
}

/** sets the number of blocks for decomposition */
void DECdecompSetNBlocks(
   DEC_DECOMP*           decdecomp,          /**< Pointer to the decdecomp instance */
   int                   nblocks             /**< number of blocks for decomposition */
   )
{
   assert(decdecomp != NULL);
   assert(nblocks >= 0);

   decdecomp->nblocks = nblocks;
}

/** gets the number of blocks for decomposition */
int DECdecompGetNBlocks(
   DEC_DECOMP*           decdecomp           /**< Pointer to the decdecomp instance */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->nblocks;
}

/** copies the input subscipvars array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetSubscipvars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_VAR***           subscipvars,        /**< Subscipvars array  */
   int*                  nsubscipvars,       /**< number of subscipvars per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   int i;
   int b;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(subscipvars != NULL);
   assert(nsubscipvars != NULL);
   assert(decdecomp->nblocks > 0);

   assert(decdecomp->subscipvars == NULL);
   assert(decdecomp->nsubscipvars == NULL);

   *valid = TRUE;

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipvars, decdecomp->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nsubscipvars, decdecomp->nblocks) );

   assert(decdecomp->subscipvars != NULL);
   assert(decdecomp->nsubscipvars != NULL);

   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      assert((subscipvars[b] == NULL) == (nsubscipvars[b] == 0));
      decdecomp->nsubscipvars[b] = nsubscipvars[b];

      if( nsubscipvars[b] < 0 )
         *valid = FALSE;
      else if( nsubscipvars[b] > 0 )
      {
         assert(subscipvars[b] != NULL);
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &(decdecomp->subscipvars[b]), subscipvars[b], nsubscipvars[b]) ); /*lint !e866*/

         for( i = 0; i < nsubscipvars[b]; ++i )
         {
            SCIP_CALL( SCIPcaptureVar(scip, decdecomp->subscipvars[b][i]) );
         }
      }
      else if( nsubscipvars[b] == 0 )
      {
         decdecomp->subscipvars[b] = NULL;
      }
   }

   return SCIP_OKAY;
}

/** returns the subscipvars array of the given decdecomp structure */
SCIP_VAR***  DECdecompGetSubscipvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->subscipvars;
}

/** returns the nsubscipvars array of the given decdecomp structure */
int*  DECdecompGetNSubscipvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->nsubscipvars;
}

/** copies the input subscipconss array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetSubscipconss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_CONS***          subscipconss,       /**< Subscipconss array  */
   int*                  nsubscipconss,      /**< number of subscipconss per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   int i;
   int b;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(subscipconss != NULL);
   assert(nsubscipconss != NULL);
   assert(valid != NULL);

   assert(decdecomp->nblocks > 0);
   assert(decdecomp->subscipconss == NULL);
   assert(decdecomp->nsubscipconss == NULL);

   *valid = TRUE;

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->subscipconss, decdecomp->nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nsubscipconss, decdecomp->nblocks) );

   assert(decdecomp->subscipconss != NULL);
   assert(decdecomp->nsubscipconss != NULL);

   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      if( nsubscipconss[b] <= 0 || subscipconss[b] == NULL )
         *valid = FALSE;

      decdecomp->nsubscipconss[b] = nsubscipconss[b];

      if( nsubscipconss[b] > 0 )
      {
         assert(subscipconss[b] != NULL);
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->subscipconss[b], subscipconss[b], nsubscipconss[b]) ); /*lint !e866*/
         for( i = 0; i < nsubscipconss[b]; ++i )
         {
            SCIP_CALL( SCIPcaptureCons(scip, decdecomp->subscipconss[b][i]) );
         }
      }
   }

   return SCIP_OKAY;
}

/** returns the subscipconss array of the given decdecomp structure */
SCIP_CONS***  DECdecompGetSubscipconss(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->subscipconss;
}

/** returns the nsubscipconss array of the given decdecomp structure */
int*  DECdecompGetNSubscipconss(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->nsubscipconss;
}

/** copies the input linkingconss array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetLinkingconss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_CONS**           linkingconss,       /**< Linkingconss array  */
   int                   nlinkingconss,      /**< number of linkingconss per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   int i;

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(linkingconss != NULL);
   assert(nlinkingconss >= 0);

   assert(decdecomp->linkingconss == NULL);
   assert(decdecomp->nlinkingconss == 0);

   decdecomp->nlinkingconss = nlinkingconss;

   if( nlinkingconss > 0 )
   {
      assert(linkingconss != NULL);

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->linkingconss, linkingconss, nlinkingconss) );

      for( i = 0; i < nlinkingconss; ++i )
      {
         SCIP_CALL( SCIPcaptureCons(scip, decdecomp->linkingconss[i]) );
      }
   }

   *valid = linkingconss != NULL || nlinkingconss == 0;

   return SCIP_OKAY;
}

/** returns the linkingconss array of the given decdecomp structure */
SCIP_CONS**  DECdecompGetLinkingconss(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->linkingconss;
}

/** returns the nlinkingconss array of the given decdecomp structure */
int  DECdecompGetNLinkingconss(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   assert(decdecomp->nlinkingconss >= 0);

   return decdecomp->nlinkingconss;
}

/** copies the input linkingvars array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetLinkingvars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_VAR**            linkingvars,        /**< Linkingvars array  */
   int                   nlinkingvars,       /**< number of linkingvars per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   int i;

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(linkingvars != NULL || nlinkingvars == 0);

   assert(decdecomp->linkingvars == NULL);
   assert(decdecomp->nlinkingvars == 0);

   decdecomp->nlinkingvars = nlinkingvars;

   if( nlinkingvars > 0 )
   {
      assert(linkingvars != NULL);

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &decdecomp->linkingvars, linkingvars, nlinkingvars) );

      for( i = 0; i < nlinkingvars; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, decdecomp->linkingvars[i]) );
      }
   }

   *valid = linkingvars != NULL || nlinkingvars == 0;

   return SCIP_OKAY;
}

/** returns the linkingvars array of the given decdecomp structure */
SCIP_VAR**  DECdecompGetLinkingvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->linkingvars;
}

/** returns the nlinkingvars array of the given decdecomp structure */
int  DECdecompGetNLinkingvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   assert(decdecomp->nlinkingvars >= 0);

   return decdecomp->nlinkingvars;
}

/** copies the input stairlinkingvars array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetStairlinkingvars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_VAR***           stairlinkingvars,   /**< Linkingvars array  */
   int*                  nstairlinkingvars,  /**< number of linkingvars per block */
   SCIP_Bool*            valid               /**< returns whether the resulting decdecomp is valid */
   )
{
   int b;
   int i;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(stairlinkingvars != NULL);
   assert(nstairlinkingvars != NULL);
   assert(valid != NULL);
   assert(decdecomp->nblocks > 0);

   assert(decdecomp->stairlinkingvars == NULL);
   assert(decdecomp->nstairlinkingvars == NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->stairlinkingvars, decdecomp->nblocks-1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &decdecomp->nstairlinkingvars, decdecomp->nblocks-1) );

   assert(decdecomp->stairlinkingvars != NULL);
   assert(decdecomp->nstairlinkingvars != NULL);

   for( b = 0; b < decdecomp->nblocks-1; ++b )
   {
      assert(nstairlinkingvars[b] > 0);
      decdecomp->nstairlinkingvars[b] = nstairlinkingvars[b];

      assert(stairlinkingvars[b] != NULL);
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(decdecomp->stairlinkingvars[b]), stairlinkingvars[b], nstairlinkingvars[b]) ); /*lint !e866 */
   }

   for( b = 0; b < decdecomp->nblocks-1; ++b )
   {
      for( i = 0; i < nstairlinkingvars[b]; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, decdecomp->stairlinkingvars[b][i]) );
      }
   }

   *valid = TRUE; /**@todo A valid check needs to be implemented */
   return SCIP_OKAY;
}

/** returns the stairlinkingvars array of the given decdecomp structure */
SCIP_VAR***  DECdecompGetStairlinkingvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->stairlinkingvars;
}

/** returns the nstairlinkingvars array of the given decdecomp structure */
int*  DECdecompGetNStairlinkingvars(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   assert(decdecomp->nstairlinkingvars != NULL );
   return decdecomp->nstairlinkingvars;
}

/** sets the vartoblock hashmap of the given decdecomp structure */
void  DECdecompSetVartoblock(
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_HASHMAP*         vartoblock,         /**< Vartoblock hashmap */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   assert(decdecomp != NULL);
   assert(vartoblock != NULL);

   *valid = TRUE;

   decdecomp->vartoblock = vartoblock;
}

/** returns the vartoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecompGetVartoblock(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->vartoblock;
}

/** sets the constoblock hashmap of the given decdecomp structure */
void  DECdecompSetConstoblock(
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_HASHMAP*         constoblock,        /**< Constoblock hashmap */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   )
{
   assert(decdecomp != NULL);
   assert(constoblock != NULL);

   *valid = TRUE;

   decdecomp->constoblock = constoblock;
}

/** returns the constoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecompGetConstoblock(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->constoblock;
}

/** sets the varindex hashmap of the given decdecomp structure */
void  DECdecompSetVarindex(
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_HASHMAP*         varindex            /**< Varindex hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(varindex != NULL);
   decdecomp->varindex = varindex;
}

/** returns the varindex hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecompGetVarindex(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->varindex;
}

/** sets the consindex hashmap of the given decdecomp structure */
void  DECdecompSetConsindex(
   DEC_DECOMP*           decdecomp,          /**< DEC_DECOMP data structure */
   SCIP_HASHMAP*         consindex           /**< Consindex hashmap */
   )
{
   assert(decdecomp != NULL);
   assert(consindex != NULL);
   decdecomp->consindex = consindex;
}

/** returns the consindex hashmap of the given decdecomp structure */
SCIP_HASHMAP*  DECdecompGetConsindex(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);
   return decdecomp->consindex;
}

/** completely initializes decdecomp from the values of the hashmaps */
SCIP_RETCODE DECfillOutDecdecompFromHashmaps(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition structure */
   SCIP_HASHMAP*         vartoblock,         /**< variable to block hashmap */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< constraint array */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            valid,              /**< pointer to indicate whether the structure is valid */
   SCIP_Bool             staircase           /**< should the decomposition be a staircase structure */
   )
{
   SCIP_HASHMAP* varindex;
   SCIP_HASHMAP* consindex;
   int* nsubscipconss;
   int* nsubscipvars;
   int* nstairlinkingvars;
   SCIP_VAR*** stairlinkingvars;
   SCIP_CONS*** subscipconss;
   SCIP_Bool success;
   int idx;
   int linkindex;
   int cindex;
   int cumindex;
   SCIP_Bool haslinking;
   int i;
   int b;
   SCIP_VAR** curvars;
   int ncurvars;
   int j;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(vartoblock != NULL);
   assert(constoblock != NULL);
   assert(nblocks >= 0);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(conss != NULL);
   assert(nconss > 0);
   assert(valid != NULL);

   DECdecompSetNBlocks(decdecomp, nblocks);
   *valid = TRUE;

   DECdecompSetType(decdecomp, DEC_DECTYPE_DIAGONAL, valid);
   SCIP_CALL( fillOutConsFromConstoblock(scip, decdecomp, constoblock, nblocks, conss, nconss, &haslinking) );

   if( haslinking )
   {
      SCIPdebugMessage("Decomposition has linking constraints and is bordered.\n");
      DECdecompSetType(decdecomp, DEC_DECTYPE_BORDERED, valid);
      assert(*valid);
   }

   SCIP_CALL( fillOutVarsFromVartoblock(scip,  decdecomp, vartoblock, nblocks, vars, nvars, &haslinking) );

   if( haslinking )
   {
      SCIPdebugMessage("Decomposition has linking variables and is arrowhead.\n");
      DECdecompSetType(decdecomp, DEC_DECTYPE_ARROWHEAD, valid);
      assert(*valid);
   }

   if( !staircase )
   {
      SCIP_CALL( DECdecompCheckConsistency(scip, decdecomp) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPhashmapCreate(&varindex, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&consindex, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &stairlinkingvars, nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nstairlinkingvars, nblocks) );

   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(stairlinkingvars[i]), nvars) ); /*lint !e866*/
      nstairlinkingvars[i] = 0;
   }

   nsubscipconss = DECdecompGetNSubscipconss(decdecomp);
   subscipconss = DECdecompGetSubscipconss(decdecomp);
   nsubscipvars = DECdecompGetNSubscipvars(decdecomp);

   idx = 0;
   cindex = 0;
   cumindex = 0;

   /* try to deduce staircase map */
   for( b = 0; b < nblocks; ++b )
   {
      cumindex += nsubscipvars[b];
      SCIPdebugMessage("block %d (%d vars):\n", b, nsubscipvars[b]);
      linkindex = 0;
      for( i = 0; i < nsubscipconss[b]; ++i )
      {
         SCIP_CONS* cons;
         cons = subscipconss[b][i];

         SCIP_CALL( SCIPhashmapInsert(consindex, cons, (void*)(size_t)(cindex+1)) );
         ++cindex;
         SCIP_CALL( SCIPgetConsNVars(scip, cons, &ncurvars, &success) );
         assert(success);

         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );

         SCIP_CALL( SCIPgetConsVars(scip, cons, curvars, ncurvars, &success) );
         assert(success);

         for( j = 0; j < ncurvars; ++j )
         {
            SCIP_VAR* probvar = SCIPvarGetProbvar(curvars[j]);

            /* if the variable is linking */
            if( (int)(size_t)SCIPhashmapGetImage(vartoblock, probvar) == nblocks+1 ) /*lint !e507*/
            {
               /* if it has not been already assigned, it links to the next block */
               if( !SCIPhashmapExists(varindex, probvar) )
               {
                  SCIPdebugMessage("assigning link var <%s> to index <%d>\n", SCIPvarGetName(probvar), cumindex+linkindex+1);
                  SCIP_CALL( SCIPhashmapInsert(varindex, probvar, (void*)(size_t)(cumindex+linkindex+1)) );
                  stairlinkingvars[b][nstairlinkingvars[b]] = probvar;
                  ++(nstairlinkingvars[b]);
                  linkindex++;
               }
            }
            else
            {
               assert(((int) (size_t) SCIPhashmapGetImage(vartoblock, probvar)) -1 == b);  /*lint !e507*/
               SCIP_CALL( SCIPhashmapInsert(varindex, probvar, (void*)(size_t)(idx+1)) );
               ++idx;
            }
         }
         SCIPfreeBufferArray(scip, &curvars);
      }
      idx += linkindex;
      cumindex += linkindex;
   }
   DECdecompSetVarindex(decdecomp, varindex);
   DECdecompSetConsindex(decdecomp, consindex);
   DECdecompSetType(decdecomp, DEC_DECTYPE_STAIRCASE, valid);
   assert(*valid);

   for( b = 0; b < nblocks; ++b )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(stairlinkingvars[b]), nstairlinkingvars[b]) ); /*lint !e866*/
   }

   SCIP_CALL( DECdecompSetStairlinkingvars(scip, decdecomp, stairlinkingvars, nstairlinkingvars, valid) );
   assert(*valid);

   for( b = 0; b < nblocks; ++b )
   {
      SCIPfreeMemoryArray(scip, &stairlinkingvars[b]);
   }
   SCIPfreeMemoryArray(scip, &stairlinkingvars);
   SCIPfreeMemoryArray(scip, &nstairlinkingvars);

   SCIP_CALL( DECdecompCheckConsistency(scip, decdecomp) );

   return SCIP_OKAY;
}

/** completely fills out detector structure from only the constraint partition */
SCIP_RETCODE DECfilloutDecdecompFromConstoblock(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decomposition structure */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< constraint array */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             staircase           /**< should the decomposition be a staircase structure */
   )
{
   SCIP_HASHMAP* vartoblock;
   int i;
   int j;

   SCIP_VAR** curvars;
   int ncurvars;
   SCIP_Bool valid;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(decdecomp != NULL);
   assert(constoblock != NULL);
   assert(nblocks >= 0);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(conss != NULL);

   assert(nconss > 0);

   SCIP_CALL( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip), nvars) );
   for( i = 0; i < nconss; ++i )
   {
      int consblock;

      consblock = (int)(size_t)SCIPhashmapGetImage(constoblock, conss[i]);  /*lint !e507*/

      assert(consblock > 0 && consblock < nblocks+2);

      SCIP_CALL( SCIPgetConsNVars(scip, conss[i], &ncurvars, &success) );
      assert(success);

      SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );

      SCIP_CALL( SCIPgetConsVars(scip, conss[i], curvars, ncurvars, &success) );
      assert(success);

      for( j = 0; j < ncurvars; ++j )
      {
         SCIP_VAR* probvar = SCIPvarGetProbvar(curvars[j]);
         assert( SCIPvarIsActive(probvar) );
         if( !SCIPhashmapExists(vartoblock, probvar) && consblock <= nblocks )
         {
            SCIP_CALL( SCIPhashmapSetImage(vartoblock, probvar, (void*) (size_t) consblock) );
         }
         else if( consblock <= nblocks )
         {
            SCIP_CALL( SCIPhashmapSetImage(vartoblock, probvar, (void*) (size_t) (nblocks+1)) );
         }

         DECdecompSetVartoblock(decdecomp, vartoblock, &valid);
         assert(valid);
         DECdecompSetConstoblock(decdecomp, constoblock, &valid);
         assert(valid);
      }

      SCIPfreeBufferArray(scip, &curvars);
   }

   for( i = 0; i < nvars; ++i )
   {
      if( !SCIPhashmapExists(vartoblock, vars[i]) )
      {
         SCIP_CALL( SCIPhashmapSetImage(vartoblock, vars[i], (void*) (size_t) (nblocks+1)) );
      }
   }

   SCIP_CALL( DECfillOutDecdecompFromHashmaps(scip, decdecomp, vartoblock, constoblock, nblocks, vars, nvars, conss, nconss, &valid, staircase) );
   assert(valid);

   return SCIP_OKAY;
}

/** sets the detector for the given decdecomp structure */
void DECdecompSetDetector(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   DEC_DETECTOR*         detector            /**< detector data structure */
   )
{
   assert(decdecomp != NULL);

   decdecomp->detector = detector;
}

/** gets the detector for the given decdecomp structure */
DEC_DETECTOR* DECdecompGetDetector(
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
   )
{
   assert(decdecomp != NULL);

   return decdecomp->detector;
}

/** transforms all constraints and variables, updating the arrays */
SCIP_RETCODE DECdecompTransform(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   )
{
   int b;
   int c;
   int v;
   SCIP_HASHMAP* newconstoblock;
   SCIP_HASHMAP* newvartoblock;
   SCIP_VAR* newvar;

   assert(SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED);

   SCIP_CALL( SCIPhashmapCreate(&newconstoblock, SCIPblkmem(scip), SCIPgetNConss(scip)) );
   SCIP_CALL( SCIPhashmapCreate(&newvartoblock, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   /* transform all constraints and put them into constoblock */
   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      for( c = 0; c < decdecomp->nsubscipconss[b]; ++c )
      {
         SCIP_CONS* newcons;
         SCIPdebugMessage("%d, %d: %s (%s)\n", b, c, SCIPconsGetName(decdecomp->subscipconss[b][c]), SCIPconsIsTransformed(decdecomp->subscipconss[b][c])?"t":"o" );
         assert(decdecomp->subscipconss[b][c] != NULL);
         newcons = SCIPfindCons(scip, SCIPconsGetName(decdecomp->subscipconss[b][c]));
         if( newcons != decdecomp->subscipconss[b][c] )
         {
            SCIP_CALL( SCIPcaptureCons(scip, newcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &(decdecomp->subscipconss[b][c])) );
            decdecomp->subscipconss[b][c] = newcons;
         }
         assert(decdecomp->subscipconss[b][c] != NULL);
         assert(!SCIPhashmapExists(newconstoblock, decdecomp->subscipconss[b][c]));
         SCIP_CALL( SCIPhashmapSetImage(newconstoblock, decdecomp->subscipconss[b][c], (void*) (size_t) (b+1)) );
      }
   }
   /* transform all variables and put them into vartoblock */
   for( b = 0; b < decdecomp->nblocks; ++b )
   {
      int idx;
      for( v = 0, idx = 0; v < decdecomp->nsubscipvars[b]; ++v )
      {
         assert(decdecomp->subscipvars[b][v] != NULL);

         SCIPdebugMessage("%d, %d: %s (%p, %s)\n", b, v, SCIPvarGetName(decdecomp->subscipvars[b][v]),
            decdecomp->subscipvars[b][v], SCIPvarIsTransformed(decdecomp->subscipvars[b][v])?"t":"o" );

         /* make sure that newvar is a transformed variable */
         SCIP_CALL( SCIPgetTransformedVar(scip, decdecomp->subscipvars[b][v], &newvar) );
         SCIP_CALL( SCIPreleaseVar(scip, &(decdecomp->subscipvars[b][v])) );

         assert(newvar != NULL);
         assert(SCIPvarIsTransformed(newvar));

         newvar = SCIPvarGetProbvar(newvar);
         assert(newvar != NULL);

         /* the probvar can also be fixed, in which case we do not need it in the block; furthermore, multiple variables
          * can resolve to the same active problem variable, so we check whether we already handled the variable
          */
         if( SCIPvarIsActive(newvar) && !SCIPhashmapExists(newvartoblock, newvar) )
         {
            decdecomp->subscipvars[b][idx] = newvar;
            SCIP_CALL( SCIPcaptureVar(scip, newvar) );
            SCIPdebugMessage("%d, %d: %s (%p, %s)\n", b, v, SCIPvarGetName(decdecomp->subscipvars[b][idx]),
               decdecomp->subscipvars[b][idx], SCIPvarIsTransformed(decdecomp->subscipvars[b][idx])?"t":"o" );

            assert(decdecomp->subscipvars[b][idx] != NULL);
            assert(!SCIPhashmapExists(newvartoblock, decdecomp->subscipvars[b][idx]));
            SCIP_CALL( SCIPhashmapSetImage(newvartoblock, decdecomp->subscipvars[b][idx], (void*) (size_t) (b+1)) );
            ++idx;
         }
      }
      decdecomp->nsubscipvars[b] = idx;
   }

   /* transform all linking constraints */
   for( c = 0; c < decdecomp->nlinkingconss; ++c )
   {
      SCIP_CONS* newcons;

      SCIPdebugMessage("m, %d: %s (%s)\n", c, SCIPconsGetName(decdecomp->linkingconss[c]), SCIPconsIsTransformed(decdecomp->linkingconss[c])?"t":"o" );
      assert(decdecomp->linkingconss[c] != NULL);
      newcons = SCIPfindCons(scip, SCIPconsGetName(decdecomp->linkingconss[c]));
      if( newcons != decdecomp->linkingconss[c] )
      {
         SCIP_CALL( SCIPcaptureCons(scip, newcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &(decdecomp->linkingconss[c])) );
         decdecomp->linkingconss[c] = newcons;
      }
      SCIP_CALL( SCIPhashmapSetImage(newconstoblock, decdecomp->linkingconss[c],(void*) (size_t) (decdecomp->nblocks+1) ) );

      assert(decdecomp->linkingconss[c] != NULL);
   }

   /* transform all linking variables */
   for( v = 0; v < decdecomp->nlinkingvars; ++v )
   {
      SCIPdebugMessage("m, %d: %s (%p, %s)\n", v, SCIPvarGetName(decdecomp->linkingvars[v]), decdecomp->linkingvars[v], SCIPvarIsTransformed(decdecomp->linkingvars[v])?"t":"o" );
      assert(decdecomp->linkingvars[v] != NULL);

      if( !SCIPvarIsTransformed(decdecomp->linkingvars[v]) )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, decdecomp->linkingvars[v], &newvar) );
         newvar = SCIPvarGetProbvar(newvar);
      }
      else
         newvar = decdecomp->linkingvars[v];
      assert(newvar != NULL);
      assert(SCIPvarIsTransformed(newvar));

      decdecomp->linkingvars[v] = newvar;
      SCIP_CALL( SCIPhashmapSetImage(newvartoblock, decdecomp->linkingvars[v], (void*) (size_t) (decdecomp->nblocks+1) ) );
      SCIPdebugMessage("m, %d: %s (%p, %s)\n", v, SCIPvarGetName(decdecomp->linkingvars[v]), decdecomp->linkingvars[v], SCIPvarIsTransformed(decdecomp->linkingvars[v])?"t":"o" );
      assert(decdecomp->linkingvars[v] != NULL);
   }

   SCIPhashmapFree(&decdecomp->constoblock);
   decdecomp->constoblock = newconstoblock;
   SCIPhashmapFree(&decdecomp->vartoblock);
   decdecomp->vartoblock = newvartoblock;

   SCIP_CALL( DECdecompCheckConsistency(scip, decdecomp) );

   return SCIP_OKAY;
}

/** prints out detailed information on the contents of decdecomp*/
void DECdecompPrintDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   )
{
   int i;
   int j;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   SCIPinfoMessage(scip, NULL, "================DEC_DECOMP===============\n");
   SCIPinfoMessage(scip, NULL, "# blocks: %i\n", decdecomp->nblocks);
   for( i = 0; i < decdecomp->nblocks; ++i )
   {
      SCIPinfoMessage(scip, NULL, "Block #%i (#vars: %i, #conss: %i):\n", i+1, decdecomp->nsubscipvars[i], decdecomp->nsubscipconss[i]);
      SCIPinfoMessage(scip, NULL, "Variables (block, index):\n");
      for( j = 0; j < decdecomp->nsubscipvars[i]; ++j )
      {
         var = decdecomp->subscipvars[i][j];
         SCIPinfoMessage(scip, NULL, "\t%s (%i, %i)\n", SCIPvarGetName(var), *(int*) SCIPhashmapGetImage(decdecomp->vartoblock, (void*) var), *(int*) SCIPhashmapGetImage(decdecomp->varindex, (void*) var));
      }
      SCIPinfoMessage(scip, NULL, "Constraints:\n");
      for( j = 0; j < decdecomp->nsubscipconss[i]; ++j )
      {
         cons = decdecomp->subscipconss[i][j];
         SCIPinfoMessage(scip, NULL, "\t%s (%i, %i)\n", SCIPconsGetName(cons), *(int*) SCIPhashmapGetImage(decdecomp->constoblock, (void*) cons), *(int*) SCIPhashmapGetImage(decdecomp->consindex, (void*) cons));
      }
      SCIPinfoMessage(scip, NULL, "========================================\n");
   }
   SCIPinfoMessage(scip, NULL, "Linking variables #%i (varindex) :\n", decdecomp->nlinkingvars);
   for( j = 0; j < decdecomp->nlinkingvars; ++j )
   {
      var = decdecomp->linkingvars[j];
      SCIPinfoMessage(scip, NULL, "\t%s (%i)\n", SCIPvarGetName(var), *(int*) SCIPhashmapGetImage(decdecomp->varindex, (void*) var));
   }
   SCIPinfoMessage(scip, NULL, "========================================\n");
   SCIPinfoMessage(scip, NULL, "Linking constraints #%i (consindex) :\n", decdecomp->nlinkingconss);
   for( j = 0; j < decdecomp->nlinkingconss; ++j )
   {
      cons = decdecomp->linkingconss[j];
      SCIPinfoMessage(scip, NULL, "\t%s (%i)\n", SCIPconsGetName(cons), *(int*) SCIPhashmapGetImage(decdecomp->consindex, (void*) cons));
   }
   SCIPinfoMessage(scip, NULL, "========================================\n");
}

/** checks the consistency of the data structure
 *
 *  In particular, it checks whether the redundant information in the structure agree and
 *  whether the variables in the structure are both existant in the arrays and in the problem
 */
SCIP_RETCODE DECdecompCheckConsistency(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decomposition data structure */
   )
{
#ifndef NDEBUG
   // SCIP_Bool* varishandled;
   // SCIP_Bool* consishandled;

   int c;
   int b;
   int v;

   //SCIP_CALL( SCIPallocMemoryArray(scip, &varishandled, SCIPgetNVars(scip)) );
   //SCIP_CALL( SCIPallocMemoryArray(scip, &consishandled, SCIPgetNConss(scip)) );
   //
   //BMSclearMemoryArray(varishandled, SCIPgetNVars(scip));
   //BMSclearMemoryArray(consishandled, SCIPgetNConss(scip));

   SCIPdebugMessage("Problem is %stransformed\n", SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED ? "": "not ");

   for( v = 0; v < SCIPgetNVars(scip); ++v )
   {
      assert(SCIPhashmapExists(DECdecompGetVartoblock(decdecomp), SCIPgetVars(scip)[v]));
   }

   for( c = 0; c < SCIPgetNConss(scip); ++c )
   {
      if( !GCGisConsGCGCons(SCIPgetConss(scip)[c]) )
      {
         assert(SCIPhashmapExists(DECdecompGetConstoblock(decdecomp), SCIPgetConss(scip)[c]));
      }
   }

   /* Check whether subscipcons are correct */
   for( b = 0; b < DECdecompGetNBlocks(decdecomp); ++b )
   {
      for( c = 0; c < DECdecompGetNSubscipconss(decdecomp)[b]; ++c )
      {
         SCIP_VAR** curvars;
         int ncurvars;
         SCIP_CONS* cons = DECdecompGetSubscipconss(decdecomp)[b][c];
         SCIPdebugMessage("Cons <%s> in block %d = %d\n", SCIPconsGetName(cons), b, ((int) (size_t) SCIPhashmapGetImage(DECdecompGetConstoblock(decdecomp), cons)) -1);  /*lint !e507*/
         assert(SCIPfindCons(scip, SCIPconsGetName(cons)) != NULL);
         assert(((int) (size_t) SCIPhashmapGetImage(DECdecompGetConstoblock(decdecomp), cons)) -1 == b); /*lint !e507*/
         ncurvars = SCIPgetNVarsXXX(scip, cons);
         SCIP_CALL( SCIPallocMemoryArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, cons, curvars, ncurvars) );

         for( v = 0; v < ncurvars; ++v )
         {
            int varblock;
            SCIP_VAR* var = SCIPvarGetProbvar(curvars[v]);
            varblock = ((int) (size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), var)) -1;  /*lint !e507*/
            SCIPdebugMessage("\tVar <%s> in block %d = %d\n", SCIPvarGetName(var), b, varblock);
            assert(SCIPfindVar(scip, SCIPvarGetName(var)) != NULL);
            assert(SCIPvarIsActive(var));
            assert(varblock == b || varblock == DECdecompGetNBlocks(decdecomp));
         }
         SCIPfreeMemoryArray(scip, &curvars);
      }

      for( v = 0; v < DECdecompGetNSubscipvars(decdecomp)[b]; ++v )
      {
         int varblock;
         SCIP_VAR* var = DECdecompGetSubscipvars(decdecomp)[b][v];
         varblock = ((int) (size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), var)) -1; /*lint !e507*/
         SCIPdebugMessage("Var <%s> in block %d = %d\n", SCIPvarGetName(var), b, varblock);
         assert(SCIPfindVar(scip, SCIPvarGetName(var)) != NULL);
         assert(SCIPvarIsActive(var));
         assert(varblock == b || varblock == DECdecompGetNBlocks(decdecomp));
      }
   }

   /* check linking constraints and variables */
   for( v = 0; v < DECdecompGetNLinkingvars(decdecomp); ++v )
   {
      assert(((int) (size_t) SCIPhashmapGetImage(DECdecompGetVartoblock(decdecomp), DECdecompGetLinkingvars(decdecomp)[v])) -1 == DECdecompGetNBlocks(decdecomp)); /*lint !e507*/
   }
   for (c = 0; c < DECdecompGetNLinkingconss(decdecomp); ++c)
   {
      assert(((int) (size_t) SCIPhashmapGetImage(DECdecompGetConstoblock(decdecomp), DECdecompGetLinkingconss(decdecomp)[c])) -1 ==  DECdecompGetNBlocks(decdecomp)); /*lint !e507*/
   }

   switch( DECdecompGetType(decdecomp) )
   {
   case DEC_DECTYPE_UNKNOWN:
         assert(FALSE);
      break;
   case DEC_DECTYPE_ARROWHEAD:
      assert(DECdecompGetNLinkingvars(decdecomp) > 0);
      break;
   case DEC_DECTYPE_BORDERED:
      assert(DECdecompGetNLinkingvars(decdecomp) == 0 && DECdecompGetNLinkingconss(decdecomp) > 0);
      break;
   case DEC_DECTYPE_DIAGONAL:
      assert(DECdecompGetNLinkingvars(decdecomp) == 0 && DECdecompGetNLinkingconss(decdecomp) == 0);
      break;
   case DEC_DECTYPE_STAIRCASE:
      assert(DECdecompGetNLinkingvars(decdecomp) > 0 && DECdecompGetNLinkingconss(decdecomp) == 0);
      break;
   default:
         assert(FALSE);
         break;
   }

#endif
   return SCIP_OKAY;
}

/** returns whether the constraint belongs to GCG or not */
SCIP_Bool GCGisConsGCGCons(
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSHDLR* conshdlr;
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   if( strcmp("origbranch", SCIPconshdlrGetName(conshdlr)) == 0 )
      return TRUE;
   else if( strcmp("masterbranch", SCIPconshdlrGetName(conshdlr)) == 0 )
      return TRUE;

   return FALSE;
}

/** creates a decomposition with all constraints in the master */
extern
SCIP_RETCODE DECcreateBasicDecomp(
   SCIP*                 scip,                /**< SCIP data structure */
   DEC_DECOMP**          decomp               /**< decomposition structure */
   )
{
   SCIP_HASHMAP* constoblock;
   SCIP_CONS** conss;
   int nconss;
   int i;
   int c;

   assert(scip != NULL);
   assert(decomp != NULL);

   SCIP_CALL( DECdecompCreate(scip, decomp) );
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss) );

   for( c = 0,i = 0; c < nconss; ++c )
   {
      if( GCGisConsGCGCons(conss[c]) )
         continue;

      SCIP_CALL( SCIPhashmapInsert(constoblock, conss[c], (void*) (size_t) 1 ) );
      ++i;
   }

   SCIP_CALL( DECfilloutDecdecompFromConstoblock(scip, *decomp, constoblock, 0, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetConss(scip), SCIPgetNConss(scip), FALSE) );

   return SCIP_OKAY;
}
