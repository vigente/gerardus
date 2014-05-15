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

/**@file   pub_decomp.h
 * @ingroup DECOMP
 * @ingroup PUBLICMETHODS
 * @brief  public methods for working with decomposition structures
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef GCG_PUB_DECOMP_H__
#define GCG_PUB_DECOMP_H__

#include "type_decomp.h"
#include "scip/type_scip.h"
#include "scip/type_retcode.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "scip/type_misc.h"
#include "type_detector.h"

#ifdef __cplusplus
extern "C" {
#endif

/** converts the DEC_DECTYPE enum to a string */
const char *DECgetStrType(
   DEC_DECTYPE           type                /**< decomposition type */
   );

/** initializes the decdecomp structure to absolutely nothing */
SCIP_RETCODE DECdecompCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP**          decdecomp           /**< decdecomp instance */
   );

/** frees the decdecomp structure */
SCIP_RETCODE DECdecompFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP**          decdecomp           /**< decdecomp instance */
   );

/** sets the type of the decomposition */
void DECdecompSetType(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   DEC_DECTYPE           type,               /**< type of the decomposition */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   );

/** gets the type of the decomposition */
DEC_DECTYPE DECdecompGetType(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** sets the presolved flag for decomposition */
void DECdecompSetPresolved(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_Bool             presolved           /**< presolved flag for decomposition */
   );

/** gets the presolved flag for decomposition */
SCIP_Bool DECdecompGetPresolved(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** sets the number of blocks for decomposition */
void DECdecompSetNBlocks(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   int                   nblocks             /**< number of blocks for decomposition */
   );

/** gets the number of blocks for decomposition */
int DECdecompGetNBlocks(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** copies the input subscipvars array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetSubscipvars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_VAR***           subscipvars,        /**< subscipvars array  */
   int*                  nsubscipvars,       /**< number of subscipvars per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   );

/** returns the subscipvars array of the given decdecomp structure */
SCIP_VAR*** DECdecompGetSubscipvars(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** returns the nsubscipvars array of the given decdecomp structure */
int* DECdecompGetNSubscipvars(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** copies the input subscipconss array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetSubscipconss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_CONS***          subscipconss,       /**< subscipconss array  */
   int*                  nsubscipconss,      /**< number of subscipconss per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   );

/** returns the subscipconss array of the given decdecomp structure */
SCIP_CONS*** DECdecompGetSubscipconss(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** returns the nsubscipconss array of the given decdecomp structure */
int*  DECdecompGetNSubscipconss(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** copies the input linkingconss array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetLinkingconss(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_CONS**           linkingconss,       /**< linkingconss array  */
   int                   nlinkingconss,      /**< number of linkingconss per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   );

/** returns the linkingconss array of the given decdecomp structure */
SCIP_CONS**  DECdecompGetLinkingconss(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** returns the nlinkingconss array of the given decdecomp structure */
int  DECdecompGetNLinkingconss(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** copies the input linkingvars array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetLinkingvars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_VAR**            linkingvars,        /**< linkingvars array  */
   int                   nlinkingvars,       /**< number of linkingvars per block */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   );

/** returns the linkingvars array of the given decdecomp structure */
SCIP_VAR** DECdecompGetLinkingvars(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** returns the nlinkingvars array of the given decdecomp structure */
int DECdecompGetNLinkingvars(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** copies the input stairlinkingvars array to the given decdecomp structure */
SCIP_RETCODE DECdecompSetStairlinkingvars(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_VAR***           stairlinkingvars,   /**< stairlinkingvars array  */
   int*                  nstairlinkingvars,  /**< number of linkingvars per block */
   SCIP_Bool*            valid               /**< returns whether the resulting decdecomp is valid */
   );

/** returns the stairlinkingvars array of the given decdecomp structure */
SCIP_VAR*** DECdecompGetStairlinkingvars(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** returns the nstairlinkingvars array of the given decdecomp structure */
int* DECdecompGetNStairlinkingvars(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** sets the vartoblock hashmap of the given decdecomp structure */
void DECdecompSetVartoblock(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_HASHMAP*         vartoblock,         /**< Vartoblock hashmap */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   );

/** returns the vartoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP* DECdecompGetVartoblock(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** sets the constoblock hashmap of the given decdecomp structure */
void DECdecompSetConstoblock(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_HASHMAP*         constoblock,        /**< Constoblock hashmap */
   SCIP_Bool*            valid               /**< pointer to indicate whether the structure is valid */
   );

/** returns the constoblock hashmap of the given decdecomp structure */
SCIP_HASHMAP* DECdecompGetConstoblock(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** sets the varindex hashmap of the given decdecomp structure */
void DECdecompSetVarindex(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_HASHMAP*         varindex            /**< Varindex hashmap */
   );

/** returns the varindex hashmap of the given decdecomp structure */
SCIP_HASHMAP* DECdecompGetVarindex(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** sets the consindex hashmap of the given decdecomp structure */
void DECdecompSetConsindex(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_HASHMAP*         consindex           /**< Consindexk hashmap */
   );

/** returns the consindex hashmap of the given decdecomp structure */
SCIP_HASHMAP* DECdecompGetConsindex(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** completely initializes decdecomp from the values of the hashmaps */
SCIP_RETCODE DECfillOutDecdecompFromHashmaps(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   SCIP_HASHMAP*         vartoblock,         /**< variable to block hashmap */
   SCIP_HASHMAP*         constoblock,        /**< constraint to block hashmap */
   int                   nblocks,            /**< number of blocks */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   SCIP_CONS**           conss,              /**< constraint array */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            valid,              /**< pointer to indicate whether the structure is valid */
   SCIP_Bool             staircase           /**< should the decomposition be a staircase structure */
   );

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
   );

/** sets the detector for the given decdecomp structure */
void DECdecompSetDetector(
   DEC_DECOMP*           decdecomp,          /**< decdecomp instance */
   DEC_DETECTOR*         detector            /**< detector data structure */
   );

/** gets the detector for the given decdecomp structure */
DEC_DETECTOR* DECdecompGetDetector(
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** transforms all constraints and variables, updating the arrays */
SCIP_RETCODE DECdecompTransform(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

/** prints detailed information on the contents of decdecomp on the command line */
void DECdecompPrintDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decdecomp instance */
   );

extern
SCIP_RETCODE DECdecompCheckConsistency(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< decomposition data structure */
   );

/** returns whether the constraint belongs to GCG or not */
extern
SCIP_Bool GCGisConsGCGCons(
   SCIP_CONS*            cons                /**< constraint to check */
   );

/** creates a decomposition with all constraints in the master */
extern
SCIP_RETCODE DECcreateBasicDecomp(
   SCIP*                 scip,                /**< SCIP data structure */
   DEC_DECOMP**          decomp               /**< decomposition structure */
   );

#ifdef __cplusplus
}
#endif
#endif
