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

/**@file   cons_origbranch.h
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 */

#ifndef GCG_CONS_ORIGBRANCH_H__
#define GCG_CONS_ORIGBRANCH_H__

#include "scip/scip.h"
#include "type_branchgcg.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns the branch orig constraint of the current node, only needs the pointer to scip */
extern
SCIP_CONS* GCGconsOrigbranchGetActiveCons(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** creates the handler for origbranch constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrOrigbranch(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a origbranch constraint */
extern
SCIP_RETCODE GCGcreateConsOrigbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_NODE*            node,               /**< the node to which this origbranch constraint belongs */
   SCIP_CONS*            parentcons,         /**< origbranch constraint associated with the father node */
   SCIP_BRANCHRULE*      branchrule,         /**< the branching rule that created the b&b node the constraint belongs to */
   GCG_BRANCHDATA*       branchdata          /**< branching data storing information about the branching restrictions at the
                                              *   corresponding node */
   );

/** returns the stack and the number of elements on it */
extern
void GCGconsOrigbranchGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   );

/** returns the branching data for a given origbranch constraint */
extern
GCG_BRANCHDATA* GCGconsOrigbranchGetBranchdata(
   SCIP_CONS*            cons                /**< origbranch constraint for which the branching data is requested */
   );

/** returns the branchrule for a given origbranch constraint */
extern
SCIP_BRANCHRULE* GCGconsOrigbranchGetBranchrule(
   SCIP_CONS*            cons                /**< origbranch constraint for which the branchrule is requested */
   );

/** returns the node in the B&B tree at which the given origbranch constraint is sticking */
extern
SCIP_NODE* GCGconsOrigbranchGetNode(
   SCIP_CONS*            cons                /**< origbranch constraint for which the corresponding node is requested */
   );

/** returns the origbranch constraint of the B&B father of the node at which the
    given origbranch constraint is sticking */
extern
SCIP_CONS* GCGconsOrigbranchGetParentcons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the origbranch constraint of
                                              *   the father node is requested */
   );

/** returns the origbranch constraint of the first child of the node at which the
    given origbranch constraint is sticking */
extern
SCIP_CONS* GCGconsOrigbranchGetChild1cons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the origbranch constraint of
                                              *   the first child node is requested */
   );

/** returns the origbranch constraint of the second child of the node at which the
    given origbranch constraint is sticking */
extern
SCIP_CONS* GCGconsOrigbranchGetChild2cons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the origbranch constraint of
                                              *   the second child node is requested */
   );

/** sets the masterbranch constraint of the node in the master program corresponding to the node
    at which the given origbranchbranch constraint is sticking */
extern
void GCGconsOrigbranchSetMastercons(
   SCIP_CONS*            cons,               /**< origbranch constraint for which the masterbranch constraint should be set */
   SCIP_CONS*            mastercons          /**< masterbranch constraint corresponding to the given origbranch constraint */
   );

/** returns the masterbranch constraint of the node in the master program corresponding to the node
    at which the given origbranchbranch constraint is sticking */
extern
SCIP_CONS* GCGconsOrigbranchGetMastercons(
   SCIP_CONS*            cons                /**< origbranch constraint for which the corresponding masterbranch
                                              *   constraint is requested */
   );

/** checks the consistency of the origbranch constraints in the problem */
extern
void GCGconsOrigbranchCheckConsistency(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds a bound change on an original variable found by propagation in the original problem
 *  to the given origbranch constraint so that is will be transferred to the master problem */
extern
SCIP_RETCODE GCGconsOrigbranchAddPropBoundChg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< origbranch constraint to which the bound change is added */
   SCIP_VAR*             var,                /**< variable on which the bound change was performed */
   SCIP_BOUNDTYPE        boundtype,          /**< bound type of the bound change */
   SCIP_Real             newbound            /**< new bound of the variable after the bound change */
   );

/** returns the array of bound changes on original variables found by propagation in the original problem
 *  at the node corresponding to the given origbranch constraint and clears the arrays */
extern
SCIP_RETCODE GCGconsOrigbranchGetPropBoundChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< origbranch constraint for which the bound changes are requested */
   SCIP_VAR***           vars,               /**< pointer to store array of variables corresponding to the bound changes */
   SCIP_BOUNDTYPE**      boundtypes,         /**< pointer to store array of the types of the bound changes */
   SCIP_Real**           newbounds,          /**< pointer to store array of the new bounds */
   int*                  npropbounds         /**< pointer to store the number of bound changes stored at the constraint */
   );

/** returns the number of bound changes on original variables found by propagation in the original problem
 *  at the node corresponding to the given origbranch constraint */
extern
int GCGconsOrigbranchGetNPropBoundChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< origbranch constraint for which the bound changes are requested */
   );

/** adds initial constraint to root node */
extern
SCIP_RETCODE SCIPconsOrigbranchAddRootCons(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
