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

/**@file   nodesel_master.c
 * @ingroup NODESELECTORS
 * @brief  node selector for coordination of master and original formulation
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "nodesel_master.h"
#include "cons_origbranch.h"
#include "cons_masterbranch.h"
#include "pricer_gcg.h"

#define NODESEL_NAME             "master"
#define NODESEL_DESC             "orig master coordination"
#define NODESEL_STDPRIORITY           0
#define NODESEL_MEMSAVEPRIORITY  100000


/** node selector data */
struct SCIP_NodeselData
{
   SCIP_Longint lastorignodenumber;
};

/*
 * Callback methods
 */

#define nodeselCopyMaster NULL

/** destructor of node selector to free user data (called when SCIP is exiting) */
static
SCIP_DECL_NODESELFREE(nodeselFreeMaster)
{
   SCIP_NODESELDATA* nodeseldata;

   assert(nodesel != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   SCIPfreeMemory(scip, &nodeseldata);

   return SCIP_OKAY;
}


/** initialization method of node selector (called after problem was transformed) */
#define nodeselInitMaster NULL


/** deinitialization method of node selector (called before transformed problem is freed) */
#define nodeselExitMaster NULL


/** solving process initialization method of node selector (called when branch and bound process is about to begin) */
#define nodeselInitsolMaster NULL


/** solving process deinitialization method of node selector (called before branch and bound process data is freed) */
#define nodeselExitsolMaster NULL


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectMaster)
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODE** nodes;
   SCIP_CONS* origcons;
   SCIP_CONS* parentorigcons;
   SCIP_CONS* parentmastercons;
   SCIP* origscip;
   int nnodes;
   SCIP_Longint orignodenumber;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   origscip = GCGpricerGetOrigprob(scip);

   *selnode = NULL;

   orignodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(origscip));

   if( orignodenumber != nodeseldata->lastorignodenumber )
   {
      nodeseldata->lastorignodenumber = orignodenumber;

      origcons = GCGconsOrigbranchGetActiveCons(origscip);
      parentorigcons = GCGconsOrigbranchGetParentcons(origcons);

      /* check whether the current node is the root node and has no parent */
      if( parentorigcons == NULL )
      {
         assert((GCGconsOrigbranchGetNode(origcons) == SCIPgetRootNode(origscip)) || ( GCGconsOrigbranchGetNode(origcons) == NULL) );
         assert(GCGconsOrigbranchGetMastercons(origcons) != NULL);
            assert((GCGconsMasterbranchGetNode(GCGconsOrigbranchGetMastercons(origcons)) == SCIPgetRootNode(scip)) || (GCGconsMasterbranchGetNode(GCGconsOrigbranchGetMastercons(origcons)) == NULL));

         *selnode = SCIPgetRootNode(scip);
         SCIPdebugMessage("selected root node in the master program\n");
      }
      else
      {
         parentmastercons = GCGconsOrigbranchGetMastercons(parentorigcons);
         assert(parentmastercons != NULL);

         assert( (GCGconsOrigbranchGetChild1cons(parentorigcons) == origcons)
            != (GCGconsOrigbranchGetChild2cons(parentorigcons) == origcons));

         /* the original cons is the left child of its parentcons,
            select the left child of the corresponding parentcons in the master*/
         if( GCGconsOrigbranchGetChild1cons(parentorigcons) == origcons )
         {
            assert(GCGconsMasterbranchGetChild1cons(parentmastercons) != NULL);
            assert(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetChild1cons(parentmastercons)) != NULL);

            *selnode = GCGconsMasterbranchGetNode(GCGconsMasterbranchGetChild1cons(parentmastercons));
            SCIPdebugMessage("Master nodeselector selected node %"SCIP_LONGINT_FORMAT" corresponding to node %"SCIP_LONGINT_FORMAT" in the original program, since the parents (%"SCIP_LONGINT_FORMAT"/o, %"SCIP_LONGINT_FORMAT"/m) are linked\n",
               SCIPnodeGetNumber(*selnode), SCIPnodeGetNumber(GCGconsOrigbranchGetNode(origcons)),
               SCIPnodeGetNumber(GCGconsOrigbranchGetNode(parentorigcons)),
               SCIPnodeGetNumber(GCGconsMasterbranchGetNode(parentmastercons)));
         }

         /* the original cons is the right child of its parentcons,
            select the right child of the corresponding parentcons in the master */
         else
         {
            assert(GCGconsOrigbranchGetChild2cons(parentorigcons) == origcons);
            assert(GCGconsMasterbranchGetChild2cons(parentmastercons) != NULL);
            assert(GCGconsMasterbranchGetNode(GCGconsMasterbranchGetChild2cons(parentmastercons)) != NULL);

            *selnode = GCGconsMasterbranchGetNode(GCGconsMasterbranchGetChild2cons(parentmastercons));
            SCIPdebugMessage("Master nodeselector selected node %"SCIP_LONGINT_FORMAT" corresponding to node %"SCIP_LONGINT_FORMAT" in the original program, since the parents (%"SCIP_LONGINT_FORMAT"/o, %"SCIP_LONGINT_FORMAT"/m) are linked\n",
               SCIPnodeGetNumber(*selnode), SCIPnodeGetNumber(GCGconsOrigbranchGetNode(origcons)),
               SCIPnodeGetNumber(GCGconsOrigbranchGetNode(parentorigcons)),
               SCIPnodeGetNumber(GCGconsMasterbranchGetNode(parentmastercons)));
         }

      }

      if( *selnode == NULL )
      {
         SCIPerrorMessage("nodesel_master could not find a node corresponding to the current original node!\n");
      }
      assert(*selnode != NULL);

      /* set the dual bound to the lower bound of the corresponding original node */
      SCIP_CALL( SCIPupdateNodeDualbound(scip, *selnode, SCIPgetNodeLowerbound(origscip, SCIPgetCurrentNode(origscip))) );
   }
   else
   {
      SCIPdebugMessage("select random node\n");

      if( SCIPgetNChildren(scip) > 0 )
      {
         SCIP_CALL( SCIPgetChildren(scip, &nodes, &nnodes) );
         *selnode = nodes[0];
      }
      else if( SCIPgetNSiblings(scip) > 0 )
      {
         SCIP_CALL( SCIPgetSiblings(scip, &nodes, &nnodes) );
         *selnode = nodes[0];
      }
      else if( SCIPgetNLeaves(scip) > 0 )
      {
         SCIP_CALL( SCIPgetLeaves(scip, &nodes, &nnodes) );
         *selnode = nodes[0];
      }
   }

#ifndef NDEBUG
   GCGconsOrigbranchCheckConsistency(origscip);
   GCGconsMasterbranchCheckConsistency(scip);
#endif

   return SCIP_OKAY;
}


/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompMaster)
{
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   if( SCIPnodeGetNumber(node1) < SCIPnodeGetNumber(node2) )
      return 1;
   else
      return -1;
}

/*
 * master specific interface methods
 */

/** creates the node selector for depth first search and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselMaster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;

   /* create master node selector data */
   SCIP_CALL( SCIPallocMemory(scip, &nodeseldata) );

   nodeseldata->lastorignodenumber = -1LL;

   /* include node selector */
   SCIP_CALL( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
         nodeselCopyMaster, nodeselFreeMaster, nodeselInitMaster, nodeselExitMaster,
         nodeselInitsolMaster, nodeselExitsolMaster, nodeselSelectMaster, nodeselCompMaster,
         nodeseldata) );

   return SCIP_OKAY;
}

