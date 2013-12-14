/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*             This file is part of the program and software framework       */
/*                  UG --- Ubquity Generator Framework                       */
/*                                                                           */
/*    Copyright (C) 2010-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  UG is distributed under the terms of the ZIB Academic Licence.           */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with UG; see the file COPYING. If not email to scip@zib.de.        */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    paraNodePth.h
 * @brief   ParaNode extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_NODE_PTH_H__
#define __PARA_NODE_PTH_H__

#include <iostream>
#include <fstream>
#include "paraComm.h"
#include "paraNode.h"

namespace UG
{

/** ParaNodeMpi class */
class ParaNodePth : public ParaNode
{
   /** create ParaNode datatype */
   ParaNodePth *createDatatype(ParaComm *comm);

public :
   /** default constructor */
   ParaNodePth(
         )
   {
   }

   /** constructor */
   ParaNodePth( NodeId inNodeId,
         NodeId inGeneratorNodeId,
         int inDepth,
         double inDualBoundValue,
         double inOriginalDualBoundValue,
         double inEstimatedValue,
         ParaDiffSubproblem *inDiffSubproblem
         )
         : ParaNode(inNodeId, inGeneratorNodeId, inDepth, inDualBoundValue, inOriginalDualBoundValue, inEstimatedValue, inDiffSubproblem)
   {
   }

   /** destructor */
   ~ParaNodePth(
         )
   {
   }

   ParaNodePth *clone(ParaComm *comm) {
      if( diffSubproblem )
      {
         return ( new
            ParaNodePth(nodeId, generatorNodeId, depth, dualBoundValue, initialDualBoundValue,
                  initialDualBoundValue,diffSubproblem->clone(comm) ) );
      }
      else
      {
         return ( new
            ParaNodePth(nodeId, generatorNodeId, depth, dualBoundValue, initialDualBoundValue,
                  initialDualBoundValue, 0 ) );
      }
   }
   int bcast( ParaComm *comm, int root );
   int send( ParaComm *comm, int destination );
   int receive( ParaComm *comm, int source );
};

}

#endif // __PARA_NODE_PTH_H__
