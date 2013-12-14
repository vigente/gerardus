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

/**@file    paraNodePth.cpp
 * @brief   ParaNode extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraCommPth.h"
#include "paraNodePth.h"

using namespace UG;

ParaNodePth *
ParaNodePth::createDatatype(
    ParaComm *comm
      )
{
   return clone(comm);
}

int
ParaNodePth::bcast(
      ParaComm *comm,
      int root
      )
{
   DEF_PARA_COMM( commPth, comm);

   if( commPth->getRank() == root )
   {
      for( int i = 0; i < commPth->getSize(); i++ )
      {
         if( i != root )
         {
            ParaNodePth *sent;
            sent = createDatatype(comm);
            assert(!(sent->mergeNodeInfo));
            sent->mergeNodeInfo = 0;
            PARA_COMM_CALL(
               commPth->uTypeSend((void *)sent, ParaNodeType, i, TagNode)
            );
         }
      }
   }
   else
   {
      ParaNodePth *received;
      PARA_COMM_CALL(
         commPth->uTypeReceive((void **)&received, ParaNodeType, root, TagNode)
      );
      nodeId = received->nodeId;
      generatorNodeId = received->generatorNodeId;
      depth = received->depth;
      dualBoundValue = received->dualBoundValue;
      initialDualBoundValue = received->initialDualBoundValue;
      estimatedValue = received->estimatedValue;
      diffSubproblemInfo = received->diffSubproblemInfo;
      if( diffSubproblemInfo )
      {
         diffSubproblem = received->diffSubproblem->clone(commPth);
      }
      basisInfo = received->basisInfo;
      mergingStatus = received->mergingStatus;
      delete received;
   }
   return 0;
}

int
ParaNodePth::send(
      ParaComm *comm,
      int destination
      )
{
    DEF_PARA_COMM( commPth, comm);

    ParaNodePth *sent;
    sent = createDatatype(comm);
    assert(!(sent->mergeNodeInfo));
    sent->mergeNodeInfo = 0;
    PARA_COMM_CALL(
       commPth->uTypeSend((void *)sent, ParaNodeType, destination, TagNode)
    );

   return 0;
}

int
ParaNodePth::receive(
      ParaComm *comm,
      int source
      )
{
   DEF_PARA_COMM( commPth, comm);

   ParaNodePth *received;
   PARA_COMM_CALL(
      commPth->uTypeReceive((void **)&received, ParaNodeType, source, TagNode)
   );
   nodeId = received->nodeId;
   generatorNodeId = received->generatorNodeId;
   depth = received->depth;
   dualBoundValue = received->dualBoundValue;
   initialDualBoundValue = received->initialDualBoundValue;
   estimatedValue = received->estimatedValue;
   diffSubproblemInfo = received->diffSubproblemInfo;
   if( diffSubproblemInfo )
   {
      diffSubproblem = received->diffSubproblem->clone(commPth);
   }
   basisInfo = received->basisInfo;
   mergingStatus = received->mergingStatus;
   delete received;

   return 0;
}
