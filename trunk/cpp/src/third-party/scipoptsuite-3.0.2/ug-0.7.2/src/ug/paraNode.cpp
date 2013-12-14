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

/**@file    paraNode.cpp
 * @brief   Base class for ParaNode.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraComm.h"
#include "paraNode.h"

using namespace UG;

void
ParaNode::write(ogzstream &out){
   out.write((char *)&nodeId.subtreeId.lcId, sizeof(int));
   out.write((char *)&nodeId.subtreeId.globalSubtreeIdInLc, sizeof(int));
   out.write((char *)&nodeId.subtreeId.solverId, sizeof(int));
   out.write((char *)&nodeId.seqNum, sizeof(long long));
   out.write((char *)&generatorNodeId.subtreeId.lcId, sizeof(int));
   out.write((char *)&generatorNodeId.subtreeId.globalSubtreeIdInLc, sizeof(int));
   out.write((char *)&generatorNodeId.subtreeId.solverId, sizeof(int));
   out.write((char *)&generatorNodeId.seqNum, sizeof(long long));
   out.write((char *)&depth, sizeof(int));
   out.write((char *)&dualBoundValue, sizeof(double));
   out.write((char *)&initialDualBoundValue, sizeof(double));
   out.write((char *)&estimatedValue, sizeof(double));
   out.write((char *)&diffSubproblemInfo, sizeof(int));
   if( !mergeNodeInfo )
   {
      if( diffSubproblemInfo ) diffSubproblem->write(out);
   }
   else
   {
      if( mergeNodeInfo->origDiffSubproblem )
      {
         mergeNodeInfo->origDiffSubproblem->write(out);
      }
      else
      {
         if( diffSubproblemInfo ) diffSubproblem->write(out);
      }
   }
   out.write((char *)&basisInfo, sizeof(int));
}

bool
ParaNode::read(ParaComm *comm, igzstream &in, bool onlyBoundChanges){
   in.read((char *)&nodeId.subtreeId.lcId, sizeof(int));
   if( in.eof() ) return false;
   in.read((char *)&nodeId.subtreeId.globalSubtreeIdInLc, sizeof(int));
   in.read((char *)&nodeId.subtreeId.solverId, sizeof(int));
   in.read((char *)&nodeId.seqNum, sizeof(long long));
   in.read((char *)&generatorNodeId.subtreeId.lcId, sizeof(int));
   in.read((char *)&generatorNodeId.subtreeId.globalSubtreeIdInLc, sizeof(int));
   in.read((char *)&generatorNodeId.subtreeId.solverId, sizeof(int));
   in.read((char *)&generatorNodeId.seqNum, sizeof(long long));
   in.read((char *)&depth, sizeof(int));
   in.read((char *)&dualBoundValue, sizeof(double));
   in.read((char *)&initialDualBoundValue, sizeof(double));
   in.read((char *)&estimatedValue, sizeof(double));
   in.read((char *)&diffSubproblemInfo, sizeof(int));
   if( diffSubproblemInfo ){
      diffSubproblem = comm->createParaDiffSubproblem();
      diffSubproblem->read(comm, in, onlyBoundChanges);
   }
   in.read((char *)&basisInfo, sizeof(int));
   return true;
}
