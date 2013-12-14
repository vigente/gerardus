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

/**@file    paraSolverTerminationState.cpp
 * @brief   This class contains solver termination state which is transferred form Solver to LC.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraDef.h"
#include "paraComm.h"
#include "paraSolverTerminationState.h"

using namespace UG;

/** stringfy ParaCalculationState */
std::string
ParaSolverTerminationState::toString(
      )
{
   std::ostringstream os;
   switch( interrupted )
   {
   case 0:
   {
      os << "######### Solver Rank = " << rank << " is terminated. #########" << std::endl;
      break;
   }
   case 1:
   {
      os << "######### Solver Rank = " << rank << " is interrupted. #########" << std::endl;
      break;
   }
   case 2:
   {
      os << "######### Solver Rank = " << rank << " is at checkpoint. #########" << std::endl;
      break;
   }
   case 3:
   {
      os << "######### Solver Rank = " << rank << " is at the end of racing process. #########" << std::endl;
      break;
   }
   default:
   {
      THROW_LOGICAL_ERROR1("invalid interrupted flag in ParaSolverTerminationState!");
   }
   }

    os << "#=== Elapsed time to terminate this Solver = " << runningTime << std::endl;
    os << "#=== Total computing time = " << runningTime - (idleTimeToFirstParaNode+idleTimeBetweenParaNodes+idleTimeAfterLastParaNode+idleTimeToWaitNotificationId+idleTimeToWaitToken ) << std::endl;
    os << "#=== Total idle time = " << (idleTimeToFirstParaNode+idleTimeBetweenParaNodes+idleTimeAfterLastParaNode+idleTimeToWaitNotificationId+idleTimeToWaitToken) << std::endl;
    os << "#=== ( Idle time to start first ParaNode = " << idleTimeToFirstParaNode
    << ", Idle time between ParaNods = " << idleTimeBetweenParaNodes
    << ", Idle Time after last ParaNode = " << idleTimeAfterLastParaNode
    << " )" << std::endl;
    os << "#=== ( Idle time to wait notification Id messages = " << idleTimeToWaitNotificationId << " )" << std::endl;
    os << "#=== ( Idle time to wait acknowledgment of completion = " << idleTimeToWaitAckCompletion << " )" << std::endl;
    os << "#=== ( Idle time to wait token = " << idleTimeToWaitToken << " )" << std::endl;
    if( nParaNodesSolved > 0 )
    {
       os << "#=== Total root node process time = " << totalRootNodeTime << " ( Mean = " << (totalRootNodeTime/nParaNodesSolved)
       << ", Min = " << minRootNodeTime << ", Max = " << maxRootNodeTime << " )" << std::endl;
    }
    else
    {
       os << "#=== Total root node process time = 0.0 ( Mean = 0.0, Min = 0.0, Max = 0.0 )" << std::endl;
    }
   os << "#=== The number of ParaNodes received in this solver = " << nParaNodesReceived << std::endl;
   os << "#=== The number of ParaNodes sent from this solver = " << totalNSent << std::endl;
   if( nParaNodesSolved > 0 )
   {
      os << "#=== The number of nodes solved in this solver = " << totalNSolved
      << " ( / Subtree : Mean = " << totalNSolved/nParaNodesSolved <<  ", Min = " << minNSolved << ", Max = " << maxNSolved << " )"<< std::endl;
   }
   else
   {
      os << "#=== The number of nodes solved in this solver = 0 ( / Subtree : Mean = 0, Min = 0, Max = 0 )" << std::endl;
   }
   os << "#=== The number of ParaNodes solved in this solver = " << nParaNodesSolved << std::endl;
   os << "#=== ( Solved at root node  =  " << nParaNodesSolvedAtRoot << ", Solved at pre-checking of root node solvability = "
   << nParaNodesSolvedAtPreCheck << " )" << std::endl;
   os << "#=== The number of improved solutions found in this solver = " << totalNImprovedIncumbent << std::endl;

   return os.str();
}

void
ParaSolverTerminationState::write(
      ogzstream &out
      )
{
   out.write((char *)&interrupted, sizeof(int));
   out.write((char *)&rank, sizeof(int));
   out.write((char *)&totalNSolved, sizeof(int));
   out.write((char *)&minNSolved, sizeof(int));
   out.write((char *)&maxNSolved, sizeof(int));
   out.write((char *)&totalNSent, sizeof(int));
   out.write((char *)&totalNImprovedIncumbent, sizeof(int));
   out.write((char *)&nParaNodesReceived, sizeof(int));
   out.write((char *)&nParaNodesSolved, sizeof(int));
   out.write((char *)&nParaNodesSolvedAtRoot, sizeof(int));
   out.write((char *)&nParaNodesSolvedAtPreCheck, sizeof(int));
   out.write((char *)&runningTime, sizeof(double));
   out.write((char *)&idleTimeToFirstParaNode, sizeof(double));
   out.write((char *)&idleTimeBetweenParaNodes, sizeof(double));
   out.write((char *)&idleTimeAfterLastParaNode, sizeof(double));
   out.write((char *)&idleTimeToWaitNotificationId, sizeof(double));
   out.write((char *)&idleTimeToWaitToken, sizeof(double));
   out.write((char *)&totalRootNodeTime, sizeof(double));
   out.write((char *)&minRootNodeTime, sizeof(double));
   out.write((char *)&maxRootNodeTime, sizeof(double));
}

bool
ParaSolverTerminationState::read(
      ParaComm *comm,
      igzstream &in
      )
{
   in.read((char *)&interrupted, sizeof(int));
   if( in.eof() ) return false;
   in.read((char *)&rank, sizeof(int));
   in.read((char *)&totalNSolved, sizeof(int));
   in.read((char *)&minNSolved, sizeof(int));
   in.read((char *)&maxNSolved, sizeof(int));
   in.read((char *)&totalNSent, sizeof(int));
   in.read((char *)&totalNImprovedIncumbent, sizeof(int));
   in.read((char *)&nParaNodesReceived, sizeof(int));
   in.read((char *)&nParaNodesSolved, sizeof(int));
   in.read((char *)&nParaNodesSolvedAtRoot, sizeof(int));
   in.read((char *)&nParaNodesSolvedAtPreCheck, sizeof(int));
   in.read((char *)&runningTime, sizeof(double));
   in.read((char *)&idleTimeToFirstParaNode, sizeof(double));
   in.read((char *)&idleTimeBetweenParaNodes, sizeof(double));
   in.read((char *)&idleTimeAfterLastParaNode, sizeof(double));
   in.read((char *)&idleTimeToWaitNotificationId, sizeof(double));
   in.read((char *)&idleTimeToWaitToken, sizeof(double));
   in.read((char *)&totalRootNodeTime, sizeof(double));
   in.read((char *)&minRootNodeTime, sizeof(double));
   in.read((char *)&maxRootNodeTime, sizeof(double));
   return true;
}

