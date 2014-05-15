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

/**@file    paraLoadCoordinatorTerminationState.cpp
 * @brief   Load coordinator termination state.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <sstream>
#include "paraLoadCoordinatorTerminationState.h"

using namespace UG;

/** stringfy ParaCalculationState */
std::string
ParaLoadCoordinatorTerminationState::toString(
      )
{
   std::ostringstream s;
   if( isCheckpointState )
   {
      s << "######### LoadCoordinator Rank = " << rank << " is at checkpoint. #########" << std::endl;
   }
   else
   {
      s << "######### LoadCoordinator Rank = " << rank << " is terminated. #########" << std::endl;
   }

   s << "#=== The number of ParaNodes received = " << nReceived << std::endl;
   s << "#=== The number of ParaNodes sent = " << nSent << std::endl;
   s << "#=== ( # sent back immediately = " << nSentBackImmediately << ", # failed to send back = " << nFailedToSendBack << " )" << std::endl;
   s << "#=== ( # sent back immediately ( another node ) = " << nSentBackImmediatelyAnotherNode
     << ", # failed to send back ( another node ) = " << nFailedToSendBackAnotherNode << " )" << std::endl;
   s << "#=== The number of ParaNodes deleted in LoadCoordinator = " << (nDeletedInLc + nDeletedByMerging)
     << " ( by merging: " << nDeletedByMerging << " )" << std::endl;
   s << "#=== Maximum usage of node pool = " << nMaxUsageOfNodePool << std::endl;
   if( nNodesInNodePool )
   {
      s << "#=== LoadCoodibator NodePool is not empty. "
         <<  nNodesInNodePool  << " nodes remained." << std::endl;
   }
   if( nNodesLeftInAllSolvers )
   {
      s << "#=== Solvers have nodes. "
        << nNodesLeftInAllSolvers << " nodes left in Solvers." << std::endl;
   }
   if( isCheckpointState )
   {
      s << "#=== Current global best dual bound value = " <<  externalGlobalBestDualBoundValue
        << "( internal value = " << globalBestDualBoundValue << " )" << std::endl;
      s << "#=== Idle time to checkpoint of this LoadCoordinator  = " << idleTime << std::endl;
      s << "#=== Elapsed time to checkpoint of this LoadCoordinator  = " << runningTime << std::endl;
   }
   else
   {
      s << "#=== Idle time to terminate this LoadCoordinator  = " << idleTime << std::endl;
      s << "#=== Elapsed time to terminate this LoadCoordinator  = " << runningTime << std::endl;
   }
   return s.str();
}

void
ParaLoadCoordinatorTerminationState::write(
      ogzstream &out
      )
{
   out.write((char *)&isCheckpointState, sizeof(bool));
   out.write((char *)&rank, sizeof(int));
   out.write((char *)&nWarmStart, sizeof(int));
   out.write((char *)&nSent, sizeof(long long));
   out.write((char *)&nSentBackImmediately, sizeof(long long));
   out.write((char *)&nSentBackImmediatelyAnotherNode, sizeof(int));
   out.write((char *)&nReceived, sizeof(long long));
   out.write((char *)&nDeletedInLc, sizeof(int));
   out.write((char *)&nFailedToSendBack, sizeof(int));
   out.write((char *)&nFailedToSendBackAnotherNode, sizeof(int));
   out.write((char *)&nMaxUsageOfNodePool, sizeof(int));
   out.write((char *)&nNodesInNodePool, sizeof(int));
   out.write((char *)&nNodesLeftInAllSolvers, sizeof(long long));
   out.write((char *)&globalBestDualBoundValue, sizeof(double));
   out.write((char *)&externalGlobalBestDualBoundValue, sizeof(double));
   out.write((char *)&idleTime, sizeof(double));
   out.write((char *)&runningTime, sizeof(double));
}

bool
ParaLoadCoordinatorTerminationState::read(
      ParaComm *comm,
      igzstream &in
      )
{
   in.read((char *)&isCheckpointState, sizeof(bool));
   if( in.eof() ) return false;
   in.read((char *)&rank, sizeof(int));
   in.read((char *)&nWarmStart, sizeof(int));
   in.read((char *)&nSent, sizeof(long long));
   in.read((char *)&nSentBackImmediately, sizeof(long long));
   in.read((char *)&nSentBackImmediatelyAnotherNode, sizeof(int));
   in.read((char *)&nReceived, sizeof(long long));
   in.read((char *)&nDeletedInLc, sizeof(int));
   in.read((char *)&nFailedToSendBack, sizeof(int));
   in.read((char *)&nFailedToSendBackAnotherNode, sizeof(int));
   in.read((char *)&nMaxUsageOfNodePool, sizeof(int));
   in.read((char *)&nNodesInNodePool, sizeof(int));
   in.read((char *)&nNodesLeftInAllSolvers, sizeof(long long));
   in.read((char *)&globalBestDualBoundValue, sizeof(double));
   in.read((char *)&externalGlobalBestDualBoundValue, sizeof(double));
   in.read((char *)&idleTime, sizeof(double));
   in.read((char *)&runningTime, sizeof(double));
   return true;
}

