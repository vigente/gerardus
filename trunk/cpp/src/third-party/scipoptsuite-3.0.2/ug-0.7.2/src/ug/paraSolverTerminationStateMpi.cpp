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

/**@file    paraSolverTerminationStateMpi.cpp
 * @brief   ParaSolverTerminationState extension for MIP communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraComm.h"
#include "paraSolverTerminationStateMpi.h"

using namespace UG;

MPI::Datatype
ParaSolverTerminationStateMpi::createDatatype(){

   const int nBlocks = 22;

   MPI::Datatype datatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   int blockLengths[nBlocks];
   MPI::Aint displacements[nBlocks];
   MPI::Datatype types[nBlocks];

   for( int i = 0; i < nBlocks; i++ ){
      blockLengths[i] = 1;
      types[i] = MPI::INT;
   }

   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( &interrupted )
   );
   displacements[0] = 0;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &rank )
   );
   displacements[1] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &totalNSolved )
   );
   displacements[2] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &minNSolved )
   );
   displacements[3] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &maxNSolved )
   );
   displacements[4] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &totalNSent )
   );
   displacements[5] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &totalNImprovedIncumbent )
   );
   displacements[6] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nParaNodesReceived )
   );
   displacements[7] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nParaNodesSolved )
   );
   displacements[8] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nParaNodesSolvedAtRoot )
   );
   displacements[9] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nParaNodesSolvedAtPreCheck )
   );
   displacements[10] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &runningTime )
   );
   displacements[11] = address - startAddress;
   types[11] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &idleTimeToFirstParaNode )
   );
   displacements[12] = address - startAddress;
   types[12] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &idleTimeBetweenParaNodes )
   );
   displacements[13] = address - startAddress;
   types[13] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &idleTimeAfterLastParaNode )
   );
   displacements[14] = address - startAddress;
   types[14] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &idleTimeToWaitNotificationId )
   );
   displacements[15] = address - startAddress;
   types[15] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &idleTimeToWaitAckCompletion )
   );
   displacements[16] = address - startAddress;
   types[16] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &idleTimeToWaitToken )
   );
   displacements[17] = address - startAddress;
   types[17] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &totalRootNodeTime )
   );
   displacements[18] = address - startAddress;
   types[18] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &minRootNodeTime )
   );
   displacements[19] = address - startAddress;
   types[19] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &maxRootNodeTime )
   );
   displacements[20] = address - startAddress;
   types[20] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &detTime )
   );
   displacements[21] = address - startAddress;
   types[21] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return datatype;

}

void
ParaSolverTerminationStateMpi::send(
      ParaComm *comm,
      int destination,
      int tag
      )
{
   DEF_PARA_COMM( commMpi, comm);

   MPI::Datatype datatype;
   datatype = createDatatype();
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->send(&interrupted, 1, datatype, destination, tag)
   );
   datatype.Free();
}

void
ParaSolverTerminationStateMpi::receive(
      ParaComm *comm,
      int source,
      int tag
      )
{
   DEF_PARA_COMM( commMpi, comm);

   MPI::Datatype datatype;
   datatype = createDatatype();
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->receive(&interrupted, 1, datatype, source, tag)
   );
   datatype.Free();
}
