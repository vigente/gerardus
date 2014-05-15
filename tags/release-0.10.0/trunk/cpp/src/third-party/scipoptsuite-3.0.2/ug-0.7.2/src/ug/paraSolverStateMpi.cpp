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

/**@file    paraSolverStateMpi.cpp
 * @brief   ParaSolverState extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraSolverStateMpi.h"

using namespace UG;

MPI::Datatype
ParaSolverStateMpi::createDatatype(
      )
{

   const int nBlocks = 9;

   MPI::Datatype datatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   int blockLengths[nBlocks];
   MPI::Aint displacements[nBlocks];
   MPI::Datatype types[nBlocks];

   for( int i = 0; i < nBlocks; i++ )
   {
      blockLengths[i] = 1;
      types[i] = MPI::INT;
   }

   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( &racingStage )
   );
   displacements[0] = 0;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &notificationId )
   );
   displacements[1] = address - startAddress;
   types[1] = MPI::UNSIGNED;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &lcId )
   );
   displacements[2] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &globalSubtreeIdInLc )
   );
   displacements[3] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nNodesSolved )
   );
   displacements[4] = address - startAddress;
#ifdef _ALIBABA
   types[4] = MPI::LONG;
#else
   types[4] = MPI::LONG_LONG;
#endif

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nNodesLeft )
   );
   displacements[5] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &bestDualBoundValue )
   );
   displacements[6] = address - startAddress;
   types[6] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &globalBestPrimalBoundValue )
   );
   displacements[7] = address - startAddress;
   types[7] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &detTime )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return datatype;
}

void
ParaSolverStateMpi::send(
      ParaComm *comm,
      int destination,
      int tag
      )
{
   assert(nNodesLeft >= 0);
   assert(bestDualBoundValue >= -1e+10);
   DEF_PARA_COMM( commMpi, comm);

   MPI::Datatype datatype;
   datatype = createDatatype();
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->send(&racingStage, 1, datatype, destination, tag)
   );
   datatype.Free();
}

void
ParaSolverStateMpi::receive(
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
      commMpi->receive(&racingStage, 1, datatype, source, tag)
   );
   datatype.Free();
}
