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

/**@file    paraCalculationStateMpi.cpp
 * @brief   CalcutationStte object extension for MPI communication
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraCalculationStateMpi.h"

using namespace UG;

MPI::Datatype
ParaCalculationStateMpi::createDatatype(
      )
{

   const int nBlocks = 14;

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
         startAddress,  MPI::Get_address( &compTime )
   );
   displacements[0] = 0;
   types[0] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &rootTime )
   );
   displacements[1] = address - startAddress;
   types[1] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nSolved )
   );
   displacements[2] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nSent )
   );
   displacements[3] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nImprovedIncumbent )
   );
   displacements[4] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &terminationState )
   );
   displacements[5] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nSolvedWithNoPreprocesses )
   );
   displacements[6] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nSimplexIterRoot )
   );
   displacements[7] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &averageSimplexIter )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nRestarts )
   );
   displacements[9] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &minIisum )
   );
   displacements[10] = address - startAddress;
   types[10] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &maxIisum )
   );
   displacements[11] = address - startAddress;
   types[11] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &minNii )
   );
   displacements[12] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &maxNii )
   );
   displacements[13] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return datatype;

}

void
ParaCalculationStateMpi::send(
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
      commMpi->send(&compTime, 1, datatype, destination, tag)
   );
   datatype.Free();
}

void
ParaCalculationStateMpi::receive(
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
      commMpi->receive(&compTime, 1, datatype, source, tag)
   );
   datatype.Free();
}
