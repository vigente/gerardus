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

/**@file    scipParaInitialStatMpi.cpp
 * @brief   ScipParaInitialStat extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


// #include "paraTagDefMpi.h"
#include "scipParaInitialStatMpi.h"

using namespace UG;
using namespace ParaSCIP;

/** create ScipDiffSubproblemPreDatatype */
MPI::Datatype
ScipParaInitialStatMpi::createDatatype1(
      )
{
   const int nBlocks = 4;

   MPI::Datatype datatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   int blockLengths[nBlocks];
   MPI::Aint displacements[nBlocks];
   MPI::Datatype types[nBlocks];

   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( &maxDepth )
   );
   blockLengths[0] = 1;
   displacements[0] = 0;
   types[0] = MPI::INT;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &maxTotalDepth )
   );
   blockLengths[1] = 1;
   displacements[1] = address - startAddress;
   types[1] = MPI::INT;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nVarBranchStatsDown )
   );
   blockLengths[2] = 1;
   displacements[2] = address - startAddress;
   types[2] = MPI::INT;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nVarBranchStatsUp )
   );
   blockLengths[3] = 1;
   displacements[3] = address - startAddress;
   types[3] = MPI::INT;

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return datatype;
}

/** create ScipDiffSubproblemDatatype */
MPI::Datatype
ScipParaInitialStatMpi::createDatatype2(
      bool memAllocNecessary
      )
{
   assert( nVarBranchStatsDown != 0 || nVarBranchStatsUp != 0 );

   int nBlocks = 0;

   MPI::Datatype datatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   int blockLengths[17];            // reserve maximum number of elements
   MPI::Aint displacements[17];     // reserve maximum number of elements
   MPI::Datatype types[17];         // reserve maximum number of elements

   if( nVarBranchStatsDown )
   {
      if( memAllocNecessary )
      {
         idxLBranchStatsVarsDown = new int[nVarBranchStatsDown];
         nVarBranchingDown = new int[nVarBranchStatsDown];
         downpscost = new SCIP_Real[nVarBranchStatsDown];
         downvsids = new SCIP_Real[nVarBranchStatsDown];
         downconflen = new SCIP_Real[nVarBranchStatsDown];
         downinfer = new SCIP_Real[nVarBranchStatsDown];
         downcutoff = new SCIP_Real[nVarBranchStatsDown];
      }

      MPI_CALL_WITH_RET_VAL(
            startAddress,  MPI::Get_address( idxLBranchStatsVarsDown )
      );
      displacements[nBlocks] = 0;
      blockLengths[nBlocks] = nVarBranchStatsDown;
      types[nBlocks] = MPI::INT;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nVarBranchingDown )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nVarBranchStatsDown;
      types[nBlocks] = MPI::INT;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downpscost )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nVarBranchStatsDown;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downvsids )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nVarBranchStatsDown;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downconflen )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nVarBranchStatsDown;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downinfer )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nVarBranchStatsDown;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downcutoff )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nVarBranchStatsDown;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;
   }

   if( nVarBranchStatsUp )
   {
      if( memAllocNecessary )
      {
         idxLBranchStatsVarsUp = new int[nVarBranchStatsUp];
         nVarBranchingUp = new int[nVarBranchStatsUp];
         uppscost = new SCIP_Real[nVarBranchStatsUp];
         upvsids = new SCIP_Real[nVarBranchStatsUp];
         upconflen = new SCIP_Real[nVarBranchStatsUp];
         upinfer = new SCIP_Real[nVarBranchStatsUp];
         upcutoff = new SCIP_Real[nVarBranchStatsUp];
      }

      if( nBlocks == 0 )
      {
         MPI_CALL_WITH_RET_VAL(
               startAddress,  MPI::Get_address( idxLBranchStatsVarsUp )
         );
         displacements[nBlocks] = 0;
      }
      else
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( idxLBranchStatsVarsUp )
         );
         displacements[nBlocks] = address - startAddress;
      }
      blockLengths[nBlocks] = nVarBranchStatsUp;
      types[nBlocks] = MPI::INT;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nVarBranchingUp )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStatsUp;
      types[nBlocks] = MPI::INT;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( uppscost )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStatsUp;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( upvsids )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStatsUp;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( upconflen )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStatsUp;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( upinfer )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStatsUp;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( upcutoff )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStatsUp;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;
   }

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
   );

   return datatype;
}

/** send solution data to the rank */
void
ScipParaInitialStatMpi::send(ParaComm *comm, int destination)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype datatype1;
   datatype1 = createDatatype1();
   datatype1.Commit();
   PARA_COMM_CALL(
      commMpi->send(&maxDepth, 1, datatype1, destination, TagInitialStat)
   );
   datatype1.Free();

   if( nVarBranchStatsDown !=0 || nVarBranchStatsUp != 0 )
   {
      MPI::Datatype datatype2;
      datatype2 = createDatatype2(false);
      datatype2.Commit();
      if( nVarBranchStatsDown )
      {
         PARA_COMM_CALL(
               commMpi->send(idxLBranchStatsVarsDown, 1, datatype2, destination, TagInitialStat1)
         );
      }
      else
      {
         if( nVarBranchStatsUp )
         {
            PARA_COMM_CALL(
                  commMpi->send(idxLBranchStatsVarsUp, 1, datatype2, destination, TagInitialStat1)
            );
         }
      }
      datatype2.Free();
   }
}

/** receive solution data from the source rank */
void
ScipParaInitialStatMpi::receive(ParaComm *comm, int source)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype datatype1;
   datatype1 = createDatatype1();
   datatype1.Commit();
   PARA_COMM_CALL(
      commMpi->receive(&maxDepth, 1, datatype1, source, TagInitialStat)
   );
   datatype1.Free();

   if( nVarBranchStatsDown !=0 || nVarBranchStatsUp != 0 )
   {
      MPI::Datatype datatype2;
      datatype2 = createDatatype2(true);
      datatype2.Commit();
      if( nVarBranchStatsDown )
      {
         PARA_COMM_CALL(
               commMpi->receive(idxLBranchStatsVarsDown, 1, datatype2, source, TagInitialStat1)
         );
      }
      else
      {
         if( nVarBranchStatsUp )
         {
            PARA_COMM_CALL(
                  commMpi->receive(idxLBranchStatsVarsUp, 1, datatype2, source, TagInitialStat1)
            );
         }
      }
      datatype2.Free();
   }
}
