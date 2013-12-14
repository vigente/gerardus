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

/**@file    scipParaDiffSubproblemMpi.cpp
 * @brief   ScipParaDiffSubproblem extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <mpi.h>
#include "scipParaCommMpi.h"
#include "scipParaDiffSubproblemMpi.h"

using namespace UG;
using namespace ParaSCIP;

/** create ScipDiffSubproblemDatatype1 */
/************************************************
 * Currently, Datatype1 is not necessary.       *
 * I create this code for the future extension. *
 ************************************************/
MPI::Datatype
ScipParaDiffSubproblemMpi::createDatatype1(
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
         startAddress,  MPI::Get_address( &localInfoIncluded )
   );
   blockLengths[0] = 1;
   displacements[0] = 0;
   types[0] = MPI::INT;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nBoundChanges )
   );
   blockLengths[1] = 1;
   displacements[1] = address - startAddress;
   types[1] = MPI::INT;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nLinearConss )
   );
   blockLengths[2] = 1;
   displacements[2] = address - startAddress;
   types[2] = MPI::INT;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nVarBranchStats )
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
ScipParaDiffSubproblemMpi::createDatatype2(
      bool memAllocNecessary
      )
{
   assert( nBoundChanges != 0 || nLinearConss != 0 );

   int nBlocks = 0;

   MPI::Datatype datatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   int blockLengths[17];            // reserve maximum number of elements
   MPI::Aint displacements[17];     // reserve maximum number of elements
   MPI::Datatype types[17];         // reserve maximum number of elements

   if( nBoundChanges )
   {
      if( memAllocNecessary )
      {
         indicesAmongSolvers = new int[nBoundChanges];
         branchBounds = new SCIP_Real[nBoundChanges];
         boundTypes = new SCIP_BOUNDTYPE[nBoundChanges];
      }

      MPI_CALL_WITH_RET_VAL(
            startAddress,  MPI::Get_address( indicesAmongSolvers )
      );
      displacements[nBlocks] = 0;
      blockLengths[nBlocks] = nBoundChanges;
      types[nBlocks] = MPI::INT;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( branchBounds )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nBoundChanges;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( boundTypes )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nBoundChanges;
      types[nBlocks] = MPI::INT;
      nBlocks++;
   }

   if( nLinearConss )
   {
      if( memAllocNecessary )
      {
         linearLhss = new SCIP_Real[nLinearConss];
         linearRhss = new SCIP_Real[nLinearConss];
         nLinearCoefs = new int[nLinearConss];
      }
      if( nBlocks == 0 )
      {
         MPI_CALL_WITH_RET_VAL(
               startAddress,  MPI::Get_address( linearLhss )
         );
         displacements[nBlocks] = 0;
         blockLengths[nBlocks] = nLinearConss;
         types[nBlocks] = MPI::DOUBLE;
         nBlocks++;
      }
      else
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( linearLhss )
         );
         displacements[nBlocks] = address - startAddress;;
         blockLengths[nBlocks] = nLinearConss;
         types[nBlocks] = MPI::DOUBLE;
         nBlocks++;
      }

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( linearRhss )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nLinearConss;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nLinearCoefs )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nLinearConss;
      types[nBlocks] = MPI::INT;
      nBlocks++;
   }

   if( nVarBranchStats )
   {
      if( memAllocNecessary )
      {
         idxLBranchStatsVars = new int[nVarBranchStats];
         downpscost = new SCIP_Real[nVarBranchStats];
         uppscost = new SCIP_Real[nVarBranchStats];
         downvsids = new SCIP_Real[nVarBranchStats];
         upvsids = new SCIP_Real[nVarBranchStats];
         downconflen = new SCIP_Real[nVarBranchStats];
         upconflen = new SCIP_Real[nVarBranchStats];
         downinfer = new SCIP_Real[nVarBranchStats];
         upinfer = new SCIP_Real[nVarBranchStats];
         downcutoff = new SCIP_Real[nVarBranchStats];
         upcutoff = new SCIP_Real[nVarBranchStats];
      }

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxLBranchStatsVars )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::INT;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downpscost )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( uppscost )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downvsids )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( upvsids )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downconflen )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( upconflen )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downinfer )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( upinfer )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( downcutoff )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( upcutoff )
      );
      displacements[nBlocks] = address - startAddress;;
      blockLengths[nBlocks] = nVarBranchStats;
      types[nBlocks] = MPI::DOUBLE;
      nBlocks++;

   }

   assert(nBlocks <= 17);

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
   );

   return datatype;
}

/** create ScipDiffSubproblemDatatype */
MPI::Datatype
ScipParaDiffSubproblemMpi::createDatatype3(
      bool memAllocNecessary
      )
{
   assert( nLinearConss != 0 );

   int nBlocks = 0;

   MPI::Datatype datatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   int blockLengths[nLinearConss*2];
   MPI::Aint displacements[nLinearConss*2];
   MPI::Datatype types[nLinearConss*2];

   if( memAllocNecessary )
   {
      linearCoefs = new SCIP_Real*[nLinearConss];
      idxLinearCoefsVars = new int*[nLinearConss];
      for(int i = 0; i < nLinearConss; i++ )
      {
         linearCoefs[i] = new SCIP_Real[nLinearCoefs[i]];
         idxLinearCoefsVars[i] = new int[nLinearCoefs[i]];
      }
   }

   for(int i = 0; i < nLinearConss; i++ )
   {
      if( i == 0 )
      {
         MPI_CALL_WITH_RET_VAL(
               startAddress,  MPI::Get_address( linearCoefs[i] )
         );
         displacements[nBlocks] = 0;
         blockLengths[nBlocks] = nLinearCoefs[i];
         types[nBlocks] = MPI::DOUBLE;
         nBlocks++;
      }
      else
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( linearCoefs[i] )
         );
         displacements[nBlocks] = address - startAddress;
         blockLengths[nBlocks] = nLinearCoefs[i];
         types[nBlocks] = MPI::DOUBLE;
         nBlocks++;
      }

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxLinearCoefsVars[i] )
      );
      displacements[nBlocks] = address - startAddress;
      blockLengths[nBlocks] = nLinearCoefs[i];
      types[nBlocks] = MPI::INT;
      nBlocks++;
   }

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
   );

   return datatype;
}

int
ScipParaDiffSubproblemMpi::bcast(ParaComm *comm, int root)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype datatype1;
   datatype1 = createDatatype1();
   datatype1.Commit();
   PARA_COMM_CALL(
      commMpi->bcast(&localInfoIncluded, 1, datatype1, root)
   );
   datatype1.Free();

   if( nBoundChanges || nLinearConss  ){
      MPI::Datatype datatype2;
      if( comm->getRank() == root )
      {
         datatype2 = createDatatype2(false);
      }
      else
      {
         datatype2 = createDatatype2(true);
      }
      datatype2.Commit();
      if( nBoundChanges )
      {
         PARA_COMM_CALL(
               commMpi->bcast(indicesAmongSolvers, 1, datatype2, root)
         );
      }
      else
      {
         PARA_COMM_CALL(
               commMpi->bcast(linearLhss, 1, datatype2, root)
         );
      }
      datatype2.Free();

      if( nLinearConss )
      {
         MPI::Datatype datatype3;
         if( comm->getRank() == root )
         {
            datatype3 = createDatatype3(false);
         }
         else
         {
            datatype3 = createDatatype3(true);
         }
         datatype3.Commit();
         PARA_COMM_CALL(
               commMpi->bcast(linearCoefs[0], 1, datatype3, root)
         );
         datatype3.Free();
      }
   }
   return 0;
}

int
ScipParaDiffSubproblemMpi::send(ParaComm *comm, int dest)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype datatype1;
   datatype1 = createDatatype1();
   datatype1.Commit();
   PARA_COMM_CALL(
      commMpi->send(&localInfoIncluded, 1, datatype1, dest, TagDiffSubproblem)
   );
   datatype1.Free();

   if( nBoundChanges || nLinearConss ){
      MPI::Datatype datatype2;
      datatype2 = createDatatype2(false);
      datatype2.Commit();
      if( nBoundChanges )
      {
         PARA_COMM_CALL(
               commMpi->send(indicesAmongSolvers, 1, datatype2, dest, TagDiffSubproblem1)
         );
      }
      else
      {
         PARA_COMM_CALL(
               commMpi->send(linearLhss, 1, datatype2, dest, TagDiffSubproblem1)
         );
      }
      datatype2.Free();
      if( nLinearConss )
      {
         MPI::Datatype datatype3;
         datatype3 = createDatatype3(false);
         datatype3.Commit();
         PARA_COMM_CALL(
               commMpi->send(linearCoefs[0], 1, datatype3, dest, TagDiffSubproblem2)
         );
         datatype3.Free();
      }
   }
   return 0;
}

int
ScipParaDiffSubproblemMpi::receive(ParaComm *comm, int source)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype datatype1;
   datatype1 = createDatatype1();
   datatype1.Commit();
   PARA_COMM_CALL(
      commMpi->receive(&localInfoIncluded, 1, datatype1, source, TagDiffSubproblem)
   );
   datatype1.Free();

   if( nBoundChanges || nLinearConss )
   {
      MPI::Datatype datatype2;
      datatype2 = createDatatype2(true);
      datatype2.Commit();
      if( nBoundChanges )
      {
         PARA_COMM_CALL(
               commMpi->receive(indicesAmongSolvers, 1, datatype2, source, TagDiffSubproblem1)
         );
      }
      else
      {
         PARA_COMM_CALL(
               commMpi->receive(linearLhss, 1, datatype2, source, TagDiffSubproblem1)
         );
      }
      datatype2.Free();
      if( nLinearConss )
      {
         MPI::Datatype datatype3;
         datatype3 = createDatatype3(true);
         datatype3.Commit();
         PARA_COMM_CALL(
               commMpi->receive(linearCoefs[0], 1, datatype3, source, TagDiffSubproblem2)
         );
         datatype3.Free();
      }
   }
   return 0;
}

/** create clone of this object */
ScipParaDiffSubproblemMpi *
ScipParaDiffSubproblemMpi::clone(ParaComm *comm)
{
   return( new ScipParaDiffSubproblemMpi(this) );

}
