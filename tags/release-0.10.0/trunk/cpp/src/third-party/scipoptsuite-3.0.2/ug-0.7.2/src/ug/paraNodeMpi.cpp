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

/**@file    paraNodeMpi.cpp
 * @brief   ParaNode extension for MIP communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <mpi.h>
#include "paraNodeMpi.h"

using namespace UG;

MPI::Datatype
ParaNodeMpi::createDatatype(
      )
{
   const int nBlocks = 15;

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
         startAddress,  MPI::Get_address( &nodeId.subtreeId.lcId )
   );
   displacements[0] = 0;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nodeId.subtreeId.globalSubtreeIdInLc )
   );
   displacements[1] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nodeId.subtreeId.solverId )
   );
   displacements[2] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nodeId.seqNum )
   );
   displacements[3] = address - startAddress;
#ifdef _ALIBABA
   types[3] = MPI::LONG;
#else
   types[3] = MPI::LONG_LONG;
#endif

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &generatorNodeId.subtreeId.lcId )
   );
   displacements[4] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &generatorNodeId.subtreeId.globalSubtreeIdInLc )
   );
   displacements[5] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &generatorNodeId.subtreeId.solverId )
   );
   displacements[6] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &generatorNodeId.seqNum )
   );
   displacements[7] = address - startAddress;
#ifdef _ALIBABA
   types[7] = MPI::LONG;
#else
   types[7] = MPI::LONG_LONG;
#endif

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &depth )
   );
   displacements[8] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &dualBoundValue )
   );
   displacements[9] = address - startAddress;
   types[9] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &initialDualBoundValue )
   );
   displacements[10] = address - startAddress;
   types[10] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &estimatedValue )
   );
   displacements[11] = address - startAddress;
   types[11] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &diffSubproblemInfo )
   );
   displacements[12] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &basisInfo )
   );
   displacements[13] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
          address,  MPI::Get_address( &mergingStatus )
   );
   displacements[14] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return datatype;
}

int
ParaNodeMpi::bcast(
      ParaComm *comm,
      int root
      )
{
    DEF_PARA_COMM( commMpi, comm);

    MPI::Datatype datatype;
    datatype = createDatatype();
    datatype.Commit();
    PARA_COMM_CALL(
       commMpi->bcast(&nodeId.subtreeId.lcId, 1, datatype, root)
    );
    datatype.Free();

   // root node does not have diffSubproblem
   if( diffSubproblemInfo ) diffSubproblem->bcast(commMpi, root);
   return 0;
}

int
ParaNodeMpi::send(
      ParaComm *comm,
      int destination
      )
{
    DEF_PARA_COMM( commMpi, comm);

    MPI::Datatype datatype;
    datatype = createDatatype();
    datatype.Commit();
    PARA_COMM_CALL(
       commMpi->send(&nodeId.subtreeId.lcId, 1, datatype, destination, TagNode)
    );
    datatype.Free();
   // root node does not have diffSubproblem
   if( diffSubproblemInfo ) diffSubproblem->send(commMpi, destination);
   return 0;
}

int
ParaNodeMpi::receive(ParaComm *comm, int source){
   DEF_PARA_COMM( commMpi, comm);

   MPI::Datatype datatype;
   datatype = createDatatype();
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->receive(&nodeId.subtreeId.lcId, 1, datatype, source, TagNode)
   );
   datatype.Free();

   if( diffSubproblemInfo )
   {
      diffSubproblem = commMpi->createParaDiffSubproblem();
      diffSubproblem->receive(commMpi, source);
   }

   return 0;
}
