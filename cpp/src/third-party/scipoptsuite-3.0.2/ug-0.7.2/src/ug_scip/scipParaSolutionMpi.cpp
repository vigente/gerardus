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

/**@file    scipParaSolutionMpi.cpp
 * @brief   ScipParaSolution extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <mpi.h>
#include "scipParaCommMpi.h"
#include "scipParaSolutionMpi.h"

using namespace UG;
using namespace ParaSCIP;

/** create clone of this object */
ScipParaSolutionMpi *
ScipParaSolutionMpi::clone(ParaComm *comm)
{
   return( new ScipParaSolutionMpi(objectiveFunctionValue, nVars, indicesAmongSolvers, values));
}

/** create ScipDiffSubproblemPreDatatype */
MPI::Datatype
ScipParaSolutionMpi::createPreDatatype(
      )
{
   const int nBlocks = 2;
   MPI::Datatype preDatatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   int blockLengths[nBlocks];
   MPI::Aint displacements[nBlocks];
   MPI::Datatype types[nBlocks];

   for( int i = 0; i < nBlocks; i++ ){
       blockLengths[i] = 1;
   }

   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( &objectiveFunctionValue )
   );
   displacements[0] = 0;
   types[0] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nVars )
   );
   displacements[1] = address - startAddress;
   types[1] = MPI::INT;

   MPI_CALL_WITH_RET_VAL(
         preDatatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return preDatatype;
}

/** create ScipDiffSubproblemDatatype */
MPI::Datatype
ScipParaSolutionMpi::createDatatype(
      bool memAllocNecessary
      )
{
   const int nBlocks = 2;

   MPI::Datatype datatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   int blockLengths[nBlocks];
   MPI::Aint displacements[nBlocks];
   MPI::Datatype types[nBlocks];

   if( nVars )
   {
      if( memAllocNecessary )
      {
         indicesAmongSolvers = new int[nVars];
         values = new SCIP_Real[nVars];
      }

      MPI_CALL_WITH_RET_VAL(
            startAddress,  MPI::Get_address( indicesAmongSolvers )
      );
      displacements[0] = 0;
      blockLengths[0] = nVars;
      types[0] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( values )
      );
      displacements[1] = address - startAddress;
      blockLengths[1] = nVars;
      types[1] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
            datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
      );
   }
   return datatype;
}

/** send solution data to the rank */
void
ScipParaSolutionMpi::bcast(ParaComm *comm, int root)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype preDatatype;
   preDatatype = createPreDatatype();
   preDatatype.Commit();
   PARA_COMM_CALL(
      commMpi->bcast(&objectiveFunctionValue, 1, preDatatype, root)
   );
   preDatatype.Free();

   if( nVars ){
      MPI::Datatype datatype;
      if( comm->getRank() == root )
      {
         datatype = createDatatype(false);
      }
      else
      {
         datatype = createDatatype(true);
      }
      datatype.Commit();
      PARA_COMM_CALL(
            commMpi->bcast(indicesAmongSolvers, 1, datatype, root)
      );
      datatype.Free();
   }
}


/** send solution data to the rank */
void
ScipParaSolutionMpi::send(ParaComm *comm, int destination)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype preDatatype;
   preDatatype = createPreDatatype();
   preDatatype.Commit();
   PARA_COMM_CALL(
      commMpi->send(&objectiveFunctionValue, 1, preDatatype, destination, TagSolution)
   );
   preDatatype.Free();

   if( nVars ){
      MPI::Datatype datatype;
      datatype = createDatatype(false);
      datatype.Commit();
      PARA_COMM_CALL(
            commMpi->send(indicesAmongSolvers, 1, datatype, destination, TagSolution1)
      );
      datatype.Free();
   }
}

/** receive solution data from the source rank */
void
ScipParaSolutionMpi::receive(ParaComm *comm, int source)
{
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype preDatatype;
   preDatatype = createPreDatatype();
   preDatatype.Commit();
   PARA_COMM_CALL(
      commMpi->receive(&objectiveFunctionValue, 1, preDatatype, source, TagSolution)
   );
   preDatatype.Free();

   if( nVars ){
      MPI::Datatype datatype;
      datatype = createDatatype(true);
      datatype.Commit();
      PARA_COMM_CALL(
            commMpi->receive(indicesAmongSolvers, 1, datatype, source, TagSolution1)
      );
      datatype.Free();
   }
}
