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

/**@file    scipDiffParamSetMpi.cpp
 * @brief   ScipDiffParamSet extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <string.h>
#include <cassert>
#include "scip/scip.h"
#include "ug/paraCommMpiWorld.h"
#include "scipDiffParamSetMpi.h"

using namespace UG;
using namespace ParaSCIP;

/** create scipDiffParamSetPreType */
MPI::Datatype
ScipDiffParamSetMpi::createDatatype1(
      )
{
   MPI::Datatype datatype;

   int blockLengthsPre[13];
   MPI::Aint displacementsPre[13];
   MPI::Datatype typesPre[13];

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   for( int i = 0; i < 13; i++ ){
       blockLengthsPre[i] = 1;
       typesPre[i] = MPI::INT;
   }

   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( &numBoolParams )
   );
   displacementsPre[0] = 0;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &boolParamNamesSize )
   );
   displacementsPre[1] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &numIntParams )
   );
   displacementsPre[2] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &intParamNamesSize )
   );
   displacementsPre[3] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &numLongintParams )
   );
   displacementsPre[4] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &longintParamNamesSize )
   );
   displacementsPre[5] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &numRealParams )
   );
   displacementsPre[6] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &realParamNamesSize )
   );
   displacementsPre[7] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &numCharParams )
   );
   displacementsPre[8] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &charParamNamesSize )
   );
   displacementsPre[9] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &numStringParams )
   );
   displacementsPre[10] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &stringParamNamesSize )
   );
   displacementsPre[11] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &stringParamValuesSize )
   );
   displacementsPre[12] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(13, blockLengthsPre, displacementsPre, typesPre)
         );

   return datatype;

}

/** create scipDiffParamSetType */
MPI::Datatype
ScipDiffParamSetMpi::createDatatype2(
      bool memAllocNecessary
      )
{
   MPI::Datatype datatype;

   int blockLengths[12];
   MPI::Aint displacements[12];
   MPI::Datatype types[12];

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   if( memAllocNecessary ){
      allocateMemoty();
   }

   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( boolParamNames )
   );
   blockLengths[0] = boolParamNamesSize;
   displacements[0] = 0;
   types[0] = MPI::CHAR;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( boolParamValues )
   );
   blockLengths[1] = numBoolParams;
   displacements[1] = address - startAddress;
   types[1] = MPI::UNSIGNED;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( intParamNames )
   );
   blockLengths[2] = intParamNamesSize;
   displacements[2] = address - startAddress;
   types[2] = MPI::CHAR;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( intParamValues )
   );
   blockLengths[3] = numIntParams;
   displacements[3] = address - startAddress;
   types[3] = MPI::INT;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( longintParamNames )
   );
   blockLengths[4] = longintParamNamesSize;
   displacements[4] = address - startAddress;
   types[4] = MPI::CHAR;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( longintParamValues )
   );
   blockLengths[5] = numLongintParams;
   displacements[5] = address - startAddress;
#ifdef _ALIBABA
   types[5] = MPI::LONG;
#else
   types[5] = MPI::LONG_LONG;
#endif
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( realParamNames )
   );
   blockLengths[6] = realParamNamesSize;
   displacements[6] = address - startAddress;
   types[6] = MPI::CHAR;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( realParamValues )
   );
   blockLengths[7] = numRealParams;
   displacements[7] = address - startAddress;
   types[7] = MPI::DOUBLE;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( charParamNames )
   );
   blockLengths[8] = charParamNamesSize;
   displacements[8] = address - startAddress;
   types[8] = MPI::CHAR;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( charParamValues )
   );
   blockLengths[9] = numCharParams;
   displacements[9] = address - startAddress;
   types[9] = MPI::CHAR;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( stringParamNames )
   );
   blockLengths[10] = stringParamNamesSize;
   displacements[10] = address - startAddress;
   types[10] = MPI::CHAR;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( stringParamValues )
   );
   blockLengths[11] = stringParamValuesSize;
   displacements[11] = address - startAddress;
   types[11] = MPI::CHAR;

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(12, blockLengths, displacements, types)
         );

   return datatype;
}

/** send solution data to the rank */
int
ScipDiffParamSetMpi::bcast(
      ParaComm *comm,
      int root
      )
{
   DEF_PARA_COMM( commMpi, comm);

   MPI::Datatype datatype = createDatatype1();
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->bcast(&numBoolParams, 1, datatype, root)
   );
   datatype.Free();

   if( comm->getRank() == root )
   {
      datatype = createDatatype2(false);
   }
   else
   {
      datatype = createDatatype2(true);
   }
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->bcast(boolParamNames, 1, datatype, root)
   );
   datatype.Free();
   return 0;
}

/** send solution data to the rank */
int
ScipDiffParamSetMpi::send(
      ParaComm *comm,
      int dest
      )
{
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype datatype = createDatatype1();
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->send(&numBoolParams, 1, datatype, dest, TagSolverDiffParamSet1)
   );
   datatype.Free();

   datatype = createDatatype2(false);
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->send(boolParamNames, 1, datatype, dest, TagSolverDiffParamSet2)
   );
   datatype.Free();
   return 0;
}

 /** receive solution data from the source rank */
int
ScipDiffParamSetMpi::receive(
       ParaComm *comm,
       int source
       )
 {
   DEF_PARA_COMM( commMpi, comm);
   MPI::Datatype datatype = createDatatype1();
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->receive(&numBoolParams, 1, datatype, source, TagSolverDiffParamSet1)
   );
   datatype.Free();

   datatype = createDatatype2(true);
   datatype.Commit();
   PARA_COMM_CALL(
       commMpi->receive(boolParamNames, 1, datatype, source, TagSolverDiffParamSet2)
   );
   datatype.Free();
   return 0;
 }
