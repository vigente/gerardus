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

/**@file    scipParaRacingRampUpParamSetMpi.cpp
 * @brief   ScipParaRacingRampUpParamSet extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scipParaCommMpi.h"
#include "scipDiffParamSetMpi.h"
#include "scipParaRacingRampUpParamSetMpi.h"
#include <cstring>

using namespace UG;
using namespace ParaSCIP;

/** create Datatype */
MPI::Datatype
ScipParaRacingRampUpParamSetMpi::createDatatype(
      )
{
   const int nBlocks = 7;

   MPI::Datatype datatype;

   int blockLengths[nBlocks];
   MPI::Aint displacements[nBlocks];
   MPI::Datatype types[nBlocks];

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   for( int i = 0; i < nBlocks; i++ )
   {
       blockLengths[i] = 1;
       types[i] = MPI::INT;
   }
   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( &terminationCriteria )
   );
   displacements[0] = 0;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nNodesLeft )
   );
   displacements[1] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &timeLimit )
   );
   displacements[2] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &scipRacingParamSeed )
   );
   displacements[3] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &permuteProbSeed )
   );
   displacements[4] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &generateBranchOrderSeed )
   );
   displacements[5] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &scipDiffParamSetInfo )
   );
   displacements[6] = address - startAddress;

   types[2] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return datatype;

}

int
ScipParaRacingRampUpParamSetMpi::send(
      ParaComm *comm,
      int dest)
{

   DEF_PARA_COMM( commMpi, comm);

   MPI::Datatype datatype;
   datatype = createDatatype();
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->send(&terminationCriteria, 1, datatype, dest, TagRacingRampUpParamSet)
   );
   datatype.Free();

   if( scipDiffParamSetInfo )
   {
      scipDiffParamSet->send(commMpi, dest);
   }

   return 0;

}

int
ScipParaRacingRampUpParamSetMpi::receive(
      ParaComm *comm,
      int source)
{

   DEF_PARA_COMM( commMpi, comm);

   MPI::Datatype datatype;
   datatype = createDatatype();
   datatype.Commit();
   PARA_COMM_CALL(
      commMpi->receive(&terminationCriteria, 1, datatype, source, TagRacingRampUpParamSet)
   );
   datatype.Free();

   if( scipDiffParamSetInfo )
   {
	  DEF_SCIP_PARA_COMM( scipParaComm, comm );
      scipDiffParamSet = scipParaComm->createScipDiffParamSet();
      scipDiffParamSet->receive(commMpi, source);
   }

   return 0;

}
