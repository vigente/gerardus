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

/**@file    scipParaInstanceMpi.cpp
 * @brief   ScipParaInstance extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scipParaCommMpi.h"
#include "scipParaInstanceMpi.h"

using namespace UG;
using namespace ParaSCIP;

/** create ScipInstancePrePreDatatype */
MPI::Datatype
ScipParaInstanceMpi::createDatatype1(
      )
{
   const int nBlocks = 19;
   MPI::Datatype prePreDatatype;

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
         startAddress,  MPI::Get_address( &lProbName )
   );
   displacements[0] = 0;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &origObjSense )
   );
   displacements[1] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &objScale )
   );
   displacements[2] = address - startAddress;
   types[2] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &objOffset )
   );
   displacements[3] = address - startAddress;
   types[3] = MPI::DOUBLE;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nVars )
   );
   displacements[4] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &lVarNames )
   );
   displacements[5] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nConss )
   );
   displacements[6] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &lConsNames )
   );
   displacements[7] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nLinearConss )
   );
   displacements[8] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nSetppcConss )
   );
   displacements[9] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nLogicorConss )
   );
   displacements[10] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nKnapsackConss )
   );
   displacements[11] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nVarboundConss )
   );
   displacements[12] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nVarBoundDisjunctionConss )
   );
   displacements[13] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nSos1Conss )
   );
   displacements[14] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nSos2Conss )
   );
   displacements[15] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nAggregatedConss )
   );
   displacements[16] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &lAggregatedVarNames )
   );
   displacements[17] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &lAggregatedConsNames )
   );
   displacements[18] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         prePreDatatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return prePreDatatype;
}

void
ScipParaInstanceMpi::allocateMemoryForDatatype2(
      )
{
   if( !lProbName )  THROW_LOGICAL_ERROR1("No problem name");
   probName = new char[lProbName+1];
   if( nVars )
   {
      varLbs = new SCIP_Real[nVars];
      varUbs = new SCIP_Real[nVars];
      objCoefs = new SCIP_Real[nVars];
      varTypes = new int[nVars];
      if(lVarNames) varNames = new char[lVarNames];
      posVarNames = new int[nVars];
   }
   if( nConss )
   {
      if( lConsNames ) consNames = new char[lConsNames];
      posConsNames = new int[nConss];
   }
   if( nLinearConss )
   {
      idxLinearConsNames = new int[nLinearConss];
      linearLhss = new SCIP_Real[nLinearConss];
      linearRhss = new SCIP_Real[nLinearConss];
      nLinearCoefs = new int[nLinearConss];
   }
   if( nSetppcConss )
   {
      idxSetppcConsNames = new int[nSetppcConss];
      nIdxSetppcVars = new int[nSetppcConss];
      setppcTypes = new int[nSetppcConss];
   }
   if( nLogicorConss )
   {
      idxLogicorConsNames = new int[nLogicorConss];
      nIdxLogicorVars = new int[nLogicorConss];
   }
   if( nKnapsackConss )
   {
      idxKnapsackConsNames = new int[nKnapsackConss];
      capacities = new SCIP_Longint[nKnapsackConss];
      nLKnapsackCoefs = new int[nKnapsackConss];
   }
   if( nVarboundConss )
   {
      idxVarboundConsNames = new int[nVarboundConss];
      varboundLhss = new SCIP_Real[nVarboundConss];
      varboundRhss = new SCIP_Real[nVarboundConss];
      idxVarboundCoefVar1s = new int[nVarboundConss];
      varboundCoef2s = new SCIP_Real[nVarboundConss];
      idxVarboundCoefVar2s = new int[nVarboundConss];
   }
   if( nVarBoundDisjunctionConss ){
      idxBoundDisjunctionConsNames = new int[nVarBoundDisjunctionConss];
      nVarsBoundDisjunction = new int[nVarBoundDisjunctionConss];

   }
   if( nSos1Conss )
   {
      idxSos1ConsNames = new int[nSos1Conss];
      nSos1Coefs = new int[nSos1Conss];
   }
   if( nSos2Conss )
   {
      idxSos2ConsNames = new int[nSos2Conss];
      nSos2Coefs = new int[nSos2Conss];
   }
   if( nAggregatedConss )
   {
      if( lAggregatedVarNames ) aggregatedVarNames = new char[lAggregatedVarNames];
      posAggregatedVarNames = new int[nAggregatedConss];
      aggregatedConsNames = new char[lAggregatedConsNames];
      posAggregatedConsNames = new int[nAggregatedConss];
      if( nAggregatedConss ) aggregatedLhsAndLhss = new SCIP_Real[nAggregatedConss];
      nAggregatedCoefs = new int [nAggregatedConss];
   }
}

/** create ScipInstancePreDatatype */
MPI::Datatype
ScipParaInstanceMpi::createDatatype2(
      bool memAllocNecessary
      )
{
   const int nBlocks = 37;
   MPI::Datatype datatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   if( memAllocNecessary )
   {
      allocateMemoryForDatatype2();
   }

   int blockLengths[nBlocks];
   MPI::Aint displacements[nBlocks];
   MPI::Datatype types[nBlocks];

   int n = 0;

   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( probName )
   );
   displacements[n] = 0;
   blockLengths[n] = lProbName + 1;
   types[n++] = MPI::CHAR;

   if( nVars )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( varLbs )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVars;
      types[n++] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( varUbs )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVars;
      types[n++] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( objCoefs )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVars;
      types[n++] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( varTypes )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVars;
      types[n++] = MPI::INT;

      if( lVarNames )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( varNames )
         );
         displacements[n] =  address - startAddress;
         blockLengths[n] = lVarNames;
         types[n++] = MPI::CHAR;
      }

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( posVarNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVars;
      types[n++] = MPI::INT;
   }

   if( nConss )
   {
      if( lConsNames )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( consNames )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = lConsNames;
         types[n++] = MPI::CHAR;
      }
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( posConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nConss;
      types[n++] = MPI::INT;
   }

   if( nLinearConss )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxLinearConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nLinearConss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( linearLhss )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nLinearConss;
      types[n++] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( linearRhss )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nLinearConss;
      types[n++] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nLinearCoefs )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nLinearConss;
      types[n++] = MPI::INT;
   }

   if( nSetppcConss )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxSetppcConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nSetppcConss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nIdxSetppcVars )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nSetppcConss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( setppcTypes )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nSetppcConss;
      types[n++] = MPI::INT;
   }

   if( nLogicorConss )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxLogicorConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nLogicorConss;
      types[n++] = MPI::INT;
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nIdxLogicorVars )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nLogicorConss;
      types[n++] = MPI::INT;
   }

   if( nKnapsackConss )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxKnapsackConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nKnapsackConss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( capacities )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nKnapsackConss;
#ifdef _ALIBABA
      types[n++] = MPI::LONG;
#else
      types[n++] = MPI::LONG_LONG;
#endif

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nLKnapsackCoefs )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nKnapsackConss;
      types[n++] = MPI::INT;
   }

   if( nVarboundConss )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxVarboundConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVarboundConss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( varboundLhss )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVarboundConss;
      types[n++] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( varboundRhss )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVarboundConss;
      types[n++] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxVarboundCoefVar1s )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVarboundConss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( varboundCoef2s )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVarboundConss;
      types[n++] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxVarboundCoefVar2s )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVarboundConss;
      types[n++] = MPI::INT;
   }

   if( nVarBoundDisjunctionConss )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxBoundDisjunctionConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVarBoundDisjunctionConss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nVarsBoundDisjunction )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nVarBoundDisjunctionConss;
      types[n++] = MPI::INT;
   }

   if( nSos1Conss )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxSos1ConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nSos1Conss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nSos1Coefs )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nSos1Conss;
      types[n++] = MPI::INT;
   }

   if( nSos2Conss )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( idxSos2ConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nSos2Conss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( nSos2Coefs )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nSos2Conss;
      types[n++] = MPI::INT;
   }

   if( nAggregatedConss )
   {
      if( lAggregatedVarNames )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( aggregatedVarNames )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = lAggregatedVarNames;
         types[n++] = MPI::CHAR;
      }

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( posAggregatedVarNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nAggregatedConss;
      types[n++] = MPI::INT;

      if( lAggregatedConsNames )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( aggregatedConsNames )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = lAggregatedConsNames;
         types[n++] = MPI::CHAR;
      }

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( posAggregatedConsNames )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nAggregatedConss;
      types[n++] = MPI::INT;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( aggregatedLhsAndLhss )
      );
      displacements[n] = address - startAddress;
      blockLengths[n] = nAggregatedConss;
      types[n++] = MPI::DOUBLE;

      MPI_CALL_WITH_RET_VAL(
             address,  MPI::Get_address( nAggregatedCoefs )
       );
       displacements[n] = address - startAddress;
       blockLengths[n] = nAggregatedConss;
       types[n++] = MPI::INT;
   }

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(n, blockLengths, displacements, types)
         );

   return datatype;
}

void
ScipParaInstanceMpi::allocateMemoryForDatatype3(
      )
{
   if( nLinearConss )
   {
      linearCoefs = new SCIP_Real*[nLinearConss];
      idxLinearCoefsVars = new int*[nLinearConss];
      for(int i = 0; i < nLinearConss; i++ )
      {
         linearCoefs[i] = new  SCIP_Real[nLinearCoefs[i]];
         idxLinearCoefsVars[i] = new int[nLinearCoefs[i]];
      }
   }
   if( nSetppcConss )
   {
      idxSetppcVars = new int*[nSetppcConss];
      for(int i = 0; i < nSetppcConss; i++ )
      {
         idxSetppcVars[i] = new int[nIdxSetppcVars[i]];
      }
   }
   if( nLogicorConss )
   {
      idxLogicorVars = new int*[nLogicorConss];
      for( int i = 0; i < nLogicorConss; i++ )
      {
         idxLogicorVars[i] = new int[nIdxLogicorVars[i]];
      }
   }
   if( nKnapsackConss )
   {
      knapsackCoefs = new SCIP_Longint*[nKnapsackConss];
      idxKnapsackCoefsVars = new int*[nKnapsackConss];
      for( int i = 0; i < nKnapsackConss; i++ )
      {
         knapsackCoefs[i] = new SCIP_Longint[nLKnapsackCoefs[i]];
         idxKnapsackCoefsVars[i] = new int[nLKnapsackCoefs[i]];
      }
   }
   if( nVarBoundDisjunctionConss )
   {
      idxVarBoundDisjunction = new int*[nVarBoundDisjunctionConss];
      boundTypesBoundDisjunction = new SCIP_BOUNDTYPE*[nVarBoundDisjunctionConss];
      boundsBoundDisjunction = new SCIP_Real*[nVarBoundDisjunctionConss];
      for( int i = 0; i < nVarBoundDisjunctionConss; i++ )
      {
         idxVarBoundDisjunction[i] = new int[nVarsBoundDisjunction[i]];
         boundTypesBoundDisjunction[i] = new SCIP_BOUNDTYPE[nVarsBoundDisjunction[i]];
         boundsBoundDisjunction[i] = new SCIP_Real[nVarsBoundDisjunction[i]];
      }
   }
   if( nSos1Conss )
   {
      sos1Coefs = new SCIP_Real*[nSos1Conss];
      idxSos1CoefsVars = new int*[nSos1Conss];
      for( int i = 0; i < nSos1Conss; i++ )
      {
         sos1Coefs[i] = new SCIP_Real[nSos1Coefs[i]];
         idxSos1CoefsVars[i] = new int[nSos1Coefs[i]];
      }
   }
   if( nSos2Conss )
   {
      sos2Coefs = new SCIP_Real*[nSos2Conss];
      idxSos2CoefsVars = new int*[nSos2Conss];
      for( int i = 0; i < nSos1Conss; i++ )
      {
         sos2Coefs[i] = new SCIP_Real[nSos2Coefs[i]];
         idxSos2CoefsVars[i] = new int[nSos2Coefs[i]];
      }
   }
   if( nAggregatedConss )
   {
      aggregatedCoefs = new SCIP_Real*[nAggregatedConss];
      idxAggregatedCoefsVars = new int*[nAggregatedConss];
      for( int i = 0; i < nAggregatedConss; i++ )
      {
         aggregatedCoefs[i] = new SCIP_Real[nAggregatedCoefs[i]];
         idxAggregatedCoefsVars[i] = new int[nAggregatedCoefs[i]];
      }
   }
}

/** create ScipInstanceDatatype */
MPI::Datatype
ScipParaInstanceMpi::createDatatype3(
      bool memAllocNecessary
      )
{
   MPI::Datatype datatype;

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   if( memAllocNecessary )
   {
      allocateMemoryForDatatype3();
   }

   int nArrays = 1 + 2*nLinearConss + nSetppcConss + nLogicorConss + 2*nKnapsackConss + 3*nVarBoundDisjunctionConss + 2*nSos1Conss + 2*nSos2Conss + 2*nAggregatedConss;
   int *blockLengths = new int[nArrays];
   MPI::Aint *displacements = new MPI::Aint[nArrays];
   MPI::Datatype *types = new MPI::Datatype[nArrays];

   int n = 0;

   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( &dummyToKeepStartPos )
   );
   displacements[n] = 0;
   blockLengths[n] = 1;
   types[n++] = MPI::INT;

   if( nLinearConss )
   {
      for(int i = 0; i <  nLinearConss; i++ )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( linearCoefs[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nLinearCoefs[i];
         types[n++] = MPI::DOUBLE;
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( idxLinearCoefsVars[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nLinearCoefs[i];
         types[n++] = MPI::INT;
      }
   }

   if( nSetppcConss )
   {
      for( int i = 0; i < nSetppcConss; i++ )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( idxSetppcVars[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nIdxSetppcVars[i];
         types[n++] = MPI::INT;
      }
   }

   if( nLogicorConss )
   {
      for( int i = 0; i < nLogicorConss; i++ )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( idxLogicorVars[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nIdxLogicorVars[i];
         types[n++] = MPI::INT;
      }
   }
   if( nKnapsackConss )
   {
      for( int i = 0; i < nKnapsackConss; i++ )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( knapsackCoefs[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nLKnapsackCoefs[i];
#ifdef _ALIBABA
         types[n++] = MPI::LONG;
#else
         types[n++] = MPI::LONG_LONG;
#endif
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( idxKnapsackCoefsVars[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nLKnapsackCoefs[i];
         types[n++] = MPI::INT;
      }
   }
   if( nVarBoundDisjunctionConss )
   {
      for( int i = 0; i < nVarBoundDisjunctionConss; i++ )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( idxVarBoundDisjunction[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nVarsBoundDisjunction[i];
         types[n++] = MPI::INT;
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( boundTypesBoundDisjunction[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nVarsBoundDisjunction[i];
         types[n++] = MPI::INT;
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( boundsBoundDisjunction[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nVarsBoundDisjunction[i];
         types[n++] = MPI::DOUBLE;
      }
   }
   if( nSos1Conss )
   {
      for( int i = 0; i < nSos1Conss; i++ )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( sos1Coefs[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nSos1Coefs[i];
         types[n++] = MPI::DOUBLE;
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( idxSos1CoefsVars[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nSos1Coefs[i];
         types[n++] = MPI::INT;
      }
   }
   if( nSos2Conss )
   {
      for( int i = 0; i < nSos1Conss; i++ )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( sos2Coefs[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nSos2Coefs[i];
         types[n++] = MPI::DOUBLE;
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( idxSos2CoefsVars[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nSos2Coefs[i];
         types[n++] = MPI::INT;
      }
   }
   if( nAggregatedConss )
   {
      for( int i = 0; i < nAggregatedConss; i++ )
      {
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( aggregatedCoefs[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nAggregatedCoefs[i];
         types[n++] = MPI::DOUBLE;
         MPI_CALL_WITH_RET_VAL(
               address,  MPI::Get_address( idxAggregatedCoefsVars[i] )
         );
         displacements[n] = address - startAddress;
         blockLengths[n] = nAggregatedCoefs[i];
         types[n++] = MPI::INT;
      }
   }

   assert(n == nArrays);

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(n, blockLengths, displacements, types)
         );

   delete [] blockLengths;
   delete [] displacements;
   delete [] types;

   return datatype;
}

int
ScipParaInstanceMpi::bcast(
      ParaComm *comm,
      int root,
      int method
      )
{
   DEF_PARA_COMM( commMpi, comm);

   switch ( method )
   {
   case 0 :
   {
      MPI::Datatype datatype = createDatatype1();
      datatype.Commit();
      PARA_COMM_CALL(
         commMpi->bcast(&lProbName, 1, datatype, root)
      );
      datatype.Free();

      if( commMpi->getRank() == root )
      {
         datatype = createDatatype2(false);
      }
      else
      {
         datatype = createDatatype2(true);
      }
      datatype.Commit();

      PARA_COMM_CALL(
         commMpi->bcast(probName, 1, datatype, root)
      );
      datatype.Free();

      if( commMpi->getRank() == root )
      {
         datatype = createDatatype3(false);
      }
      else
      {
         datatype = createDatatype3(true);
      }
      datatype.Commit();
      PARA_COMM_CALL(
         commMpi->bcast(&dummyToKeepStartPos, 1, datatype, root)
      );
      datatype.Free();
      break;
   }
   case 1:
   case 2:
   {
      MPI::Datatype datatype = createDatatype1();
      datatype.Commit();
      PARA_COMM_CALL(
         commMpi->bcast(&lProbName, 1, datatype, root)
      );
      datatype.Free();
      break;
   }
   default:
      THROW_LOGICAL_ERROR1("Undefined instance transfer method");
   }

   return 0;

}
