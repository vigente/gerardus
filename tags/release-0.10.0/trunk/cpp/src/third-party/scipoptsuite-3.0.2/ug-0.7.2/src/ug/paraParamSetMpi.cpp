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

/**@file    paraParamSetMpi.cpp
 * @brief   ParaParamSet extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "paraCommMpiWorld.h"
#include "paraParamSetMpi.h"
#include <cstring>

using namespace UG;

/** allocate memory for transfer */
void
ParaParamSetMpi::allocateMemory(
      )
{
   if( ParaParamsBoolN )
   {
      boolParams = new int[ParaParamsBoolN];
      boolParamValues = new char[ParaParamsBoolN];
   }
   if( ParaParamsIntN )
   {
      intParams = new int[ParaParamsIntN];
      intParamValues = new int[ParaParamsIntN];
   }
   if( ParaParamsLongintN )
   {
      longintParams = new int[ParaParamsLongintN];
      longintParamValues = new long long[ParaParamsLongintN];
   }
   if( ParaParamsRealN )
   {
      realParams = new int[ParaParamsRealN];
      realParamValues = new double[ParaParamsRealN];
   }
   if( ParaParamsCharN )
   {
      charParams = new int[ParaParamsCharN];
      charParamValues = new char[ParaParamsCharN];
   }
   if( ParaParamsStringN )
   {
      int allStringParamValusesSize = 0;
      stringParams = new int[ParaParamsStringN];
      for( int i = ParaParamsStringFirst; i <= ParaParamsStringLast; i++ )
      {
         ParaParamString *paraParamString = dynamic_cast< ParaParamString * >(paraParams[i]);
         if( !paraParamString->isDefaultValue() )
         {
            allStringParamValusesSize += (std::strlen(paraParamString->getValue()) + 1);
         }
      }
      stringParamValues = new char[allStringParamValusesSize];
   }
}

/** free memory for transfer */
void
ParaParamSetMpi::freeMemory(
      )
{
   if( ParaParamsBoolN )
   {
      delete [] boolParams;
      delete [] boolParamValues;
   }
   if( ParaParamsIntN )
   {
      delete [] intParams;
      delete [] intParamValues;
   }
   if( ParaParamsLongintN )
   {
      delete [] longintParams;
      delete [] longintParamValues;
   }
   if( ParaParamsRealN )
   {
      delete [] realParams;
      delete [] realParamValues;
   }
   if( ParaParamsCharN )
   {
      delete [] charParams;
      delete [] charParamValues;
   }
   if( ParaParamsStringN )
   {
      delete [] stringParams;
      delete [] stringParamValues;
   }
}


/** constructor with scip */
void
ParaParamSetMpi::createDiffParams(
      )
{
   /** compute param name size for each type */
   if ( ParaParamsBoolN )
   {
      nBoolParams = 0;
      for( int i = ParaParamsBoolFirst; i <= ParaParamsBoolLast; i++ )
      {
         ParaParamBool *paraParamBool = dynamic_cast< ParaParamBool * >(paraParams[i]);
         if( !paraParamBool->isDefaultValue() )
         {
            boolParams[nBoolParams] = i;
            if( paraParamBool->getValue() )
            {
               boolParamValues[nBoolParams] = 'T';
            }
            else
            {
               boolParamValues[nBoolParams] = 'F';
            }
            nBoolParams++;
         }
      }
   }
   if ( ParaParamsIntN )
   {
      nIntParams = 0;
      for( int i = ParaParamsIntFirst; i <= ParaParamsIntLast; i++ )
      {
         ParaParamInt *paraParamInt = dynamic_cast< ParaParamInt * >(paraParams[i]);
         if( !paraParamInt->isDefaultValue() )
         {
            intParams[nIntParams] = i;
            intParamValues[nIntParams] = paraParamInt->getValue();
            nIntParams++;
         }
      }
   }
   if ( ParaParamsLongintN )
   {
      nLongintParams = 0;
      for( int i = ParaParamsLongintFirst; i <= ParaParamsLongintLast; i++ )
      {
         ParaParamLongint *paraParamLongint = dynamic_cast< ParaParamLongint * >(paraParams[i]);
         if( !paraParamLongint->isDefaultValue() )
         {
            longintParams[nLongintParams] = i;
            intParamValues[nLongintParams] = paraParamLongint->getValue();
            nLongintParams++;
         }
      }
   }
   if ( ParaParamsRealN )
   {
      nRealParams = 0;
      for( int i = ParaParamsRealFirst; i <= ParaParamsRealLast; i++ )
      {
         ParaParamReal *paraParamReal = dynamic_cast< ParaParamReal * >(paraParams[i]);
         if( !paraParamReal->isDefaultValue() )
         {
            realParams[nRealParams] = i;
            realParamValues[nRealParams] = paraParamReal->getValue();
            nRealParams++;
         }
      }
   }
   if ( ParaParamsCharN )
   {
      nCharParams = 0;
      for( int i = ParaParamsCharFirst; i <= ParaParamsCharLast; i++ )
      {
         ParaParamChar *paraParamChar = dynamic_cast< ParaParamChar * >(paraParams[i]);
         if( !paraParamChar->isDefaultValue() )
         {
            charParams[nCharParams] = i;
            charParamValues[nCharParams] = paraParamChar->getValue();
            nCharParams++;
         }
      }
   }
   if ( ParaParamsStringN )
   {
      nStringParams = 0;
      stringParamValuesSize = 0;
      int pos = 0;
      for( int i = ParaParamsStringFirst; i <= ParaParamsStringLast; i++ )
      {
         ParaParamString *paraParamString = dynamic_cast< ParaParamString * >(paraParams[i]);
         if( !paraParamString->isDefaultValue() )
         {
            stringParams[nStringParams] = i;
            std::strcpy( &stringParamValues[pos], paraParamString->getValue() );
            int len = (std::strlen(paraParamString->getValue()) + 1);
            pos += len;
            stringParamValuesSize += len;
            nStringParams++;
         }
      }
   }
}

/** set these parameter values in scip environment */
void
ParaParamSetMpi::setDiffParams(
      )
{
   if ( nBoolParams )
   {
      for( int i = 0; i < nBoolParams; i++ )
      {
         ParaParamBool *paraParamBool = dynamic_cast< ParaParamBool * >(paraParams[boolParams[i]]);
         if( boolParamValues[i] == 'T')
         {
            paraParamBool->setValue(true);
         }
         else
         {
            paraParamBool->setValue(false);
         }
      }
   }
   if ( nIntParams )
   {
      for( int i = 0; i < nIntParams; i++ )
      {
         ParaParamInt *paraParamInt = dynamic_cast< ParaParamInt * >(paraParams[intParams[i]]);
         paraParamInt->setValue(intParamValues[i]);
      }
   }
   if ( nLongintParams )
   {
      for( int i = 0; i < nLongintParams; i++ )
      {
         ParaParamLongint *paraParamLongint = dynamic_cast< ParaParamLongint * >(paraParams[longintParams[i]]);
         paraParamLongint->setValue(longintParamValues[i]);
      }
   }
   if ( nRealParams )
   {
      for( int i = 0; i < nRealParams; i++ )
      {
         ParaParamReal *paraParamReal = dynamic_cast< ParaParamReal * >(paraParams[realParams[i]]);
         paraParamReal->setValue(realParamValues[i]);
      }
   }
   if ( nIntParams )
   {
      for( int i = 0; i < nCharParams; i++ )
      {
         ParaParamChar *paraParamChar = dynamic_cast< ParaParamChar * >(paraParams[charParams[i]]);
         paraParamChar->setValue(charParamValues[i]);
      }
   }
   int pos = 0;
   if ( nStringParams )
   {
      for( int i = 0; i < nStringParams; i++ )
      {
         ParaParamString *paraParamString = dynamic_cast< ParaParamString * >(paraParams[stringParams[i]]);
         char *value = new char[std::strlen(&stringParamValues[pos]) + 1];
         std::strcpy(value, &stringParamValues[pos]);
         paraParamString->setValue(value);
         pos += ( std::strlen(&stringParamValues[pos]) + 1 );
      }
   }
}

/** create Datatype1 */
MPI::Datatype
ParaParamSetMpi::createDatatype1(
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
         startAddress,  MPI::Get_address( &nBoolParams )
   );
   displacements[0] = 0;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nIntParams )
   );
   displacements[1] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nLongintParams )
   );
   displacements[2] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nRealParams )
   );
   displacements[3] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nCharParams )
   );
   displacements[4] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &nStringParams )
   );
   displacements[5] = address - startAddress;
   MPI_CALL_WITH_RET_VAL(
         address,  MPI::Get_address( &stringParamValuesSize )
   );
   displacements[6] = address - startAddress;

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return datatype;

}

/** create Datatype2 */
MPI::Datatype
ParaParamSetMpi::createDatatype2(
      bool reallocateStringParamValues
      )
{
   const int nMaxBlocks = 13;

   MPI::Datatype datatype;

   int blockLengths[nMaxBlocks];
   MPI::Aint displacements[nMaxBlocks];
   MPI::Datatype types[nMaxBlocks];

   MPI::Aint startAddress = 0;
   MPI::Aint address = 0;

   int nBlocks = 0;

   MPI_CALL_WITH_RET_VAL(
         startAddress,  MPI::Get_address( &nBoolParams )
   );
   blockLengths[nBlocks] = 1;
   displacements[nBlocks] = 0;
   types[nBlocks] = MPI::INT;
   nBlocks++;

   if( nBoolParams )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( boolParams )
      );
      blockLengths[nBlocks] = nBoolParams;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::INT;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( boolParamValues )
      );
      blockLengths[nBlocks] = nBoolParams;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::CHAR;
      nBlocks++;
   }

   if( nIntParams )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( intParams )
      );
      blockLengths[nBlocks] = nIntParams;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::INT;
      nBlocks++;

      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( intParamValues )
      );
      blockLengths[nBlocks] = nIntParams;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::INT;
      nBlocks++;
   }

   if( nLongintParams )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( longintParams )
      );
      blockLengths[nBlocks] = nLongintParams;
      displacements[nBlocks] = address - startAddress;
#ifdef _ALIBABA
      types[nBlocks] = MPI::LONG;
#else
      types[nBlocks] = MPI::LONG_LONG;
#endif
      nBlocks++;
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( longintParamValues )
      );
      blockLengths[nBlocks] = nLongintParams;
      displacements[nBlocks] = address - startAddress;
#ifdef _ALIBABA
      types[nBlocks] = MPI::LONG;
#else
      types[nBlocks] = MPI::LONG_LONG;
#endif
      nBlocks++;
   }

   if( nRealParams )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( realParams )
      );
      blockLengths[nBlocks] = nRealParams;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::INT;
      nBlocks++;
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( realParamValues )
      );
      blockLengths[nBlocks] = nRealParams;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::DOUBLE;

      nBlocks++;
   }

   if( nCharParams )
   {
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( charParams )
      );
      blockLengths[nBlocks] = nCharParams;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::INT;
      nBlocks++;
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( charParamValues )
      );
      blockLengths[nBlocks] = nCharParams;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::CHAR;
      nBlocks++;
   }

   if( nStringParams )
   {
      if( reallocateStringParamValues )
      {
         delete[] stringParamValues;
         stringParamValues = new char[stringParamValuesSize];
      }
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( stringParams )
      );
      blockLengths[nBlocks] = nStringParams;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::INT;
      nBlocks++;
      MPI_CALL_WITH_RET_VAL(
            address,  MPI::Get_address( stringParamValues )
      );
      blockLengths[nBlocks] = stringParamValuesSize;
      displacements[nBlocks] = address - startAddress;
      types[nBlocks] = MPI::CHAR;
      nBlocks++;
   }

   MPI_CALL_WITH_RET_VAL(
         datatype, MPI::Datatype::Create_struct(nBlocks, blockLengths, displacements, types)
         );

   return datatype;
}

int
ParaParamSetMpi::bcast(ParaComm *comm, int root)
{

   DEF_PARA_COMM( commMpi, comm);

   allocateMemory();
   if( commMpi->getRank() == root )
   {
      createDiffParams();
   }

   MPI::Datatype datatype1;
   datatype1 = createDatatype1();
   datatype1.Commit();
   PARA_COMM_CALL(
      commMpi->bcast(&nBoolParams, 1, datatype1, root)
   );
   datatype1.Free();

   MPI::Datatype datatype2;
   if( commMpi->getRank() == root )
   {
      datatype2 = createDatatype2(false);
   }
   else
   {
      datatype2 = createDatatype2(true);
   }
   datatype2.Commit();
   PARA_COMM_CALL(
      commMpi->bcast(&nBoolParams, 1, datatype2, root)  // nBoolParams sending twice, to fix start point
   );
   datatype2.Free();
   if( commMpi->getRank() != root )
   {
      setDiffParams();
   }
   freeMemory();

   return 0;

}
