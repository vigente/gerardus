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

/**@file    scipDiffParamSetPth.cpp
 * @brief   ScipDiffParamSet extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cstring>
#include <cassert>
#include "scip/scip.h"
#include "scipDiffParamSetPth.h"

using namespace UG;
using namespace ParaSCIP;

/** create scipDiffParamSetPreType */
ScipDiffParamSetPth *
ScipDiffParamSetPth::createDatatype(
      )
{
   return clone();
}

/** create clone */
ScipDiffParamSetPth *
ScipDiffParamSetPth::clone(
      )
{
   ScipDiffParamSetPth *newParam = new ScipDiffParamSetPth();

   newParam->numBoolParams = numBoolParams;
   newParam->boolParamNamesSize = boolParamNamesSize;
   newParam->numIntParams = numIntParams;
   newParam->intParamNamesSize = intParamNamesSize;
   newParam->numLongintParams = numLongintParams;
   newParam->longintParamNamesSize = longintParamNamesSize;
   newParam->numRealParams = numRealParams;
   newParam->realParamNamesSize = realParamNamesSize;
   newParam->numCharParams = numCharParams;
   newParam->charParamNamesSize = charParamNamesSize;
   newParam->numStringParams = numStringParams;
   newParam->stringParamNamesSize = stringParamNamesSize;
   newParam->stringParamValuesSize = stringParamValuesSize;

   newParam->allocateMemoty();

   char *cloneParamName;
   char *paramName;

   /** copy boolean parameter values */
   cloneParamName = newParam->boolParamNames;
   paramName = boolParamNames;
   for(int i = 0; i < newParam->numBoolParams; i++ )
   {
      newParam->boolParamValues[i] = boolParamValues[i];
      strcpy(cloneParamName, paramName);
      cloneParamName += strlen(cloneParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy int parameter values */
   cloneParamName = newParam->intParamNames;
   paramName = intParamNames;
   for(int i = 0; i < newParam->numIntParams; i++ )
   {
      newParam->intParamValues[i] = intParamValues[i];
      strcpy(cloneParamName, paramName);
      cloneParamName += strlen(cloneParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy longint parameter values */
   cloneParamName = newParam->longintParamNames;
   paramName = longintParamNames;
   for(int i = 0; i < newParam->numLongintParams; i++ )
   {
      newParam->longintParamValues[i] = longintParamValues[i];
      strcpy(cloneParamName, paramName);
      cloneParamName += strlen(cloneParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy real parameter values */
   cloneParamName = newParam->realParamNames;
   paramName = realParamNames;
   for(int i = 0; i < newParam->numRealParams; i++ )
   {
      newParam->realParamValues[i] = realParamValues[i];
      strcpy(cloneParamName, paramName);
      cloneParamName += strlen(cloneParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy char parameter values */
   cloneParamName = newParam->charParamNames;
   paramName = charParamNames;
   for(int i = 0; i < newParam->numCharParams; i++ )
   {
      newParam->charParamValues[i] = charParamValues[i];
      strcpy(cloneParamName, paramName);
      cloneParamName += strlen(cloneParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy string parameter values */
   char *cloneParamValue = newParam->stringParamValues;
   char *paramValue = stringParamValues;
   cloneParamName = newParam->stringParamNames;
   paramName = stringParamNames;
   for(int i = 0; i < newParam->numStringParams; i++ )
   {
      strcpy(cloneParamValue, paramValue);
      cloneParamValue += strlen(cloneParamValue) + 1;
      paramValue += strlen(paramValue) + 1;
      strcpy(cloneParamName, paramName);
      cloneParamName += strlen(cloneParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   return newParam;
}

void
ScipDiffParamSetPth::setValues(
      ScipDiffParamSetPth *from
      )
{

   numBoolParams = from->numBoolParams;
   boolParamNamesSize = from->boolParamNamesSize;
   numIntParams = from->numIntParams;
   intParamNamesSize = from->intParamNamesSize;
   numLongintParams = from->numLongintParams;
   longintParamNamesSize = from->longintParamNamesSize;
   numRealParams = from->numRealParams;
   realParamNamesSize = from->realParamNamesSize;
   numCharParams = from->numCharParams;
   charParamNamesSize = from->charParamNamesSize;
   numStringParams = from->numStringParams;
   stringParamNamesSize = from->stringParamNamesSize;
   stringParamValuesSize = from->stringParamValuesSize;

   allocateMemoty();

   char *paramName;
   char *fromParamName;

   /** copy boolean parameter values */
   paramName = boolParamNames;
   fromParamName = from->boolParamNames;
   for(int i = -0; i < from->numBoolParams; i++ )
   {
      boolParamValues[i] = from->boolParamValues[i];
      strcpy(paramName, fromParamName);
      fromParamName += strlen(fromParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy int parameter values */
   paramName = intParamNames;
   fromParamName = from->intParamNames;
   for(int i = -0; i < from->numIntParams; i++ )
   {
      intParamValues[i] = from->intParamValues[i];
      strcpy(paramName, fromParamName);
      fromParamName += strlen(fromParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy longint parameter values */
   paramName = longintParamNames;
   fromParamName = from->longintParamNames;
   for(int i = -0; i < from->numLongintParams; i++ )
   {
      longintParamValues[i] = from->longintParamValues[i];
      strcpy(paramName, fromParamName);
      fromParamName += strlen(fromParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy real parameter values */
   paramName = realParamNames;
   fromParamName = from->realParamNames;
   for(int i = -0; i < from->numRealParams; i++ )
   {
      realParamValues[i] = from->realParamValues[i];
      strcpy(paramName, fromParamName);
      fromParamName += strlen(fromParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy char parameter values */
   paramName = charParamNames;
   fromParamName = from->charParamNames;
   for(int i = -0; i < from->numCharParams; i++ )
   {
      charParamValues[i] = from->charParamValues[i];
      strcpy(paramName, fromParamName);
      fromParamName += strlen(fromParamName) + 1;
      paramName += strlen(paramName) + 1;
   }

   /** copy string parameter values */
   char *paramValue = stringParamValues;
   char *fromParamValue = from->stringParamValues;
   paramName = stringParamNames;
   fromParamName = from->stringParamNames;
   for(int i = -0; i < from->numStringParams; i++ )
   {
      strcpy(paramValue, fromParamValue);
      fromParamValue += strlen(fromParamValue) + 1;
      paramValue += strlen(paramValue) + 1;
      strcpy(paramName, fromParamName);
      fromParamName += strlen(fromParamName) + 1;
      paramName += strlen(paramName) + 1;
   }
}

/** send solution data to the rank */
int
ScipDiffParamSetPth::bcast(
      ParaComm *comm,
      int root
      )
{
   DEF_PARA_COMM( commPth, comm);

   if( commPth->getRank() == root )
   {
      for( int i = 0; i < commPth->getSize(); i++ )
      {
         if( i != root )
         {
            ScipDiffParamSetPth *sent;
            sent = createDatatype();
            PARA_COMM_CALL(
               commPth->uTypeSend((void *)sent, ParaSolverDiffParamType, i, TagSolverDiffParamSet)
            );
         }
      }
   }
   else
   {
      ScipDiffParamSetPth *received;
      PARA_COMM_CALL(
         commPth->uTypeReceive((void **)&received, ParaSolverDiffParamType, root, TagSolverDiffParamSet)
      );
      setValues(received);
      delete received;
   }

   return 0;
}

/** send solution data to the rank */
int
ScipDiffParamSetPth::send(
      ParaComm *comm,
      int dest
      )
{
   DEF_PARA_COMM( commPth, comm);

   PARA_COMM_CALL(
      commPth->uTypeSend((void *)createDatatype(), ParaSolverDiffParamType, dest, TagSolverDiffParamSet)
   );

   return 0;
}

 /** receive solution data from the source rank */
int
ScipDiffParamSetPth::receive(
       ParaComm *comm,
       int source
       )
 {
   DEF_PARA_COMM( commPth, comm);

   ScipDiffParamSetPth *received;
   PARA_COMM_CALL(
      commPth->uTypeReceive((void **)&received, ParaSolverDiffParamType, source, TagSolverDiffParamSet)
   );
   setValues(received);
   delete received;

   return 0;
 }
