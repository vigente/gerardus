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

/**@file    scipDiffParamSet.cpp
 * @brief   SCIP parameter set to be transferred ( Only keep difference between default settings ).
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <string.h>
#include <cassert>
#include "scipDiffParamSet.h"

using namespace UG;
using namespace ParaSCIP;

/** allocate memory for names and values */
void
ScipDiffParamSet::allocateMemoty(
      )
{
   boolParamNames = new char[boolParamNamesSize];
   boolParamValues = new unsigned int[numBoolParams];
   intParamNames = new char[intParamNamesSize];
   intParamValues = new int[numIntParams];
   longintParamNames = new char[longintParamNamesSize];
   longintParamValues = new long long[numLongintParams];
   realParamNames = new char[realParamNamesSize];
   realParamValues = new double[numRealParams];
   charParamNames = new char[charParamNamesSize];
   charParamValues = new char[numCharParams];
   stringParamNames = new char[stringParamNamesSize];
   stringParamValues = new char[stringParamValuesSize];
}

/** constructor with scip */
ScipDiffParamSet::ScipDiffParamSet(
      SCIP *scip
      )
      : numBoolParams(0), boolParamNamesSize(0), boolParamNames(0), boolParamValues(0),
      numIntParams(0), intParamNamesSize(0), intParamNames(0), intParamValues(0),
      numLongintParams(0), longintParamNamesSize(0), longintParamNames(0), longintParamValues(0),
      numRealParams(0), realParamNamesSize(0), realParamNames(0), realParamValues(0),
      numCharParams(0), charParamNamesSize(0), charParamNames(0), charParamValues(0),
      numStringParams(0), stringParamNamesSize(0), stringParamNames(0), stringParamValuesSize(0), stringParamValues(0)
{
   /** get Parameters */
   int nParams = SCIPgetNParams (scip);
   SCIP_PARAM **params = SCIPgetParams (scip);

   /** count the number of parameters for each type */
   for( int i = 0; i < nParams; i++ )
   {
      if( !SCIPparamIsDefault (params[i]) )
      {
         switch (SCIPparamGetType (params[i]))
         {
         case SCIP_PARAMTYPE_BOOL:
            numBoolParams++;
            boolParamNamesSize += strlen(SCIPparamGetName(params[i]) ) + 1;
            break;
         case SCIP_PARAMTYPE_INT:
            numIntParams++;
            intParamNamesSize += strlen(SCIPparamGetName(params[i]) ) + 1;
            break;
         case SCIP_PARAMTYPE_LONGINT:
            numLongintParams++;
            longintParamNamesSize += strlen(SCIPparamGetName(params[i]) ) + 1;
            break;
         case SCIP_PARAMTYPE_REAL:
            numRealParams++;
            realParamNamesSize += strlen(SCIPparamGetName(params[i]) ) + 1;
            break;
         case SCIP_PARAMTYPE_CHAR:
            numCharParams++;
            charParamNamesSize += strlen(SCIPparamGetName(params[i]) ) + 1;
            break;
         case SCIP_PARAMTYPE_STRING:
            numStringParams++;
            stringParamNamesSize += strlen(SCIPparamGetName(params[i]) ) + 1;
            stringParamValuesSize += strlen( SCIPparamGetString(params[i]) ) + 1;
            break;
         default:
            THROW_LOGICAL_ERROR2("Invalid SCIP parameter type: Type = ", params[i]);
         }
      }
   }

   /** allocate memory */
   allocateMemoty();

   /** set names and values for each type */
   int iBool = 0;
   int iInt = 0;
   int iLongint = 0;
   int iReal = 0;
   int iChar = 0;
   int iString = 0;
   char *boolName = boolParamNames;
   char *intName = intParamNames;
   char *longintName = longintParamNames;
   char *realName = realParamNames;
   char *charName = charParamNames;
   char *stringName = stringParamNames;
   char *stringValue = stringParamValues;
   for( int i = 0; i < nParams; i++ )
   {
      if( !SCIPparamIsDefault (params[i]) )
      {
         switch (SCIPparamGetType (params[i]))
         {
         case SCIP_PARAMTYPE_BOOL:
            strcpy(boolName, SCIPparamGetName(params[i]));
            boolParamValues[iBool] = SCIPparamGetBool(params[i]);
            boolName += strlen(SCIPparamGetName(params[i]) ) + 1;
            iBool++;
            break;
         case SCIP_PARAMTYPE_INT:
            strcpy(intName, SCIPparamGetName(params[i]));
            intParamValues[iInt] = SCIPparamGetInt(params[i]);
            intName += strlen( SCIPparamGetName(params[i]) ) + 1;
            iInt++;
            break;
         case SCIP_PARAMTYPE_LONGINT:
            strcpy(longintName, SCIPparamGetName(params[i]));
            longintParamValues[iLongint] = SCIPparamGetLongint(params[i]);
            longintName += strlen( SCIPparamGetName(params[i]) ) + 1;
            iLongint++;
            break;
         case SCIP_PARAMTYPE_REAL:
            strcpy(realName, SCIPparamGetName(params[i]));
            realParamValues[iReal] = SCIPparamGetReal(params[i]);
            realName += strlen( SCIPparamGetName(params[i]) ) + 1;
            iReal++;
            break;
         case SCIP_PARAMTYPE_CHAR:
            strcpy(charName, SCIPparamGetName(params[i]));
            charParamValues[iChar] = SCIPparamGetChar(params[i]);
            charName += strlen( SCIPparamGetName(params[i]) ) + 1;
            iChar++;
            break;
         case SCIP_PARAMTYPE_STRING:
            strcpy( stringName, SCIPparamGetName(params[i]) );
            strcpy( stringValue, SCIPparamGetString(params[i]) );
            stringName += strlen( SCIPparamGetName(params[i]) )  + 1;
            stringValue += strlen( SCIPparamGetString(params[i]) ) + 1;
            iString++;
            break;
         default:
            THROW_LOGICAL_ERROR2("Invalid SCIP parameter type: Type = ", params[i]);
         }
      }
   }

   assert(
         iBool == numBoolParams && iInt == numIntParams && iLongint == numLongintParams &&
         iReal == numRealParams && iChar == numCharParams && iString == numStringParams &&
         boolParamNamesSize == ( boolName - boolParamNames ) &&
         intParamNamesSize == ( intName - intParamNames ) &&
         longintParamNamesSize == ( longintName - longintParamNames ) &&
         realParamNamesSize == ( realName - realParamNames ) &&
         charParamNamesSize == ( charName - charParamNames ) &&
         stringParamNamesSize == ( stringName - stringParamNames ) &&
         stringParamValuesSize == ( stringValue - stringParamValues )
         );
}

/** set these parameter values in scip environment */
void
ScipDiffParamSet::setParametersInScip(
      SCIP *scip
      )
{
   char *paramName;

   /** set boolean parameter values in scip environment */
   paramName = boolParamNames;
   for( int i = 0; i < numBoolParams; i++ )
   {
      SCIP_CALL_ABORT(
            SCIPsetBoolParam(scip, paramName, boolParamValues[i])
            );
      paramName += strlen(paramName) + 1;
   }
   assert( boolParamNamesSize == ( paramName - boolParamNames ) );

   /** set int parameter values in scip environment */
   paramName = intParamNames;
   for( int i = 0; i < numIntParams; i++ )
   {
      SCIP_CALL_ABORT(
            SCIPsetIntParam(scip, paramName, intParamValues[i])
            );
      paramName += strlen(paramName) + 1;
   }
   assert( intParamNamesSize == ( paramName - intParamNames ) );

   /** set longint parameter values in scip environment */
   paramName = longintParamNames;
   for( int i = 0; i < numLongintParams; i++ )
   {
      SCIP_CALL_ABORT(
            SCIPsetLongintParam(scip, paramName, longintParamValues[i])
            );
      paramName += strlen(paramName) + 1;
   }
   assert( longintParamNamesSize == ( paramName - longintParamNames ) );

   /** set real parameter values in scip environment */
   paramName = realParamNames;
   for( int i = 0; i < numRealParams; i++ )
   {
      SCIP_CALL_ABORT(
            SCIPsetRealParam(scip, paramName, realParamValues[i])
            );
      paramName += strlen(paramName) + 1;
   }
   assert( realParamNamesSize == ( paramName - realParamNames ) );

   /** set char parameter values in scip environment */
   paramName = charParamNames;
   for( int i = 0; i < numCharParams; i++ )
   {
      SCIP_CALL_ABORT(
            SCIPsetCharParam(scip, paramName, charParamValues[i])
            );
      paramName += strlen(paramName) + 1;
   }
   assert( charParamNamesSize == ( paramName - charParamNames ) );

   /** set string parameter values in scip environment */
   paramName = stringParamNames;
   char *paramValue = stringParamValues;
   for( int i = 0; i < numStringParams; i++ )
   {
      SCIP_CALL_ABORT(
            SCIPsetStringParam(scip, paramName, paramValue)
            );
      paramName += strlen(paramName) + 1;
      paramValue += strlen(paramValue) + 1;
   }
   assert( stringParamNamesSize == ( paramName - stringParamNames ) );
   assert( stringParamValuesSize == ( paramValue - stringParamValues ) );

}

/** stringfy DiffParamSet */
std::string
ScipDiffParamSet::toString(
      )
{
    std::ostringstream s;

    char *paramName;

    /** stringfy boolean parameter values */
    paramName = boolParamNames;
    for( int i = 0; i < numBoolParams; i++ )
    {
       s << paramName << " = ";
       if( boolParamValues[i] )
       {
          s << "TRUE";
       }
       else
       {
          s << "FALSE";
       }
       s << std::endl;
       paramName += strlen(paramName) + 1;
    }
    assert( boolParamNamesSize == ( paramName - boolParamNames ) );

    /** stringfy int parameter values */
    paramName = intParamNames;
    for( int i = 0; i < numIntParams; i++ )
    {
       s << paramName << " = " << intParamValues[i] << std::endl;
       paramName += strlen(paramName) + 1;
    }
    assert( intParamNamesSize == ( paramName - intParamNames ) );

    /** stringfy longint parameter values */
    paramName = longintParamNames;
    for( int i = 0; i < numLongintParams; i++ )
    {
       s << paramName << " = " <<  longintParamValues[i] << std::endl;
       paramName += strlen(paramName) + 1;
    }
    assert( longintParamNamesSize == ( paramName - longintParamNames ) );

    /** stringfy real parameter values */
    paramName = realParamNames;
    for( int i = 0; i < numRealParams; i++ )
    {
       s << paramName << " = " << realParamValues[i] << std::endl;
       paramName += strlen(paramName) + 1;
    }
    assert( realParamNamesSize == ( paramName - realParamNames ) );

    /** stringfy char parameter values */
    paramName = charParamNames;
    for( int i = 0; i < numCharParams; i++ )
    {
       s << paramName << " = " << charParamValues[i] << std::endl;
       paramName += strlen(paramName) + 1;
    }
    assert( charParamNamesSize == ( paramName - charParamNames ) );

    /** stringfy string parameter values */
    paramName = stringParamNames;
    char *paramValue = stringParamValues;
    for( int i = 0; i < numStringParams; i++ )
    {
       s << paramName << " = \"" << paramValue << "\"" << std::endl;
       paramName += strlen(paramName) + 1;
       paramValue += strlen(paramValue) + 1;
    }
    return s.str();

}

/** write ScipDiffParamSet */
void
ScipDiffParamSet::write(
      ogzstream &out
      )
{
   /** write boolean parameter names and values */
   out.write((char *)&numBoolParams, sizeof(int));
   out.write((char *)&boolParamNamesSize, sizeof(int));
   if( numBoolParams > 0 )
   {
      out.write(boolParamNames, boolParamNamesSize);
      for( int i = 0; i < numBoolParams; i++ )
      {
         out.write((char *)&boolParamValues[i], sizeof(unsigned int));
      }
   }


   /** write int parameter names and values */
   out.write((char *)&numIntParams, sizeof(int));
   out.write((char *)&intParamNamesSize, sizeof(int));
   if( numIntParams > 0 )
   {
      out.write(intParamNames, intParamNamesSize);
      for( int i = 0; i < numIntParams; i++ )
      {
         out.write((char *)&intParamValues[i], sizeof(int));
      }
   }


   /** write longint parameter names and values */
   out.write((char *)&numLongintParams, sizeof(int));
   out.write((char *)&longintParamNamesSize, sizeof(int));
   if( numLongintParams > 0 )
   {
      out.write(longintParamNames, longintParamNamesSize);
      for( int i = 0; i < numLongintParams; i++ )
      {
         out.write((char *)&longintParamValues[i], sizeof(long long));
      }
   }

   /** write real parameter names and values */
   out.write((char *)&numRealParams, sizeof(int));
   out.write((char *)&realParamNamesSize, sizeof(int));
   if( numRealParams > 0 )
   {
      out.write(realParamNames, realParamNamesSize);
      for( int i = 0; i < numRealParams; i++ )
      {
         out.write((char *)&realParamValues[i], sizeof(double));
      }
   }


   /** write char parameter values */
   out.write((char *)&numCharParams, sizeof(int));
   out.write((char *)&charParamNamesSize, sizeof(int));
   if( numCharParams > 0 )
   {
      out.write(charParamNames, charParamNamesSize);
      for( int i = 0; i < numCharParams; i++ )
      {
         out.write((char *)&charParamValues[i], sizeof(char));
      }
   }

   /** write string parameter values */
   out.write((char *)&numStringParams, sizeof(int));
   out.write((char *)&stringParamNamesSize, sizeof(int));
   out.write((char *)&stringParamValuesSize, sizeof(int));
   if( numStringParams > 0 )
   {
      out.write(stringParamNames, stringParamNamesSize);
      out.write(stringParamValues, stringParamValuesSize);
   }
}

/** read ScipDiffParamSet */
bool
ScipDiffParamSet::read(
      ParaComm *comm,
      igzstream &in
      )
{
   /** read boolean parameter names and values */
   in.read((char *)&numBoolParams, sizeof(int));
   if( in.eof() ) return false;
   in.read((char *)&boolParamNamesSize, sizeof(int));
   if( numBoolParams > 0 )
   {
      boolParamNames = new char[boolParamNamesSize];
      boolParamValues = new unsigned int[numBoolParams];
      in.read(boolParamNames, boolParamNamesSize);
      for( int i = 0; i < numBoolParams; i++ )
      {
         in.read((char *)&boolParamValues[i], sizeof(unsigned int));
      }
   }

   /** read int parameter names and values */
   in.read((char *)&numIntParams, sizeof(int));
   in.read((char *)&intParamNamesSize, sizeof(int));
   if( numIntParams > 0 )
   {
      intParamNames = new char[intParamNamesSize];
      intParamValues = new int[numIntParams];
      in.read(intParamNames, intParamNamesSize);
      for( int i = 0; i < numIntParams; i++ )
      {
         in.read((char *)&intParamValues[i], sizeof(int));
      }
   }

   /** read longint parameter names and values */
   in.read((char *)&numLongintParams, sizeof(int));
   in.read((char *)&longintParamNamesSize, sizeof(int));
   if( numLongintParams > 0 )
   {
      longintParamNames = new char[longintParamNamesSize];
      longintParamValues = new long long[numLongintParams];
      in.read(longintParamNames, longintParamNamesSize);
      for( int i = 0; i < numLongintParams; i++ )
      {
         in.read((char *)&longintParamValues[i], sizeof(long long));
      }
   }

   /** read real parameter names and values */
   in.read((char *)&numRealParams, sizeof(int));
   in.read((char *)&realParamNamesSize, sizeof(int));
   if( numRealParams > 0 )
   {
      realParamNames = new char[realParamNamesSize];
      realParamValues = new double[numRealParams];
      in.read(realParamNames, realParamNamesSize);
      for( int i = 0; i < numRealParams; i++ )
      {
         in.read((char *)&realParamValues[i], sizeof(double));
      }
   }

   /** read char parameter values */
   in.read((char *)&numCharParams, sizeof(int));
   in.read((char *)&charParamNamesSize, sizeof(int));
   if( numCharParams > 0 )
   {
      charParamNames = new char[charParamNamesSize];
      charParamValues = new char[numCharParams];
      in.read(charParamNames, charParamNamesSize);
      for( int i = 0; i < numCharParams; i++ )
      {
         in.read((char *)&charParamValues[i], sizeof(char));
      }
   }

   /** read string parameter values */
   in.read((char *)&numStringParams, sizeof(int));
   in.read((char *)&stringParamNamesSize, sizeof(int));
   in.read((char *)&stringParamValuesSize, sizeof(int));
   if( numStringParams > 0 )
   {
      stringParamNames = new char[stringParamNamesSize];
      stringParamValues = new char[stringParamValuesSize];
      in.read(stringParamNames, stringParamNamesSize);
      in.read(stringParamValues, stringParamValuesSize);
   }

   return true;
}
