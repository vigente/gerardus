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

/**@file    paraParamSetMpi.h
 * @brief   ParaParamSet extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_PARAM_SET_MPI_H__
#define __PARA_PARAM_SET_MPI_H__
#include <mpi.h>
#include "paraCommMpiWorld.h"
#include "paraParamSet.h"

namespace UG
{

class ParaParamSetMpi : public ParaParamSet {
   /** only transfer NOT default values */
   int          nBoolParams;              /**< the number of bool parameters */
   int          *boolParams;              /**< boolean parameter ids */
   char         *boolParamValues;         /**< boolean parameter values */

   int          nIntParams;               /**< the number of int parameters */
   int          *intParams;               /**< int parameter ids */
   int          *intParamValues;          /**< int parameter values */

   int          nLongintParams;           /**< the number of longint parameters */
   int          *longintParams;           /**< longint parameter ids */
   long long    *longintParamValues;      /**< longint parameter values */

   int          nRealParams;              /**< the number of real parameters */
   int          *realParams;              /**< real parameter ids */
   double       *realParamValues;         /**< real parameter values */

   int          nCharParams;              /**< the number of char parameters */
   int          *charParams;              /**< char parameter ids */
   char         *charParamValues;         /**< char parameter values */

   int          nStringParams;            /**< the number of string parameters */
   int          *stringParams;            /**< string parameter ids */
   int          stringParamValuesSize;    /**< size of stringParameterValues area */
   char         *stringParamValues;       /**< string parameter values: values are concatenated */

   void allocateMemory();
   void freeMemory();
   void createDiffParams();
   void setDiffParams();

   /** create ParaParamSetDatatype1 */
   MPI::Datatype createDatatype1();
   /** create ParaParamSetDatatype2 */
   MPI::Datatype createDatatype2(bool reallocateStringPramsValue);
public:
    ParaParamSetMpi() :
       nBoolParams(0), boolParams(0), boolParamValues(0),
       nIntParams(0), intParams(0), intParamValues(0),
       nLongintParams(0), longintParams(0), longintParamValues(0),
       nRealParams(0), realParams(0), realParamValues(0),
       nCharParams(0), charParams(0), charParamValues(0),
       nStringParams(0), stringParams(0), stringParamValuesSize(0), stringParamValues(0) {}
    int bcast(ParaComm *comm, int root);
    ~ParaParamSetMpi(){}
};

}

#endif  // __PARA_PARAM_SET_MPI_H__
