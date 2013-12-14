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

/**@file    scipParaSolutionPth.cpp
 * @brief   ScipParaSolution extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


// #include "paraTagDefPth.h"
#include "scipParaSolutionPth.h"

using namespace UG;
using namespace ParaSCIP;

/** create clone of this object */
ScipParaSolutionPth *
ScipParaSolutionPth::clone(UG::ParaComm *comm)
{
   return( new ScipParaSolutionPth(objectiveFunctionValue, nVars, indicesAmongSolvers, values));
}

/** create ScipDiffSubproblemPreDatatype */
ScipParaSolutionPth *
ScipParaSolutionPth::createDatatype(
      UG::ParaComm *comm
      )
{
   return clone(comm);
}

/** send solution data to the rank */
void
ScipParaSolutionPth::bcast(ParaComm *comm, int root)
{
   DEF_PARA_COMM( commPth, comm);

   if( commPth->getRank() == root )
   {
      for( int i = 0; i < commPth->getSize(); i++ )
      {
         if( i != root )
         {
            PARA_COMM_CALL(
               commPth->uTypeSend((void *)createDatatype(comm), ParaSolutionType, i, TagSolution)
            );
         }
      }
   }
   else
   {
      ScipParaSolutionPth *received;
      PARA_COMM_CALL(
         commPth->uTypeReceive((void **)&received, ParaSolutionType, root, TagSolution)
      );

      objectiveFunctionValue = received->objectiveFunctionValue;
      nVars = received->nVars;
      indicesAmongSolvers = new int[nVars];
      values = new SCIP_Real[nVars];
      for( int i = 0; i < nVars; i++ )
      {
         indicesAmongSolvers[i] = received->indicesAmongSolvers[i];
         values[i] = received->values[i];
      }
      delete received;
   }
}

/** send solution data to the rank */
void
ScipParaSolutionPth::send(ParaComm *comm, int destination)
{
   DEF_PARA_COMM( commPth, comm);
   PARA_COMM_CALL(
      commPth->uTypeSend((void *)createDatatype(comm), ParaSolutionType, destination, TagSolution)
   );
}

/** receive solution data from the source rank */
void
ScipParaSolutionPth::receive(ParaComm *comm, int source)
{
   DEF_PARA_COMM( commPth, comm);

   ScipParaSolutionPth *received;
   PARA_COMM_CALL(
      commPth->uTypeReceive((void **)&received, ParaSolutionType, source, TagSolution)
   );

   objectiveFunctionValue = received->objectiveFunctionValue;
   nVars = received->nVars;
   indicesAmongSolvers = new int[nVars];
   values = new SCIP_Real[nVars];
   for( int i = 0; i < nVars; i++ )
   {
      indicesAmongSolvers[i] = received->indicesAmongSolvers[i];
      values[i] = received->values[i];
   }
   delete received;

}
