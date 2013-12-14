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

/**@file    scipParaSolutionMpi.h
 * @brief   ScipParaSolution extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_SOLUTION_MPI_H__
#define __SCIP_PARA_SOLUTION_MPI_H__

#include <mpi.h>
#include "ug/paraTagDefMpi.h"
#include "scipParaSolution.h"

namespace ParaSCIP
{

/** ScipSolution class */
class ScipParaSolutionMpi : public ScipParaSolution
{
   /** create scipSolutionPreDatatype */
   MPI::Datatype createPreDatatype();
   /** create scipSolutionDatatype */
   MPI::Datatype createDatatype(bool memAllocNecessary);

public:

   /** default constructor */
   ScipParaSolutionMpi(
	      )
   {
   }

   /** constructor */
   ScipParaSolutionMpi(
         SCIP_Real      objval,
         int            inNvars,
         SCIP_VAR **    vars,
         SCIP_Real *    vals
         )
        : ScipParaSolution(objval, inNvars, vars, vals){}

   /** constructor */
   ScipParaSolutionMpi(
         double inObjectiveFunctionValue,
         int inNVars,                       /**< number of variables */
         int *inIndicesAmongSolvers,        /**< array of variable indices ( probindex )  */
         SCIP_Real *inValues                /**< array of bounds which the branchings     */
         ): ScipParaSolution(inObjectiveFunctionValue, inNVars, inIndicesAmongSolvers, inValues) {}

   /** destructor */
   ~ScipParaSolutionMpi(
         )
   {
   }

   /** create clone of this object */
   ScipParaSolutionMpi *clone(UG::ParaComm *comm);

   /** broadcast solution data to from the root rank */
   void bcast(UG::ParaComm *comm, int root);

   /** send solution data to the rank */
   void send(UG::ParaComm *comm, int destination);

   /** receive solution data from the source rank */
   void receive(UG::ParaComm *comm, int source);

};

typedef ScipParaSolutionMpi *ScipParaSolutionMpiPtr;

}

#endif // __SCIP_PARA_SOLUTION_MPI_H__

