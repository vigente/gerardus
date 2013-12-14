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

/**@file    scipParaDiffSubproblemMpi.h
 * @brief   ScipParaDiffSubproblem extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_DIFF_SUBPROBLEM_MPI_H__
#define __SCIP_PARA_DIFF_SUBPROBLEM_MPI_H__

#include <mpi.h>
#include "ug/paraComm.h"
#include "scipParaDiffSubproblem.h"

namespace ParaSCIP
{

/** The difference between instance and subproblem: this is base class */
class ScipParaDiffSubproblemMpi : public ScipParaDiffSubproblem
{

   /** create scipDiffSubproblem datatype1 */
   MPI::Datatype createDatatype1();
   /** create scipDiffSubproblem datatype2 */
   MPI::Datatype createDatatype2(bool memAllocNecessary);
   /** create scipDiffSubproblem datatype2 */
   MPI::Datatype createDatatype3(bool memAllocNecessary);

public:
   /** default constructor */
   ScipParaDiffSubproblemMpi()
   {
      assert( localInfoIncluded == 0 && nBoundChanges == 0 && nLinearConss == 0 );
   }

   /** Constructor */
   ScipParaDiffSubproblemMpi(
         SCIP *inScip,
         ScipParaSolver *inScipParaSolver,
         int inNNewBranchVars,
         SCIP_VAR **inNewBranchVars,
         SCIP_Real *inNewBranchBounds,
         SCIP_BOUNDTYPE *inNewBoundTypes
         ) : ScipParaDiffSubproblem(inScip, inScipParaSolver,
               inNNewBranchVars, inNewBranchVars, inNewBranchBounds,inNewBoundTypes)
   {
   }

   /** Constructor */
   ScipParaDiffSubproblemMpi(
         ScipParaDiffSubproblem *paraDiffSubproblem
         ) : ScipParaDiffSubproblem(paraDiffSubproblem)
   {
   }


   /** destractor */
   ~ScipParaDiffSubproblemMpi()
   {
   }

   /** create clone of this object */
   ScipParaDiffSubproblemMpi *clone(UG::ParaComm *comm);

   int bcast(UG::ParaComm *comm, int root);

   int send(UG::ParaComm *comm, int dest);

   int receive(UG::ParaComm *comm, int source);
};

}

#endif    // __SCIP_PARA_DIFF_SUBPROBLEM_MPI_H__

