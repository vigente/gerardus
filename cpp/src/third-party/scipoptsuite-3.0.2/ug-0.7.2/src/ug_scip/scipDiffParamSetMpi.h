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

/**@file    scipDiffParamSetMpi.h
 * @brief   ScipDiffParamSet extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_DIFF_PARAM_SET_MPI_H__
#define __SCIP_DIFF_PARAM_SET_MPI_H__

#include <mpi.h>
#include "ug/paraComm.h"
#include "scipDiffParamSet.h"

namespace ParaSCIP
{

/** ScipDiffParamSet class */
class ScipDiffParamSetMpi: public ScipDiffParamSet
{

   /** create scipDiffParamSetPreType */
   MPI::Datatype createDatatype1();
   /** create scipDiffParamSetType */
   MPI::Datatype createDatatype2(bool memAllocNecessary);

public:
   /** constructor */
   ScipDiffParamSetMpi(
         )
	  {
   }

   /** constructor with scip */
   ScipDiffParamSetMpi(
         SCIP *scip
         )
         : ScipDiffParamSet(scip)
   {
   }

   /** destructor */
   ~ScipDiffParamSetMpi(
		   )
   {
   }

   /** broadcast scipDiffParamSet */
   int bcast(UG::ParaComm *comm, int root);

   /** send scipDiffParamSet to the rank */
   int send(UG::ParaComm *comm, int destination);

   /** receive scipDiffParamSet from the source rank */
   int receive(UG::ParaComm *comm, int source);

};

typedef ScipDiffParamSet *ScipDiffParamSetPtr;

}

#endif // _SCIP_DIFF_PARAM_SET_MPI_H__

