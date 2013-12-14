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

/**@file    scipParaInstanceMpi.h
 * @brief   ScipParaInstance extension for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_INSTANCE_MPI_H__
#define __SCIP_PARA_INSTANCE_MPI_H__

#include <mpi.h>
#include "ug/paraDef.h"
#include "ug/paraComm.h"
#include "scipParaInstance.h"

namespace ParaSCIP
{

/** ScipInstanceMpi */
class ScipParaInstanceMpi : public ScipParaInstance
{

   int dummyToKeepStartPos;
   char *fileName;

   /** create scipDiffParamSetPreType */
   MPI::Datatype createDatatype1();
   /** create scipDiffParamSetPreType */
   MPI::Datatype createDatatype2(bool memAllocNecessary);
   /** create scipDiffParamSetType */
   MPI::Datatype createDatatype3(bool memAllocNecessary);

   void allocateMemoryForDatatype2();
   void allocateMemoryForDatatype3();

   char *getFileName(){ return fileName; }

public:
   /** constructor */
   ScipParaInstanceMpi(
         )
         : dummyToKeepStartPos(0), fileName(0)
   {
   }

   /** constructor : only called from ScipInitiator */
   ScipParaInstanceMpi(
         SCIP *scip,
         int  method
         ) : ScipParaInstance(scip, method), dummyToKeepStartPos(0)
   {
   }

   /** destractor */
   ~ScipParaInstanceMpi(
	        )
   {
   }

  /** create presolved problem instance that is solved by ParaSCIP form scip environment in this object */
  void copyScipEnvironment(
        SCIP **scip
        )
  {
     /** this routine for Pthread version. So, this should not be used **/
     abort();
  }

  SCIP *getScip(
        )
  {
     /** this routine for Pthread version. So, this should not be used **/
     abort();
  }

  void setFileName(char *inFileName)
  {
     fileName = inFileName;
  }

   /** broadcasts instance to all solvers */
   int bcast(UG::ParaComm *comm, int rank, int method);
};

typedef ScipParaInstanceMpi *ScipParaInstanceMpiPtr;

}

#endif  // __SCIP_PARA_INSTANCE_MPI_H__

