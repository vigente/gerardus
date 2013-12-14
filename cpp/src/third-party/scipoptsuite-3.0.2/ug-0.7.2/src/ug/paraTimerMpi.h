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

/**@file    paraTimerMpi.h
 * @brief   ParaTimer extension for MPI timer.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_TIMER_MPI_H__
#define __PARA_TIMER_MPI_H__
#include <mpi.h>
#include "paraComm.h"
#include "paraTimer.h"

namespace UG
{

class ParaTimerMpi : public ParaTimer
{
   double startTime;
public:
   ParaTimerMpi()
   {
      startTime = MPI::Wtime();
   }
   ~ParaTimerMpi() {}
   void init(ParaComm *comm);
   double getElapsedTime()
   {
      return ( MPI::Wtime() - startTime );
   }
};

}

#endif // __PARA_TIMER_MPI_H__
