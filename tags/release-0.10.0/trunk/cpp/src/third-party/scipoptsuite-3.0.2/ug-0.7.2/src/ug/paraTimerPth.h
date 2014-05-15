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

/**@file    paraTimerPth.h
 * @brief   ParaTimer extension for Pthreads.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_TIMER_PTH_H__
#define __PARA_TIMER_PTH_H__
#include "paraTimer.h"
#include "paraSysTimer.h"

namespace UG
{

class ParaTimerPth : public ParaTimer
{
   ParaSysTimer timer;
public:
   ParaTimerPth()
   {
   }
   ~ParaTimerPth() {}
   void init(ParaComm *comm)
   {
      timer.start();
   }
   double getElapsedTime()
   {
      return timer.getRTimeInterval();
   }
};

}

#endif // __PARA_TIMER_PTH_H__
