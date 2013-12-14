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

/**@file    paraTimeLimitMonitorPth.h
 * @brief   Time limit monitor thread class.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_TIME_LIMIT_MONITOR_PTH_H__
#define __PARA_TIME_LIMIT_MONITOR_PTH_H__

#include <algorithm>
#include "paraDef.h"
#include "paraCommPth.h"

namespace UG
{

class ParaTimeLimitMonitorPth
{
protected:
   ParaCommPth  *paraComm;                  /**< ParaCommunicator object */
   double       hardTimeLimit;              /**< hard time limit */
public:
   ParaTimeLimitMonitorPth(){ THROW_LOGICAL_ERROR1("Default constructor of ParaTimeLimitMonitor is called"); }
   ParaTimeLimitMonitorPth(
       ParaCommPth *comm,
       double timelimit
       ) : paraComm(comm)
   {
      hardTimeLimit = timelimit + std::min(60.0, 0.1 * timelimit);   // set hard time limit + 60 seconds longer than time limit
   }

   virtual ~ParaTimeLimitMonitorPth(){}
   void run(){
      sleep(static_cast<unsigned int>(hardTimeLimit));
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, ParaBYTE, 0, TagHardTimeLimit)
      );
   }
};

}

#endif // __PARA_TIME_LIMIT_MONITOR_PTH_H__

