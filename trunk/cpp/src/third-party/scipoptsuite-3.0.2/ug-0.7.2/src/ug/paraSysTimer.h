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

/**@file    paraSysTimer.h
 * @brief   System timer.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_SYS_TIMER_H__
#define __PARA_SYS_TIMER_H__

//-------------------
// include files
//-------------------
#include <cstdlib>
#include "paraTimer.h"

#ifdef __APPLE__
#define BSD
#endif /* SUN_OS */

#ifdef SUN_OS
#define BSD
#endif /* SUN_OS */

#ifdef SOLARIS
#define SYSV
#endif /* SOLARIS */

#ifdef linux
#define SYSV
#endif /* linux */

#ifdef BlueGene
#define BSD
#endif

#if !(defined WIN32 || defined SYSV || defined BSD )
#error cannot detect timer type!
#endif

#ifdef BSD
#include <sys/time.h>
#include <sys/resource.h>
#endif /* BSD */

#ifdef WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#include <time.h>
#include <windows.h>
#endif /* WIN32 */

#ifdef SYSV
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#define TICKS HZ
#endif /* SYSV */

namespace UG
{
//-------------------------------------------
// Class definition
//-------------------------------------------

class ParaSysTimer : public ParaTimer {
public:
   ParaSysTimer(
         )
   {
   }

   ~ParaSysTimer(
         )
   {
   }

   void init(
         ParaComm* paraComm
         )
   {
      start();
   }

   double getElapsedTime(
         )
   {
      return getRTimeInterval();
   }

   void    start(void);
   void    stop(void);
   double  getStartTime(void);
   double  getRTimeInterval(void);
   double  getRTime(void);
   double  getUTime(void);
   double  getSTime(void);

private:

#ifdef BSD
   struct timeval   stTvTimeStart, stTvTimeStop;
   struct rusage    stRuStart, stRuStop;
#  endif /* BSD */

#ifdef WIN32
   struct _timeb timebStart, timebStop;
   FILETIME ftCreationTime, ftExitTime,
	        ftKernelTimeStart, ftUserTimeStart,
			ftKernelTimeStop,  ftUserTimeStop;
   HANDLE   hCurrentProcess;
#  endif /* WIN32 */

#  ifdef SYSV
   long lTimeStart, lTimeStop;
   struct tms stTmsStart, stTmsStop;
/*   long times(); */
#  endif /* SYSV */

};

}

#endif  // __PARA_SYS_TIMER_H__
