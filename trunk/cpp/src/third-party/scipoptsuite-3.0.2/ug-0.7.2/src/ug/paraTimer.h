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

/**@file    paraTimer.h
 * @brief   Base class for Timer.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_TIMER_H__
#define __PARA_TIMER_H__

namespace UG
{

class ParaComm;

class ParaTimer
{
public:
   ParaTimer() {}
   virtual ~ParaTimer() {}
   /**********************************************
    * if you want to set original initial time,  *
    * you can do it init()                       *
    **********************************************/
   virtual void init(ParaComm* paraComm) = 0;
   virtual double getElapsedTime() = 0;
};

}

#endif  // __PARA_TIMER_H__
