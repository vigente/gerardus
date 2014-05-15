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

/**@file    paraDeterministicTimer.h
 * @brief   Base class for deterministic timer
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_DETERMINISTIC_TIMER_H__
#define __PARA_DETERMINISTIC_TIMER_H__

#include "paraComm.h"

namespace UG
{

class ParaDeterministicTimer
{
public:
   ParaDeterministicTimer() {}
   virtual ~ParaDeterministicTimer() {}
   /**********************************************
    * if you want to set original initial time,  *
    * you can do it init()                       *
    **********************************************/
   // virtual void init() = 0;   // arguments shuld be different depending on MIP solver parallelized
                                 // So, add init() function in derived classes and use it in constuctor of
                                 // xxxParaSolver class
   virtual void normalize(ParaComm *comm){}
   virtual void update(double value) = 0;
   virtual double getElapsedTime() = 0;
};

}

#endif  // __PARA_TIMER_H__
