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

/**@file    scipParaDeterministicTimer.h
 * @brief   ParaDeterministicTimer extension for SCIP.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_DETERMINISITC_TIMER_H__
#define __SCIP_PARA_DETERMINISITC_TIMER_H__

#include "ug/paraDeterministicTimer.h"

namespace UG
{

class ScipParaDeterministicTimer : public ParaDeterministicTimer
{
   double current; 
   int    normalizeFactor;
public:
   ScipParaDeterministicTimer() : current(0.0), normalizeFactor(1) {}
   virtual ~ScipParaDeterministicTimer() {}
   /**********************************************
    * if you want to set original initial time,  *
    * you can do it init()                       *
    **********************************************/
   void normalize(ParaComm *comm){ normalizeFactor = comm->getSize() - 1; }
   void update(double value) { current += value; }
   double getElapsedTime() { return current/normalizeFactor; }
};

}

#endif  // __SCIP_PARA_DETERMINISTIC_TIMER_H__
