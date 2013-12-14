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

/**@file    paraInitialStat.h
 * @brief   Base class for initial statistics collecting class
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_INITIAL_STAT_H__
#define __PARA_INITIAL_STAT_H__

#include <iostream>
#include "paraComm.h"

namespace UG
{

/** The initial statistics collecting  */
class ParaInitialStat {
   /*****************************
    * DO NOT HAVE DATA MEMBER!! *
    ****************************/
public:
   /** default constructor */
   ParaInitialStat()
   {
   }

   /** destractor */
   virtual ~ParaInitialStat()
   {
   }

   /** create clone of this object */
   virtual ParaInitialStat *clone(ParaComm *comm) = 0;

   /** user should implement send method */
   virtual void send(ParaComm *comm, int dest) = 0;

   /** user should implement receive method */
   virtual void receive(ParaComm *comm, int source) = 0;

   /** get maximum depth */
   virtual int getMaxDepth() = 0;

   /** stringfy subproblem ( for debugging ) */
   virtual const std::string toString() = 0;
};

}

#endif    // __PARA_INITIAL_STAT_H__

