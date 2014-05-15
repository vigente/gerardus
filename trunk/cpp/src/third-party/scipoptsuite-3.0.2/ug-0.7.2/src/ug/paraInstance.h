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

/**@file    paraInstance.h
 * @brief   Base class for instance data.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_INSTANCE_H__
#define __PARA_INSTANCE_H__

#include "paraComm.h"

namespace UG
{

/** Instance (base class) */
class ParaInstance
{
   /*****************************
    * DO NOT HAVE DATA MEMBER!! *
    ****************************/
public:
   /** constructor */
   ParaInstance(
         )
   {
   }

   /** destractor */
   virtual ~ParaInstance(
        )
   {
   }

   /** get problem name */
   virtual const char* getProbName() = 0;

   /** get number of variables */
   virtual int getNVars() = 0;

   /** broadcasts instance to all solvers */
   virtual int bcast(ParaComm *comm, int rank, int method) = 0;

   /** for debgu */
   virtual const std::string toString() = 0;

};

}

#endif  // __PARA_INSTANCE_H__
