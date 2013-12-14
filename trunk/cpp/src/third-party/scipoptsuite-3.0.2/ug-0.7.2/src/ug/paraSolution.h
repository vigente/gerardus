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

/**@file    paraSolution.h
 * @brief   Base class for solution.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_SOLUTION_H__
#define __PARA_SOLUTION_H__

#include "paraComm.h"
#include "gzstream.h"

namespace UG
{

/** Solution class ( Base class ) */
class ParaSolution
{
   /*****************************
    * DO NOT HAVE DATA MEMBER!! *
    ****************************/
public:

   /** constructor */
   ParaSolution(
      )
   {
   }

   /** destructor */
   virtual ~ParaSolution(
       )
   {
   }

   /** get objective function value */
   virtual double getObjectiveFuntionValue() = 0;

   /** create clone of this object */
   virtual ParaSolution *clone(ParaComm *comm) = 0;

   /** broadcast solution data to from the root rank */
   virtual void bcast(ParaComm *comm, int root) = 0;

   /** send solution data to the rank */
   virtual void send(ParaComm *comm, int destination) = 0;

   /** receive solution data from the source rank */
   virtual void receive(ParaComm *comm, int source) = 0;

   /** user should implement write method */
   virtual void write(ogzstream &out) = 0;

   /** user should implement read method */
   virtual bool read(ParaComm *comm, igzstream &in) = 0;

   virtual const std::string toString(){ return std::string(""); }

};

typedef ParaSolution *ParaSolutionPtr;

}

#endif // __PARA_SOLUTION_H__
