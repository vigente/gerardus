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

/**@file    paraDiffSubproblem.h
 * @brief   Base class for a container which has difference between instance and subproblem.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_DIFF_SUBPROBLEM_H__
#define __PARA_DIFF_SUBPROBLEM_H__

#include <iostream>
#include <fstream>
#include "gzstream.h"
#include "paraComm.h"
#include "paraInstance.h"
#include "paraMergeNodesStructs.h"

namespace UG
{

/** The difference between instance and subproblem: this is base class */
class ParaDiffSubproblem {
   /*****************************
    * DO NOT HAVE DATA MEMBER!! *
    ****************************/
public:
   /** default constructor */
   ParaDiffSubproblem()
   {
   }

   /** destractor */
   virtual ~ParaDiffSubproblem()
   {
   }

   /** create clone of this object */
   virtual ParaDiffSubproblem *clone(ParaComm *comm) = 0;

   /** user should implement send method */
   virtual int bcast(ParaComm *comm, int root) = 0;

   /** user should implement send method */
   virtual int send(ParaComm *comm, int dest) = 0;

   /** user should implement receive method */
   virtual int receive(ParaComm *comm, int source) = 0;

   /** user should implement write method */
   virtual void write(ogzstream &out) = 0;

   /** user should implement read method */
   virtual void read(ParaComm *comm, igzstream &in, bool onlyBoundChanges) = 0;

   /** get the number of bound changes */
   virtual int getNBoundChanges() = 0;

   /** get fixed variables **/
   virtual int getFixedVariables(ParaInstance *instance, ParaFixedVariable **fixedVars ) = 0;

   /** create new ParaDiffSubproblem using fixed variables information */
   virtual ParaDiffSubproblem* createDiffSubproblem(ParaComm *comm, int n, ParaFixedVariable *fixedVars ) = 0;

   /** stringfy subproblem ( for debugging ) */
   virtual const std::string toString() = 0;
};

}

#endif    // __PARA_DIFF_SUBPROBLEM_H__

