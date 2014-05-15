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

/**@file    scipParaDiffSubproblemPth.h
 * @brief   ScipParaDiffSubproblem extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_DIFF_SUBPROBLEM_PTH_H__
#define __SCIP_PARA_DIFF_SUBPROBLEM_PTH_H__

#include "ug/paraDef.h"
#include "ug/paraCommPth.h"
#include "scipParaDiffSubproblem.h"
#include "scip/scip.h"

namespace ParaSCIP
{

class ScipParaSolver;

/** The difference between instance and subproblem: this is base class */
class ScipParaDiffSubproblemPth : public ScipParaDiffSubproblem
{
public:
   /** default constructor */
   ScipParaDiffSubproblemPth()
   {
   }

   /** Constructor */
   ScipParaDiffSubproblemPth(
         SCIP *inScip,
         ScipParaSolver *inScipParaSolver,
         int inNNewBranchVars,
         SCIP_VAR **inNewBranchVars,
         SCIP_Real *inNewBranchBounds,
         SCIP_BOUNDTYPE *inNewBoundTypes
         ) : ScipParaDiffSubproblem(inScip, inScipParaSolver,
               inNNewBranchVars, inNewBranchVars, inNewBranchBounds,inNewBoundTypes)
   {
   }

   /** Constructor */
   ScipParaDiffSubproblemPth(
         ScipParaDiffSubproblem *paraDiffSubproblem
         ) : ScipParaDiffSubproblem(paraDiffSubproblem)
   {
   }


   /** destractor */
   ~ScipParaDiffSubproblemPth()
   {
   }

   /** create clone of this object */
   ScipParaDiffSubproblemPth *clone(
         UG::ParaComm *comm
         )
   {
      return(
            new ScipParaDiffSubproblemPth(this)
      );
   }

   int bcast(
         UG::ParaComm *comm,
         int root
         )
   {
      THROW_LOGICAL_ERROR1("bcast is issued in ScipParaDiffSubproblemPth");
   }

   int send(
         UG::ParaComm *comm,
         int dest
         )
   {
      THROW_LOGICAL_ERROR1("send is issued in ScipParaDiffSubproblemPth");
   }

   int receive(
         UG::ParaComm *comm,
         int source
         )
   {
      THROW_LOGICAL_ERROR1("receive is issued in ScipParaDiffSubproblemPth");
   }
};

}

#endif    // __SCIP_PARA_DIFF_SUBPROBLEM_PTH_H__

