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

/**@file    scipParaObjBranchRule.h
 * @brief   Branching rule plug-in for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PARA_OBJ_BRANCHRULE_H__
#define __SCIP_PARA_OBJ_BRANCHRULE_H__

#include <cassert>
#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"
#include "objscip/objscip.h"
#include "scipParaSolver.h"

namespace ParaSCIP
{

/** C++ wrapper object for branching rules */
class ScipParaObjBranchRule : public scip::ObjBranchrule
{
public:
   /*lint --e{1540}*/

   /** SCIP ParaSolver */
   ScipParaSolver *scipParaSolver;

   /** parasolver constructor */
   ScipParaObjBranchRule(
      ScipParaSolver *solver       /**< SCIP ParaSolver */
      )
   : ObjBranchrule(solver->getScip(), "ScipParaObjBranchRule", "Branch rule plug-in for ParaSCIP",
         999999, -1, 1.0 ), scipParaSolver(solver)
   {
   }

   /** destructor */
   virtual ~ScipParaObjBranchRule()
   {
   }

   /** destructor of branching rule to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** initialization method of branching rule (called after problem was transformed) */
   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** deinitialization method of branching rule (called before transformed problem is freed) */
   virtual SCIP_RETCODE scip_exit(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
   virtual SCIP_RETCODE scip_initsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
   virtual SCIP_RETCODE scip_exitsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_BRANCHRULE*   branchrule          /**< the branching rule itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** branching execution method for fractional LP solutions
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the current node was detected to be infeasible
    *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict clause) was generated; this result code must not be
    *                      returned, if allowaddcons is FALSE
    *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current LP solution infeasible
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_BRANCHED   : branching was applied
    *  - SCIP_DIDNOTRUN  : the branching rule was skipped
    */
   virtual SCIP_RETCODE scip_execlp(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_BRANCHRULE*   branchrule,         /**< the branching rule itself */
      SCIP_Bool          allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
      SCIP_RESULT*       result              /**< pointer to store the result of the branching call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
      {
          /** set integer infeaibility statistics */
          int ncands = 0;
          SCIP_Real* fracs = 0;
          SCIP_Real iisum = 0.0;
          SCIP_CALL_ABORT( SCIPgetLPBranchCands(scip, NULL,NULL, &fracs, &ncands, NULL) );
          for( int i = 0; i < ncands; ++i )
            iisum += fracs[i];
          scipParaSolver->setII(iisum, ncands);
      }
      return SCIP_OKAY;
   }
   
   /** branching execution method for external candidates
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the current node was detected to be infeasible
    *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict clause) was generated; this result code must not be
    *                      returned, if allowaddcons is FALSE
    *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current pseudo solution infeasible
    *  - SCIP_BRANCHED   : branching was applied
    *  - SCIP_DIDNOTRUN  : the branching rule was skipped
    */
   virtual SCIP_RETCODE scip_execext(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_BRANCHRULE*   branchrule,         /**< the branching rule itself */
      SCIP_Bool          allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
      SCIP_RESULT*       result              /**< pointer to store the result of the branching call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** branching execution method for not completely fixed pseudo solutions
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the current node was detected to be infeasible
    *  - SCIP_CONSADDED  : an additional constraint (e.g. a conflict clause) was generated; this result code must not be
    *                      returned, if allowaddcons is FALSE
    *  - SCIP_REDUCEDDOM : a domain was reduced that rendered the current pseudo solution infeasible
    *  - SCIP_BRANCHED   : branching was applied
    *  - SCIP_DIDNOTRUN  : the branching rule was skipped
    */
   virtual SCIP_RETCODE scip_execps(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_BRANCHRULE*   branchrule,         /**< the branching rule itself */
      SCIP_Bool          allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
      SCIP_RESULT*       result              /**< pointer to store the result of the branching call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
};

} /* namespace ParaSCIP */
#endif
