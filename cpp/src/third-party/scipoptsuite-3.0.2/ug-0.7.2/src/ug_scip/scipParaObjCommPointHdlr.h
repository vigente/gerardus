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

/**@file    scipParaObjCommPointHdlr.h
 * @brief   Event handlr for communication point.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PARA_COMM_POINT_HDLR_H__
#define __SCIP_PARA_COMM_POINT_HDLR_H__

#include <cstring>
#include "scipParaComm.h"
#include "scipParaInstance.h"
#include "scipParaSolver.h"
#include "objscip/objeventhdlr.h"
#include "scip/scipdefplugins.h"

namespace ParaSCIP
{

/** C++ wrapper object for event handlers */
class ScipParaObjCommPointHdlr : public scip::ObjEventhdlr
{
   UG::ParaComm   *paraComm;
   ScipParaSolver *scipParaSolver;
   SCIP           *scipToCheckRootSolvability;
   SCIP           *originalScip;
   bool           needToSendNode;
   bool           originalSelectionStrategy;
   int            originalPriority;
   SCIP_Longint   previousNNodesSolved;
   SCIP_Longint   previousLpIter;
   const char     *changeNodeSelName;
   // double         previousCommTime;
   bool           cloned;                  // indicate that this hander is cloned or not
   bool           interrupting;            // indicate that this handler called interrupt or not
   void           processNewSolution(SCIP *scip, SCIP_EVENT* event);
   void           checkRootNodeSolvabilityAndSendParaNode(SCIP *scip);
   void           sendNode( SCIP *scip, SCIP_NODE* node, int depth, int nBranchVars, SCIP_VAR **branchVars, SCIP_Real *branchBounds, SCIP_BOUNDTYPE *boundTypes );
   void           changeSearchStrategy(SCIP *scip);
public:
   ScipParaObjCommPointHdlr(
		    UG::ParaComm   *comm,
            ScipParaSolver *solver
         )
         : scip::ObjEventhdlr::ObjEventhdlr(solver->getScip(), "ScipParaObjCommPointHdlr", "Event handler to communicate with LC"),
           paraComm(comm), scipParaSolver(solver), scipToCheckRootSolvability(0), originalScip(0), needToSendNode(false),
           originalSelectionStrategy(true), previousNNodesSolved(0), previousLpIter(0),
           cloned(false), interrupting(false)
   {
      changeNodeSelName = scipParaSolver->getChangeNodeSelName();
      if( !cloned && scipParaSolver->getParaParamSet()->getBoolParamValue(UG::RootNodeSolvabilityCheck) )
      {
         /* initialize SCIP to check root solvability */
         SCIP_CALL_ABORT( SCIPcreate(&scipToCheckRootSolvability) );
         /* include default SCIP plugins */
         SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scipToCheckRootSolvability) );
         ScipParaInstance* scipParaInstance = dynamic_cast< ScipParaInstance* >(scipParaSolver->getParaInstance());
         scipParaInstance->createProblem(scipToCheckRootSolvability,
               solver->getParaParamSet()->getIntParamValue(UG::InstanceTransferMethod),
               solver->getParaParamSet()->getBoolParamValue(UG::NoPreprocessingInLC),
               NULL
         );   // LC presolving setting file should not be set, when it does this check!
      }
   }

   ScipParaObjCommPointHdlr(
		    UG::ParaComm   *comm,
            ScipParaSolver *solver,
            SCIP           *subScip,
            SCIP           *inOriginalScip,
            bool inCloned
         )          : scip::ObjEventhdlr::ObjEventhdlr(subScip, "ScipParaObjCommPointHdlr", "Event handler to communicate with LC"),
         paraComm(comm), scipParaSolver(solver), scipToCheckRootSolvability(0), originalScip(inOriginalScip), needToSendNode(false),
         originalSelectionStrategy(true), previousNNodesSolved(0), previousLpIter(0),
         cloned(inCloned), interrupting(false)
   {
   }

   /** destructor */
   ~ScipParaObjCommPointHdlr(
         )
   {
      if( scipToCheckRootSolvability )
      {
         SCIP_CALL_ABORT( SCIPfree(&scipToCheckRootSolvability) );
      }
   }

   /** clone method, used to copy plugins which are not constraint handlers or variable pricer plugins */
   ObjCloneable* clone(
      SCIP*           scip                /**< SCIP data structure */
      ) const
   {
      return new ScipParaObjCommPointHdlr(paraComm, scipParaSolver, scip, scipParaSolver->getScip(), true);
   }
   /** returns whether the objective plugin is copyable */
   SCIP_Bool iscloneable(
      void
      ) const
   {
      return true;
   }

   // SCIP *getParentScip(){ return parentScip; }
   bool isColne(){return cloned;}

   /** destructor of event handler to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of event handler (called after problem was transformed) */
   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {  /*lint --e{715}*/
      SCIP_CALL( SCIPcatchEvent( scip,
            ( SCIP_EVENTTYPE_GBDCHANGED |
                  SCIP_EVENTTYPE_BOUNDTIGHTENED |
                  SCIP_EVENTTYPE_LPEVENT |
                  SCIP_EVENTTYPE_ROWEVENT |
                  SCIP_EVENTTYPE_NODEFOCUSED |
                  SCIP_EVENTTYPE_BESTSOLFOUND
                  )
            , eventhdlr, NULL, NULL) );
      interrupting = false;
      return SCIP_OKAY;
   }

   /** deinitialization method of event handler (called before transformed problem is freed) */
   virtual SCIP_RETCODE scip_exit(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of event handler (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The event handler may use this call to initialize its branch and bound specific data.
    *
    */
   virtual SCIP_RETCODE scip_initsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of event handler (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The event handler should use this call to clean up its branch and bound data.
    */
   virtual SCIP_RETCODE scip_exitsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr           /**< the event handler itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** frees specific constraint data */
   virtual SCIP_RETCODE scip_delete(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr,          /**< the event handler itself */
      SCIP_EVENTDATA**   eventdata           /**< pointer to the event data to free */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** execution method of event handler
    *
    *  Processes the event. The method is called every time an event occurs, for which the event handler
    *  is responsible. Event handlers may declare themselves resposible for events by calling the
    *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
    *  given event handler and event data.
    */
   virtual SCIP_RETCODE scip_exec(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr,          /**< the event handler itself */
      SCIP_EVENT*        event,              /**< event to process */
      SCIP_EVENTDATA*    eventdata           /**< user data for the event */
      );
};

}; /* namespace ParaSCIP */

#endif // __SCIP_PARA_COMM_POINT_HDLR_H__
