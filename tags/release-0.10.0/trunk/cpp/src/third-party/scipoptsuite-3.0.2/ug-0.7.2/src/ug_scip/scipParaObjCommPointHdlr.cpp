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

/**@file    scipParaObjCommPointHdlr.cpp
 * @brief   Event handlr for communication point.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <string>
#include <sstream>
#include <iostream>

#include "scipParaObjCommPointHdlr.h"
#include "scip/struct_nodesel.h"

#include "scipParaSolution.h"
#include "scipParaDiffSubproblem.h"

using namespace ParaSCIP;

void
ScipParaObjCommPointHdlr::processNewSolution(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENT*        event               /**< event to process */
      )
{
   SCIP_SOL *sol = SCIPeventGetSol(event);
   int nVars = SCIPgetNOrigVars(scip);
   SCIP_VAR **vars = SCIPgetOrigVars(scip);
   SCIP_Real *vals = new SCIP_Real[nVars];
   SCIP_CALL_ABORT( SCIPgetSolVals(scip, sol, nVars, vars, vals) );

   DEF_SCIP_PARA_COMM( scipParaComm, paraComm);
   scipParaSolver->saveIfImprovedSolutionWasFound(
         scipParaComm->createScipParaSolution(
               SCIPgetSolOrigObj(scip, sol),
               nVars,
               vars,
               vals
               )
         );
   SCIP_CALL_ABORT( SCIPprintSol(scip,sol,NULL, FALSE) );
   SCIP_Bool feasible;
   SCIP_CALL_ABORT( SCIPcheckSolOrig(scip,sol,&feasible,1,1 ) );

   delete [] vals;
}

SCIP_RETCODE
ScipParaObjCommPointHdlr::scip_exec(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_EVENTHDLR*    eventhdlr,          /**< the event handler itself */
      SCIP_EVENT*        event,              /**< event to process */
      SCIP_EVENTDATA*    eventdata           /**< user data for the event */
      )
{
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), "ScipParaObjCommPointHdlr") == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) &
         ( SCIP_EVENTTYPE_PRESOLVEROUND |
         SCIP_EVENTTYPE_NODEFOCUSED |
         SCIP_EVENTTYPE_LPEVENT |
         SCIP_EVENTTYPE_ROWEVENT |
         SCIP_EVENTTYPE_BESTSOLFOUND
         )
         );

#if SCIP_VERSION == 211 && SCIP_SUBVERSION == 0
   if( SCIPgetStage(scip) >= SCIP_STAGE_FREESOLVE )
#else
   if( SCIPgetStage(scip) >= SCIP_STAGE_EXITSOLVE )
#endif
   {
      return  SCIP_OKAY;
   }

   double detTime = -1.0;
   SCIP *subScip = 0;
   if( cloned )
   {
      subScip = scip;
      scip = originalScip;
   }

   if( scipParaSolver->getParaParamSet()->getBoolParamValue(UG::Deterministic) )
   {
      assert(scipParaSolver->getDeterministicTimer());
      if( scipParaSolver->getParaParamSet()->getBoolParamValue(UG::EventWeightedDeterministic) )
      {
         if( SCIPeventGetType(event) == SCIP_EVENTTYPE_FIRSTLPSOLVED )
         {
            SCIP_Longint totalLpIter = 0;
            if( cloned )
            {
               totalLpIter = SCIPgetNLPIterations(subScip);
            }
            else
            {
               totalLpIter = SCIPgetNLPIterations(scip);
            }
            SCIP_Longint lpIter = totalLpIter - previousLpIter;
            if( lpIter < 0 )
            {
               lpIter = totalLpIter;
            }
            scipParaSolver->getDeterministicTimer()->update(lpIter);
            previousLpIter = totalLpIter;

         }
         else if ( SCIPeventGetType(event) == SCIP_EVENTTYPE_LPSOLVED )
         {
            SCIP_Longint totalLpIter = 0;
            if( cloned )
            {
               totalLpIter = SCIPgetNLPIterations(subScip);
            }
            else
            {
               totalLpIter = SCIPgetNLPIterations(scip);
            }
            SCIP_Longint lpIter = totalLpIter - previousLpIter;
            if( lpIter < 0 )
            {
               lpIter = totalLpIter;
            }
            scipParaSolver->getDeterministicTimer()->update(lpIter);
            previousLpIter = totalLpIter;
         }
         else if ( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND
               || SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEFOCUSED )
         {
         }
         else
         {
            scipParaSolver->getDeterministicTimer()->update(1.0);
         }
      }
      else
      {
         if ( SCIPeventGetType(event) != SCIP_EVENTTYPE_BESTSOLFOUND
               && SCIPeventGetType(event) != SCIP_EVENTTYPE_NODEFOCUSED )
         {
            scipParaSolver->getDeterministicTimer()->update(1.0);
         }
      }
      detTime = scipParaSolver->getDeterministicTimer()->getElapsedTime();
      if( ( detTime - scipParaSolver->getPreviousCommTime() )
            < ( scipParaSolver->getParaParamSet()->getRealParamValue(UG::NotificationInterval) ) )
      {
         if( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND )
         {
            if( cloned )
            {
               return SCIP_OKAY;
            }
            else
            {
               processNewSolution(scip, event);
            }
         }
      }
      do
      {
         scipParaSolver->iReceiveMessages();
      } while( !scipParaSolver->waitToken(scipParaSolver->getRank()) );
   }
   else
   {
      if( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND )
      {
         if( cloned )
         {
            return  SCIP_OKAY;
         }
         else
         {
            processNewSolution(scip, event);
         }
      }
   }

   //
   //  the following interrupt routine moved from the inside of the following if clause
   //  and  scipParaSolver->iReceiveMessages() before solution sending is also removed (comment out)
   //
   /*********************************************************************************************
     * check if there are messages: this is the first thing to do for communication *
     ********************************************************************************************/
   scipParaSolver->iReceiveMessages();

   /*********************************************************************************************
     * update solution. This have to do right after receiving message *
     ********************************************************************************************/
   scipParaSolver->sendLocalSolution();   // if local solution exists, it should be sent right after iReceiveMessages
   if( !cloned )
   {
      scipParaSolver->updatePendingSolution();
   }

   if( scipParaSolver->newParaNodeExists() ||
         scipParaSolver->isInterrupting() ||
         scipParaSolver->isRacingWinnerParamReceived() ||
         scipParaSolver->isTerminationRequested() )
   {
      if( !interrupting  )
      {
         if( cloned )
         {
            SCIP_CALL_ABORT( SCIPinterruptSolve(subScip) );
         }
         else
         {
            SCIP_CALL_ABORT( SCIPinterruptSolve(scip) );
         }
         interrupting = true;
      }
      return SCIP_OKAY;
   }

   SCIP_Longint nNodes = 0;
   if( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING
         && SCIPgetStage(scip) != SCIP_STAGE_INITSOLVE )
   {
      nNodes = SCIPgetNNodes(scip);
   }

   if( previousNNodesSolved == nNodes )
   {
      if( scipParaSolver->getParaParamSet()->getBoolParamValue(UG::Deterministic)
            && ( !scipParaSolver->isInCollectingMode() )
            && ( !scipParaSolver->newParaNodeExists() )
            && ( !scipParaSolver->isInterrupting() )
            && ( !scipParaSolver->isRacingWinnerParamReceived() )
            && ( !scipParaSolver->isTerminationRequested() )
            )
      {
         scipParaSolver->passToken(scipParaSolver->getRank());
         scipParaSolver->setPreviousCommTime(detTime);
      }

      return SCIP_OKAY;
   }
   previousNNodesSolved = nNodes;

   /** if root node is solved, set root node time */
   if( nNodes == 2 )
   {
      /** when a problem is solved at root,
       * its root node process time is set on paraSolver main loop */
      scipParaSolver->setRootNodeTime();
      scipParaSolver->setRootNodeSimplexIter(scipParaSolver->getSimplexIter());
      if( scipParaSolver->getCurrentSolivingNodeMergingStatus() == 0 &&
            SCIPgetDualbound(scip) <
            scipParaSolver->getCurrentSolvingNodeInitialDualBound() * ( 1.0 - scipParaSolver->getParaParamSet()->getRealParamValue(UG::AllowableRegressionRatioInMerging)) )
      {
         if( !interrupting  )
         {
            if( cloned )
            {
               SCIP_CALL_ABORT( SCIPinterruptSolve(subScip) );
            }
            else
            {
               scipParaSolver->getCurrentNode()->setMergingStatus(3);  // cannot be merged
               SCIP_CALL_ABORT( SCIPinterruptSolve(scip) );
            }
            interrupting = true;
         }
         return SCIP_OKAY;
      }
   }

   /*****************************************************************************
    * sends solver state message if it is necessary, that is,                   *
    * notification interval time has been passed from the previous notification *
    *****************************************************************************/
   double bestDualBoundValue = -SCIPinfinity(scip);
   if( SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMED ||
         SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING ||
         SCIPgetStage(scip) == SCIP_STAGE_INITSOLVE
         )
   {
      bestDualBoundValue = scipParaSolver->getCurrentNode()->getDualBoundValue();
   }
   else
   {
      bestDualBoundValue = SCIPgetDualbound(scip);
   }

   /************************************************************
    *  check if global best incumbent value is updated or not. *
    *  if it is updated, set cutoff value                      *
    ************************************************************/
   if( scipParaSolver->isGlobalIncumbentUpdated() &&
         SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING
                  && SCIPgetStage(scip) != SCIP_STAGE_INITSOLVE )
   {
      if( !cloned
            && SCIPeventGetType(event) != SCIP_EVENTTYPE_BOUNDTIGHTENED
            )
      {
         /** set cutoff value */
         // if( scipParaSolver->getGlobalBestIncumbentValue() > bestDualBoundValue )
         // {
         SCIP_CALL_ABORT( SCIPsetObjlimit(scip, scipParaSolver->getGlobalBestIncumbentValue()) );
         // }
         scipParaSolver->globalIncumbnetValueIsReflected();
      }
   }

   if( scipParaSolver->getGlobalBestIncumbentValue() < bestDualBoundValue )
   {
      if( !interrupting  )
      {
         if( cloned )
         {
            if( SCIPgetStage(subScip) != SCIP_STAGE_INITSOLVE )
            {
               SCIP_CALL_ABORT( SCIPinterruptSolve(subScip) );
               interrupting = true;
            }

         }
         else
         {
            if( SCIPgetStage(scip) != SCIP_STAGE_INITSOLVE )
            {
               SCIP_CALL_ABORT( SCIPinterruptSolve(scip) );
               interrupting = true;
            }
         }
      }
      return SCIP_OKAY;
   }

   if ( needToSendNode
         || scipParaSolver->notificationIsNecessary()
         || nNodes == 2
         || scipParaSolver->getParaParamSet()->getBoolParamValue(UG::Deterministic) )  // always notify
   {
      if( bestDualBoundValue >= -1e+10 )  // only after some dual bound value is computed, notify its status
      {
         int nNodesLeft = 0;
         if( SCIPgetStage(scip) != SCIP_STAGE_TRANSFORMED &&
               SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING &&
               SCIPgetStage(scip) != SCIP_STAGE_INITSOLVE )
         {
            nNodesLeft = SCIPgetNNodesLeft(scip);
         }
         scipParaSolver->sendSolverState(nNodes, nNodesLeft, bestDualBoundValue, detTime);
         scipParaSolver->waitMessageIfNecessary();
         if( scipParaSolver->getParaParamSet()->getBoolParamValue(UG::Deterministic)
               && ( !needToSendNode )
               && ( !scipParaSolver->isManyNodesCollectionRequested() )
               && !( !scipParaSolver->isRacingStage()
                     &&  scipParaSolver->isAggressivePresolvingSpecified()
                     &&  scipParaSolver->getSubMipDepth() <
                     ( scipParaSolver->getOffsetDepth() + scipParaSolver->getAggresivePresolvingStopDepth() ) &&
                       SCIPgetDepth(scip) > scipParaSolver->getAggresivePresolvingDepth() )
           )
         {
            scipParaSolver->passToken(scipParaSolver->getRank());
            scipParaSolver->setPreviousCommTime(detTime);
         }
      }
   }

   if( scipParaSolver->getNStopSolvingMode() > 0 &&
         !scipParaSolver->isRacingStage() &&
         scipParaSolver->getCurrentSolivingNodeMergingStatus() != 0 &&
         ( nNodes < scipParaSolver->getNStopSolvingMode() ) &&
         ( bestDualBoundValue - scipParaSolver->getLcBestDualBoundValue() ) > 0.0 &&
          REALABS( ( bestDualBoundValue - scipParaSolver->getLcBestDualBoundValue() )
         / scipParaSolver->getLcBestDualBoundValue() ) > scipParaSolver->getBoundGapForStopSolving() )
   {
      if( scipParaSolver->getBigDualGapSubtreeHandlingStrategy() == 0 )
      {
         scipParaSolver->sendAnotherNodeRequest(bestDualBoundValue);   // throw away nodes
      }
      else if( scipParaSolver->getBigDualGapSubtreeHandlingStrategy() == 1 )
      {
         scipParaSolver->setSendBackAllNodes();
      }
      else
      {
         THROW_LOGICAL_ERROR2("Invalid big dual gap handling strategy. startegy = ",
               scipParaSolver->getBigDualGapSubtreeHandlingStrategy() );
      }
   }
   /*******************************************************************
    * if new ParaNode was received or Interrupt message was received, *
    * stop solving the current ParaNode                               *
    ******************************************************************/
   if( cloned ||
         ( SCIPeventGetType(event) != SCIP_EVENTTYPE_NODEFOCUSED
         && !scipParaSolver->isInCollectingMode() ) )
   {
      return SCIP_OKAY; // process until SCIP_EVENTTYPE_NODEFOCUSED event
   }

   /******************************************************************
    * if node depth is smaller than that specified by parameter,
    * not to send nodes and keep them.
    *****************************************************************/
   if( SCIPnodeGetDepth( SCIPgetCurrentNode( scip ) )
         < scipParaSolver->getParaParamSet()->getIntParamValue(UG::KeepNodesDepth) &&
         !scipParaSolver->isCollectingAllNodes() )
   {
      return SCIP_OKAY;
   }

   /*******************************************************************
    * if ParaNode transfer has been requested, send a ParaNode        *
    ******************************************************************/
   if( ( needToSendNode &&
           (
                SCIPgetNNodesLeft( scip ) > scipParaSolver->getThresholdValue(SCIPgetNNodes(scip)) ||
                ( SCIPgetNNodesLeft( scip ) > 1 &&
                      ( scipParaSolver->isAggressiveCollecting() ||
                            ( !scipParaSolver->isRampUp() &&
                                  ( bestDualBoundValue - scipParaSolver->getLcBestDualBoundValue() < 0.0 ||
                                         ( ( bestDualBoundValue - scipParaSolver->getLcBestDualBoundValue() ) > 0.0 &&
                                                          ( REALABS( ( bestDualBoundValue - scipParaSolver->getLcBestDualBoundValue() )
                                                                / scipParaSolver->getLcBestDualBoundValue() ) < scipParaSolver->getBoundGapForCollectingMode() )
                                          )
                                   )
                            )
                      )
                )
           )
       )
       || scipParaSolver->isManyNodesCollectionRequested()
       || ( !scipParaSolver->isRacingStage() &&
             scipParaSolver->isAggressivePresolvingSpecified() &&
             scipParaSolver->getSubMipDepth() <
             ( scipParaSolver->getOffsetDepth() + scipParaSolver->getAggresivePresolvingStopDepth() ) &&
             SCIPgetDepth(scip) > scipParaSolver->getAggresivePresolvingDepth() )
   )
   {
      checkRootNodeSolvabilityAndSendParaNode(scip);
      if( ( scipParaSolver->isManyNodesCollectionRequested() && (!scipParaSolver->isCollectingAllNodes() ) ) &&
            scipParaSolver->isBreaking() &&
            ( scipParaSolver->isTransferLimitReached() ||
                  ( bestDualBoundValue > scipParaSolver->getTargetBound() ) ) )
      {
         scipParaSolver->resetBreakingInfo();
      }
   }

   if( !scipParaSolver->isManyNodesCollectionRequested() )
   {
      if( !scipParaSolver->isRampUp() )
      {
         if( ( !scipParaSolver->isWarmStarted() ) && scipParaSolver->isRacingRampUp() )
         {
            if( scipParaSolver->isRacingWinner() || !scipParaSolver->isRacingStage() )
            {
               if( originalSelectionStrategy )
               {
                  changeSearchStrategy(scip);
               }
               needToSendNode = true;
            }
         }
         else
         {
            if( ( scipParaSolver->getParaParamSet()->getBoolParamValue(UG::BreakFirstSubtree)
                  && paraComm->getRank() == 1 ) ||
                  bestDualBoundValue - scipParaSolver->getLcBestDualBoundValue() < 0.0 ||
                  ( ( bestDualBoundValue - scipParaSolver->getLcBestDualBoundValue() ) > 0.0 &&
                        ( REALABS( ( bestDualBoundValue - scipParaSolver->getLcBestDualBoundValue() )
                              / scipParaSolver->getLcBestDualBoundValue() ) < scipParaSolver->getBoundGapForCollectingMode() ) ) )
            {
               if( originalSelectionStrategy )
               {
                  changeSearchStrategy(scip);
               }
               if( needToSendNode ) needToSendNode = false;
               else needToSendNode = true;
            }
            else
            {
               needToSendNode = false;
            }
         }
      }
      else
      {
         /*********************************************
          * check if solver is in collecting mode     *
          ********************************************/
         if( scipParaSolver->isInCollectingMode() )
         {
            if( originalSelectionStrategy )
            {
               changeSearchStrategy(scip);
            }
            if( scipParaSolver->getNSendInCollectingMode() > 0 || scipParaSolver->isAggressiveCollecting() )
            {
               needToSendNode = true;
            }
            else
            {
               if( needToSendNode ) needToSendNode = false;
               else needToSendNode = true;
            }
         }
         else
         {
            if( !originalSelectionStrategy )
            {
               SCIP_NODESEL *nodesel = SCIPgetNodesel(scip);
               SCIP_CALL_ABORT( SCIPsetNodeselStdPriority(scip, nodesel, originalPriority ) );
               originalSelectionStrategy = true;
            }
            needToSendNode = false;   // In iReceive, more than on collecting mode message may be received
         }
      }
   }

   return SCIP_OKAY;
}

void
ScipParaObjCommPointHdlr::checkRootNodeSolvabilityAndSendParaNode(
      SCIP*  scip
      )
{
   SCIP_NODE* node = SCIPgetCurrentNode( scip );
   int depth = SCIPnodeGetDepth( node );
   SCIP_VAR **branchVars = new SCIP_VAR*[depth];
   SCIP_Real *branchBounds = new SCIP_Real[depth];
   SCIP_BOUNDTYPE *boundTypes = new  SCIP_BOUNDTYPE[depth];
   int nBranchVars;
   SCIPnodeGetAncestorBranchings( node, branchVars, branchBounds, boundTypes, &nBranchVars, depth );
   if( nBranchVars > depth )  // did not have enough memory, then reallocate
   {
      delete [] branchVars;
      delete [] branchBounds;
      delete [] boundTypes;
      branchVars = new SCIP_VAR*[nBranchVars];
      branchBounds = new SCIP_Real[nBranchVars];
      boundTypes = new  SCIP_BOUNDTYPE[nBranchVars];
      SCIPnodeGetAncestorBranchings( node, branchVars, branchBounds, boundTypes, &nBranchVars, nBranchVars );
   }

   if( scipParaSolver->getParaParamSet()->getBoolParamValue(UG::AllBoundChangesTransfer) &&
         !( scipParaSolver->getParaParamSet()->getBoolParamValue(UG::NoAllBoundChangesTransferInRacing) &&
               scipParaSolver->isRacingStage() ) )
   {
      int nVars = SCIPgetNVars(scip);
      SCIP_VAR **vars = SCIPgetVars(scip);
      int *iBranchVars = new int[nBranchVars];
      /* create the variable mapping hash map */
      SCIP_HASHMAP* varmapLb;
      SCIP_HASHMAP* varmapUb;
      SCIP_CALL_ABORT( SCIPhashmapCreate(&varmapLb, SCIPblkmem(scip), nVars) );
      SCIP_CALL_ABORT( SCIPhashmapCreate(&varmapUb, SCIPblkmem(scip), nVars) );
      for( int i = 0; i < nBranchVars; i++ )
      {
         iBranchVars[i] = i;
         if( boundTypes[i] == SCIP_BOUNDTYPE_LOWER )
         {
            SCIP_CALL_ABORT( SCIPhashmapInsert(varmapLb, branchVars[i], &iBranchVars[i] ) );
         }
         else
         {
            SCIP_CALL_ABORT( SCIPhashmapInsert(varmapUb, branchVars[i], &iBranchVars[i] ) );
         }
      }
      SCIP_VAR **preBranchVars = branchVars;
      SCIP_Real *preBranchBounds = branchBounds;
      SCIP_BOUNDTYPE *preBboundTypes = boundTypes;
      branchVars = new SCIP_VAR*[nBranchVars+nVars*2];
      branchBounds = new SCIP_Real[nBranchVars+nVars*2];
      boundTypes = new  SCIP_BOUNDTYPE[nBranchVars+nVars*2];
      for( int i = 0; i < nBranchVars; i++ )
      {
         branchVars[i] = preBranchVars[i];
         branchBounds[i] = preBranchBounds[i];
         boundTypes[i] = preBboundTypes[i];
      }
      int *iBranchVar = NULL;
      for( int i = 0; i < nVars; i++ )
      {
         iBranchVar =  (int *)SCIPhashmapGetImage(varmapLb, vars[i]);
         if( iBranchVar )
         {
            // assert( EPSGE(preBranchBounds[*iBranchVar], SCIPvarGetLbLocal(vars[i]) ,DEFAULT_NUM_EPSILON ) );
            if( EPSLT(preBranchBounds[*iBranchVar], SCIPvarGetLbLocal(vars[i]) ,DEFAULT_NUM_EPSILON ) )
            {
               branchBounds[*iBranchVar] = SCIPvarGetLbLocal(vars[i]);  // node is current node
            }
         }
         else
         {
            if( EPSGT( SCIPvarGetLbLocal(vars[i]), SCIPvarGetLbGlobal(vars[i]), MINEPSILON ) )
            {
               branchVars[nBranchVars] = vars[i];
               branchBounds[nBranchVars] = SCIPvarGetLbLocal(vars[i]);
               boundTypes[nBranchVars] = SCIP_BOUNDTYPE_LOWER;
               nBranchVars++;
            }
         }
         iBranchVar = (int *)SCIPhashmapGetImage(varmapUb, vars[i]);
         if( iBranchVar )
         {
            // assert( EPSLE(preBranchBounds[*iBranchVar], SCIPvarGetUbLocal(vars[i]) ,DEFAULT_NUM_EPSILON ) );
            if( EPSGT(preBranchBounds[*iBranchVar], SCIPvarGetUbLocal(vars[i]) ,DEFAULT_NUM_EPSILON ) )
            {
               branchBounds[*iBranchVar] = SCIPvarGetUbLocal(vars[i]); // node is current node
            }
         }
         else
         {
            if( EPSLT( SCIPvarGetUbLocal(vars[i]), SCIPvarGetUbGlobal(vars[i]), MINEPSILON ) )
            {
               branchVars[nBranchVars] = vars[i];
               branchBounds[nBranchVars] = SCIPvarGetUbLocal(vars[i]);
               boundTypes[nBranchVars] = SCIP_BOUNDTYPE_UPPER;
               nBranchVars++;
            }
         }
      }
      SCIPhashmapFree(&varmapLb);
      SCIPhashmapFree(&varmapUb);
      delete [] preBranchVars;
      delete [] preBranchBounds;
      delete [] preBboundTypes;
      delete [] iBranchVars;
   }

   /** check root node solvability */
   if( scipParaSolver->getParaParamSet()->getBoolParamValue(UG::RootNodeSolvabilityCheck) )
   {
      SCIP_CALL_ABORT( SCIPtransformProb(scipToCheckRootSolvability) );
      SCIP_VAR **copyVars = SCIPgetVars(scipToCheckRootSolvability);
      for(int v = 0; v < nBranchVars; v++)
      {
         int index = SCIPvarGetIndex(branchVars[v]);
         assert(index != -1);
         assert(std::string(SCIPvarGetName(branchVars[v]))==std::string(SCIPvarGetName(copyVars[index])));
         if( boundTypes[v] == SCIP_BOUNDTYPE_LOWER )
         {
            SCIP_CALL_ABORT(SCIPchgVarLbGlobal(scipToCheckRootSolvability,copyVars[index], branchBounds[v]));
         }
         else if (boundTypes[v] == SCIP_BOUNDTYPE_UPPER)
         {
            SCIP_CALL_ABORT(SCIPchgVarUbGlobal(scipToCheckRootSolvability,copyVars[index], branchBounds[v]));
         }
         else
         {
            THROW_LOGICAL_ERROR2("Invalid bound type: type = ", static_cast<int>(boundTypes[v]));
         }
      }
      SCIP_CALL_ABORT(SCIPsetLongintParam(scipToCheckRootSolvability,"limits/nodes", 1));
      SCIP_CALL_ABORT(SCIPsetObjlimit(scipToCheckRootSolvability, scipParaSolver->getGlobalBestIncumbentValue()));
      SCIP_CALL_ABORT(SCIPsolve(scipToCheckRootSolvability));

      SCIP_STATUS status = SCIPgetStatus(scipToCheckRootSolvability);

      switch(status)
      {
      case SCIP_STATUS_OPTIMAL :
      {
         SCIP_SOL *bestSol = SCIPgetBestSol( scip );
         int nVars = SCIPgetNOrigVars(scip);
         SCIP_VAR **vars = SCIPgetOrigVars(scip);
         SCIP_Real *vals = new SCIP_Real[nVars];
         SCIP_CALL_ABORT( SCIPgetSolVals(scip, bestSol, nVars, vars, vals) );
         DEF_SCIP_PARA_COMM( scipParaComm, paraComm);
         scipParaSolver->saveIfImprovedSolutionWasFound(
               scipParaComm->createScipParaSolution(
                     SCIPgetSolOrigObj(scip, bestSol),
                     nVars,
                     vars,
                     vals
                  )
               );
         delete [] vals;
         /** remove the node sent from SCIP environment */
         SCIP_CALL_ABORT( SCIPcutoffNode( scip, node) );
         needToSendNode = false;   // This means that a node is sent in the next time again,
                                   // because this flag is flipped when collecting mode is checked
         scipParaSolver->countInPrecheckSolvedParaNodes();
         break;
      }
      case SCIP_STATUS_INFEASIBLE :
      case SCIP_STATUS_INFORUNBD :
      {
         /** remove the node sent from SCIP environment */
         SCIP_CALL_ABORT( SCIPcutoffNode( scip, node ) );
         needToSendNode = false;   // This means that a node is sent in the next time again,
                                    // because this flag is flipped when collecting mode is checked
         scipParaSolver->countInPrecheckSolvedParaNodes();
         break;
      }
      case SCIP_STATUS_NODELIMIT :
      {
         sendNode(scip, node, depth, nBranchVars, branchVars, branchBounds, boundTypes);
         break;
      }
      default:
         THROW_LOGICAL_ERROR2("Invalid status after root solvability check: status = ", static_cast<int>(status));
      }

      SCIP_CALL_ABORT( SCIPfreeTransform(scipToCheckRootSolvability) );
   }
   else
   {
      sendNode(scip, node, depth, nBranchVars, branchVars, branchBounds, boundTypes);
   }

   delete [] branchVars;
   delete [] branchBounds;
   delete [] boundTypes;
}

void
ScipParaObjCommPointHdlr::sendNode(
      SCIP *scip,
      SCIP_NODE* node,
      int depth,
      int nNewBranchVars,
      SCIP_VAR **newBranchVars,
      SCIP_Real *newBranchBounds,
      SCIP_BOUNDTYPE *newBoundTypes
      )
{

   DEF_SCIP_PARA_COMM( scipParaComm, paraComm);
   ScipParaDiffSubproblem *diffSubproblem = scipParaComm->createScipParaDiffSubproblem(
         scip,
         scipParaSolver,
         nNewBranchVars,
         newBranchVars,
         newBranchBounds,
         newBoundTypes
         );

   long long n = SCIPnodeGetNumber( node );
   double dualBound = SCIPgetDualbound(scip);
   assert(SCIPisFeasGE(scip, SCIPnodeGetLowerbound(node) , SCIPgetLowerbound(scip)));
   double estimateValue = SCIPnodeGetEstimate( node );
   assert(SCIPisFeasGE(scip, estimateValue, SCIPnodeGetLowerbound(node) ));
   scipParaSolver->sendParaNode(n, depth, dualBound, estimateValue, diffSubproblem);

   /** remove the node sent from SCIP environment */
   SCIP_CALL_ABORT( SCIPcutoffNode( scip, node) );
}

void
ScipParaObjCommPointHdlr::changeSearchStrategy(
      SCIP*  scip
      )
{
   int numnodesels = SCIPgetNNodesels( scip );
   SCIP_NODESEL** nodesels = SCIPgetNodesels( scip );
   int i;
   for( i = 0; i < numnodesels; ++i )
   {
      std::string nodeselname(SCIPnodeselGetName(nodesels[i]));
      if( std::string(nodeselname) == std::string(changeNodeSelName) )
      {
         originalPriority = SCIPnodeselGetStdPriority(nodesels[i]);
         SCIP_CALL_ABORT( SCIPsetNodeselStdPriority(scip, nodesels[i], 536870911 ) );
         originalSelectionStrategy = false;
         break;
      }
   }
   assert( i != numnodesels );
}
