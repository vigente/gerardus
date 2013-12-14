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

/**@file    scipParaSolver.h
 * @brief   ParaSolver extension for SCIP: Parallelized solver implementation for SCIP.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_SOLVER_H__
#define __SCIP_PARA_SOLVER_H__

#include <list>
#include "ug/paraSolver.h"
#include "scipUserPlugins.h"
#include "scipParaDiffSubproblem.h"
#include "scipDiffParamSet.h"

#define ENFORCED_THRESHOLD 5

namespace ParaSCIP
{
typedef struct LocalNodeInfo_t
{
   /********************************
    * for local cuts and conflicts *
    * *****************************/
   SCIP_Real      linearLhs;           /**< array of lhs */
   SCIP_Real      linearRhs;           /**< array of rhs */
   int            nLinearCoefs;         /**< array of number of coefficient values for linear constrains */
   SCIP_Real      *linearCoefs;         /**< array of non-zero coefficient values of linear constrains */
   int            *idxLinearCoefsVars;  /**< array of indices of no-zero coefficient values of linear constrains */
} LocalNodeInfo;
typedef LocalNodeInfo * LocalNodeInfoPtr;

class ScipParaSolver : public UG::ParaSolver
{
   SCIP                      *scip;
   SCIP                      *scipToCheckEffectOfRootNodeProcesses;
   ScipDiffParamSet          *scipDiffParamSetRoot;
   ScipDiffParamSet          *scipDiffParamSet;
   SCIP_MESSAGEHDLR          *messagehdlr;
   FILE                      *logfile;
   ScipDiffParamSet          *originalParamSet;
   std::list<LocalNodeInfoPtr> *conflictConsList;
   ScipUserPlugins           *userPlugins;

   int           nOrgVars;                   /**< number of original variables */
   SCIP_Real     *orgVarLbs;                 /**< array of original lower bound of variable */
   SCIP_Real     *orgVarUbs;                 /**< array of original upper bound of variable */

   SCIP_CONS** addedConss;
   SCIP_CONS* addedDualCons;

   void setRacingParams(UG::ParaRacingRampUpParamSet *inRacingParams, bool winnerParam);
   void setWinnerRacingParams(UG::ParaRacingRampUpParamSet *inRacingParams);
   void createSubproblem();
   void freeSubproblem();
   void solve();
   long long getNNodesSolved();
   int getNNodesLeft();
   double getDualBoundValue();

   void solveToCheckEffectOfRootNodePreprocesses();
   void saveOrgProblemBounds();

public:
   ScipParaSolver(int argc, char **argv,
         UG::ParaComm *comm, UG::ParaParamSet *paraParamSet, UG::ParaInstance *paraInstance, UG::ParaDeterministicTimer *detTimer);
   ScipParaSolver(int argc, char **argv,
         UG::ParaComm *comm, UG::ParaParamSet *paraParamSet, UG::ParaInstance *paraInstance, UG::ParaDeterministicTimer *detTimer, bool thread);
   ~ScipParaSolver();
   const char *getChangeNodeSelName(){
      if( paraParams->getIntParamValue(UG::NodeTransferMode) == 0 ) return "estimate";
      else return "bfs";
   }

   ScipParaDiffSubproblem *getParentDiffSubproblem(){ return dynamic_cast< ScipParaDiffSubproblem * >(currentNode->getDiffSubproblem() ); }
   void writeCurrentNodeProblem(const std::string& filename);
   void tryNewSolution(UG::ParaSolution *sol);
   void setLightWeightRootNodeProcess();
   void setOriginalRootNodeProcess();

   int getOffsetDepth()
   {
      if( currentNode->getDiffSubproblem() )
      {
         ScipParaDiffSubproblem *scipDiffSubproblem = dynamic_cast<ScipParaDiffSubproblem *>(currentNode->getDiffSubproblem());
         return scipDiffSubproblem->getOffset();
      }
      else
      {
         return 0;
      }
   }

   SCIP *getScip(){
      return scip;
   }

   std::list<LocalNodeInfoPtr> *getConflictConsList(
         )
   {
      return conflictConsList;
   }

   void writeSubproblem();

   long long getSimplexIter(
         )
   {
      SCIP_STAGE stage = SCIPgetStage(scip);
      if( stage == SCIP_STAGE_PRESOLVED || stage == SCIP_STAGE_SOLVING || stage == SCIP_STAGE_SOLVED )
      {
         return SCIPgetNLPIterations(scip);
      }
      else
      {
         return 0;
      }
   }

   int getNRestarts()
   {
      SCIP_STAGE stage = SCIPgetStage(scip);
      if( stage != SCIP_STAGE_INIT )
      {
         return SCIPgetNRuns(scip);
      }
      else
      {
         return 0;
      }
   }

   /** set user plugins */
   void setUserPlugins(ScipUserPlugins *inUi) { userPlugins = inUi; }

   /** include user plugins */
   void includeUserPlugins(SCIP *inScip)
   {
      if( userPlugins )
      {
         (*userPlugins)(inScip);
      }
   }

   bool wasTerminatedNormally(
         )
   {
      return true;
   }

};

}

#endif // __SCIP_PARA_SOLVER_H__
