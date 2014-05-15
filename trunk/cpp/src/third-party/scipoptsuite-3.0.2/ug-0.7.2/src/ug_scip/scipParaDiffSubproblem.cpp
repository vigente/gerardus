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

/**@file    scipParaDiffSubproblem.cpp
 * @brief   ParaDiffSubproblem extension for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "ug/paraComm.h"
#include "scipParaSolver.h"
#include "scipParaInstance.h"

using namespace UG;
using namespace ParaSCIP;

ScipParaDiffSubproblem::ScipParaDiffSubproblem(
      SCIP *scip,
      ScipParaSolver *scipParaSolver,
      int nNewBranchVars,
      SCIP_VAR **newBranchVars,
      SCIP_Real *newBranchBounds,
      SCIP_BOUNDTYPE *newBoundTypes
      ) : localInfoIncluded(0),
      nLinearConss(0), linearLhss(0), linearRhss(0),
      nLinearCoefs(0), linearCoefs(0), idxLinearCoefsVars(0),
      offset(0), nVarBranchStats(0), idxLBranchStatsVars(0), downpscost(0), uppscost(0),
      downvsids(0), upvsids(0), downconflen(0), upconflen(0), downinfer(0), upinfer(0),
      downcutoff(0), upcutoff(0)
{
   nBoundChanges = nNewBranchVars;
   int nParentBranchVars = 0;
   int i = 0;

   ScipParaDiffSubproblem *parentDiffSubproblem = scipParaSolver->getParentDiffSubproblem();

   if( parentDiffSubproblem )
   {
      nParentBranchVars = parentDiffSubproblem->getNBoundChanges();
      nBoundChanges += nParentBranchVars;
   }
   indicesAmongSolvers = new int[nBoundChanges];
   branchBounds = new SCIP_Real[nBoundChanges];
   boundTypes = new SCIP_BOUNDTYPE[nBoundChanges];
   if( parentDiffSubproblem )
   {
      for( i = 0; i < nParentBranchVars; i ++ )
      {
         indicesAmongSolvers[i] = parentDiffSubproblem->getIndex(i);
         branchBounds[i] = parentDiffSubproblem->getBranchBound(i);
         boundTypes[i] = parentDiffSubproblem->getBoundType(i);
      }
   }

   for( int v = nNewBranchVars -1 ; v >= 0; --v )
   {
      SCIP_VAR *transformVar = newBranchVars[v];
      SCIP_Real scalar = 1.0;
      SCIP_Real constant = 0.0;
      SCIP_CALL_ABORT( SCIPvarGetOrigvarSum(&transformVar, &scalar, &constant ) );
      if( transformVar == NULL ) continue;
      assert( (!scipParaSolver->getParaParamSet()->getBoolParamValue(NoSolverPresolvingAtRoot) ) ||
              (!scipParaSolver->getCurrentNode()->isRootNode() ) ||
              ( scipParaSolver->getParaParamSet()->getBoolParamValue(NoSolverPresolvingAtRoot)  &&
                  scipParaSolver->getCurrentNode()->isRootNode() &&
                  scalar == 1.0 && constant == 0.0 )
            );
      indicesAmongSolvers[i] = SCIPvarGetIndex(transformVar);
      branchBounds[i] = ( newBranchBounds[v] - constant ) / scalar;
      if( scalar > 0.0 )
      {
         boundTypes[i] = newBoundTypes[v];
      }
      else
      {
         boundTypes[i] = SCIP_BoundType(1 - newBoundTypes[v]);
      }
      if( SCIPvarGetType(transformVar) != SCIP_VARTYPE_CONTINUOUS )
      {
         if( boundTypes[i] == SCIP_BOUNDTYPE_LOWER )
         {
            branchBounds[i] = SCIPfeasCeil(scip, branchBounds[i]);
         }
         else
         {
            branchBounds[i] = SCIPfeasFloor(scip, branchBounds[i]);
         }
      }
      i++;
   }
   nBoundChanges = i;

   if( scipParaSolver->getParaParamSet()->getBoolParamValue(TransferLocalCuts) ||
         scipParaSolver->getParaParamSet()->getBoolParamValue(TransferConflicts)
         )
   {
      addLocalNodeInfo(scip, scipParaSolver);
   }

   if( scipParaSolver->getParaParamSet()->getBoolParamValue(TransferBranchStats) &&
         (!( scipParaSolver->isRacingStage()
               && scipParaSolver->getParaParamSet()->getBoolParamValue(RacingStatBranching) ) )
         )
   {
      addBranchVarStats(scip, scipParaSolver);
   }
}

void
ScipParaDiffSubproblem::addLocalNodeInfo(
      SCIP *scip,
      ScipParaSolver *scipParaSolver
      )
{
   std::list<LocalNodeInfoPtr> localCutsList;

   if( scipParaSolver->getParaParamSet()->getBoolParamValue(TransferLocalCuts) )
   {
      SCIP_CUT** cuts;
      int ncuts;

      cuts = SCIPgetPoolCuts(scip);
      ncuts = SCIPgetNPoolCuts(scip);

      for( int c = 0; c < ncuts; ++c )
      {
         SCIP_ROW* row;

         row = SCIPcutGetRow(cuts[c]);
         assert(!SCIProwIsLocal(row));
         assert(!SCIProwIsModifiable(row));
         if( SCIPcutGetAge(cuts[c]) == 0 && SCIProwIsInLP(row) )
         {
            SCIP_COL** cols;
            int ncols;
            int i;

            /* create a linear constraint out of the cut */
            cols = SCIProwGetCols(row);
            ncols = SCIProwGetNNonz(row);

            LocalNodeInfo *localNodeInfo = new LocalNodeInfo;
            localNodeInfo->nLinearCoefs = ncols;
            localNodeInfo->idxLinearCoefsVars = new int[ncols];
            localNodeInfo->linearCoefs = new double[ncols];
            double lhs, rhs;
            if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
            {
               lhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
            }
            else
            {
               lhs = -SCIPinfinity(scip);
            }
            if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) )
            {
               rhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);
            }
            else
            {
               rhs = SCIPinfinity(scip);
            }

            SCIP_Real *vals = SCIProwGetVals(row);
            for( i = 0; i < ncols; ++i )
            {
               SCIP_VAR *transformVar = SCIPcolGetVar(cols[i]);
               SCIP_Real scalar = vals[i];
               SCIP_Real constant = 0.0;
               if( SCIPvarGetOrigvarSum(&transformVar, &scalar, &constant ) ==  SCIP_INVALIDDATA )
                  break;
               assert(transformVar != NULL);
               if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
               {
                  lhs -= constant;
               }
               if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) )
               {
                  rhs -= constant;
               }
               localNodeInfo->idxLinearCoefsVars[i] = SCIPvarGetIndex(transformVar);
               localNodeInfo->linearCoefs[i] = scalar;
            }
            if( i == ncols )
            {
               assert( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) == !SCIPisInfinity(scip, -lhs)  );
               assert( !SCIPisInfinity(scip, SCIProwGetRhs(row)) == !SCIPisInfinity(scip, rhs)  );
               localNodeInfo->linearLhs = lhs;
               localNodeInfo->linearRhs = rhs;
               localCutsList.push_back(localNodeInfo);
            }
            else
            {
               delete [] localNodeInfo->idxLinearCoefsVars;
               delete [] localNodeInfo->linearCoefs;
               delete localNodeInfo;
            }
         }
      }
   }

   ScipParaDiffSubproblem *parentDiffSubproblem = scipParaSolver->getParentDiffSubproblem();
   std::list<LocalNodeInfoPtr> *conflictConsList = scipParaSolver->getConflictConsList();

   int nParentConss = 0;
   if( parentDiffSubproblem ) nParentConss +=  parentDiffSubproblem->nLinearConss;
   int nConflicts = 0;
   if( conflictConsList ) nConflicts += conflictConsList->size();
   nLinearConss = nParentConss + localCutsList.size() + nConflicts;
   if( nLinearConss > 0 )
   {
      linearLhss = new SCIP_Real[nLinearConss];
      linearRhss = new SCIP_Real[nLinearConss];
      nLinearCoefs = new int[nLinearConss];
      linearCoefs = new SCIP_Real*[nLinearConss];
      idxLinearCoefsVars = new int*[nLinearConss];

      int i = 0;
      for(; i < nParentConss; i++ )
      {
         linearLhss[i] = parentDiffSubproblem->linearLhss[i];
         linearRhss[i] = parentDiffSubproblem->linearRhss[i];
         nLinearCoefs[i] = parentDiffSubproblem->nLinearCoefs[i];
         linearCoefs[i] = new SCIP_Real[nLinearCoefs[i]];
         idxLinearCoefsVars[i] = new int[nLinearCoefs[i]];
         for( int j = 0; j < nLinearCoefs[i]; j++ )
         {
            linearCoefs[i][j] = parentDiffSubproblem->linearCoefs[i][j];
            idxLinearCoefsVars[i][j] = parentDiffSubproblem->idxLinearCoefsVars[i][j];
         }
      }

      int nLocalCuts = localCutsList.size();
      for(; i < ( nParentConss + nLocalCuts ); i++ )
      {
         assert(!localCutsList.empty());
         LocalNodeInfo *cutInfo = localCutsList.front();
         localCutsList.pop_front();
         linearLhss[i] = cutInfo->linearLhs;
         linearRhss[i] = cutInfo->linearRhs;
         nLinearCoefs[i] = cutInfo->nLinearCoefs;
         linearCoefs[i] = cutInfo->linearCoefs;
         idxLinearCoefsVars[i] = cutInfo->idxLinearCoefsVars;
         delete cutInfo;
      }

      if( i < nLinearConss )
      {
         assert( conflictConsList );
         std::list<LocalNodeInfoPtr>::iterator pos;
         pos = conflictConsList->begin();
         for(; i < nLinearConss; i++ )
         {
            assert( pos != conflictConsList->end() );
            linearLhss[i] = (*pos)->linearLhs;
            linearRhss[i] = (*pos)->linearRhs;
            nLinearCoefs[i] = (*pos)->nLinearCoefs;
            linearCoefs[i] = new SCIP_Real[nLinearCoefs[i]];
            idxLinearCoefsVars[i] = new int[nLinearCoefs[i]];
            for( int j = 0; j < nLinearCoefs[i]; j++ )
            {
               linearCoefs[i][j] = (*pos)->linearCoefs[j];
               idxLinearCoefsVars[i][j] = (*pos)->idxLinearCoefsVars[j];
            }
            pos++;
         }
      }
   }
}

void
ScipParaDiffSubproblem::addBranchVarStats(
      SCIP *scip,
      ScipParaSolver *scipParaSolver
      )
{
   int nvars;                                /* number of variables                           */
   int nbinvars;                             /* number of binary variables                    */
   int nintvars;                             /* number of integer variables                   */
   SCIP_VAR** vars;                          /* transformed problem's variables               */
   SCIP_CALL_ABORT( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   int ngenvars = nbinvars+nintvars;
   nVarBranchStats = ngenvars;
   idxLBranchStatsVars = new int[ngenvars];
   downpscost = new SCIP_Real[ngenvars];
   uppscost = new SCIP_Real[ngenvars];
   downvsids = new SCIP_Real[ngenvars];
   upvsids = new SCIP_Real[ngenvars];
   downconflen = new SCIP_Real[ngenvars];
   upconflen = new SCIP_Real[ngenvars];
   downinfer = new SCIP_Real[ngenvars];
   upinfer = new SCIP_Real[ngenvars];
   downcutoff = new SCIP_Real[ngenvars];
   upcutoff = new SCIP_Real[ngenvars];
   for( int i = 0; i < ngenvars; ++i )
   {
      assert( SCIPvarGetType(vars[i]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_INTEGER );
      SCIP_VAR *transformVar = vars[i];
      SCIP_Real scalar = 1.0;
      SCIP_Real constant = 0.0;
      SCIP_CALL_ABORT( SCIPvarGetOrigvarSum(&transformVar, &scalar, &constant ) );
      assert(transformVar != NULL);
      idxLBranchStatsVars[i] = SCIPvarGetIndex(transformVar);

      SCIP_BRANCHDIR branchdir1;
      SCIP_BRANCHDIR branchdir2;
      if( scalar > 0.0 )
      {
         branchdir1 = SCIP_BRANCHDIR_DOWNWARDS;
         branchdir2 = SCIP_BRANCHDIR_UPWARDS;
      }
      else
      {
         branchdir1 = SCIP_BRANCHDIR_UPWARDS;
         branchdir2 = SCIP_BRANCHDIR_DOWNWARDS;
      }
      downpscost[i] = SCIPgetVarPseudocost(scip, vars[i], branchdir1);
      uppscost[i] = SCIPgetVarPseudocost(scip, vars[i], branchdir2);
      downvsids[i] = SCIPgetVarVSIDS(scip, vars[i], branchdir1);
      upvsids[i] = SCIPgetVarVSIDS(scip, vars[i], branchdir2);
      downconflen[i] = SCIPgetVarAvgConflictlength(scip, vars[i], branchdir1);
      upconflen[i] = SCIPgetVarAvgConflictlength(scip, vars[i], branchdir2);
      downinfer[i] = SCIPgetVarAvgInferences(scip, vars[i], branchdir1);
      upinfer[i] = SCIPgetVarAvgInferences(scip, vars[i], branchdir2);
      downcutoff[i] = SCIPgetVarAvgCutoffs(scip, vars[i], branchdir1);
      upcutoff[i] = SCIPgetVarAvgCutoffs(scip, vars[i], branchdir2);
   }
}

void
ScipParaDiffSubproblem::addInitialBranchVarStats(
      int  inMinDepth,
      int  inMaxDepth,
      SCIP *scip
      )
{
   if( nVarBranchStats != 0 ) return;       /* if this is not initial transfer to the other solvers, do nothing */

   offset = inMinDepth;

   int nvars;                                /* number of variables                           */
   int nbinvars;                             /* number of binary variables                    */
   int nintvars;                             /* number of integer variables                   */
   SCIP_VAR** vars;                          /* transformed problem's variables               */
   SCIP_CALL_ABORT( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   int ngenvars = nbinvars+nintvars;
   nVarBranchStats = ngenvars;
   idxLBranchStatsVars = new int[ngenvars];
   downpscost = new SCIP_Real[ngenvars];
   uppscost = new SCIP_Real[ngenvars];
   downvsids = new SCIP_Real[ngenvars];
   upvsids = new SCIP_Real[ngenvars];
   downconflen = new SCIP_Real[ngenvars];
   upconflen = new SCIP_Real[ngenvars];
   downinfer = new SCIP_Real[ngenvars];
   upinfer = new SCIP_Real[ngenvars];
   downcutoff = new SCIP_Real[ngenvars];
   upcutoff = new SCIP_Real[ngenvars];
   for( int i = 0; i < ngenvars; ++i )
   {
      assert( SCIPvarGetType(vars[i]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_INTEGER );
      idxLBranchStatsVars[i] = SCIPvarGetProbindex(vars[i]);
      downpscost[i] = SCIPgetVarPseudocost(scip, vars[i], SCIP_BRANCHDIR_DOWNWARDS);
      uppscost[i] = SCIPgetVarPseudocost(scip, vars[i], SCIP_BRANCHDIR_UPWARDS);
      downvsids[i] = SCIPgetVarVSIDS(scip, vars[i], SCIP_BRANCHDIR_DOWNWARDS);
      upvsids[i] = SCIPgetVarVSIDS(scip, vars[i], SCIP_BRANCHDIR_UPWARDS);
      downconflen[i] = SCIPgetVarAvgConflictlength(scip, vars[i], SCIP_BRANCHDIR_DOWNWARDS);
      upconflen[i] = SCIPgetVarAvgConflictlength(scip, vars[i], SCIP_BRANCHDIR_UPWARDS);
      downinfer[i] = SCIPgetVarAvgInferences(scip, vars[i], SCIP_BRANCHDIR_DOWNWARDS);
      upinfer[i] = SCIPgetVarAvgInferences(scip, vars[i], SCIP_BRANCHDIR_UPWARDS);
      downcutoff[i] = SCIPgetVarAvgCutoffs(scip, vars[i], SCIP_BRANCHDIR_DOWNWARDS);
      upcutoff[i] = SCIPgetVarAvgCutoffs(scip, vars[i], SCIP_BRANCHDIR_UPWARDS);
   }
}

void
ScipParaDiffSubproblem::write(
      ogzstream &out
      )
{
   out.write((char *)&localInfoIncluded, sizeof(int));
   out.write((char *)&nBoundChanges, sizeof(int));
   for(int i = 0; i < nBoundChanges; i++ )
   {
      out.write((char *)&indicesAmongSolvers[i], sizeof(int));
      out.write((char *)&branchBounds[i], sizeof(SCIP_Real));
      out.write((char *)&boundTypes[i], sizeof(int));
   }
   out.write((char *)&nLinearConss, sizeof(int));
   for(int i = 0; i < nLinearConss; i++ )
   {
      out.write((char *)&linearLhss[i], sizeof(SCIP_Real));
      out.write((char *)&linearRhss[i], sizeof(SCIP_Real));
      out.write((char *)&nLinearCoefs[i], sizeof(int));
      for(int j = 0; j < nLinearCoefs[i]; j++ )
      {
         out.write((char *)&linearCoefs[i][j], sizeof(SCIP_Real));
         out.write((char *)&idxLinearCoefsVars[i][j], sizeof(int));
      }
   }
   out.write((char *)&offset, sizeof(int));
   out.write((char *)&nVarBranchStats, sizeof(int));
   for(int i = 0; i < nVarBranchStats; i++ )
   {
      out.write((char *)&idxLBranchStatsVars[i], sizeof(int));
      out.write((char *)&downpscost[i], sizeof(SCIP_Real));
      out.write((char *)&uppscost[i], sizeof(SCIP_Real));
      out.write((char *)&downvsids[i], sizeof(SCIP_Real));
      out.write((char *)&upvsids[i], sizeof(SCIP_Real));
      out.write((char *)&downconflen[i], sizeof(SCIP_Real));
      out.write((char *)&upconflen[i], sizeof(SCIP_Real));
      out.write((char *)&downinfer[i], sizeof(SCIP_Real));
      out.write((char *)&upinfer[i], sizeof(SCIP_Real));
      out.write((char *)&downcutoff[i], sizeof(SCIP_Real));
      out.write((char *)&upcutoff[i], sizeof(SCIP_Real));
   }

}

void
ScipParaDiffSubproblem::read(
      ParaComm *comm,
      igzstream &in,
      bool onlyBoundChanges
      )
{
   in.read((char *)&localInfoIncluded, sizeof(int));
   in.read((char *)&nBoundChanges, sizeof(int));
   indicesAmongSolvers = new int[nBoundChanges];
   branchBounds = new SCIP_Real[nBoundChanges];
   boundTypes = new SCIP_BOUNDTYPE[nBoundChanges];
   for(int i = 0; i < nBoundChanges; i++ )
   {
      in.read((char *)&indicesAmongSolvers[i], sizeof(int));
      in.read((char *)&branchBounds[i], sizeof(SCIP_Real));
      in.read((char *)&boundTypes[i], sizeof(int));
   }
   if( !onlyBoundChanges )
   {
      in.read((char *)&nLinearConss, sizeof(int));
      if( nLinearConss > 0 )
      {
         linearLhss = new SCIP_Real[nLinearConss];
         linearRhss = new SCIP_Real[nLinearConss];
         nLinearCoefs = new int[nLinearConss];
         linearCoefs = new SCIP_Real*[nLinearConss];
         idxLinearCoefsVars = new int*[nLinearConss];
         for(int i = 0; i < nLinearConss; i++ )
         {
            in.read((char *)&linearLhss[i], sizeof(SCIP_Real));
            in.read((char *)&linearRhss[i], sizeof(SCIP_Real));
            in.read((char *)&nLinearCoefs[i], sizeof(int));
            assert(nLinearCoefs[i] > 0);
            linearCoefs[i] = new SCIP_Real[nLinearCoefs[i]];
            idxLinearCoefsVars[i] = new int[nLinearCoefs[i]];
            for(int j = 0; j < nLinearCoefs[i]; j++ )
            {
               in.read((char *)&linearCoefs[i][j], sizeof(SCIP_Real));
               in.read((char *)&idxLinearCoefsVars[i][j], sizeof(int));
            }
         }
      }
#if !( SCIP_VERSION == 211 && SCIP_SUBVERSION == 0 )
      in.read((char *)&offset, sizeof(int));
      in.read((char *)&nVarBranchStats, sizeof(int));
      if( nVarBranchStats > 0 )
      {
         idxLBranchStatsVars = new int[nVarBranchStats];
         downpscost = new SCIP_Real[nVarBranchStats];
         uppscost = new SCIP_Real[nVarBranchStats];
         downvsids = new SCIP_Real[nVarBranchStats];
         upvsids = new SCIP_Real[nVarBranchStats];
         downconflen = new SCIP_Real[nVarBranchStats];
         upconflen = new SCIP_Real[nVarBranchStats];
         downinfer = new SCIP_Real[nVarBranchStats];
         upinfer = new SCIP_Real[nVarBranchStats];
         downcutoff = new SCIP_Real[nVarBranchStats];
         upcutoff = new SCIP_Real[nVarBranchStats];
      }
      for(int i = 0; i < nVarBranchStats; i++ )
      {
         in.read((char *)&idxLBranchStatsVars[i], sizeof(int));
         in.read((char *)&downpscost[i], sizeof(SCIP_Real));
         in.read((char *)&uppscost[i], sizeof(SCIP_Real));
         in.read((char *)&downvsids[i], sizeof(SCIP_Real));
         in.read((char *)&upvsids[i], sizeof(SCIP_Real));
         in.read((char *)&downconflen[i], sizeof(SCIP_Real));
         in.read((char *)&upconflen[i], sizeof(SCIP_Real));
         in.read((char *)&downinfer[i], sizeof(SCIP_Real));
         in.read((char *)&upinfer[i], sizeof(SCIP_Real));
         in.read((char *)&downcutoff[i], sizeof(SCIP_Real));
         in.read((char *)&upcutoff[i], sizeof(SCIP_Real));
      }
#endif
   }
   else
   {
      assert( nLinearConss == 0);
      SCIP_Real dummyReal;
      int       dummyInt;
      int       nConss;
      in.read((char *)&nConss, sizeof(int));
      if( nConss > 0 )
      {
         for(int i = 0; i < nConss; i++ )
         {
            int nCoefs;
            in.read((char *)&dummyReal, sizeof(SCIP_Real));
            in.read((char *)&dummyReal, sizeof(SCIP_Real));
            in.read((char *)&nCoefs, sizeof(int));
            assert(nCoefs > 0);
            for(int j = 0; j < nCoefs; j++ )
            {
               in.read((char *)&dummyReal, sizeof(SCIP_Real));
               in.read((char *)&dummyInt, sizeof(int));
            }
         }
      }
#if !( SCIP_VERSION == 211 && SCIP_SUBVERSION == 0 )
      assert( nVarBranchStats == 0 );
      int nStats;
      in.read((char *)&dummyInt, sizeof(int));
      in.read((char *)&nStats, sizeof(int));

      for(int i = 0; i < nStats; i++ )
      {
         in.read((char *)&dummyInt, sizeof(int));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
         in.read((char *)&dummyReal, sizeof(SCIP_Real));
      }
#endif
   }
}

/** get fixed variables **/
int
ScipParaDiffSubproblem::getFixedVariables(
      ParaInstance *instance,
      ParaFixedVariable **fixedVars
      )
{
   ScipParaInstance *scipInstance = dynamic_cast< ScipParaInstance* >(instance);
   *fixedVars = new ParaFixedVariable[nBoundChanges];   //  allocate space for the maximum possible numbers
   int n = 0;
   int k = 0;
   for( int i = 0; i < nBoundChanges; i++ )
   {
      switch( static_cast<SCIP_VARTYPE>(scipInstance->getVarType( indicesAmongSolvers[i] ) ) )
      {
      case SCIP_VARTYPE_BINARY:
         if( boundTypes[i] == SCIP_BOUNDTYPE_LOWER )
         {
            assert( EPSEQ(branchBounds[i], scipInstance->getVarUb(indicesAmongSolvers[i]), DEFAULT_NUM_EPSILON ) );
         }
         else
         {
            assert( EPSEQ(branchBounds[i], scipInstance->getVarLb(indicesAmongSolvers[i]), DEFAULT_NUM_EPSILON ) );
         }
         for( k = 0 ; k < n; k++ )
         {
            if( indicesAmongSolvers[i] == (*fixedVars)[k].index )   // when I checked, this case happened!
               break;
         }
         if( k != n ) break;
         (*fixedVars)[n].nSameValue = 0;
         (*fixedVars)[n].index = indicesAmongSolvers[i];
         (*fixedVars)[n].value = branchBounds[i];
         (*fixedVars)[n].mnode = 0;
         (*fixedVars)[n].next = 0;
         (*fixedVars)[n].prev = 0;
         n++;
         break;
      case SCIP_VARTYPE_INTEGER:
      case SCIP_VARTYPE_IMPLINT:
      case SCIP_VARTYPE_CONTINUOUS:
         if( boundTypes[i] == SCIP_BOUNDTYPE_LOWER )
         {
            if( !EPSEQ(branchBounds[i], scipInstance->getVarUb(indicesAmongSolvers[i]), DEFAULT_NUM_EPSILON ) )
            {
               int j = i + 1;
               for( ; j < nBoundChanges; j++ )
               {
                  if( indicesAmongSolvers[i] == indicesAmongSolvers[j] &&
                        boundTypes[j] == SCIP_BOUNDTYPE_UPPER &&
                        EPSEQ( branchBounds[i], branchBounds[j], DEFAULT_NUM_EPSILON )
                        )
                  {
                     break;
                  }
               }
               if( j >= nBoundChanges ) break;
            }
         }
         else
         {
            if( !EPSEQ(branchBounds[i], scipInstance->getVarLb(indicesAmongSolvers[i]), DEFAULT_NUM_EPSILON ) )
            {
               int j = i + 1;
               for( ; j < nBoundChanges; j++ )
               {
                  if( indicesAmongSolvers[i] == indicesAmongSolvers[j] &&
                        boundTypes[j] == SCIP_BOUNDTYPE_LOWER &&
                        EPSEQ( branchBounds[i], branchBounds[j], DEFAULT_NUM_EPSILON )
                        )
                  {
                     break;
                  }
               }
               if( j >= nBoundChanges ) break;
            }
         }
         for(k = 0; k < n; k++ )
         {
            if( indicesAmongSolvers[i] == (*fixedVars)[k].index )   // when I checked, this case happened!
               break;
         }
         if( k != n ) break;
         (*fixedVars)[n].nSameValue = 0;
         (*fixedVars)[n].index = indicesAmongSolvers[i];
         (*fixedVars)[n].value = branchBounds[i];
         (*fixedVars)[n].mnode = 0;
         (*fixedVars)[n].next = 0;
         (*fixedVars)[n].prev = 0;
         n++;

         break;
      default:
         THROW_LOGICAL_ERROR2("Invalid Variable Type = ", static_cast<int>(scipInstance->getVarType( indicesAmongSolvers[i] ) ) );
      }
   }
   if( n == 0 )
   {
      delete [] *fixedVars;
      *fixedVars = 0;
   }
   return n;
}

/** create new ParaDiffSubproblem using fixed variables information */
ParaDiffSubproblem*
ScipParaDiffSubproblem::createDiffSubproblem(
      ParaComm          *comm,
      int                n,
      ParaFixedVariable *fixedVars
      )
{
   ScipParaDiffSubproblem *diffSubproblem = dynamic_cast<ScipParaDiffSubproblem*>( comm->createParaDiffSubproblem() );
   int nNewBranches = 0;
   diffSubproblem->indicesAmongSolvers = new int[nBoundChanges];
   diffSubproblem->branchBounds = new SCIP_Real[nBoundChanges];
   diffSubproblem->boundTypes = new SCIP_BOUNDTYPE[nBoundChanges];
   for( int i = 0; i < nBoundChanges; i++ )
   {
      for( int j = 0; j < n; j++ )
      {
         if( indicesAmongSolvers[i] == fixedVars[j].index )
         {
            diffSubproblem->indicesAmongSolvers[nNewBranches] = indicesAmongSolvers[i];
            diffSubproblem->branchBounds[nNewBranches] = branchBounds[i];
            diffSubproblem->boundTypes[nNewBranches] = boundTypes[i];
            nNewBranches++;
            break;
         }
      }
   }
   diffSubproblem->nBoundChanges = nNewBranches;
   return diffSubproblem;
}

/** stringfy ParaCalculationState */
const std::string
ScipParaDiffSubproblem::toString(
      )
{
   std::ostringstream s;

   s << "localInfoIncluded = " << localInfoIncluded << std::endl;
   s << "nBranches = " <<  nBoundChanges << std::endl;
   for(int i = 0; i < nBoundChanges; i++ )
   {
      s << "indicesAmongSolvers[" << i << "] = " << indicesAmongSolvers[i] << std::endl;
      s << "branchBounds[" << i << "] = " << branchBounds[i] << std::endl;
      s << "boudTypes[" << i << "] = " << static_cast<int>(boundTypes[i]) << std::endl;
   }
   s << "nLinearConss = " << nLinearConss << std::endl;
   if( nLinearConss > 0 )
   {
      for(int i = 0; i < nLinearConss; i++ )
      {
         s << "linearLhss[" << i << "] = " << linearLhss[i] << std::endl;
         s << "linearRhss[" << i << "] = " << linearRhss[i] << std::endl;
         s << "nLinearCoefs[" << i << "] = " << nLinearCoefs[i] << std::endl;
         for( int j = 0; j < nLinearCoefs[i]; j++ )
         {
            s << "linearCoefs[" << i << "][" << j << "] = " << linearCoefs[i][j] << std::endl;
            s << "idxLinearCoefsVars[" << i << "][" << j << "] = " << idxLinearCoefsVars[i][j] << std::endl;
         }
      }
   }
   s << "offset = " << offset << std::endl;
   if( nVarBranchStats > 0 )
   {
      for(int i = 0; i < nVarBranchStats; i++ )
      {
         s << "idxLBranchStatsVars[" << i << "] = " << idxLBranchStatsVars[i] << std::endl;
         s << "downpscost[" << i << "] = " << downpscost[i] << std::endl;
         s << "uppscost[" << i << "] = " << uppscost[i] << std::endl;
         s << "downvsids[" << i << "] = " << downvsids[i] << std::endl;
         s << "upvsids[" << i << "] = " << upvsids[i] << std::endl;
         s << "downconflen[" << i << "] = " << downconflen[i] << std::endl;
         s << "upconflen[" << i << "] = " << upconflen[i] << std::endl;
         s << "downinfer[" << i << "] = " << downinfer[i] << std::endl;
         s << "upinfer[" << i << "] = " << upinfer[i] << std::endl;
         s << "downcutoff[" << i << "] = " << downcutoff[i] << std::endl;
         s << "upcutoff[" << i << "] = " << upcutoff[i] << std::endl;
      }

   }
   return s.str();
}


