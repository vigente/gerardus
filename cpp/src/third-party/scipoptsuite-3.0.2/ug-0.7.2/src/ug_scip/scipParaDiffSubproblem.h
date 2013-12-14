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

/**@file    scipParaDiffSubproblem.h
 * @brief   ParaInitialStat extension for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_DIFF_SUBPROBLEM_H__
#define __SCIP_PARA_DIFF_SUBPROBLEM_H__

#include <iostream>
#include <fstream>
#include "ug/paraComm.h"
#include "ug/gzstream.h"
#include "ug/paraInstance.h"
#include "ug/paraDiffSubproblem.h"
#include "scip/scip.h"

namespace ParaSCIP
{

class ScipParaSolver;

/** The difference between instance and subproblem: this is base class */
class ScipParaDiffSubproblem : public UG::ParaDiffSubproblem
{
protected:
   int localInfoIncluded;                /**< 0 (0000 0000): not included
                                              1 (0000 0001): if local cuts are included
                                              2 (0000 0010): if conflicts are included
                                              3 (0000 0011): if local cuts and conflicts are included */
   /******************************
    * for variable bound changes *
    * ***************************/
   int             nBoundChanges;            /**< number of branching variables */
   int             *indicesAmongSolvers; /**< array of variable indices ( unique index )  */
   SCIP_Real       *branchBounds;        /**< array of bounds which the branchings     */
   SCIP_BOUNDTYPE  *boundTypes;          /**< array of boundtypes which the branchings */
   /********************************
    * for local cuts and conflicts *
    * *****************************/
   int            nLinearConss;          /**< number of linear constrains */
   SCIP_Real      *linearLhss;           /**< array of lhs */
   SCIP_Real      *linearRhss;           /**< array of rhs */
   int            *nLinearCoefs;         /**< array of number of coefficient values for linear constrains */
   SCIP_Real      **linearCoefs;         /**< array of non-zero coefficient values of linear constrains */
   int            **idxLinearCoefsVars;  /**< array of indices of no-zero coefficient values of linear constrains */
   /********************************
    * for local var brnach stats   *
    * *****************************/
   int           offset;                /**< root node depth of the collected nodes */
   int           nVarBranchStats;       /**< number of branch stats */
   int           *idxLBranchStatsVars;  /**< indecies of indices of branch stats vars */
   SCIP_Real     *downpscost;           /**< values to which pseudocosts for downwards branching */
   SCIP_Real     *uppscost;             /**< values to which pseudocosts for upwards branching */
   SCIP_Real     *downvsids;            /**< values to which VSIDS score for downwards branching */
   SCIP_Real     *upvsids;              /**< values to which VSIDS score for upwards branching */
   SCIP_Real     *downconflen;          /**< values to which conflict length score for downwards branching */
   SCIP_Real     *upconflen;            /**< values to which conflict length score for upwards branching */
   SCIP_Real     *downinfer;            /**< values to which inference counter for downwards branching */
   SCIP_Real     *upinfer;              /**< values to which inference counter for upwards branching */
   SCIP_Real     *downcutoff;           /**< values to which cutoff counter for downwards branching */
   SCIP_Real     *upcutoff;             /**< values to which cutoff counter for upwards branching */
public:
   /** default constructor */
   ScipParaDiffSubproblem(
         )
         : localInfoIncluded(0),
           nBoundChanges(0), indicesAmongSolvers(0), branchBounds(0), boundTypes(0),
           nLinearConss(0), linearLhss(0), linearRhss(0),
           nLinearCoefs(0), linearCoefs(0), idxLinearCoefsVars(0),
           offset(0), nVarBranchStats(0), idxLBranchStatsVars(0), downpscost(0), uppscost(0),
           downvsids(0), upvsids(0), downconflen(0), upconflen(0), downinfer(0), upinfer(0),
           downcutoff(0), upcutoff(0)
   {
   }

   ScipParaDiffSubproblem(
         SCIP *scip,
         ScipParaSolver *scipParaSolver,
         int nNewBranchVars,
         SCIP_VAR **newBranchVars,
         SCIP_Real *newBranchBounds,
         SCIP_BOUNDTYPE *newBoundTypes
         );

   ScipParaDiffSubproblem(
         ScipParaDiffSubproblem *diffSubproblem
         ) : localInfoIncluded(0),
         nBoundChanges(0), indicesAmongSolvers(0), branchBounds(0), boundTypes(0),
         nLinearConss(0), linearLhss(0), linearRhss(0),
         nLinearCoefs(0), linearCoefs(0), idxLinearCoefsVars(0),
         offset(0), nVarBranchStats(0), idxLBranchStatsVars(0), downpscost(0), uppscost(0),
         downvsids(0), upvsids(0), downconflen(0), upconflen(0), downinfer(0), upinfer(0),
         downcutoff(0), upcutoff(0)
   {
      if( !diffSubproblem ) return;

      localInfoIncluded = diffSubproblem->localInfoIncluded;
      nBoundChanges = diffSubproblem->nBoundChanges;
      if( nBoundChanges )
      {
         indicesAmongSolvers = new int[nBoundChanges];
         branchBounds = new SCIP_Real[nBoundChanges];
         boundTypes = new SCIP_BOUNDTYPE[nBoundChanges];
         for( int i = 0; i < nBoundChanges; i++ )
         {
            indicesAmongSolvers[i] = diffSubproblem->indicesAmongSolvers[i];
            branchBounds[i] = diffSubproblem->branchBounds[i];
            boundTypes[i] = diffSubproblem->boundTypes[i];
         }
      }
      nLinearConss = diffSubproblem->nLinearConss;
      if( nLinearConss > 0 )
      {
         linearLhss = new SCIP_Real[nLinearConss];
         linearRhss = new SCIP_Real[nLinearConss];
         nLinearCoefs = new int[nLinearConss];
         linearCoefs = new SCIP_Real*[nLinearConss];
         idxLinearCoefsVars = new int*[nLinearConss];
         for( int c = 0; c < nLinearConss; c++ )
         {
            linearLhss[c] = diffSubproblem->linearLhss[c];
            linearRhss[c] = diffSubproblem->linearRhss[c];
            nLinearCoefs[c] = diffSubproblem->nLinearCoefs[c];
            linearCoefs[c] = new SCIP_Real[nLinearCoefs[c]];
            idxLinearCoefsVars[c] = new int[nLinearCoefs[c]];
            for( int v = 0; v < nLinearCoefs[c]; v++ )
            {
               linearCoefs[c][v] = diffSubproblem->linearCoefs[c][v];
               idxLinearCoefsVars[c][v] = diffSubproblem->idxLinearCoefsVars[c][v];
            }
         }
      }

      offset = diffSubproblem->offset;
      nVarBranchStats = diffSubproblem->nVarBranchStats;
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
         for( int i = 0; i < nVarBranchStats; ++i )
         {
            idxLBranchStatsVars[i] = diffSubproblem->idxLBranchStatsVars[i];
            downpscost[i] = diffSubproblem->downpscost[i];
            uppscost[i] = diffSubproblem->uppscost[i];
            downvsids[i] = diffSubproblem->downvsids[i];
            upvsids[i] = diffSubproblem->upvsids[i];
            downconflen[i] = diffSubproblem->downconflen[i];
            upconflen[i] = diffSubproblem->upconflen[i];
            downinfer[i] = diffSubproblem->downinfer[i];
            upinfer[i] = diffSubproblem->upinfer[i];
            downcutoff[i] = diffSubproblem->downcutoff[i];
            upcutoff[i] = diffSubproblem->upcutoff[i];
         }
      }
   }

   /** destractor */
   virtual ~ScipParaDiffSubproblem()
   {
      if( indicesAmongSolvers ) delete[] indicesAmongSolvers;
      if( branchBounds ) delete[] branchBounds;
      if( boundTypes ) delete[] boundTypes;
      if( nLinearConss > 0 )
      {
         if( linearLhss ) delete[] linearLhss;
         if( linearRhss ) delete[] linearRhss;
         for( int c = 0; c < nLinearConss; c++ )
         {
            if( nLinearCoefs[c] > 0 )
            {
               if( linearCoefs[c] ) delete[] linearCoefs[c];
               if( idxLinearCoefsVars[c] ) delete[] idxLinearCoefsVars[c];
            }
         }
         if( linearCoefs ) delete[] linearCoefs;
         if( idxLinearCoefsVars ) delete[] idxLinearCoefsVars;
      }
      if( nVarBranchStats > 0 )
      {
         if( idxLBranchStatsVars ) delete[] idxLBranchStatsVars;
         if( downpscost ) delete[] downpscost;
         if( uppscost ) delete[] uppscost;
         if( downvsids ) delete[] downvsids;
         if( upvsids ) delete[] upvsids;
         if( downconflen ) delete[] downconflen;
         if( upconflen ) delete[] upconflen;
         if( downinfer ) delete[] downinfer;
         if( upinfer ) delete[] upinfer;
         if( downcutoff ) delete[] downcutoff;
         if( upcutoff ) delete[] upcutoff;
      }
   }

   int getNBoundChanges(){ return nBoundChanges; }
   int getIndex(int i){ return indicesAmongSolvers[i]; }
   SCIP_Real getBranchBound(int i){ return branchBounds[i]; }
   SCIP_BOUNDTYPE getBoundType(int i){ return boundTypes[i]; }

   int getNLinearConss(){ return nLinearConss; }
   SCIP_Real getLinearLhs(int i){ return linearLhss[i]; }
   SCIP_Real getLinearRhs(int i){ return linearRhss[i]; }
   int getNLinearCoefs(int i){ return nLinearCoefs[i]; }
   SCIP_Real getLinearCoefs(int i, int j){ return linearCoefs[i][j]; }
   int getIdxLinearCoefsVars(int i, int j){ return idxLinearCoefsVars[i][j]; }

   int getNVarBranchStats(){ return nVarBranchStats; };
   int getIdxLBranchStatsVars(int i){ return idxLBranchStatsVars[i]; }
   SCIP_Real getDownpscost(int i){ return downpscost[i]; }
   SCIP_Real getUppscost(int i){ return uppscost[i]; }
   SCIP_Real getDownvsids(int i){ return downvsids[i]; }
   SCIP_Real getUpvsids(int i){ return upvsids[i]; }
   SCIP_Real getDownconflen(int i){ return downconflen[i]; }
   SCIP_Real getUpconflen(int i){ return upconflen[i]; }
   SCIP_Real getDowninfer(int i){ return downinfer[i]; }
   SCIP_Real getUpinfer(int i){ return upinfer[i]; }
   SCIP_Real getDowncutoff(int i){ return downcutoff[i]; }
   SCIP_Real getUpcutoff(int i){ return upcutoff[i]; }

   void addLocalNodeInfo(
         SCIP *scip,
         ScipParaSolver *scipParaSolver
         );
   void addBranchVarStats(
         SCIP *scip,
         ScipParaSolver *scipParaSolver
         );
   void addInitialBranchVarStats(
         int minDepth,
         int maxDepth,
         SCIP *scip
         );
   int getOffset(){ return offset; }

   void write(ogzstream &out);
   void read(UG::ParaComm *comm, igzstream &in, bool onlyBoundChanges);

   /** get fixed variables **/
   int getFixedVariables(UG::ParaInstance *instance, UG::ParaFixedVariable **fixedVars );

   /** create new ParaDiffSubproblem using fixed variables information */
   ParaDiffSubproblem* createDiffSubproblem(
		   UG::ParaComm *comm,
		   int n,
		   UG::ParaFixedVariable *fixedVars );

   /** stringfy ParaCalculationState */
   const std::string toString(
         );
};

}

#endif    // __SCIP_PARA_DIFF_SUBPROBLEM_H__

