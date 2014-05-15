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

/**@file    paraSolverPool.h
 * @brief   Solver pool.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_SOLVER_POOL_H__
#define __PARA_SOLVER_POOL_H__
#include <map>
#include "paraNode.h"
#include "paraTimer.h"
#include "paraNodePool.h"
#include "paraRacingRampUpParamSet.h"
#include "paraSolverTerminationState.h"
#include "paraDeterministicTimer.h"

namespace UG
{

class ParaSolverPoolElement;
typedef ParaSolverPoolElement * ParaSolverPoolElementPtr;

enum SolverStatus {Inactive, Racing, RacingEvaluation, Active, Dead};

#define SOLVER_POOL_INDEX( rank )   ( rank - originRank )

/** This class includes information about a Solver status */
class ParaSolverPoolElement
{
   int                   rank;                       /**< rank of the solver */
   SolverStatus          status;                     /**< status of the solver */
   bool                  collectingMode;             /**< indicate if current solver is in collecting mode or not */
   bool                  candidateOfCollecting;      /**< indicate that this solver is a candidate of collecting mode solver */
   bool                  generator;                  /**< this solver can generate subproblems */
   ParaNode              *currentNode;               /**< solving node */
   ParaSolverPoolElementPtr *selectionHeapElement;   /**< pointer to selection heap element */
   ParaSolverPoolElementPtr *collectingModeSolverHeapElement; /**< pointer to collecting mode heap element */
   /** the following values are used to make decision for load balancing */
   long long            numOfNodesSolved;            /**< number of nodes solved. -1 is the special value which means never updated in racing */
   int                  numOfDiffNodesSolved;        /**< number of nodes solved difference between current number and
                                                       * that in the previous notification time */
   int                  numOfNodesLeft;              /**< number of nodes left */
   int                  numOfDiffNodesLeft;          /**< number of nodes left difference between current number
                                                       * and that in the previous notification time */
   double               bestDualBoundValue;          /**< best dual bound value of the Solver */
   ParaRacingRampUpParamSet   *racingRampUpParamSet; /**< parameter set for racing ramp-up */
   ParaSolverTerminationState *termState;            /**< solver termination statistics: for checkpoint */
public :
   /** constructor */
   ParaSolverPoolElement(
         int inRank
         )
         : rank(inRank), status(Inactive), collectingMode(false), candidateOfCollecting(false), generator(true), currentNode(0),
         selectionHeapElement(0), collectingModeSolverHeapElement(0), numOfNodesSolved(0), numOfNodesLeft(0), numOfDiffNodesLeft(0), bestDualBoundValue(0.0),
         racingRampUpParamSet(0), termState(0)
   {
   }

   /** destractor */
   ~ParaSolverPoolElement(
         )
   {
      if( currentNode ) delete currentNode;
      if( termState ) delete termState;
   }

   /** get rank of the solver */
   int getRank(
         )
   {
      return rank;
   }

   /** get current solving node */
   ParaNode *getCurrentNode(
         )
   {
      return currentNode;
   }

   /** extract current solving node */
   ParaNode *extractCurrentNode(
         )
   {
      ParaNode *node = currentNode;
      currentNode = 0;
      return node;
   }

   /** get selection heap element */
   ParaSolverPoolElementPtr *getSelectionHeapElement(
         )
   {
      return selectionHeapElement;
   }

   /** set selection heap element */
   void setSelectionHeapElement(
         ParaSolverPoolElementPtr *inSelectionHeapElement
         )
   {
       selectionHeapElement = inSelectionHeapElement;
   }

   /** get collecting mode solver heap element */
   ParaSolverPoolElementPtr *getCollectingModeSolverHeapElement(
         )
   {
      return collectingModeSolverHeapElement;
   }

   /** set collecting mode solver heap lelemtn */
   void setCollectingModeSolverHeapElement(
         ParaSolverPoolElementPtr *inCollectingModeSolverHeapElement
         )
   {
      collectingModeSolverHeapElement = inCollectingModeSolverHeapElement;
   }

   /** get number of nodes solved */
   long long getNumOfNodesSolved(
         )
   {
      return numOfNodesSolved;
   }

   /** get number of nodes left */
   int getNumOfNodesLeft(
         )
   {
      return numOfNodesLeft;
   }

   /** set number of nodes solved */
   void setNumOfNodesSolved(
         long long inNumOfNodesSolved
         )
   {
      numOfNodesSolved = inNumOfNodesSolved;
   }

   /** get number of nodes left difference between current number and that in the previous notification time */
   int getNumOfDiffNodesSolved(){ return numOfDiffNodesSolved; }

   /** set number of nodes left difference between current number and that in the previous notification time */
   void setNumOfDiffNodesSolved(
         int inNumOfDiff
         )
   {
       numOfDiffNodesSolved = inNumOfDiff;
   }

   /** set number of nodes left */
   void setNumOfNodesLeft(
         int inNumOfNodesLeft
         )
   {
      numOfNodesLeft = inNumOfNodesLeft;
   }

   /** set dual bound value on paraNode */
   void setDualBoundValue(
         double dualBoundValue
         )
   {
      currentNode->setDualBoundValue(dualBoundValue);
   }

   /** get number of nodes left difference between current number and that in the previous notification time */
   int getNumOfDiffNodesLeft(){ return numOfDiffNodesLeft; }

   /** set number of nodes left difference between current number and that in the previous notification time */
   void setNumOfDiffNodesLeft(
         int inNumOfDiff
         )
   {
       numOfDiffNodesLeft = inNumOfDiff;
   }

   /** get best dual bound value */
   double getBestDualBoundValue(
         )
   {
      return bestDualBoundValue;
   }

   /** set best dual bound value */
   void setBestDualBoundValue(
         double inBestDualBoundValue
         )
   {
      bestDualBoundValue = inBestDualBoundValue;
   }

   /** activate this solver */
   void activate(
         ParaNode *inNode
         )
   {
      status = Active;
      currentNode = inNode;
      if( currentNode )
         bestDualBoundValue = inNode->getDualBoundValue();
      else
         bestDualBoundValue = -DBL_MAX;
      numOfNodesSolved = 0;
      numOfNodesLeft = 1;  // at least 1 node should be left
   }

   /** activate for racing */
   void racingActivate(
         )
   {
      status = Racing;
      currentNode = 0;
      numOfNodesSolved = -1;
      bestDualBoundValue = -DBL_MAX;
      numOfNodesLeft = 1;  // at least 1 node should be left
   }

   /** inactivate this solver */
   void inactivate(
         )
   {
      status = Inactive;
      if( currentNode )
      {
         delete currentNode;
         currentNode = 0;
      }
      collectingMode = false;
      candidateOfCollecting = false;
      /** do not touch "generator" **/
      numOfNodesSolved = 0;
      numOfNodesLeft = 0;
      numOfDiffNodesLeft = 0;
      bestDualBoundValue = 0.0;
   }

   /** this solver died */
   ParaNode *died(
         )
   {
      status = Dead;
      collectingMode = false;
      candidateOfCollecting = false;
      numOfNodesSolved = 0;
      numOfNodesLeft = 0;
      numOfDiffNodesLeft = 0;
      bestDualBoundValue = 0.0;
      return currentNode;
   }

   /** get Solver status */
   SolverStatus getStatus(
         )
   {
      return status;
   }

   /** check if this solver is active or not */
   bool isActive(
         )
   {
      return ( status == Active );
   }

   /** check if this solver is out of collecting mode or not */
   bool isOutCollectingMode(
         )
   {
      return (!collectingMode);
   }

   /** check if this solver is in collecting mode or not */
   bool isInCollectingMode(
         )
   {
      return collectingMode;
   }

   /** set collecting mode */
   void setCollectingMode(
         bool b
        )
   {
      collectingMode = b;
   }

   /** check if this solver is candidate of collecting mode solver */
   bool isCandidateOfCollecting(
         )
   {
      return candidateOfCollecting;
   }

   /** set candidate of collecting mode solver */
   void setCandidateOfCollecting(
         bool b
         )
   {
      candidateOfCollecting = b;
   }

   /** set termination state */
   void setTermState(
         ParaSolverTerminationState *inTermState
         )
   {
      if( termState ) delete termState;
      termState = inTermState;
   }

   /** get termination state */
   ParaSolverTerminationState *getTermState(
         )
   {
      return termState;
   }

   /** switch into evaluation stage */
   void switchIntoEvaluation(
         )
   {
      assert( status == Racing );
      status = RacingEvaluation;
   }

   /** switch into evaluation stage */
   void switchOutEvaluation(
         )
   {
      assert( status == RacingEvaluation );
      status = Racing;
   }

   /** check if the solver is in racing stage */
   bool isRacingStage(
         )
   {
      return ( status == Racing );
   }

   /** check if the solver is in evaluation stage */
   bool isEvaluationStage(
         )
   {
      return ( status == RacingEvaluation );
   }

   /** make this solver No generator */
   void setNoGenerator(
         )
   {
      generator = false;
   }

   /** check if this solver is generator or not */
   bool isGenerator(
         )
   {
      return generator;
   }

};

/** Selection Heap class */
class SelectionHeap
{
public:
   enum ResultOfInsert
   {
      SUCCEEDED,
      FAILED_BY_FULL
   };
   /** constructor */
   SelectionHeap(
         int size
         );

   /** destructor */
   virtual ~SelectionHeap();

   /** insert ParaSolverPoolElementPtr to Selection Heap */
   ResultOfInsert insert(
         ParaSolverPoolElementPtr solver
         );



   /** obtain top priority ParaSolverPoolElementPtr */
   ParaSolverPoolElementPtr top(
         ) const
   {
      return heap[1];
   }

   /** remove top priority ParaSolverPoolElementPtr from Selection Heap */
   ParaSolverPoolElementPtr remove(
         );

   /** resize Selection Heap */
   void resize(
         int size
         );

   /** get current used heap size */
   inline int getHeapSize(
         ) const
   {
      return heapSize;
   }

   /** get max heap size */
   inline int getMaxHeapSize(
         ) const
   {
      return maxHeapSize;
   }

   /** update selection heap by a new dual bound value of the solver */
   virtual void updateDualBoundValue( ParaSolverPoolElementPtr, double ) = 0;

   /** delete ParaSolverPoolElementPtr from Selection Heap */
   virtual void deleteElement( ParaSolverPoolElementPtr solver ) = 0;

   /** up heap */
   virtual void upHeap(int pos) = 0;

   /** down heap */
   virtual void downHeap(int pos) = 0;

   //------------
   // for debug
   //------------
   const std::string toString();

protected:
   int                  maxHeapSize;             /**< maximum size of this heap */
   int                  heapSize;                /**< current used heap size */
   // long long            currentSequenceNumber;   /**< current sequence number */
   ParaSolverPoolElementPtr *heap;                   /**< heap : contents are ParaSolverPoolElementPtr */
};

class DescendingSelectionHeap : public SelectionHeap {
public:
   DescendingSelectionHeap(int size);
   virtual ~DescendingSelectionHeap(){}
   void updateDualBoundValue(ParaSolverPoolElementPtr, double);
   void deleteElement(ParaSolverPoolElementPtr solver);
   void upHeap(int pos);
   void downHeap(int pos);
};

class AscendingSelectionHeap : public SelectionHeap {
public:
   AscendingSelectionHeap(int size);
   virtual ~AscendingSelectionHeap(){}
   void updateDualBoundValue(ParaSolverPoolElementPtr, double);
   void deleteElement(ParaSolverPoolElementPtr solver);
   void upHeap(int pos);
   void downHeap(int pos);
};

/** Collecting Mode Solver Heap class */
class CollectingModeSolverHeap
{
public:
   enum ResultOfInsert
   {
      SUCCEEDED,
      FAILED_BY_FULL
   };
   /** constructor */
   CollectingModeSolverHeap(
         int size
         );

   /** destructor */
   virtual ~CollectingModeSolverHeap();

   /** insert ParaSolverPoolElementPtr to CollectingModeSolver Heap */
   ResultOfInsert insert(
         ParaSolverPoolElementPtr solver
         );



   /** obtain top priority ParaSolverPoolElementPtr */
   ParaSolverPoolElementPtr top(
         ) const
   {
      return heap[1];
   }

   /** remove top priority ParaSolverPoolElementPtr from CollectingModeSolver Heap */
   ParaSolverPoolElementPtr remove(
         );

   /** resize CollectingModeSolver Heap */
   void resize(
         int size
         );

   /** get current used heap size */
   inline int getHeapSize(
         ) const
   {
      return heapSize;
   }

   /** get max heap size */
   inline int getMaxHeapSize(
         ) const
   {
      return maxHeapSize;
   }

   /** update selection heap by a new dual bound value of the solver */
   virtual void updateDualBoundValue( ParaSolverPoolElementPtr, double ) = 0;

   /** delete ParaSolverPoolElementPtr from CollectingModeSolver Heap */
   virtual void deleteElement( ParaSolverPoolElementPtr solver ) = 0;

   /** up heap */
   virtual void upHeap(int pos) = 0;

   /** down heap */
   virtual void downHeap(int pos) = 0;

   //------------
   // for debug
   //------------
   const std::string toString();

protected:
   int                  maxHeapSize;             /**< maximum size of this heap */
   int                  heapSize;                /**< current used heap size */
   // long long            currentSequenceNumber;   /**< current sequence number */
   ParaSolverPoolElementPtr *heap;                   /**< heap : contents are ParaSolverPoolElementPtr */
};

class DescendingCollectingModeSolverHeap : public CollectingModeSolverHeap {
public:
   DescendingCollectingModeSolverHeap(int size);
   virtual ~DescendingCollectingModeSolverHeap(){}
   void updateDualBoundValue(ParaSolverPoolElementPtr, double);
   void deleteElement(ParaSolverPoolElementPtr solver);
   void upHeap(int pos);
   void downHeap(int pos);
};

class AscendingCollectingModeSolverHeap : public CollectingModeSolverHeap {
public:
   AscendingCollectingModeSolverHeap(int size);
   virtual ~AscendingCollectingModeSolverHeap(){}
   void updateDualBoundValue(ParaSolverPoolElementPtr, double);
   void deleteElement(ParaSolverPoolElementPtr solver);
   void upHeap(int pos);
   void downHeap(int pos);
};

class ParaRacingSolverPool;

/** Solver Pool base class */
class ParaSolverPool {
protected:
   bool                                 active;            /**< indicate if this pool is active or not */
   double                               bgap;              /**< threshold value of gap */
   double                               mp;                /**< multiplier of the threshold value p */
   double                               mBgap;             /**< multiplier of the bgap value */
   double                               absoluteGap;       /**< allowable absolute dual bound gap to the best solver */
   int                                  originRank;        /**< origin rank of solvers managed by this solver pool */
   int                                  nSolvers;          /**< number of solvers */
   int                                  nGenerator;        /**< number of generators */
   int                                  nCollectingModeSolvers; /**< number of collecting mode solvers */
   int                                  nMaxCollectingModeSolvers; /**< maximum number of solvers that can be in collecting mode */
   int                                  nLimitCollectingModeSolvers; /**< limit number of solvers that can be in collecting mode */
   long long                            nNodesSolvedInSolvers; /**< number of nodes solved in current running solvers */
   long long                            nTotalNodesSolved; /**< number of nodes solved : updated at termination of subtree comp. */
   long long                            nNodesInSolvers;   /**< number of nodes in all solvers */
   // bool                                 rampUpPhase;       /**< indicate that this system is ramp-up phase or not */
   bool                                 collectingMode;    /**< indicate that this system is in collecting mode or not */
   bool                                 breakingFirstSubtree; /**< breaking the first subtree */
   std::map< int, ParaSolverPoolElementPtr > inactiveSolvers;   /**< pointers to inactive solvers */
   std::map< int, ParaSolverPoolElementPtr > activeSolvers;     /**< pointers to active solvers */
   std::map< int, ParaSolverPoolElementPtr > deadSolvers;       /**< pointers to dead solvers */
   std::multimap< double, ParaSolverPoolElementPtr > candidatesOfCollectingModeSolvers;
                                                                /**< pointers to candidates of collecting mode solvers */
   ParaSolverPoolElementPtr             *pool;             /**< solver pool indexed by solver's rank */
   SelectionHeap                        *selectionHeap;    /**< pointers to active solvers in ascending or descending order */
   CollectingModeSolverHeap             *collectingModeSolverHeap; /**< pointers to collecting mode solvers in ascending or descending order */
   ParaComm                             *paraComm;         /**< communicator */
   ParaParamSet                         *paraParams;       /**< runtime parameters for parallelization */

protected:
   void switchInCollectingToSolver(
         int rank,
         ParaNodePool *paraNodePool
         );
public:
   /** constructor */
   ParaSolverPool(
         double inMp,
         double inBgap,
         double inMBgap,
         int inOriginRank,
         ParaComm *inParaComm,
         ParaParamSet *inParaParams
         )
         : active(false) ,bgap(inBgap), mp(inMp), mBgap(inMBgap), originRank(inOriginRank), nCollectingModeSolvers(0),
           nNodesSolvedInSolvers(0), nTotalNodesSolved(0), nNodesInSolvers(0),  // rampUpPhase(false), 
           collectingMode(false),  breakingFirstSubtree(false), paraComm(inParaComm), paraParams(inParaParams)
   {
      nSolvers = paraComm->getSize() - inOriginRank;
      nGenerator = nSolvers;
      if( paraParams->getIntParamValue(UgCplexRunningMode) == 1 )  // some solver becomes no-generator
      {
         if( paraParams->getRealParamValue(GeneratorRatio) )
         {
            nGenerator = std::max((int)(nSolvers * paraParams->getRealParamValue(GeneratorRatio)), 1);
         }
      }
      else if ( paraParams->getIntParamValue(UgCplexRunningMode) == 2 ) // No generator
      {
         nGenerator = 0;
      }
      pool = new ParaSolverPoolElementPtr[nSolvers];
      for( int i = 0; i < nSolvers; i++ )
      {
         pool[i] = new ParaSolverPoolElement(originRank+i);
         if( i >= nGenerator ) pool[i]->setNoGenerator();
         inactiveSolvers.insert(std::make_pair((originRank+i),pool[i]));
      }
      if( paraParams->getIntParamValue(MaxNumberOfCollectingModeSolvers) > 0 )
      {
         nMaxCollectingModeSolvers = std::min(
               paraParams->getIntParamValue(MaxNumberOfCollectingModeSolvers), nSolvers );
      }
      else if (paraParams->getIntParamValue(MaxNumberOfCollectingModeSolvers) == 0 )
      {
         nMaxCollectingModeSolvers = std::max(nSolvers / 2, 1);  // at least one solver has to be in collecting mode
      }
      else
      {
         nMaxCollectingModeSolvers = nSolvers;
      }
      nLimitCollectingModeSolvers = paraParams->getIntParamValue(MinNumberOfCollectingModeSolvers);
      absoluteGap = paraParams->getRealParamValue(ABgapForSwitchingToBestSolver);
      if( paraParams->getBoolParamValue(BreakFirstSubtree) )
      {
         breakingFirstSubtree = true;
      }
   }

   /** destructor */
   virtual ~ParaSolverPool(
         )
   {
      for( int i = 0; i < nSolvers; i++ )
      {
         delete pool[i];
      }
      if( pool ) delete[] pool;
   }

   /** activate this pool */
   void activate(
         )
   {
      active = true;
   }

   /** check if this pool is active or not */
   bool isActive(
         )
   {
      return active;
   }

   /** get number of solvers */
   int getNSolvers(
         )
   {
      return nSolvers;
   }

   /** get number of nodes solved in all solvers */
   long long getNnodesSolvedInSolvers(
         )
   {
      return nNodesSolvedInSolvers;
   }

   /** get number of nodes solved in all solvers */
   long long getTotalNodesSolved(
         )
   {
      return nTotalNodesSolved;
   }

   /** get number of nodes solved in all solvers */
   void addTotalNodesSolved(
         long long num
         )
   {
      nTotalNodesSolved += num;
   }

   /** get number of nodes in all solvers */
   long long getNnodesInSolvers(
         )
   {
      return nNodesInSolvers;
   }

   /** get number of active solvers */
   int getNumActiveSolvers(
         )
   {
      return activeSolvers.size();
   }

   /** get number of incactive solvers */
   int getNumInactiveSolvers(
         )
   {
      return inactiveSolvers.size();
   }

   /** get Max numb of solvers that can be in collecting mode */
   bool canIncreaseLimitNLimitCollectingModeSolvers(
         )
   {
      return (nLimitCollectingModeSolvers < nMaxCollectingModeSolvers);
   }

   /** get limit number of solvers that can be in collecting mode */
   int getNLimitCollectingModeSolvers(
         )
   {
      return nLimitCollectingModeSolvers;
   }

   // Node *getCurrentNode(int rank){ return  pool[rank-1]->getCurrentNode(); }

   /** check if this system is in collecting mode or not*/
   bool isInCollectingMode(
         )
   {
      return collectingMode;
   }

   /** check if the solver specified by rank is active or not */
   bool isSolverActive(
         int rank
         )
   {
      return pool[SOLVER_POOL_INDEX(rank)]->isActive();
   }

   /** get collecting mode of the solver specified by rank */
   bool isSolverInCollectingMode(
         int rank
         )
   {
      return pool[SOLVER_POOL_INDEX(rank)]->isInCollectingMode();
   }

   /** get current solving ParaNode by rank of solver */
   ParaNode *getCurrentNode(
         int rank
         )
   {
      return  pool[SOLVER_POOL_INDEX(rank)]->getCurrentNode();
   }

   /** extract current solving ParaNode by rank of solver */
   ParaNode *extractCurrentNodeAndInactivate(
         int rank,                                  /**< rank of the solver to be inactivated */
         ParaNodePool *paraNodePool                 /**< pointer to ParaNodePool to pass it to inactivateSolver routine */
         )
   {
      ParaNode *node = pool[SOLVER_POOL_INDEX(rank)]->extractCurrentNode();
      // for checkpointing, it should be updated always.
      // node->setDualBoundValue( getDualBoundValue(rank) );
      inactivateSolver(rank, 0, paraNodePool);  // number of nodes solved should be added when LC recives termination message */
      return node;
   }

   /** check if solving ParaNode has descendant or not */
   bool currentSolvingNodehaeDescendant(
         int rank
         )
   {
      return pool[SOLVER_POOL_INDEX(rank)]->getCurrentNode()->hasDescendant();
   }

   /** add number of nodes solved */
   void addNumNodesSolved(
         long long numOfNodesSolved                /**< number of nodes solved */
         )
   {
      nNodesSolvedInSolvers += numOfNodesSolved;
   }

   /** get the number of nodes left by rank of solver */
   int getNumOfNodesLeft(
         int rank
         )
   {
      assert(isSolverActive(rank));
      return pool[SOLVER_POOL_INDEX(rank)]->getNumOfNodesLeft();
   }

   /** get the number of nodes left in the solver which has the best dual bound value */
   int getNumOfNodesLeftInBestSolver(
         )
   {
      if( selectionHeap->getHeapSize() > 0 )
      {
         return selectionHeap->top()->getNumOfNodesLeft();
      }
      else
      {
         return 0;   // no nodes exist
      }
   }

   /** get Best Solver rank */
   int getBestSolver(
         )
   {
      if( selectionHeap->getHeapSize() > 0 )
      {
         return selectionHeap->top()->getRank();
      }
      else
      {
         return -1;   // no nodes exist
      }
   }

   /** get dual bound value of solving ParaNode by rank of solver */
   double getDualBoundValue(
         int rank
         )
   {
      assert(isSolverActive(rank));
      return pool[SOLVER_POOL_INDEX(rank)]->getBestDualBoundValue();
   }

   /** set SolverTerminationState */
   void setTermState(
         int rank,
         ParaSolverTerminationState *inTermState
         )
   {
      pool[SOLVER_POOL_INDEX(rank)]->setTermState(inTermState);
   }

   /** get SolverTermination state */
   ParaSolverTerminationState *getTermState(
         int rank
         )
   {
      return pool[SOLVER_POOL_INDEX(rank)]->getTermState();
   }

   /** update dual bound values of saving nodes to their dual bound values for subtrees */
   void updateDualBoundsForSavingNodes(
         )
   {
      for( int i = 1; i < paraComm->getSize(); i++ )
      {
         if( getCurrentNode(i) && getCurrentNode(i)->getAncestor() == 0 )
         {
            getCurrentNode(i)->updateInitialDualBoundToSubtreeDualBound();
         }
      }
   }

   /** write ParaNodes to checkpoint file */
   int writeParaNodesToCheckpointFile(
         ogzstream &out
         )
   {
      int n = 0;
      for( int i = 1; i < paraComm->getSize(); i++ )
      {
         if( getCurrentNode(i) && getCurrentNode(i)->getAncestor() == 0 )
         {
            getCurrentNode(i)->write(out);
            n++;
         }
      }
      return n;
   }

   /** write Solver statistics to checkpoint file */
   int writeSolverStatisticsToCheckpointFile(
         ogzstream &out
         )
   {
      int n = 0;
      for( int i = 1; i < paraComm->getSize(); i++ )
      {
         if( pool[SOLVER_POOL_INDEX(i)]->getTermState() )
         {
            pool[SOLVER_POOL_INDEX(i)]->getTermState()->write(out);
            n++;
         }
      }
      return n;
   }

   /** increase the limit number of solvers getting into collecting mode */
   void incNLimitCollectingModeSolvers(
         )
   {
      nLimitCollectingModeSolvers = std::min(nLimitCollectingModeSolvers+1,
            std::min(nMaxCollectingModeSolvers, nSolvers));
   }

   /** activate the solver specified by rank with specified node which has been sent */
   void activateSolver(
         int rank,
         ParaNode *node
         );


   /** activate the solver specified by rank with specified node which has been sent and nNodesLeft
    *  This method is for the racing winner */
   void activateSolver(
         int rank,
         ParaNode *node,
         int nNodesLeft
         );

   /** activate a solver with specified node which is sent within this SolverPool */
   int activateSolver(
         ParaNode *node,
         ParaRacingSolverPool *paraRacingSolverPool,
         bool rampUpPhase
         );

   /** inactivate the solver specified by rank */
   void inactivateSolver(
         int rank,                                 /**< rank of the solver to be inactivated */
         long long numOfNodesSolved,               /**< number of nodes solved */
         ParaNodePool *paraNodePool                /**< pointer to ParaNodePool to change into collecting mode */
         );

   /** the solver specified by rank died */
   ParaNode *solverDied(
         int rank
         );

   /** send incumbent value */
   void sendIncumbentValue(
         double incumbentValue
         );

   /** switch out collecting mode */
   void switchOutCollectingMode(
         );

   /** switch out collecting mode to the specified rank if it is necessary */
   void enforcedSwitchOutCollectingMode(
         int rank             /**< rank of the solver */
         );

   /** switch out collecting mode to the specified rank if it is necessary */
   void sendSwitchOutCollectingModeIfNecessary(
         int rank             /**< rank of the solver */
         );

   /** get global best dual bound value */
   virtual double getGlobalBestDualBoundValue(
         ) = 0;

   /** switch in collecting mode */
   virtual void switchInCollectingMode(
         ParaNodePool *paraNodePool
         ) = 0;

   /** update solver status */
   virtual void updateSolverStatus(
         int rank,
         long long numNodesSolved,
         int numNodesLeft,
         double solverLocalBestBound,
         ParaNodePool *paraNodePool
         ) = 0;

};

/** Solver Pool for Minimization Problems */
class ParaSolverPoolForMinimization : virtual public ParaSolverPool {
public:
   /** constructor */
   ParaSolverPoolForMinimization(
         double inMp,
         double inBgap,
         double inMBgap,
         int inOriginRank,
         ParaComm *inParaComm,
         ParaParamSet *inParaParams
         )
         : ParaSolverPool(inMp, inBgap, inMBgap, inOriginRank, inParaComm, inParaParams)
   {
      selectionHeap = new AscendingSelectionHeap(nSolvers);
      if( paraParams->getIntParamValue(SolverOrderInCollectingMode) == 0) // ordered by dual bound value
      {
         collectingModeSolverHeap = new AscendingCollectingModeSolverHeap(nMaxCollectingModeSolvers);
      }
   }

   /** destructor */
   ~ParaSolverPoolForMinimization(
         )
   {
      if( selectionHeap ) delete selectionHeap;
      if( collectingModeSolverHeap ) delete collectingModeSolverHeap;
   }

   /** get global best dual bound value */
   double getGlobalBestDualBoundValue(
         )
   {
      if( selectionHeap->getHeapSize() > 0 )
      {
         return selectionHeap->top()->getBestDualBoundValue();
      }
      else
      {
         return DBL_MAX;   // no nodes exist
      }
   }

   /** switch in collecting mode */
   void switchInCollectingMode(
         ParaNodePool *paraNodePool
         );

   /** update solver status */
   void updateSolverStatus(
         int rank,
         long long numNodesSolved,
         int numNodesLeft,
         double solverLocalBestBound,
         ParaNodePool *paraNodePool
         );
};

/** Racing Solver Pool */
class ParaRacingSolverPool
{
   int                                  winnerRank;        /**< winner rank of racing ramp-up, -1: not decided yet */
   int                                  originRank;        /**< origin rank of solvers managed by this solver pool */
   int                                  nSolvers;          /**< number of solvers */
   int                                  nEvaluationStage;  /**< number of solvers that are in evaluation stage */
   long long                            nNodesSolvedInBestSolver; /**< number of nodes solved in the best solver */
   long long                            nNodesInBestSolver;/**< number of nodes in the best solver */
   int                                  nActiveSolvers;    /**< number of active solvers */
   int                                  nInactiveSolvers;  /**< number of inactive solvers */
   double                               bestDualBound;
   double                               bestDualBoundInSolvers;
   ParaSolverPoolElementPtr             *pool;             /**< solver pool indexed by solver's rank */
   SelectionHeap                        *selectionHeap;    /**< pointers to active solvers in ascending or descending order */
   ParaComm                             *paraComm;         /**< communicator */
   ParaParamSet                         *paraParams;       /**< runtime parameters for parallelization */
   ParaNode                             *rootNode;         /**< root node */
   ParaTimer                            *paraTimer;        /**< para timer */
   ParaDeterministicTimer               *paraDetTimer;     /**< para deterministic timer */
public:
   ParaRacingSolverPool(
         int inOriginRank,
         ParaComm *inParaComm,
         ParaParamSet *inParaParams,
         ParaTimer    *inParaTimer,
         ParaDeterministicTimer *inParaDetTimer
         )
   : winnerRank(-1), originRank(inOriginRank),  nEvaluationStage(0), nNodesSolvedInBestSolver(0), nNodesInBestSolver(0), nActiveSolvers(0), nInactiveSolvers(0),
     bestDualBound(-DBL_MAX), bestDualBoundInSolvers(-DBL_MAX),paraComm(inParaComm), paraParams(inParaParams), paraTimer(inParaTimer), paraDetTimer(inParaDetTimer)
   {
       nSolvers = paraComm->getSize() - inOriginRank;
       pool = new ParaSolverPoolElementPtr[nSolvers];
       for( int i = 0; i < nSolvers; i++ )
       {
           pool[i] = new ParaSolverPoolElement(originRank+i);
       }
       selectionHeap = new DescendingSelectionHeap(nSolvers);
   }

   /** destructor */
   ~ParaRacingSolverPool(
         )
   {
      if( selectionHeap ) delete selectionHeap;
      for( int i = 0; i < nSolvers; i++ )
      {
         delete pool[i];
      }
      if( pool ) delete[] pool;
      if( rootNode ) delete rootNode;
   }

   /** extract racing root node */
   ParaNode *extractNode(
         )
   {
      rootNode->setDualBoundValue(getGlobalBestDualBoundValue());
      rootNode->setInitialDualBoundValue(getGlobalBestDualBoundValue());
      return rootNode->clone(paraComm);
   }

   ParaNode *getCurrentNode(
         int rank
         )
   {
      return rootNode;
   }

   /** get dual bound value of solving ParaNode by rank of solver */
   double getDualBoundValue(
         int rank
         )
   {
      return pool[SOLVER_POOL_INDEX(rank)]->getBestDualBoundValue();
   }

   /** get number of nodes solved */
   long long getNumOfNodesSolved(
         int rank
         )
   {
      return pool[SOLVER_POOL_INDEX(rank)]->getNumOfNodesSolved();
   }

   /** get number of nodes left in the solver specified by rank */
   int getNumNodesLeft(
         int rank
         )
   {
      return pool[SOLVER_POOL_INDEX(rank)]->getNumOfNodesLeft();
   }

   /** get global best dual bound value */
   double getGlobalBestDualBoundValue(
         )
   {
      if( selectionHeap->getHeapSize() > 0 )
      {
         return bestDualBound;
      }
      else
      {
         return -DBL_MAX;   // no nodes exist
      }
   }

   /** get winner rank */
   int getWinner(
         )
   {
      // assert( winnerRank > 0 );
      return winnerRank;   // -1 means that winner is not decided
   }

   /** get number of nodes in the best solver */
   long long getNnodesSolvedInBestSolver(
         )
   {
      return nNodesSolvedInBestSolver;
   }

   /** get number of nodes in the best solver */
   long long getNnodesLeftInBestSolver(
         )
   {
      return nNodesInBestSolver;
   }

   /** get dual bound value in inactivated solvers */
   double getBestDualBoundInInactivatedSolvers(
         )
   {
      return bestDualBoundInSolvers;
   }

   /** activate racing ramp-up solver pool */
   void activate(
         ParaNode *node
         )
   {
      rootNode = node;
      for( int rank = originRank; rank < (originRank + nSolvers); rank++ )
      {
          pool[SOLVER_POOL_INDEX( rank )]->racingActivate();
          if( rank != originRank ||
                ( paraParams->getIntParamValue(UgCplexRunningMode) != 1 ) )
          {
             selectionHeap->insert(pool[SOLVER_POOL_INDEX( rank )]);  // this should be called after activate: dual bound value need to be set
          }
          nActiveSolvers++;
      }
      nNodesInBestSolver = 1;
   }

   /** check if solver is active or not */
   bool isActive(
         int rank
         )
   {
      if(pool[SOLVER_POOL_INDEX( rank )]->isRacingStage() ||
            pool[SOLVER_POOL_INDEX( rank )]->isEvaluationStage() )
      {
         return true;
      }
      else
      {
         return false;
      }
   }

   /** update solver status */
   void updateSolverStatus(
         int rank,
         long long numNodesSolved,
         int numNodesLeft,
         double solverLocalBestBound
         );

   /** check racing termination criteria */
   bool isWinnerDecided(bool feasibleSol);

   /** inactivate the solver specified by rank */
   void inactivateSolver(
         int rank
         )
   {
      if(pool[SOLVER_POOL_INDEX( rank )]->isRacingStage() ||
            pool[SOLVER_POOL_INDEX( rank )]->isEvaluationStage() )
      {
         if( pool[SOLVER_POOL_INDEX( rank )]->getBestDualBoundValue() > bestDualBoundInSolvers )
         {
            bestDualBoundInSolvers = pool[SOLVER_POOL_INDEX( rank )]->getBestDualBoundValue();
         }
         pool[SOLVER_POOL_INDEX( rank )]->inactivate();
         nInactiveSolvers++;
         nActiveSolvers--;
      }
      else
      {
          THROW_LOGICAL_ERROR2("Rank = ", rank);
      }
   }

   /** get number of active solvers */
   int getNumActiveSolvers(
         )
   {
      return nActiveSolvers;
   }

   /** get number of inactive solvers */
   int getNumInactiveSolvers(
         )
   {
      return nInactiveSolvers;
   }

   /** check if the solver specified in an argument is evaluation stage or not */
   bool isEvaluationStage(
         int rank
         )
   {
      return ( pool[SOLVER_POOL_INDEX( rank )]->isEvaluationStage() );
   }

   /** get active solver number string */
   std::string getStrActiveSolerNumbers(
         )
   {
      std::ostringstream oss;
      oss << "   ";
      for( int rank = originRank; rank < (originRank + nSolvers); rank++ )
      {
          if( pool[SOLVER_POOL_INDEX( rank )]->isRacingStage() )
          {
             oss << rank << " ";

          }
      }
      return oss.str();
   }
};

}

#endif // __PARA_SOLVER_POOL_H__

