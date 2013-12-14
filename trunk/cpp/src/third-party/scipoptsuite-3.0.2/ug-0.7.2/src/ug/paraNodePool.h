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

/**@file    paraNodePool.h
 * @brief   ParaNode Pool.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_NODE_POOL_H__
#define __PARA_NODE_POOL_H__
#include <map>
#include <queue>
#include <cfloat>
#include <cmath>
#include "paraNode.h"
#include "paraInstance.h"
#include "paraMergeNodesStructs.h"

namespace UG
{

static const double eps = 1.0e-12;

class ParaNodeSortCriterion
{
public:
   bool operator()(const ParaNodePtr& n1, const ParaNodePtr& n2) const
   {

      return EPSLT(n1->getDualBoundValue(),n2->getDualBoundValue(),eps) ||
              ( EPSEQ(n1->getDualBoundValue(),n2->getDualBoundValue(),eps) &&
                  n1->getDiffSubproblem() && n2->getDiffSubproblem() &&
                  n1->getDiffSubproblem()->getNBoundChanges() < n2->getDiffSubproblem()->getNBoundChanges() );

   }
};

class ParaNodePool
{
protected:
   double        bgap;
   unsigned int  maxUsageOfPool;
   std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion > ascendingPool;
public:
   ParaNodePool(double inBgap) : bgap(inBgap), maxUsageOfPool(0)  {}
   virtual ~ParaNodePool(){}
   virtual void insert(ParaNodePtr node) = 0;
   virtual bool isEmpty() = 0;
   virtual ParaNodePtr extractNode() = 0;
   virtual double getBestDualBoundValue() = 0;
   virtual int getNumOfGoodNodes(double globalBestBound) = 0;
   virtual int getNumOfNodes() = 0;
   virtual int removeBoundedNodes(double incumbentValue) = 0;
   virtual void updateDualBoundsForSavingNodes() = 0;
   virtual int writeParaNodesToCheckpointFile(ogzstream &out) = 0;
   virtual const std::string toString() = 0;
   int getMaxUsageOfPool(){ return maxUsageOfPool; }
   int removeMergedNodes(
         ParaMergedNodeListElement *head
         )
   {
      int nDeleted = 0;
      assert( ascendingPool.size() > 0 );
      std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion >::iterator p;
      for( p = ascendingPool.begin(); p != ascendingPool.end() && head; )
      {
         assert( p->second );
         if( p->second->getMergingStatus() == 2 )
         {
            ParaMergedNodeListElement *prev = head;
            ParaMergedNodeListElement *cur = head;
            for( ; cur; cur=cur->next )
            {
               if( p->second == cur->node )
               {
                  break;
               }
               prev = cur;
            }
            assert(cur);
            if( prev == head )
            {
               if( cur == prev )
               {
                  head = head->next;
               }
               else
               {
                  assert( cur == prev->next );
                  prev->next = prev->next->next;
               }
            }
            else
            {
               assert( cur == prev->next );
               prev->next = prev->next->next;
            }
            assert(  p->second == cur->node );
            assert( cur->node->getMergeNodeInfo() );
            delete cur;
            nDeleted++;
            delete p->second;
            ascendingPool.erase(p++);
         }
         else
         {
            p++;
         }
      }
      assert(!head);
      return nDeleted;
   }
};

class ParaNodePoolForMinimization : virtual public ParaNodePool
{

public:
   ParaNodePoolForMinimization(double inBgap) : ParaNodePool(inBgap)  {}
   ~ParaNodePoolForMinimization(
         )
   {
      if( ascendingPool.size() > 0 )
      {
         std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion >::iterator p;
         for( p = ascendingPool.begin(); p != ascendingPool.end(); )
         {
            if( p->second ) delete p->second;
            ascendingPool.erase(p++);
         }
      }
   }

   void insert(
         ParaNodePtr paraNode
         )
   {
      ascendingPool.insert(std::make_pair(paraNode,paraNode));
      if( maxUsageOfPool < ascendingPool.size() )
      {
         maxUsageOfPool = ascendingPool.size();
      }
   }

   bool isEmpty(){ return ( ascendingPool.size() == 0 ); }

   ParaNodePtr extractNode(
         )
   {
      ParaNodePtr extracted = 0;
      std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion >::iterator p;
      p = ascendingPool.begin();
      while( p != ascendingPool.end() )
      {
         if( p->second->getMergeNodeInfo() )
         {
            assert( p->second->getMergeNodeInfo()->status != ParaMergeNodeInfo::PARA_MERGING );
            if( p->second->getMergingStatus() == 4 )  // merging representative was already deleted
            {
               delete p->second;              // this delete is not counted in statistics
               ascendingPool.erase(p++);
            }
            else
            {
               if( p->second->getMergeNodeInfo()->status == ParaMergeNodeInfo::PARA_MERGE_CHECKING_TO_OTHER_NODE )
               {
                  assert(dynamic_cast<ParaNodePtr>(p->second)->getMergeNodeInfo()->mergedTo->status == ParaMergeNodeInfo::PARA_MERGED_RPRESENTATIVE );
                  p++;
               }
               else
               {
                  extracted = p->second;
                  ascendingPool.erase(p);
                  break;
               }
            }
         }
         else
         {
            extracted = p->second;
            ascendingPool.erase(p);
            break;
         }
      }
      assert( ( p == ascendingPool.end() && (!extracted) ) || ( p != ascendingPool.end() && extracted ) );
      return extracted;
   }

   void updateDualBoundsForSavingNodes(
         )
   {
      std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion >::iterator p;
      for( p = ascendingPool.begin(); p != ascendingPool.end(); ++p )
      {
         if( p->second->getAncestor() == 0 )
         {
            p->second->updateInitialDualBoundToSubtreeDualBound();
         }
      }
   }

   int writeParaNodesToCheckpointFile(
         ogzstream &out
         )
   {
      int n = 0;
      std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion >::iterator p;
      for( p = ascendingPool.begin(); p != ascendingPool.end(); ++p )
      {
         if( p->second->getAncestor() == 0 )
         {
            p->second->write(out);
            n++;
         }
      }
      return n;
   }

   double getBestDualBoundValue(
         )
   {
      std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion >::iterator p;
      p = ascendingPool.begin();
      if( p != ascendingPool.end() )
      {
         return p->second->getDualBoundValue();
      }
      else
      {
         return DBL_MAX;  // no nodes exist
      }
   }

   int getNumOfGoodNodes(
         double globalBestBound
         )
   {
      /** The following code is not a good idea,
       * because only a node is received from a solver, LC can switch out
      if( globalBestBound > getBestDualBoundValue() )
         globalBestBound = getBestDualBoundValue();
      */
      int num = 0;
      std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion >::iterator p;
      for( p = ascendingPool.begin(); p != ascendingPool.end() &&
              ( ( ( p->second->getDualBoundValue() ) - globalBestBound ) /
                    std::max( fabs(globalBestBound) , 1.0 ) ) < bgap;
            ++p )
      {
         num++;
      }
      return num;
   }

   int getNumOfNodes(){ return ascendingPool.size(); }

   int removeBoundedNodes(
         double incumbentValue
         )
   {
      int nDeleted = 0;
      if( ascendingPool.size() > 0 )
      {
         std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion >::iterator p;
         for( p = ascendingPool.begin(); p != ascendingPool.end(); )
         {
            assert( p->second );
            if( !p->second->getMergeNodeInfo() )
            {
               if( p->second->getDualBoundValue() > incumbentValue || p->second->getMergingStatus() == 4 )
               {
                  nDeleted++;
                  delete p->second;
                  ascendingPool.erase(p++);
               }
               else
               {
                  p++;
               }
            }
            else
            {
               if( p->second->getMergeNodeInfo()->status == ParaMergeNodeInfo::PARA_MERGE_CHECKING_TO_OTHER_NODE )
               {
                  if( p->second->getDualBoundValue() > incumbentValue || p->second->getMergingStatus() == 4 )
                   {
                      nDeleted++;
                      delete p->second;
                      ascendingPool.erase(p++);
                   }
                   else
                   {
                      p++;
                   }
               }
               else
               {
                  p++;
               }
            }
         }
      }
      return nDeleted;
   }

   const std::string toString(
         )
   {
      std::ostringstream s;
      std::multimap<ParaNodePtr, ParaNodePtr, ParaNodeSortCriterion >::iterator p;
      for( p = ascendingPool.begin(); p != ascendingPool.end(); ++p )
      {
         s << p->second->toString();
      }
      return s.str();
   }
};

}

#endif // __PARA_NODE_POOL_H__
