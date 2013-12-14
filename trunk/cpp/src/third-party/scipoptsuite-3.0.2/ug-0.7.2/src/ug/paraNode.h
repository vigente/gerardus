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

/**@file    paraNode.h
 * @brief   Base class for ParaNode.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_NODE_H__
#define __PARA_NODE_H__

#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include "ug/paraDef.h"
#include "ug/paraComm.h"
#include "ug/gzstream.h"
#include "ug/paraMergeNodesStructs.h"
#include "ug/paraDiffSubproblem.h"

namespace UG
{

static const int ParaNodeLocalPtr  = 0;
static const int ParaNodeRemotePtr = 1;

/** SubtreeId class */
class SubtreeId
{
public:
   int lcId;                 /**< LoadCoordinator ID */
   int globalSubtreeIdInLc;  /**< Global Subtree ID in Solvers managed by LoadCoordinator */
   int solverId;             /**< Solver ID */

   /** default constructor */
   SubtreeId(
         )
         : lcId(-1), globalSubtreeIdInLc(-1), solverId(-1)
   {
   }

   /** constructor */
   SubtreeId(
         int inLcId, int inGlobalSubtreeIdInLc, int inSolverId
         )
        : lcId(inLcId), globalSubtreeIdInLc(inGlobalSubtreeIdInLc), solverId(inSolverId)
   {
   }

   /** destructor */
   ~SubtreeId(
         )
   {
   }

   /** == operator definition */
   bool operator == (
         const SubtreeId& inSid
         ) const
   {
      if( lcId == inSid.lcId &&
            globalSubtreeIdInLc == inSid.globalSubtreeIdInLc &&
            solverId == inSid.solverId )
         return true;
      else
         return false;
   }

   /** != operator definition */
   bool operator != (
         const SubtreeId& inSid
         ) const
   {
      if( lcId != inSid.lcId ||
            globalSubtreeIdInLc != inSid.globalSubtreeIdInLc ||
            solverId != inSid.solverId )
         return true;
      else
         return false;
   }

   /** < operator definition */
   bool operator < (
         const SubtreeId& inSid
         ) const
   {
      if( lcId < inSid.lcId )
         return false;
      if( lcId == inSid.lcId && globalSubtreeIdInLc < inSid.globalSubtreeIdInLc )
         return true;
      if( lcId == inSid.lcId && globalSubtreeIdInLc == inSid.globalSubtreeIdInLc && solverId < inSid.solverId )
         return true;
      else
         return false;
   }

   /** getter of LoadCoordinator ID */
   int getLcId(
         ) const
   {
      return lcId;
   }

   /** getter of global Subtree ID in Solvers managed by the LoadCoordinator */
   int getGlobalSubtreeIdInLc(
         ) const
   {
      return globalSubtreeIdInLc;
   }

   /** getter of Solver ID */
   int getSolverId(
         ) const
   {
      return solverId;
   }

   /** Stringfy SubtreeId */
   std::string toString(
         )
   {
      std::ostringstream s;
      s << "(" << lcId << "," << globalSubtreeIdInLc << "," << solverId << ")";
      return s.str();
   }
};

/** NodeId class */
class NodeId
{
public:
    SubtreeId subtreeId; /**< Subtree ID */
    long long seqNum;    /**< Sequential Number in the Subtree */

    /** default constructor */
    NodeId(
          )
          : subtreeId(SubtreeId()), seqNum(-1)
    {
    }

    /** constructor */
    NodeId(
          SubtreeId inSubtreeId, int inSeqNum
          )
          : subtreeId(inSubtreeId), seqNum(inSeqNum)
    {
    }

    /** destructor */
    ~NodeId(
          )
    {
    }

    /** == operator definition */
    bool operator == (
          const NodeId& inNid
          ) const
   {
      if( subtreeId == inNid.subtreeId &&
             seqNum == inNid.seqNum )
          return true;
      else
           return false;
   }

   /** != operator definition */
   bool operator != (
         const NodeId& inNid
         ) const
   {
      if( subtreeId != inNid.subtreeId ||
            seqNum != inNid.seqNum )
         return true;
      else
         return false;
   }

   /** < operator definition */
   bool operator < (
         const NodeId& inNid
         ) const
   {
      if( subtreeId < inNid.subtreeId ) return true;
      if( subtreeId == inNid.subtreeId && seqNum < inNid.seqNum )
         return true;
      else
         return false;
   }

   /** getter of Subtree ID */
   SubtreeId getSubtreeId(
         ) const
   {
      return subtreeId;
   }

   /** getter of Sequence Number */
   int getSeqNum(
         ) const
   {
      return seqNum;
   }

   /** Stringfy NodeId */
   std::string toString(
         )
   {
      std::ostringstream s;
      // for debug
      s << "[" << (subtreeId.toString()) << ":" <<  seqNum << "]";
      return s.str();
   }
};

/** pointer to indicate a ParaNode genealogical relation */
class ParaNode;
class ParaNodeGenealogicalPtr
{
   NodeId   genealogicalNodeId;        /**< descendant NodeId or ascendant NodeId */
public:
   ParaNodeGenealogicalPtr(NodeId nodeId) : genealogicalNodeId(nodeId) {}
   virtual ~ParaNodeGenealogicalPtr() {}
   virtual int getType() = 0;
   NodeId   getNodeId() { return genealogicalNodeId; }
};
class ParaNodeGenealogicalLocalPtr : public ParaNodeGenealogicalPtr
{
   ParaNode     *paraNodePtr;
public:
   ParaNodeGenealogicalLocalPtr() : ParaNodeGenealogicalPtr(NodeId()), paraNodePtr(0) {}
   ParaNodeGenealogicalLocalPtr(NodeId nodeId, ParaNode *ptr) : ParaNodeGenealogicalPtr(nodeId), paraNodePtr(ptr) {}
   ~ParaNodeGenealogicalLocalPtr(){}
   int getType(){ return ParaNodeLocalPtr; }
   ParaNode *getPointerValue() { return paraNodePtr; }
};

class ParaNodeGenealogicalRemotePtr : public ParaNodeGenealogicalPtr
{
   int       transferringLcId;          /**< LoadCoordinator Id transfer to or trasferred from */
public:
   ParaNodeGenealogicalRemotePtr() : ParaNodeGenealogicalPtr(NodeId()), transferringLcId(-1){}
   ParaNodeGenealogicalRemotePtr(NodeId nodeId, int lcId) : ParaNodeGenealogicalPtr(nodeId), transferringLcId(lcId) {}
   ~ParaNodeGenealogicalRemotePtr(){}
   int getType(){ return ParaNodeRemotePtr; }
   int getPointerValue() { return transferringLcId; }
};

typedef ParaNodeGenealogicalPtr *ParaNodeGenealogicalPtrPtr;

/** ParaNode class */
class ParaNode
{
public:
   // solving node information
   NodeId          nodeId;                 /**< solving node ID */
   NodeId          generatorNodeId;        /**< subtree root node ID of generator */
   ParaNodeGenealogicalPtr *ancestor;      /**< pointer to ancestor ParaNode : This field is not transferred */
   std::map< NodeId, ParaNodeGenealogicalPtrPtr > descendants; /**< collection of pointers to descendants : This filed is not transferred  */
protected:
   int             depth;                  /**< depth from the root node of original tree */
   double          dualBoundValue;         /**< dual bound value */
   double          initialDualBoundValue;  /**< dual bound value when this node is created.
                                                This value is updated to precise one when there is guarantee */
   double          estimatedValue;         /**< estimate value */
   int             diffSubproblemInfo;     /**< 1: with diffSubproblem, 0: no diffSubproblem */
   ParaDiffSubproblem *diffSubproblem;     /**< difference between solving instance data and subproblem data */
   int             basisInfo;
   int             mergingStatus;          /**< merging status:
                                                  -1 - no merging node,
                                                   0 - checking,
                                                   1 - merged (representative)
                                                   2 - merged to the other node
                                                   3 - cannot be merged
                                                   4 - merging representative was deleted
                                                   */
   ParaMergeNodeInfo *mergeNodeInfo;       /**< pointer to mergeNodeInfo. Not zero means merging */

public :
   /** default constructor */
   ParaNode(
         )
         : nodeId(NodeId()), generatorNodeId(NodeId()), ancestor(0),
         depth(-1), dualBoundValue(-DBL_MAX), initialDualBoundValue(0.0),
         estimatedValue(0.0), diffSubproblemInfo(0), diffSubproblem(0), mergingStatus(-1), mergeNodeInfo(0)
   {
   }

   /** constructor */
   ParaNode( NodeId inNodeId,
         NodeId inGeneratorNodeId,
         int inDepth,
         double inDualBoundValue,
         double inOriginalDualBoundValue,
         double inEstimatedValue,
         ParaDiffSubproblem *inDiffSubproblem
         )
         : nodeId(inNodeId), generatorNodeId(inGeneratorNodeId), ancestor(0),
         depth(inDepth), dualBoundValue(inDualBoundValue), initialDualBoundValue(inOriginalDualBoundValue),
         estimatedValue(inEstimatedValue), diffSubproblem(inDiffSubproblem), mergingStatus(-1), mergeNodeInfo(0)
   {
      if( diffSubproblem ) diffSubproblemInfo = 1;
      else diffSubproblemInfo = 0;
      basisInfo = 0;
   }

   /** destructor */
   virtual ~ParaNode(
         )
   {
      assert((mergingStatus != -1) || (mergingStatus == -1 && mergeNodeInfo == 0) );
      if( diffSubproblem ) delete diffSubproblem;
      if( ancestor )
      {
         ParaNodeGenealogicalLocalPtr *localPtrAncestor = dynamic_cast< ParaNodeGenealogicalLocalPtr * >(ancestor);
         if( !descendants.empty() )
         {
            std::map< NodeId, ParaNodeGenealogicalPtrPtr >::iterator pos;
            for( pos = descendants.begin(); pos != descendants.end(); )
            {
               if( pos->second->getType() == ParaNodeLocalPtr )
               {
                  ParaNodeGenealogicalLocalPtr *localPtrDescendant = dynamic_cast< ParaNodeGenealogicalLocalPtr * >(pos->second);
                  assert( localPtrDescendant->getPointerValue()->ancestor->getNodeId() == nodeId );
                  assert( localPtrAncestor->getNodeId() == localPtrAncestor->getPointerValue()->nodeId );
                  assert( localPtrDescendant->getNodeId() == localPtrDescendant->getPointerValue()->nodeId );
                  localPtrDescendant->getPointerValue()->setAncestor(
                        new ParaNodeGenealogicalLocalPtr( localPtrAncestor->getNodeId(), localPtrAncestor->getPointerValue() ) );
                  localPtrAncestor->getPointerValue()->addDescendant(
                        new ParaNodeGenealogicalLocalPtr( localPtrDescendant->getNodeId(), localPtrDescendant->getPointerValue() ) );
               }
               else
               {  /** not implemented yet **/
                  THROW_LOGICAL_ERROR1("remote pointer is not implemented yet, but it is called!");
               }
               delete pos->second;
               descendants.erase(pos++);
            }
         }
         if( ancestor->getType() == ParaNodeLocalPtr )
         {
             assert( localPtrAncestor->getNodeId() == localPtrAncestor->getPointerValue()->nodeId );
             localPtrAncestor->getPointerValue()->removeDescendant(nodeId);
         }
         else
         {   /** not implemented yet **/
             THROW_LOGICAL_ERROR1("remote pointer is not implemented yet, but it is called!");
         }
         delete ancestor;
      }
      else
      {
         if( !descendants.empty() )
         {
            std::map< NodeId, ParaNodeGenealogicalPtrPtr >::iterator pos;
            for( pos = descendants.begin(); pos != descendants.end(); )
            {
               if( pos->second->getType() == ParaNodeLocalPtr )
               {
                  ParaNodeGenealogicalLocalPtr *localPtrDescendant = dynamic_cast< ParaNodeGenealogicalLocalPtr * >(pos->second);
                  localPtrDescendant->getPointerValue()->setAncestor(0);
               }
               else
               {  /** not implemented yet **/
                  THROW_LOGICAL_ERROR1("remote pointer is not implemented yet, but it is called!");
               }
               delete pos->second;
               descendants.erase(pos++);
            }
         }
      }
      if( mergeNodeInfo )
      {
         if( mergeNodeInfo->nMergedNodes == 0 && mergeNodeInfo->mergedTo )
         {
            assert(mergeNodeInfo->status == ParaMergeNodeInfo::PARA_MERGE_CHECKING_TO_OTHER_NODE);
            assert(mergeNodeInfo->mergedTo->status ==  ParaMergeNodeInfo::PARA_MERGED_RPRESENTATIVE);
            mergeNodeInfo->mergedTo->nMergedNodes--;
         }

         assert( mergeNodeInfo->status != ParaMergeNodeInfo::PARA_MERGING);

         if( mergeNodeInfo->status == ParaMergeNodeInfo::PARA_MERGED_RPRESENTATIVE )
         {
            for( ParaFixedVariable *traverse = mergeNodeInfo->fixedVariables[mergeNodeInfo->keyIndex].next;
                  traverse;
                  traverse = traverse->next )
            {
               if( traverse->mnode->nMergedNodes == 0 && mergeNodeInfo == traverse->mnode->mergedTo )
               {
                  traverse->mnode->mergedTo->nMergedNodes--;
                  traverse->mnode->mergedTo = 0;
                  if( traverse->mnode->paraNode->getDualBoundValue() < mergeNodeInfo->paraNode->getDualBoundValue() )
                  {
                     traverse->mnode->paraNode->setMergingStatus(0);
                     traverse->mnode->status = ParaMergeNodeInfo::PARA_MERGING;
                     traverse->mnode->nMergedNodes = -1;
                     traverse->mnode->nSameValueVariables = -1;
                     traverse->mnode->keyIndex = -1;
                  }
                  else
                  {
                     traverse->mnode->paraNode->setMergingStatus(4);  // merging representative was deleted -> this node should be deleted
                     traverse->mnode->status = ParaMergeNodeInfo::PARA_DELETED;
                     traverse->mnode->nMergedNodes = -1;
                     traverse->mnode->nSameValueVariables = -1;
                     traverse->mnode->keyIndex = -1;
                  }
               }
            }
         }

         if( mergeNodeInfo->fixedVariables )
         {
            for( int i = 0; i < mergeNodeInfo->nFixedVariables; i++ )
            {
               for( ParaFixedVariable *traverse = mergeNodeInfo->fixedVariables[i].prev;
                     traverse;
                     traverse = traverse->prev
                     )
               {
                  traverse->nSameValue--;
               }
               if( mergeNodeInfo->fixedVariables[i].prev )
               {
                  mergeNodeInfo->fixedVariables[i].prev->next = mergeNodeInfo->fixedVariables[i].next;
                  if( mergeNodeInfo->fixedVariables[i].next )
                  {
                     mergeNodeInfo->fixedVariables[i].next->prev = mergeNodeInfo->fixedVariables[i].prev;
                  }
                  else
                  {
                     mergeNodeInfo->fixedVariables[i].prev->next = 0;
                  }
               }
               else
               {
                  if( mergeNodeInfo->fixedVariables[i].next )
                  {
                     mergeNodeInfo->fixedVariables[i].next->prev = 0;
                  }
               }
            }
            delete [] mergeNodeInfo->fixedVariables;
         }
         delete mergeNodeInfo;
      }
   }

   /** check if root node or not */
   bool isRootNode(){
      // we want to know on which solver on which LC is managed is to generate the root
      if( nodeId.subtreeId.lcId == -1 &&
            nodeId.subtreeId.globalSubtreeIdInLc == -1 &&
            nodeId.subtreeId.solverId == -1 &&
            nodeId.seqNum == -1 ) return true;
      else return false;
   }

   /** check if this node ID is the same as argument node ID */
   bool isSameNodeIdAs(
         const ParaNode& inNode
         )
   {
      if( nodeId == inNode.nodeId ) return true;
      else return false;
   }

   /** check if this node's parent ID is the same as that of argument node ID */
   bool isSameParetntNodeIdAs(
         const ParaNode& inNode
         )
   {
      if( generatorNodeId == inNode.generatorNodeId ) return true;
      else return false;
   }

   /** check if this node's parent ID is the same as that of argument node ID */
   bool isSameParetntNodeSubtreeIdAs(
         const NodeId& inNodeId
         )
   {
      if( generatorNodeId.subtreeId == inNodeId.subtreeId ) return true;
      else return false;
   }

   /** check if this node's subtree ID is the same as that of argument node ID */
   bool isSameSubtreeIdAs(
         const ParaNode& inNode
         )
   {
      if( nodeId.subtreeId == inNode.nodeId.subtreeId )
         return true;
      else return false;
   }

   /** check if this node's Global subtree ID in LC is the same as that of argument node ID */
   bool isSameLcIdAs(
         const ParaNode& inNode
         )
   {
      if( nodeId.subtreeId.lcId
            == inNode.nodeId.subtreeId.lcId )
         return true;
         else return false;
   }

   /** check if this node's Global subtree ID in LC is the same as that of argument */
   bool isSameLcIdAs(
         const int lcId
         )
   {
      if( nodeId.subtreeId.lcId == lcId )
         return true;
      else return false;
   }

   /** check if this node's Global subtree ID in LC is the same as that of argument node ID */
   bool isSameGlobalSubtreeIdInLcAs(
         const ParaNode& inNode
         )
   {
      if( nodeId.subtreeId.globalSubtreeIdInLc
            == inNode.nodeId.subtreeId.globalSubtreeIdInLc )
         return true;
      else return false;
   }

   /** check if this node's Global subtree ID in LC is the same as that of argument */
   bool isSameGlobalSubtreeIdInLcAs(
         const int globalSubtreeIdInLc
         )
   {
      if( nodeId.subtreeId.globalSubtreeIdInLc == globalSubtreeIdInLc )
         return true;
      else return false;
   }

   /** getter of  LoadCoordinator ID */
   int getLcId(
         )
   {
      return nodeId.subtreeId.lcId;
   }

   /** getter of  global subtree ID in Solvers managed by LoadCoordinator */
   int getGlobalSubtreeIdInLc(
         )
   {
      return nodeId.subtreeId.globalSubtreeIdInLc;
   }

   /** setter of global subtree ID */
   void setGlobalSubtreeId(
         int lcId,
         int subtreeId
         )
   {
      nodeId.subtreeId.lcId = lcId;
      nodeId.subtreeId.globalSubtreeIdInLc = subtreeId;
   }

   /** getter of Solver ID */
   int getSolverId(
         )
   {
      return nodeId.subtreeId.solverId;
   }

   /** setter of Solver ID */
   void setSolverId(
         int id
         )
   {
      nodeId.subtreeId.solverId = id;
   }

   /** getter of node ID */
   NodeId getNodeId(
         )
   {
      return nodeId;
   }

   /** setter of node ID */
   void setNodeId(
         NodeId inNodeId
         )
   {
      nodeId = inNodeId;
   }

   /** getter of generator node ID */
   NodeId getGeneratorNodeId(
         )
   {
      return generatorNodeId;
   }

   /** setter of generator node ID */
   void setGeneratorNodeId(
         NodeId inGeneratorNodeId
         )
   {
      generatorNodeId = inGeneratorNodeId;
   }

   /** getter of depth */
   int getDepth(
         )
   {
      return depth;
   }

   /** setter of depth */
   void setDepth(
         int inDepth
         )
   {
      depth = inDepth;
   }

   /** getter of dual bound */
   double getDualBoundValue(
         )
   {
      return dualBoundValue;
   }

   /** getter of true dual bound */
   double getInitialDualBoundValue(
         )
   {
      return initialDualBoundValue;
   }

   /** setter of dual bound */
   void setDualBoundValue(
         double inDualBoundValue
         )
   {
      dualBoundValue = inDualBoundValue;
   }

   /** setter of true dual bound */
   void setInitialDualBoundValue(
         double inTrueDualBoundValue
         )
   {
      initialDualBoundValue = inTrueDualBoundValue;
   }

   /** reset dual bound value */
   void resetDualBoundValue()
   {
      dualBoundValue = initialDualBoundValue;
   }

   /** getter of estimated value */
   double getEstimatedValue(
         )
   {
      return estimatedValue;
   }

   /** setter of estimated value */
   void setEstimatedValue(
         double inEstimatedValue
         )
   {
      estimatedValue = inEstimatedValue;
   }

   /** getter of diffSubproblem */
   ParaDiffSubproblem *getDiffSubproblem(
         )
   {
      return diffSubproblem;
   }

   /** setter of diffSubproblem */
   void setDiffSubproblem(
         ParaDiffSubproblem *inDiffSubproblem
         )
   {
      diffSubproblem = inDiffSubproblem;
   }

   /** getter of ancestor */
   ParaNodeGenealogicalPtr *getAncestor(
         )
   {
      return ancestor;
   }

   /** setter of ancestor */
   void setAncestor(
         ParaNodeGenealogicalPtr *inAncestor
         )
   {
      if( ancestor ) delete ancestor;
      ancestor = inAncestor;
   }

   /** remove a descendant */
   void removeDescendant(
         NodeId removeNodeId
         )
   {
      std::map< NodeId, ParaNodeGenealogicalPtrPtr >::iterator pos;
      pos = descendants.find(removeNodeId);
      if( pos != descendants.end() )
      {
         delete pos->second;
         descendants.erase(pos);
      }
      else
      {
         for( pos = descendants.begin(); pos != descendants.end(); )
         {
            if( pos->second->getType() == ParaNodeLocalPtr )
            {
               ParaNodeGenealogicalLocalPtr *localPtrDescendant = dynamic_cast< ParaNodeGenealogicalLocalPtr * >(pos->second);
               std::cout << "Descendant NodeId = " << localPtrDescendant->getNodeId().toString() << std::endl;
            }
            else
            {
               /** not implemented yet */
            }
            pos++;
         }
         THROW_LOGICAL_ERROR1("invalid NodeId removed!");
      }
   }

   /** check if this node has descendant or not */
   bool hasDescendant(
         )
   {
      return !(descendants.empty());
   }

   /** add a descendant */
   void addDescendant(
         ParaNodeGenealogicalPtr *inDescendant
         )
   {
      descendants.insert(std::make_pair(inDescendant->getNodeId(),inDescendant));
   }

   void updateInitialDualBoundToSubtreeDualBound(
         )
   {
      // dualBoundValue = getMinimumDualBoundInDesendants(dualBoundValue);
      /** More accurate dual bound of this node is obtained */
      initialDualBoundValue = getMinimumDualBoundInDesendants(dualBoundValue);
   }

   double getMinimumDualBoundInDesendants(
         double value
         )
   {
      if( descendants.empty() ) return value;
      std::map< NodeId, ParaNodeGenealogicalPtrPtr >::iterator pos;
      for( pos = descendants.begin(); pos != descendants.end(); )
      {
         if( pos->second->getType() == ParaNodeLocalPtr )
         {
            ParaNodeGenealogicalLocalPtr *localPtrDescendant = dynamic_cast< ParaNodeGenealogicalLocalPtr * >(pos->second);
            value = std::min( value, localPtrDescendant->getPointerValue()->getDualBoundValue() );
            value = std::min(
                  value,
                  localPtrDescendant->getPointerValue()->getMinimumDualBoundInDesendants(value));
         }
         else
         {
            /** not implemented yet */
         }
         pos++;
      }
      return value;
   }

   virtual ParaNode* clone(ParaComm *comm) = 0;
   virtual int bcast( ParaComm *comm, int root) = 0;
   virtual int send( ParaComm *comm, int destination ) = 0;
   virtual int receive( ParaComm *comm, int source ) = 0;

   void write(ogzstream &out);
   bool read(ParaComm *comm, igzstream &in, bool onlyBoundChanges);

   const std::string toString(
         )
   {
      std::ostringstream s;
      s << "ParaNodeId = " << (nodeId.toString()) << ", GeneratorNodeId = " << (generatorNodeId.toString())
      << ", depth = " << depth << ", dual bound value = " << dualBoundValue
      << ", initialDualBoundValue = " << initialDualBoundValue
      << ", estimated value = " << estimatedValue << std::endl;
      if( diffSubproblem )
      {
         // s << diffSubproblem->toString();
      }
      return s.str();
   }

   const std::string toSimpleString(
         )
   {
      std::ostringstream s;
      s << nodeId.toString()
            << ", "
            << generatorNodeId.toString()
            << ", "
            << depth
            << ", "
            << initialDualBoundValue
            << ", "
            << dualBoundValue;
      return s.str();
   }

   void setMergeNodeInfo( ParaMergeNodeInfo *mNode)
   {
      assert(mergingStatus != -1);
      mergeNodeInfo = mNode;
   }
   ParaMergeNodeInfo *getMergeNodeInfo() { return mergeNodeInfo; }

   void setMergingStatus(int status){ mergingStatus = status; }
   int getMergingStatus(){ return mergingStatus; }

};

typedef ParaNode *ParaNodePtr;

}

#endif // __PARA_NODE_H__

