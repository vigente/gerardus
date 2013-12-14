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

/**@file    paraMergeNodesStructs.h
 * @brief   Structs used for merging nodes.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_MERGE_NODES_STRUCTS_H__
#define __PARA_MERGE_NODES_STRUCTS_H__

namespace UG
{

class ParaNode;
class ParaDiffSubproblem;

typedef struct ParaFixedValue_            ParaFixedValue;
typedef struct ParaMergeNodeInfo_         ParaMergeNodeInfo;
typedef struct ParaFixedVariable_         ParaFixedVariable;
typedef struct ParaSortedVariable_        ParaSortedVariable;
typedef struct ParaFixedValue_ *          ParaFixedValuePtr;
typedef struct ParaMergedNodeListElement_ ParaMergedNodeListElement;


struct ParaFixedValue_ {
   double            value;            /**< value for a fixed variable */
   ParaFixedVariable *head;            /**< point the head of the ParaFixedVariable */
   ParaFixedVariable *tail;            /**< point the tail of the ParaFixedVarialbe */
   ParaFixedValue    *next;            /**< point next ParaFixedValue struct */
};

struct ParaMergeNodeInfo_ {
   enum {
      PARA_MERGING,
      PARA_MERGED_RPRESENTATIVE,
      PARA_MERGE_CHECKING_TO_OTHER_NODE,
      PARA_MERGED_TO_OTHER_NODE,
      PARA_CANNOT_MERGE,
      PARA_DELETED
   }   status;                        /**< status of this ParaMargeNodeInfo */
   int nSameValueVariables;           /**< the number of fixed values which are the same as those of the merged node
                                              - This value < 0 means that this node is not merging target    */
   int nMergedNodes;                  /**< the number of merged nodes with this node.
                                              - This value > 0 : head
                                              - This value = 0 : merging to the other node
                                              - This value < 0 : no merging node  */
   int keyIndex;                      /**< The fixedVar of this index can reach all merging nodes */
   int nFixedVariables;               /**< the number of fixed variables */
   ParaFixedVariable *fixedVariables; /**< array of fixed variable info */
   ParaMergeNodeInfo *mergedTo;       /**< pointer to merge node info to which this node is merged */
   ParaNode *paraNode;                /**< ParaNode corresponding to this ParaMergeModeInfo */
   ParaDiffSubproblem *origDiffSubproblem;   /**< original DiffSubproblem */
   ParaDiffSubproblem *mergedDiffSubproblem; /**< merged DiffSubproblem, in case this node is merged and this is the head */
   ParaMergeNodeInfo *next;           /**< pointer to the next ParaMergeNodeInfo */
};

struct ParaFixedVariable_ {
   int    nSameValue;                /**< the number of same value fixed variables in the following nodes */
   int    index;                     /**< index of the variable among all solvers                 */
   double value;                     /**< fixed value                                             */
   ParaMergeNodeInfo *mnode;         /**< pointer to merge node info struct to which this info is belonging */
   ParaFixedVariable *next;          /**< pointer to the next node which has the same fixed value */        
   ParaFixedVariable *prev;          /**< pointer to the previous node which has the same fixed value */  
};

struct ParaSortedVariable_ {
   int   idxInFixedVariabes;          /**< index in the fixedVariables array */
   ParaFixedVariable *fixedVariable;  /**< pointer to the fixedVariable */
};

struct ParaMergedNodeListElement_ {
   ParaNode  *node;
   ParaMergedNodeListElement *next;
};

};

#endif // __PARA_MERGE_NODES_STRUCTS_H__
