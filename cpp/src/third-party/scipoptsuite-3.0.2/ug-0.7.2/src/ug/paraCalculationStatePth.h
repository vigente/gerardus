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

/**@file    paraCalculationStatePth.h
 * @brief   CalcutationStte object extension for Pthreads communication
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_CALCULATION_STATE_PTH_H__
#define __PARA_CALCULATION_STATE_PTH_H__

#include "paraCalculationState.h"

namespace UG
{

/** Calculation state in a ParaSolver */
class ParaCalculationStatePth : public ParaCalculationState
{
   /** create ParaNode datatype */
   ParaCalculationStatePth* createDatatype();
public:
   /** default constructor */
   ParaCalculationStatePth(){}
   /** constructor */
   ParaCalculationStatePth(
         double inCompTime,                   /**< computation time of this ParaNode */
         double inRootTime,                   /**< computation time of the root node */
         int    inNSolved,                    /**< the number of nodes solved   */
         int    inNSent,                      /**< the number of ParaNodes sent */
         int    inNImprovedIncumbent,         /**< the number of improved solution generated in this ParaSolver */
         int    inTerminationState,           /**< indicate whether if this computation is terminationState or not. 0: no, 1: terminationState */
         int    inNSolvedWithNoPreprocesses,  /**< number of nodes solved when it is solved with no preprocesses */
         int    inNSimplexIterRoot,           /**< number of simplex iteration at root node */
         double inAverageSimplexIter,         /**< average number of simplex iteration except root node */
         int    inNRestarts,                  /**< number of restarts */
         double inMinIisum,                   /**< minimum sum of integer infeasibility */
         double inMaxIisum,                   /**< maximum sum of integer infeasibility */
         int    inMinNii,                     /**< minimum number of integer infeasibility */
         int    inMaxNii                      /**< maximum number of integer infeasibility */
         )
         : ParaCalculationState(inCompTime,inRootTime, inNSolved, inNSent,inNImprovedIncumbent,inTerminationState,inNSolvedWithNoPreprocesses,
               inNSimplexIterRoot, inAverageSimplexIter, inNRestarts, inMinIisum, inMaxIisum, inMinNii, inMaxNii)
   {}
   /** destructor */
   ~ParaCalculationStatePth(){}
   void send(ParaComm *comm, int destination, int tag);
   void receive(ParaComm *comm, int source, int tag);

};

}

#endif // __PARA_CALCULATION_STATE_PTH_H__

