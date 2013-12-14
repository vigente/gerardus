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

/**@file    paraLoadCoordinatorTerminationState.h
 * @brief   Load coordinator termination state.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_LOADCOORDINATOR_TERMINATION_STATE_H__
#define __PARA_LOADCOORDINATOR_TERMINATION_STATE_H__

#include <string>
#include <cfloat>
#include "ug/paraComm.h"
#include "gzstream.h"

namespace UG
{

/** Calculation state in a ParaLoadCoordinator */
class ParaLoadCoordinatorTerminationState
{
public:

   bool         isCheckpointState;          /**< This state is at checkpoint or not */
   int          rank;                       /**< rank of this ParaLoadCoordinator */
   /** Counters related to this ParaLoadCoordinator */
   int          nWarmStart;                 /**< number of warm starts */
   long long    nSent;
   long long    nSentBackImmediately;
   int          nSentBackImmediatelyAnotherNode;
   long long    nReceived;
   int          nDeletedInLc;
   int          nDeletedByMerging;
   int          nFailedToSendBack;
   int          nFailedToSendBackAnotherNode;
   int          nMaxUsageOfNodePool;
   int          nNodesInNodePool;
   long long    nNodesLeftInAllSolvers;
   /** current dual bound value */
   double       globalBestDualBoundValue;
   double       externalGlobalBestDualBoundValue;
   /** times of this solver */
   double       idleTime;                   /**< idle time of this LoadCoordinator */
   double       runningTime;                /**< this ParaLoadCoordinator running time */


   /** default constructor */
   ParaLoadCoordinatorTerminationState(
         )
         : isCheckpointState(true), rank(0), nWarmStart(0), nSent(0), nSentBackImmediately(0),
         nSentBackImmediatelyAnotherNode(0), nReceived(0),
         nDeletedInLc(0), nDeletedByMerging(0), nFailedToSendBack(0), nFailedToSendBackAnotherNode(0),
         nMaxUsageOfNodePool(0), nNodesInNodePool(0), nNodesLeftInAllSolvers(0),
         globalBestDualBoundValue(-DBL_MAX), externalGlobalBestDualBoundValue(-DBL_MAX),
         idleTime(0.0), runningTime(0.0)
   {
   }
   /** destructor */
   virtual ~ParaLoadCoordinatorTerminationState(
	        )
   {
   }

   /** stringfy ParaCalculationState */
   std::string toString();

   /** write to checkpoint file */
   void write(ogzstream &out);

   /** read from checkpoint file */
   bool read(ParaComm *comm, igzstream &in);

};

}

#endif // __PARA_LOADCOORDINATOR_TERMINATION_STATE_H__

