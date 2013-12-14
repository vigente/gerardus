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

/**@file    paraSolverState.h
 * @brief   This class has solver state to be transferred.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_SOLVER_STATE_H__
#define __PARA_SOLVER_STATE_H__

#include <cfloat>
#include "paraComm.h"

namespace UG
{

/** ParaSolver state object for notification message */
class ParaSolverState
{
protected:
   int    racingStage;          /**< if this value is 1, solver is in racing stage */
   unsigned int notificationId; /**< id for this notification */
   int    lcId;                 /**< lc id of current ParaNode */
   int    globalSubtreeIdInLc;  /**< global subtree id of current ParaNode */
   long long nNodesSolved;      /**< number of nodes solved */
   int    nNodesLeft;           /**< number of remaining nodes  */
   double bestDualBoundValue;   /**< best dual bound value in that of remaining nodes */
   double globalBestPrimalBoundValue; /**< global best primal bound value */
   double detTime;             /**< deterministic time, -1: should be non-deterministic */
public:
   /** default constructor */
   ParaSolverState(
         )
         : racingStage(0), notificationId(0), lcId(-1), globalSubtreeIdInLc(-1), nNodesLeft(-1),  bestDualBoundValue(0.0),
           globalBestPrimalBoundValue(DBL_MAX), detTime(-1.0)
   {
   }

   /** constructor */
   ParaSolverState(
         int inRacingStage,
         unsigned int inNotificationId,
         int inLcId,
         int inGlobalSubtreeId,
         long long inNodesSolved,
         int inNodesLeft,
         double inBestDualBoundValue,
         double inGlobalBestPrimalBoundValue,
         double inDetTime
         ) : racingStage(inRacingStage), notificationId(inNotificationId), lcId(inLcId), globalSubtreeIdInLc(inGlobalSubtreeId),
         nNodesSolved(inNodesSolved), nNodesLeft(inNodesLeft), bestDualBoundValue(inBestDualBoundValue),
         globalBestPrimalBoundValue(inGlobalBestPrimalBoundValue), detTime(inDetTime)
   {
   }

   /** destractor */
   virtual ~ParaSolverState(
         )
   {
   }

   /** getter of isRacingStage */
   bool isRacingStage(
         )
   {
      return (racingStage == 1);
   }

   /** getter of notificationId */
   unsigned int getNotificaionId(
         )
   {
      return notificationId;
   }

   /** getter of lcId */
   int getLcId(
         )
   {
      return lcId;
   }

   /** getter of global subtree Id */
   int getGlobalSubtreeId(
         )
   {
      return globalSubtreeIdInLc;
   }

   /** gettter of bestDualBoundValue */
   double getSolverLocalBestDualBoundValue(
         )
   {
      return bestDualBoundValue;
   }

   /** get global best primal bound value that the notification solver has */
   double getGlobalBestPrimalBoundValue(
         )
   {
      return globalBestPrimalBoundValue;
   }

   /** getter of nodes solved */
   long long getNNodesSolved(
         )
   {
      return nNodesSolved;
   }


   /** getter of nodes left */
   int getNNodesLeft(
         )
   {
      return nNodesLeft;
   }

   /** getter of deterministic time */
   double getDeterministicTime(
         )
   {
      return detTime;
   }

   /** stringfy ParaSolverState */
   std::string toString(
         )
   {
      std::ostringstream s;
      s << "racingStage = " << racingStage << ", notificationId = " << notificationId << ": ";
      s << "[" << lcId << ":" << globalSubtreeIdInLc << "]"
      << " Best dual bound value = " << bestDualBoundValue
      << " number of nodes solved = " << nNodesSolved
      << ", number of nodes left = " << nNodesLeft;
      return s.str();
   }

   virtual void send(ParaComm *comm, int destination, int tag) = 0;
   virtual void receive(ParaComm *comm, int source, int tag) = 0;
};

}

#endif // __PARA_SOLVER_STATE_H__
