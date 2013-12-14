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

/**@file    paraInitiator.h
 * @brief   Base class of initiator that maintains original problem and incumbent solution.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_INITIATOR_H__
#define __PARA_INITIATOR_H__

#include <string>
#include "paraComm.h"
#include "paraParamSet.h"
#include "paraInstance.h"
#include "paraDiffSubproblem.h"
#include "paraSolution.h"
#include "paraNode.h"
#include "paraRacingRampUpParamSet.h"
#include "paraInitialStat.h"

namespace UG
{

enum FinalSolverState {
   InitialNodesGenerated,
   Aborted,
   HardTimeLimitIsReached,
   ComputingWasInterrupted,
   ProblemWasSolved
};

/** Initiator class */
class ParaInitiator
{
protected:
   ParaComm   *paraComm;
   char       *prefixWarm;            /**< prefix of warm start files */
   igzstream  checkpointNodesStream;  /**< stream for checkpoint nodes file */
   bool       solvedAtInit;           /**< solved at init */
public:
   /** constructor */
   ParaInitiator(
         ParaComm *inComm
         ) : paraComm(inComm), prefixWarm(0), solvedAtInit(false)
   {
   }

   /** destructor */
   virtual ~ParaInitiator(
         )
   {
   }

   /** check if warm started or not */
   bool isWarmStarted() { return prefixWarm != 0; }

   /** read ParaNode form checkpoint file */
   ParaNode *readParaNodeFromCheckpointFile(
         bool onlyBoundChanges
         )
   {
      ParaNode *paraNode = paraComm->createParaNode();
      if( paraNode->read(paraComm, checkpointNodesStream, onlyBoundChanges) )
         return paraNode;
      else
      {
         delete paraNode;
         checkpointNodesStream.close();
         return 0;
      }
   }

   /** get prefix of warm start files */
   const char *getPrefixWarm(){ return prefixWarm; }

   /** check if problem is solved at init or not */
   bool isSolvedAtInit(){ return solvedAtInit; }

   virtual int init(ParaParamSet *params, int argc, char **argv) = 0;

   /** get instance */
   virtual ParaInstance *getParaInstance() = 0;

   /** make DiffSubproblem object for root node */
   virtual ParaDiffSubproblem *makeRootNodeDiffSubproblem() = 0;

   virtual void sendSolverInitializationMessage() = 0;

   // virtual void sendRacingRampUpPhaseParameters() = 0;
   // virtual int getWinnerRank() = 0;
   virtual void generateRacingRampUpParameterSets(int nParamSets, ParaRacingRampUpParamSet **racingRampUpParamSets) = 0;

   /** convert an internal value to external value */
   virtual double convertToExternalValue(double internalValue) = 0;

   virtual ParaSolution *getGlobalBestIncumbentSolution() = 0;

   /** try to set incumbent solution */
   virtual bool tryToSetIncumbentSolution(ParaSolution *sol) = 0;

   /** get gap */
   virtual double getGap(double dualBoundValue) = 0;

   /** get epsilon */
   virtual double getEpsilon() = 0;

   /** write solution */
   virtual void writeSolution(const std::string& message) = 0;

   /** write ParaInstance */
   virtual void writeParaInstance(const std::string& filename) = 0;

   /** write checkpoint solution */
   virtual void writeCheckpointSolution(const std::string& filename) = 0;

   /** read solution from checkpoint file */
   virtual double readSolutionFromCheckpointFile(char *afterCheckpointingSolutionFileName) = 0;

   /** write solver runtime parameters */
   virtual void writeSolverParameters(std::ostream *os) = 0;

   /** set final solver status */
   virtual void setFinalSolverStatus(FinalSolverState status) = 0;

   /** set number of nodes solved */
   virtual void setNumberOfNodesSolved(long long n) = 0;

   /** set final dual bound  */
   virtual void setDualBound(double bound) = 0;

   /** output solution status */
   virtual void outputFinalSolverStatistics(std::ostream *os, double time) = 0;

   /** get solving status string */
   virtual std::string getStatus() = 0;

   /** print solver version **/
   virtual void printSolverVersion(std::ostream *os) = 0;    /**< output file (or NULL for standard output) */

   /** check if feasilbe soltuion exists or not */
   virtual bool isFeasibleSolution() = 0;

   /** set initial stat on initiator */
   virtual void accumulateInitialStat(ParaInitialStat *initialStat) = 0;

   /** set initial stat on DiffSubproblem */
   virtual void setInitialStatOnDiffSubproblem(int minDepth, int maxDepth, ParaDiffSubproblem *diffSubproblem) = 0;

};

typedef ParaInitiator *ParaInitiatorPtr;

}

#endif // __PARA_INITIATOR_HPP__
