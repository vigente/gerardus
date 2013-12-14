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

/**@file    scipParaInitiator.h
 * @brief   ParaInitiator extension for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_INITIATOR_H__
#define __SCIP_PARA_INITIATOR_H__

#include <string>
#include "ug/paraDef.h"
#include "scipParaComm.h"
#include "ug/paraInitiator.h"
#include "scipUserPlugins.h"
#include "scipDiffParamSet.h"
#include "scipParaInstance.h"
#include "scipParaSolution.h"
#include "scipParaDiffSubproblem.h"
#include "objscip/objscip.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

namespace ParaSCIP
{

/** Initiator class */
class ScipParaInitiator : public UG::ParaInitiator
{
   UG::ParaParamSet     *paraParams;
   ScipParaInstance     *instance;
   ScipParaSolution     *solution;
   ScipDiffParamSet     *scipDiffParamSetRoot;
   ScipDiffParamSet     *scipDiffParamSet;
   SCIP_MESSAGEHDLR     *messagehdlr;
   FILE                 *logfile;
   FILE                 *solutionFile;
   FILE                 *transSolutionFile;
   SCIP                 *scip;
   char                 *probname;
   char                 *settingsNameLC;
   char                 *settingsNameRoot;
   char                 *settingsName;
   char                 *racingSettingsName;
   char                 *logname;
   char                 *isolname;
   char                 *solutionFileName;
   ScipUserPlugins      *userPlugins;
   SCIP_Real            finalDualBound;
   UG::FinalSolverState finalState;
   long long            nSolved;

   bool addRootNodeCuts();
   void outputProblemInfo(int *nNonLinearConsHdlrs);
   bool onlyLinearConsHandler();

public:
   /** constructor */
   ScipParaInitiator(
         UG::ParaComm *inComm
         )
         :  ParaInitiator(inComm), paraParams(0), instance(0), solution(0), scipDiffParamSetRoot(0), scipDiffParamSet(0), messagehdlr(0), logfile(0),
            solutionFile(0), transSolutionFile(0), scip(0), probname(0), settingsNameLC(0), settingsNameRoot(0), settingsName(0), racingSettingsName(0),
            logname(0), isolname(0), solutionFileName(0), userPlugins(0), finalDualBound(-DBL_MAX), finalState(UG::Aborted), nSolved(0)
   {
   }

   /** destructor */
   ~ScipParaInitiator(
         )
   {

      if( instance ) delete instance;
      if( solution ) delete solution;
      if( scipDiffParamSetRoot ) delete scipDiffParamSetRoot;
      if( scipDiffParamSet ) delete scipDiffParamSet;
      if( userPlugins ) delete userPlugins;

      // message handler is mangaed within scip. It is freed at SCIPfree
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
      if( messagehdlr )
      {
         SCIP_CALL_ABORT( SCIPsetDefaultMessagehdlr() );
         SCIP_CALL_ABORT( SCIPfreeObjMessagehdlr(&messagehdlr) );
      }
#endif
      /******************
       * Close files *
       ******************/
      fclose(solutionFile);
      if( transSolutionFile )
      {
         fclose(transSolutionFile);
      }

      /********************
       * Deinitialization *
       ********************/
      if( !paraParams->getBoolParamValue(UG::Quiet) )
      {
         SCIP_CALL_ABORT( SCIPprintStatistics(scip, NULL) );    // output statistics (only for problem info)
      }
      if( scip )
      {
    	  SCIP_CALL_ABORT( SCIPfree(&scip) );
      }

      if( logfile != NULL )
         fclose(logfile);

      BMScheckEmptyMemory();
   }

   /** init function */
   int init(
         UG::ParaParamSet *paraParams,
         int          argc,
         char**       argv
         );



   /** get instance */
   UG::ParaInstance *getParaInstance(
         )
   {
      return instance;
   }

   /** make DiffSubproblem object for root node */
   ScipParaDiffSubproblem *makeRootNodeDiffSubproblem(
         )
   {
      return 0;
   }

   /** try to set incumbent solution */
   bool tryToSetIncumbentSolution(UG::ParaSolution *sol);

   /** send solver initialization message */
   void sendSolverInitializationMessage();

   /** generate racing ramp-up parameter sets */
   void generateRacingRampUpParameterSets(int nParamSets, UG::ParaRacingRampUpParamSet **racingRampUpParamSets);

   UG::ParaSolution *getGlobalBestIncumbentSolution()
   {
      return solution;
   }

   /** convert an internal value to external value */
   double convertToExternalValue(
         double internalValue
         )
   {
      return instance->convertToExternalValue(internalValue);
   }

   /** convert an external value to internal value */
   double convertToInternalValue(
         double externalValue
         )
   {
      return instance->convertToInternalValue(externalValue);
   }

   /** get solution file name */
   char *getSolutionFileName(
         )
   {
      return solutionFileName;
   }

   /** get gap */
   double getGap(double dualBoundValue);

   /** get epsilon */
   double getEpsilon();

   /** write solution */
   void writeSolution(const std::string& message);

   /** write ParaInstance */
   void writeParaInstance(const std::string& filename);

   /** write solver runtime parameters */
   void writeSolverParameters(std::ostream *os);

   /** write checkpoint solution */
   void writeCheckpointSolution(const std::string& filename);

   /** read solution from checkpoint file */
   double readSolutionFromCheckpointFile(char *afterCheckpointingSolutionFileName);

   /** get solving status string */
   std::string getStatus();

   /** print solver version **/
   void printSolverVersion(std::ostream *os);   /**< output file (or NULL for standard output) */

   /** check if feasilbe soltuion exists or not */
   bool isFeasibleSolution()
   {
      return ( SCIPgetBestSol(scip) != NULL );
   }

   /** set initial stat on initiator */
   void accumulateInitialStat(UG::ParaInitialStat *initialStat);

   /** set initial stat on DiffSubproblem */
   void setInitialStatOnDiffSubproblem(int minDepth, int maxDepth, UG::ParaDiffSubproblem *diffSubproblem);

   /** set final solver status */
   void setFinalSolverStatus(UG::FinalSolverState status);

   /** set number of nodes solved */
   void setNumberOfNodesSolved(long long n);

   /** set final dual bound  */
   void setDualBound(double bound);

   /** output solution status */
   void outputFinalSolverStatistics(std::ostream *os, double time);

   /** set user plugins */
   void setUserPlugins(ScipUserPlugins *inUi) { userPlugins = inUi; }

   /** include user plugins */
   void includeUserPlugins(SCIP *inScip)
   {
      if( userPlugins )
      {
         (*userPlugins)(inScip);
      }
   }

};

typedef ScipParaInitiator *ScipParaInitiatorPtr;

}

#endif // __SCIP_PARA_INITIATOR_H__

