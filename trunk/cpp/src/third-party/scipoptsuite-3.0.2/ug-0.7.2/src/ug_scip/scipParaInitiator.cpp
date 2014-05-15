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

/**@file    scipParaInitiator.cpp
 * @brief   ParaInitiator extension for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cctype>
#include <sstream>
#include "scipParaInstance.h"
#include "scipParaInitiator.h"
#include "scipParaObjMessageHdlr.h"
#include "scipParaInitialStat.h"

using namespace UG;
using namespace ParaSCIP;

bool
ScipParaInitiator::addRootNodeCuts(
      )
{
   SCIP_Longint originalLimitsNodes;
   SCIP_CALL_ABORT( SCIPgetLongintParam(scip, "limits/nodes", &originalLimitsNodes) );
   SCIP_CALL_ABORT( SCIPsetLongintParam(scip, "limits/nodes", 1) );
   if( scipDiffParamSetRoot ) scipDiffParamSetRoot->setParametersInScip(scip);
   SCIP_RETCODE ret = SCIPsolve(scip);
   if( ret != SCIP_OKAY )
   {
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
      SCIPprintError(ret, NULL);
#else
      SCIPprintError(ret);
#endif
      abort();
   }
   // Then, solver status should be checked
   SCIP_STATUS status = SCIPgetStatus(scip);
   if( status == SCIP_STATUS_OPTIMAL )   // when sub-MIP is solved at root node, the solution may not be saved
   {
      return false;
   }
   else
   {
      if( status == SCIP_STATUS_MEMLIMIT  )
      {
         std::cout << "Warning: SCIP was interrupted because the memory limit was reached" << std::endl;
         return false;
      }
   }

   SCIP_CUT** cuts;
   int ncuts;
   int ncutsadded;

   ncutsadded = 0;
   cuts = SCIPgetPoolCuts(scip);
   ncuts = SCIPgetNPoolCuts(scip);
   for( int c = 0; c < ncuts; ++c )
   {
      SCIP_ROW* row;

      row = SCIPcutGetRow(cuts[c]);
      assert(!SCIProwIsLocal(row));
      assert(!SCIProwIsModifiable(row));
      if( SCIPcutGetAge(cuts[c]) == 0 && SCIProwIsInLP(row) )
      {
         char name[SCIP_MAXSTRLEN];
         SCIP_CONS* cons;
         SCIP_COL** cols;
         SCIP_VAR** vars;
         int ncols;
         int i;

         /* create a linear constraint out of the cut */
         cols = SCIProwGetCols(row);
         ncols = SCIProwGetNNonz(row);

         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vars, ncols) );
         for( i = 0; i < ncols; ++i )
            vars[i] = SCIPcolGetVar(cols[i]);

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d", SCIProwGetName(row), SCIPgetNRuns(scip));
         SCIP_CALL_ABORT( SCIPcreateConsLinear(scip, &cons, name, ncols, vars, SCIProwGetVals(row),
               SCIProwGetLhs(row) - SCIProwGetConstant(row), SCIProwGetRhs(row) - SCIProwGetConstant(row),
               TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
         SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
         SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );

         SCIPfreeBufferArray(scip, &vars);

         ncutsadded++;
      }
   }
   SCIP_CALL_ABORT( SCIPsetLongintParam(scip, "limits/nodes", originalLimitsNodes) );
   return true;
}

/** init function */
int
ScipParaInitiator::init(
      ParaParamSet *inParaParams,
      int     argc,
      char**  argv
      )
{
   int i;
   bool quiet = false;
   paraParams = inParaParams;

   probname = argv[2];

   /*********
    * Setup *
    *********/
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   /** user include plugins */
   includeUserPlugins(scip);
   finalDualBound = SCIPinfinity(scip);
   printSolverVersion(NULL);
   /********************
    * Parse parameters *
    ********************/
   if( std::string(paraParams->getStringParamValue(SolverSettingsForInitialPresolving)) != "" )
   {
      settingsNameLC = const_cast<char*> (paraParams->getStringParamValue(SolverSettingsForInitialPresolving));
   }
   if( std::string(paraParams->getStringParamValue(SolverSettingsAtRootNode)) != "" )
   {
      settingsNameRoot = const_cast<char*> (paraParams->getStringParamValue(SolverSettingsAtRootNode));
   }
   if( std::string(paraParams->getStringParamValue(SolverSettingsExceptRootNode)) != "" )
   {
      settingsName = const_cast<char*> (paraParams->getStringParamValue(SolverSettingsExceptRootNode));
   }
   if( std::string(paraParams->getStringParamValue(SolverSettingsAtRacing)) != "" )
   {
      racingSettingsName = const_cast<char*> (paraParams->getStringParamValue(SolverSettingsAtRacing));
   }
   for( i = 3; i < argc; ++i )   /** the first argument is runtime parameter file for ParaSCIP */
                                 /** the second argument is problem file name */
   {
      if( strcmp(argv[i], "-l") == 0 )
      {
         i++;
         if( i < argc )
            logname = argv[i];
         else
         {
            std::cerr << "missing log filename after parameter '-l'" << std::endl;
            exit(1);
         }
      }
      else if( strcmp(argv[i], "-q") == 0 )
         quiet = true;
      else if( strcmp(argv[i], "-s") == 0 )
      {
         i++;
         if( i < argc )
            settingsName = argv[i];
         else
         {
            std::cerr << "missing settings filename after parameter '-s'" << std::endl;
            exit(1);
         }
      }
      else if( strcmp(argv[i], "-sr") == 0 )
      {
         i++;
         if( i < argc )
            settingsNameRoot = argv[i];
         else
         {
            std::cerr << "missing settings filename after parameter '-sr'" << std::endl;
            exit(1);
         }
      }
      else if( strcmp(argv[i], "-sl") == 0 )
      {
         i++;
         if( i < argc )
            settingsNameLC = argv[i];
         else
         {
            std::cerr << "missing settings filename after parameter '-sl'" << std::endl;
            exit(1);
         }
      }
      else if( strcmp(argv[i], "-w") == 0)
      {
         i++;
         if( i < argc )
         {
            prefixWarm = argv[i];
            char nodesFileName[256];
            sprintf(nodesFileName,"%s_nodes_LC0.gz",prefixWarm);
            checkpointNodesStream.open(nodesFileName, std::ios::in | std::ios::binary);
            if( !checkpointNodesStream.good() ){
                std::cerr << "ERROR: Opening file `" << nodesFileName << "' failed.\n";
                exit(1);
            }
         }
         else
         {
            std::cerr << "missing settings filename after parameter '-w'" << std::endl;
            exit(1);
         }
      }
      else if( strcmp(argv[i], "-racing") == 0 )
      {
         i++;
         if( i < argc )
         {
            racingSettingsName = argv[i];
         }
         else
         {
            std::cerr << "missing settings filename after parameter '-racing'" << std::endl;
            exit(1);
         }
      }
      else if ( strcmp(argv[i], "-isol") == 0 )
      {
         i++;
         if( i < argc )
         {
            isolname = argv[i];
          }
          else
          {
             std::cerr << "missing settings filename after parameter '-isol'" << std::endl;
             exit(1);
          }
      }
      else if ( strcmp(argv[i], "-sth") == 0 )
      {
         i++;  // just omit this parameter and the following number.
      }
      else if ( strcmp(argv[i], "-fsol" ) == 0 )
      {
         i++;
         if( i < argc )
         {
            solutionFileName = argv[i];
          }
          else
          {
             std::cerr << "missing solution filename after parameter '-fsol'" << std::endl;
             exit(1);
          }
      }
      else
      {
         THROW_LOGICAL_ERROR3("invalid parameter <", argv[i], ">");
      }
   }

   /*************************************************
    * set quiet message handler, if it is necessary *
    *************************************************/
   messagehdlr = NULL;
   if( paraParams->getBoolParamValue(Quiet) )
   {
      SCIP_CALL_ABORT( SCIPcreateObjMessagehdlr(&messagehdlr, new ScipParaObjMessageHdlr(paraComm, NULL, TRUE, FALSE), TRUE) );
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
      SCIP_CALL_ABORT( SCIPsetMessagehdlr(messagehdlr) );
#else
      SCIP_CALL_ABORT( SCIPsetMessagehdlr(scip, messagehdlr) );
      SCIP_CALL_ABORT( SCIPmessagehdlrRelease(&messagehdlr));
#endif
   }
   else
   {
      if( logname != NULL || quiet  )
      {
         if( logname != NULL )
         {
            std::ostringstream os;
            os << logname << paraComm->getRank();
            logfile = fopen(os.str().c_str(), "a"); // append to log file */
            if( logfile == NULL )
            {
               THROW_LOGICAL_ERROR3("cannot open log file <", logname, "> for writing");
            }
         }
         SCIP_CALL_ABORT( SCIPcreateObjMessagehdlr(&messagehdlr, new ScipParaObjMessageHdlr(paraComm, logfile, quiet, FALSE), TRUE) );
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
         SCIP_CALL_ABORT( SCIPsetMessagehdlr(messagehdlr) );
#else
         SCIP_CALL_ABORT( SCIPsetMessagehdlr(scip, messagehdlr) );
         SCIP_CALL_ABORT( SCIPmessagehdlrRelease(&messagehdlr));
#endif
      }
   }

   if( probname != NULL )
   {
      /***********************
       * Version information *
       ***********************/
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
      SCIPprintVersion(NULL);
#else
      SCIPprintVersion(scip, NULL);
#endif
      SCIPinfoMessage(scip, NULL, "\n");
      /*****************
       * Load settings *
       *****************/
      DEF_SCIP_PARA_COMM( scipParaComm, paraComm );
      if( settingsName != NULL )
      {
         SCIP_CALL( SCIPreadParams(scip, settingsName) );
         scipDiffParamSet =  scipParaComm->createScipDiffParamSet(scip);
         SCIP_CALL( SCIPresetParams(scip) );
      }
      else
      {
         scipDiffParamSet = scipParaComm->createScipDiffParamSet(scip);
      }
      if( settingsNameRoot != NULL )
      {
         SCIP_CALL( SCIPreadParams(scip, settingsNameRoot) );
         scipDiffParamSetRoot = scipParaComm->createScipDiffParamSet(scip);
         SCIP_CALL( SCIPresetParams(scip) );
      }
      else
      {
         scipDiffParamSetRoot = scipParaComm->createScipDiffParamSet(scip);
      }
      if( settingsNameLC != NULL )
      {
         SCIP_CALL( SCIPreadParams(scip, settingsNameLC) );
      }

      /**************
       * Start SCIP *
       **************/
      // Problem Creation
      int retcode = SCIPreadProb(scip, probname, NULL);
      if( retcode != SCIP_OKAY )
      {
         std::cout << "error reading file <" << probname << ">" << std::endl;
         SCIP_CALL( SCIPfreeProb(scip) );
         exit(1);
      }

      /* read initail solution, if it is specified */
      if( isolname )
      {
         SCIP_CALL_ABORT( SCIPtransformProb(scip));
         // NOTE:
         // When CPLEX license file cannot find, SCIPtransformProb(scip) may fail
         SCIP_CALL_ABORT( SCIPreadSol(scip, isolname) );
      }

      /* change problem name */
      char *probNameFromFileName;
      SCIPsplitFilename(probname, NULL, &probNameFromFileName, NULL, NULL);
      SCIP_CALL_ABORT( SCIPsetProbName(scip, probNameFromFileName));

      /* presolve problem */
      if( paraParams->getBoolParamValue(NoPreprocessingInLC) )
      {
         SCIP_CALL_ABORT( SCIPsetIntParam(scip, "presolving/maxrounds", 0));
         std::cout << "No LC presolving is specified." << std::endl;
      }
      else
      {
         if( settingsNameLC )
         {
            SCIP_CALL_ABORT( SCIPreadParams(scip, settingsNameLC) );
            std::cout << "LC presolving settings file is specified." << std::endl;
         }
         else
         {
            // SCIP_CALL_ABORT( SCIPsetPresolving(scip, SCIP_PARAMSETTING_FAST, TRUE) );
            std::cout << "Default LC presolving (default)." << std::endl;
         }
      }
      // SCIP_CALL_ABORT( SCIPsetIntParam(scip, "constraints/quadratic/replacebinaryprod", 0));
      // SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "constraints/nonlinear/reformulate", FALSE));

#ifdef _COMM_PTH
      SCIP_CALL( SCIPsetBoolParam(scip, "misc/catchctrlc", FALSE) );
#endif

#ifndef _COMM_PTH
      SCIP_Bool originalUpgradeKnapsack;
      SCIP_Bool originalUpgradeLogicor;
      SCIP_Bool originalUpgradeSetppc;
      SCIP_Bool originalUpgradeVarbound;
      bool onlyLinearConss = onlyLinearConsHandler();
      if( onlyLinearConss )
      {
         SCIP_CALL_ABORT( SCIPgetBoolParam(scip, "constraints/linear/upgrade/knapsack", &originalUpgradeKnapsack));
         SCIP_CALL_ABORT( SCIPgetBoolParam(scip, "constraints/linear/upgrade/logicor", &originalUpgradeLogicor));
         SCIP_CALL_ABORT( SCIPgetBoolParam(scip, "constraints/linear/upgrade/setppc", &originalUpgradeSetppc));
         SCIP_CALL_ABORT( SCIPgetBoolParam(scip, "constraints/linear/upgrade/varbound", &originalUpgradeVarbound));
         SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "constraints/linear/upgrade/knapsack", FALSE));
         SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "constraints/linear/upgrade/logicor", FALSE));
         SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "constraints/linear/upgrade/setppc", FALSE));
         SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "constraints/linear/upgrade/varbound", FALSE));
      }
#endif

      SCIP_CALL( SCIPpresolve(scip) );
      SCIP_STATUS scipStatus = SCIPgetStatus(scip);

      // output presolved instance information
      int nNonLinearConsHdlrs = 0;
      outputProblemInfo(&nNonLinearConsHdlrs);
#ifndef _COMM_PTH
      PARA_COMM_CALL(
            paraComm->bcast( &nNonLinearConsHdlrs, 1, ParaINT, 0 )
      );
      if( nNonLinearConsHdlrs > 0 )
      {
         paraParams->setIntParamValue(InstanceTransferMethod,2);
      }
#endif

      if( scipStatus == SCIP_STATUS_OPTIMAL )   // when sub-MIP is solved at root node, the solution may not be saved
      {
         solvedAtInit = true;
         setFinalSolverStatus(ProblemWasSolved);
         finalDualBound = SCIPgetDualbound(scip);
         std::cout << "=== solved at Init ===" << std::endl;
      }
      else
      {
         /* adding root node cuts, if necessary */
         if( paraParams->getBoolParamValue(UseRootNodeCuts) )
         {
            if( !addRootNodeCuts() )
            {
               solvedAtInit = true;
               setFinalSolverStatus(ProblemWasSolved);
               finalDualBound = SCIPgetDualbound(scip);
               std::cout << "=== solved at Init ===" << std::endl;
            }
         }
      }

      if( SCIPgetNVars(scip) == 0 )  // all variables were fixed in presolve
      {
         SCIP_CALL( SCIPsolve(scip) );
         solvedAtInit = true;
         setFinalSolverStatus(ProblemWasSolved);
         finalDualBound = SCIPgetDualbound(scip);
         std::cout << "=== solved at Init ===" << std::endl;
      }

      /** check if feasible solution is found or not. If it was found, then generate paraSolution */
      SCIP_SOL *sol = SCIPgetBestSol(scip);
      if( sol )
      {
         int nVars = SCIPgetNVars(scip);
         SCIP_VAR **vars = SCIPgetVars(scip);
         SCIP_Real *vals = new SCIP_Real[nVars];
         SCIP_CALL_ABORT( SCIPgetSolVals(scip, sol, nVars, vars, vals) );
         solution = scipParaComm->createScipParaSolution(
                        SCIPgetSolTransObj(scip, sol),  // Only this value may be used
                        nVars,
                        vars,
                        vals
                        );
         delete [] vals;
      }

      // instance = new PARA_INSTANCE_TYPE(scip, paraParams->getIntParamValue(InstanceTransferMethod));
      instance = scipParaComm->createScipParaInstance(scip, paraParams->getIntParamValue(InstanceTransferMethod));

#ifndef _COMM_PTH
      if( onlyLinearConss )
      {
         // restore original parameters
         SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "constraints/linear/upgrade/knapsack", originalUpgradeKnapsack));
         SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "constraints/linear/upgrade/logicor", originalUpgradeLogicor));
         SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "constraints/linear/upgrade/setppc", originalUpgradeSetppc));
         SCIP_CALL_ABORT( SCIPsetBoolParam(scip, "constraints/linear/upgrade/varbound", originalUpgradeVarbound));
      }
#endif

      std::ostringstream os;
      if( solutionFileName )
      {
         os << solutionFileName;
      }
      else
      {
         os << paraParams->getStringParamValue(SolutionFilePath);
         os << instance->getProbName() << ".sol";
      }
      solutionFile = fopen(os.str().c_str(), "a");  // if solution file exists, append
      if( solutionFile == NULL )
      {
         THROW_LOGICAL_ERROR3("cannot open solution file <", os.str(), "> for writing");
      }
      int maxrounds = 0;
      SCIP_CALL_ABORT( SCIPgetIntParam(scip, "presolving/maxrounds", &maxrounds));
      if( !paraParams->getBoolParamValue(Quiet) && maxrounds != 0 )
      {
         os << ".trans";
         transSolutionFile = fopen(os.str().c_str(), "a");  // if trans. solution file exists, append
         if( transSolutionFile == NULL )
         {
            THROW_LOGICAL_ERROR3("cannot open solution file <", os.str(), "> for writing");
         }
      }
   }
   else
   {
      std::cout << std::endl;
      std::cout << "syntax: " << argv[0] << "#solvers ppscip_param_file problem_file_name "
                << "[-l <logfile>] [-q] [-sl <settings>] [-s <settings>] [-sr <root_settings>] [-w <prefix_warm>] [-sth <number>]" << std::endl;
      std::cout << "  -l <logfile>        : copy output into log file" << std::endl;
      std::cout << "  -q                  : suppress screen messages" << std::endl;
      std::cout << "  -sl <settings>      : load parameter settings (.set) file for LC presolving" << std::endl;
      std::cout << "  -s <settings>       : load parameter settings (.set) file for solvers" << std::endl;
      std::cout << "  -sr <root_settings> : load parameter settings (.set) file for root" << std::endl;
      std::cout << "  -w <prefix_warm>    : warm start file prefix ( prefix_warm_nodes.gz and prefix_warm_solution.txt are read )" << std::endl;
      std::cout << "  -sth <number>       : the number of solver threads used(FiberSCIP)" << std::endl;
      THROW_LOGICAL_ERROR1("invalid parameter");
   }

   if( solution )
   {
      if( !paraParams->getBoolParamValue(Quiet) )
      {
         writeSolution("");
      }
   }

   if( solvedAtInit )
   {
      writeSolution("Final Solution");
      return 1;
   }
   else
   {
      return 0;
   }
}

bool
ScipParaInitiator::tryToSetIncumbentSolution(
      ParaSolution *sol
      )
{
   ScipParaSolution *tempSol = dynamic_cast< ScipParaSolution * >(sol);
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   paraComm->lockApp();  /* lock is necessary, if Solver runs as thread */

   /* the solution from the working node should have at least as many variables as we have in the load coordinator scip
    * it may be more if inactive variable had to be copied, i.e.,
    * SCIPgetNVars(scip) is the number of active variables in the load coordinator scip
    * tempSol->getNVars() is the number of original variables in the working node scip (scipParaSolver)
    */
   SCIP_VAR** vars = 0;
   int maxrounds = 0;
   SCIP_CALL_ABORT( SCIPgetIntParam(scip, "presolving/maxrounds", &maxrounds));
   // if( paraParams->getIntParamValue(InstanceTransferMethod) == 2    // original file read
   if( maxrounds == 0 )                                                 // nopreprocessing
   {
      //assert(SCIPgetNVars(scip) <= tempSol->getNVars());
      /* create new solution for the original problem */
      SCIP_CALL_ABORT( SCIPcreateOrigSol(scip, &newsol, 0) );
      vars = SCIPgetOrigVars(scip);
   }
   else
   {
      assert(SCIPgetNVars(scip) <= tempSol->getNVars());
      /* create new solution for the original problem */
      SCIP_CALL_ABORT( SCIPcreateSol(scip, &newsol, 0) );
      vars = SCIPgetVars(scip);
   }

   int i;
   for( i = 0; i < tempSol->getNVars(); i++ )
   {
      /* if index is larger-equal than number of active vars, then it's probably an inactive variable which had to be copied via SCIPcopy
       * so we just ignore its value
       */
      if( tempSol->indexAmongSolvers(i) >= SCIPgetNVars(scip) ) continue;
      if( tempSol->indexAmongSolvers(i) > ( tempSol->getNVars() - 1 ) ) break;
      SCIP_CALL_ABORT( SCIPsetSolVal(scip, newsol, vars[tempSol->indexAmongSolvers(i)], tempSol->getValues()[i]) );
   }
   if( i != tempSol->getNVars() )
   {
      /** the given solution should be generated in original space,
       * therefore the solution values cannot use for ParaSCIP
       */
      SCIP_CALL_ABORT( SCIPfreeSol(scip, &newsol) );
      delete tempSol;
      std::cout << "solution size mismatch! Call Yuji!" << std::endl;
      return false;
   }
   SCIP_Bool success;

#ifndef UG_SCIP_SOL_FEASIBILITY_CHECK_IN_LC

   SCIP_CALL_ABORT( SCIPaddSolFree(scip, &newsol, &success) );

   paraComm->unlockApp();

   if( success )
   {
      if( solution )
      {
         delete solution;
     }
      solution = tempSol;
      if( !paraParams->getBoolParamValue(Quiet) )
      {
         writeSolution("");
      }
      return true;
   }
   else
   {
      delete tempSol;
      return false;
   }

#else

   SCIP_CALL_ABORT( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &success) );

   paraComm->unlockApp();

   if( success )
   {
      SCIP_CALL_ABORT( SCIPfreeSol(scip, &newsol) );
      if( solution )
      {
         delete solution;
      }
      solution = tempSol;
      if( !paraParams->getBoolParamValue(Quiet) )
      {
         writeSolution("");
      }
      return true;
   }
   else
   {
      SCIP_VAR* var = 0;
      for( i = 0; i < tempSol->getNVars(); i++ )
      {
         if( tempSol->indexAmongSolvers(i) >= SCIPgetNVars(scip) ) continue;
         if( tempSol->indexAmongSolvers(i) > ( tempSol->getNVars() - 1 ) ) break;
         var = vars[tempSol->indexAmongSolvers(i)];
         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS ) continue;
         if( tempSol->getValues()[i] > 0.0 )
         {
            tempSol->setValue(i, SCIPfeasFloor(scip,tempSol->getValues()[i]));
         }
         else
         {
            tempSol->setValue(i, SCIPfeasCeil(scip, tempSol->getValues()[i]) );
         }
         SCIP_CALL_ABORT( SCIPsetSolVal(scip, newsol, var, tempSol->getValues()[i]) );
      }

      SCIP_CALL_ABORT( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, &success) );
      if( success )
      {
         tempSol->setObjectiveFuntionValue(convertToInternalValue(SCIPsolGetOrigObj(newsol)));
         SCIP_CALL_ABORT( SCIPfreeSol(scip, &newsol) );
         if( solution )
         {
            delete solution;
         }
         solution = tempSol;
         if( !paraParams->getBoolParamValue(Quiet) )
         {
            writeSolution("");
         }
         return true;
      }
      else
      {
         fprintf(solutionFile, "*** Rejected Solution ***\n");
         SCIP_CALL_ABORT(SCIPprintSol(scip, newsol, solutionFile, FALSE));
         if( transSolutionFile )
         {
            fprintf(transSolutionFile, "*** Rejected Solution ***\n");
            SCIP_CALL_ABORT(SCIPprintTransSol(scip, newsol, transSolutionFile, FALSE));
         }

         if( SCIPisFeasLT( scip, convertToExternalValue(tempSol->getObjectiveFuntionValue()), SCIPgetPrimalbound(scip) ) )
         {
            paraComm->lockApp();
            std::cout << "Current scip primal value = " << SCIPgetPrimalbound(scip) << std::endl;
            std::cout << "Objective value = " << convertToExternalValue(tempSol->getObjectiveFuntionValue()) << std::endl;
            std::cout << "Initiator did not accept solution!" << std::endl;
            paraComm->unlockApp();
         }
         SCIP_CALL_ABORT( SCIPfreeSol(scip, &newsol) );
         delete tempSol;
         return false;
      }
   }

#endif

}

void
ScipParaInitiator::sendSolverInitializationMessage(
      )
{
   assert(scipDiffParamSetRoot && scipDiffParamSet);
   scipDiffParamSetRoot->bcast(paraComm, 0);
   scipDiffParamSet->bcast(paraComm, 0);
   int warmStarted = 0;
   if( isWarmStarted() )
   {
      warmStarted = 1;
   }
   paraComm->bcast(&warmStarted,1, ParaINT, 0);
   double incumbentValue;
   if( solution )
   {
      incumbentValue = solution->getObjectiveFuntionValue();
   }
   else
   {
      SCIP_SOL *sol = SCIPgetBestSol(scip);
	   if ( sol )
	   {
	      int nVars = SCIPgetNVars(scip);
	      SCIP_VAR **vars = SCIPgetVars(scip);
	      SCIP_Real *vals = new SCIP_Real[nVars];
	      SCIP_CALL_ABORT( SCIPgetSolVals(scip, sol, nVars, vars, vals) );
	      DEF_SCIP_PARA_COMM( scipParaComm, paraComm);
	      solution = scipParaComm->createScipParaSolution(
	            SCIPgetSolTransObj(scip,sol),
	            nVars,
	             vars,
	             vals
	             );
	      delete [] vals;
	      incumbentValue = solution->getObjectiveFuntionValue();
	   }
	   else
	   {
	      incumbentValue = DBL_MAX;
	   }
   }
   paraComm->bcast(&incumbentValue, 1, ParaDOUBLE, 0);

   if( paraParams->getBoolParamValue(NoUpperBoundTransferInRacing) )
   {
      int solutionExists = 0;
      paraComm->bcast(&solutionExists, 1, ParaINT, 0);
   }
   else
   {
      /** if a feasible solution exists, broadcast the solution */
       if( paraParams->getBoolParamValue(DistributeBestPrimalSolution) )
       {
          /* bcast solution if it is necessary */
          int solutionExists = 0;
          if( solution )
          {
             solutionExists = 1;
             paraComm->bcast(&solutionExists, 1, ParaINT, 0);
             solution->bcast(paraComm, 0);
          }
          else
          {
             paraComm->bcast(&solutionExists, 1, ParaINT, 0);
          }
       }
   }
}

/** get gap */
double
ScipParaInitiator::getGap(
      double dualBoundValue
      )
{
   if( !solution ) return SCIPinfinity(scip);
   SCIP_Real primalbound = instance->convertToExternalValue(solution->getObjectiveFuntionValue());
   SCIP_Real dualbound = instance->convertToExternalValue(dualBoundValue);

   if( SCIPisEQ(scip, primalbound, dualbound) )
      return 0.0;
   else if( SCIPisZero(scip, dualbound)
      || SCIPisInfinity(scip, REALABS(primalbound))
      || primalbound * dualbound < 0.0 )
      return SCIPinfinity(scip);
   else
      return REALABS((primalbound - dualbound)/MIN(REALABS(dualbound),REALABS(primalbound)));
}

/** get epsilon */
double
ScipParaInitiator::getEpsilon(
      )
{
   SCIP_Real epsilon;
   SCIP_CALL_ABORT( SCIPgetRealParam(scip, "numerics/epsilon", &epsilon));
   return epsilon;
}

void
ScipParaInitiator::writeSolution(
      const std::string& message
      )
{
   if( message == "Final Solution" )
   {
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
      if( paraParams->getBoolParamValue(Quiet) )
      {
         SCIPmessageSetDefaultHandler();    // If no message handler is set, it cannot write solution,too.
      }
#endif
      fprintf(solutionFile, "[ Final Solution ]\n");
      if( transSolutionFile )
      {
         fprintf(transSolutionFile, "[ Final Solution ]\n");
      }
   }
   else
   {
      fprintf(solutionFile,"%s\n",message.c_str());
      if( transSolutionFile )
      {
         fprintf(transSolutionFile,"%s\n", message.c_str());
      }
   }
   SCIP_SOL* sol = SCIPgetBestSol(scip);
   if( sol )
   {
      SCIP_CALL_ABORT( SCIPprintBestSol( scip, solutionFile, FALSE) );
      if( transSolutionFile )
      {
         if( SCIPsolGetOrigin(sol) != SCIP_SOLORIGIN_ORIGINAL )
         {
            SCIP_CALL_ABORT( SCIPprintBestTransSol( scip, transSolutionFile, FALSE) );
         }
         else
         {
            fprintf(transSolutionFile, "best solution is defined in original space - cannot print it as transformed solution\n");
         }
      }
   }
   else
   {
      fprintf(solutionFile, "No Solution\n");
      if( transSolutionFile )
      {
         fprintf(transSolutionFile, "No Solution\n");
      }
   }
}

void
ScipParaInitiator::writeParaInstance(
      const std::string& filename
      )
{
   FILE *file = fopen(filename.c_str(),"a");
   if( !file )
   {
      std::cout << "file : " << filename << "cannot open." << std::endl;
      exit(1);
   }
   SCIP_CALL_ABORT( SCIPprintTransProblem(scip, file, "lp", FALSE) );
}

/** write solver runtime parameters */
void
ScipParaInitiator::writeSolverParameters(
      std::ostream *os
      )
{
   if( scipDiffParamSetRoot->nDiffParams() == 0 )
   {
      *os << "[ SCIP parameters for root Solver are all default values ]" << std::endl;
   }
   else
   {
      *os << "[ Not default SCIP parameters for root Solver are as follows ]" << std::endl;
      *os << scipDiffParamSetRoot->toString();
   }

   if( scipDiffParamSet->nDiffParams() == 0 )
   {
       *os << "[ SCIP parameters for NOT root Solvers are all default values ]" << std::endl;
   }
   else
   {
      *os << "[ Not default SCIP parameters for NOT root Solvers are as follows ]" << std::endl;
      *os << scipDiffParamSet->toString();
   }

}

/** write checkpoint solution */
void
ScipParaInitiator::writeCheckpointSolution(
      const std::string& filename
      )
{
   ogzstream checkpointSolutionStream;
   checkpointSolutionStream.open(filename.c_str(), std::ios::out | std::ios::binary);
   if( !checkpointSolutionStream )
   {
      std::cout << "Checkpoint file for solution cannot open. file name = " << filename << std::endl;
      exit(1);
   }
   if( solution )
      solution->write(checkpointSolutionStream);
   checkpointSolutionStream.close();     /** empty solution file is necessary,
                                          * because it is removed next at the next checkpoint */
}

/** read solution from checkpoint file */
double
ScipParaInitiator::readSolutionFromCheckpointFile(
      char *afterCheckpointingSolutionFileName
      )
{
   char tempSolutionFileName[256];
   sprintf(tempSolutionFileName,"%s_solution.gz", prefixWarm );
   igzstream  checkpointSolutionStream;
   checkpointSolutionStream.open(tempSolutionFileName, std::ios::in | std::ios::binary);
   if( !checkpointSolutionStream )
   {
      std::cout << "checkpoint solution file cannot open: file name = " <<  tempSolutionFileName << std::endl;
      exit(1);
   }
   if( solution )
   {
       ScipParaSolution *sol = dynamic_cast<ScipParaSolution*>(paraComm->createParaSolution());
       if( sol->read(paraComm, checkpointSolutionStream) )
       {
          if( solution->getObjectiveFuntionValue() > sol->getObjectiveFuntionValue() )
          {
             delete solution;
             solution = sol;
          }
       }
       else
       {
          delete sol;
       }
   }
   else
   {
      solution = dynamic_cast<ScipParaSolution*>(paraComm->createParaSolution());
      if( !solution->read(paraComm, checkpointSolutionStream) )
      {
         delete solution;
         solution = 0;
         checkpointSolutionStream.close();
      }
   }
   checkpointSolutionStream.close();
   if( solution )
   {
      if( !tryToSetIncumbentSolution(solution->clone(paraComm)) )
      {
         std::cout << "***** Given solution is wrong! ***************************" << std::endl;
         std::cout << "***** If the solution was given from checkpoint file,  ***" << std::endl;
         std::cout << "***** it might be generated in original problem space   **" << std::endl;
         std::cout << "***** Only primal value is used. *************************" << std::endl;
         std::cout << "***** You should better to use -isol option.  ************" << std::endl;
         std::cout << "***** Or, better to use no distribute solution option. ***" << std::endl;
      }
   }

   /** check if after checkpoing solution file exists or not */
   checkpointSolutionStream.open(afterCheckpointingSolutionFileName, std::ios::in | std::ios::binary);
   if( checkpointSolutionStream )
   {
      /** set up from after checkpointing solution file */
      ScipParaSolution *sol = dynamic_cast<ScipParaSolution*>(paraComm->createParaSolution());
      if( sol->read(paraComm, checkpointSolutionStream) )
      {
         if( !solution )
         {
            solution = sol;
            if( tryToSetIncumbentSolution(solution->clone(paraComm)) )
            {
               std::cout << "***** After checkpoint solution is RIGHT! ****************" << std::endl;
            }
         }
         else
         {
            if( solution->getObjectiveFuntionValue() > sol->getObjectiveFuntionValue() )
            {
               delete solution;
               solution = sol;
               if( tryToSetIncumbentSolution(solution->clone(paraComm)) )
               {
                  std::cout << "***** After checkpoint solution is RIGHT! ****************" << std::endl;
               }
            }
         }
      }
      else
      {
         delete sol;
      }
      checkpointSolutionStream.close();
   }

   if( solution )
   {
      return solution->getObjectiveFuntionValue();
   }
   else
   {
      return DBL_MAX;
   }
}

/** generate racing ramp-up parameter sets */
void
ScipParaInitiator::generateRacingRampUpParameterSets(
      int nParamSets,
      ParaRacingRampUpParamSet **racingRampUpParamSets
      )
{
   ScipDiffParamSet *racingScipDiffParamSet;

   if( racingSettingsName )
   {
      SCIP_CALL_ABORT( SCIPresetParams(scip) );
      SCIP_CALL_ABORT( SCIPreadParams(scip, racingSettingsName) );
      DEF_SCIP_PARA_COMM( scipParaComm, paraComm );
      racingScipDiffParamSet = scipParaComm->createScipDiffParamSet(scip);;
      SCIP_CALL_ABORT( SCIPresetParams(scip) );
   }
   else
   {
      racingScipDiffParamSet = 0;  //  all default
   }


   int n = 0;     /**< keep the number of generated params */
   int npm = -1;  /**< keep the number of variable permutation seed; start from default: -1 */
   int nbo = 0;   /**< keep the number of branching order seed */

   DEF_SCIP_PARA_COMM( scipParaComm, paraComm);

   for(;;)
   {
      for( int i = 0; i < paraParams->getIntParamValue(MaxNRacingParamSetSeed); i++ )
      {
         if( npm > ( paraParams->getIntParamValue(TryNVariablegOrderInRacing) - 1 ) ) npm = -1;
         if( nbo > paraParams->getIntParamValue(TryNBranchingOrderInRacing) ) nbo = 0;
         racingRampUpParamSets[n] = scipParaComm->createScipParaRacingRampUpParamSet(
               paraParams->getIntParamValue(RacingRampUpTerminationCriteria),
               paraParams->getIntParamValue(StopRacingNumberOfNodesLeft),
               paraParams->getRealParamValue(StopRacingTimeLimit),
               i,
               npm,
               nbo,
               racingScipDiffParamSet
               );
         npm++;
         nbo++;
         n++;
         if( n >= nParamSets ) return;
      }
   }
}

/** get solving status string */
std::string
ScipParaInitiator::getStatus(
      )
{
   SCIP_SOL* sol = SCIPgetBestSol(scip);
   if( sol )
   {
      return std::string("solution found exist");
   }
   else
   {
      return std::string("no solution");
   }
}

/** print solver version **/
void
ScipParaInitiator::printSolverVersion(
      std::ostream *os           /**< output file (or NULL for standard output) */
      )
{
#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
   SCIPprintVersion( NULL );
#else
   SCIPprintVersion( scip, NULL );
#endif

   SCIPprintExternalCodes(scip, NULL);
}

/** set initial stat on initiator */
void
ScipParaInitiator::accumulateInitialStat(
      ParaInitialStat *initialStat
      )
{
   ScipParaInitialStat *scipInitialStat = dynamic_cast<ScipParaInitialStat *>(initialStat);
   scipInitialStat->accumulateOn(scip);
}

/** set initial stat on DiffSubproblem */
void
ScipParaInitiator::setInitialStatOnDiffSubproblem(
      int                minDepth,
      int                maxDepth,
      ParaDiffSubproblem *diffSubproblem
      )
{
   ScipParaDiffSubproblem *scipDiffSubproblem = dynamic_cast<ScipParaDiffSubproblem *>(diffSubproblem);
   scipDiffSubproblem->addInitialBranchVarStats(minDepth, maxDepth, scip);
}

/** set final solver status */
void
ScipParaInitiator::setFinalSolverStatus(
      FinalSolverState state
      )
{
   finalState = state;
}

/** set number of nodes solved */
void
ScipParaInitiator::setNumberOfNodesSolved(
      long long n
      )
{
   nSolved = n;
}

/** set final dual bound  */
void
ScipParaInitiator::setDualBound(
      double bound
      )
{
   finalDualBound = bound;
}

/** output solution status */
void
ScipParaInitiator::outputFinalSolverStatistics(
      std::ostream *os,
      double time
      )
{
   if( os == 0 )
   {
      os = &std::cout;
   }

   if( finalState !=  Aborted )
   {
      *os << "SCIP Status        : ";
   }

   switch ( finalState )
   {
   case InitialNodesGenerated:
      *os << "initial nodes were generated" << std::endl;
      break;
   case Aborted:
      *os << std::endl;
      break;
   case HardTimeLimitIsReached:
      *os << "solving was interrupted [ hard time limit reached ]" << std::endl;
      break;
   case ComputingWasInterrupted:
      *os << "solving was interrupted" << std::endl;
      break;
   case ProblemWasSolved:
      *os << "problem is solved" << std::endl;
      break;
   default:
      THROW_LOGICAL_ERROR1("invalid final state");
   }

   *os << "Total Time         : " << time << std::endl;
   *os << "  solving          : " << time << std::endl;
   *os << "  presolving       : " << SCIPgetPresolvingTime(scip) << " (included in solving)" << std::endl;
   *os << "B&B Tree           :" << std::endl;
   *os << "  nodes (total)    : " << nSolved << std::endl;
   *os << "Solution           :" << std::endl;
   *os << "  Solutions found  : " << SCIPgetNSols(scip) << std::endl;
   SCIP_Real primalbound = SCIPinfinity(scip);
   if( SCIPgetNSols(scip) != 0 )
   {
      primalbound = SCIPgetPrimalbound(scip);
   }
   else
   {
      if( solution )
      {
         primalbound = solution->getObjectiveFuntionValue();   // This solution should be in original problem space.
                                                               // So, do not have to convert to original space.
      }
   }
   *os << "  Primal Bound     : ";
   if( SCIPisInfinity(scip, REALABS(primalbound) ) )
   {
      *os << "infeasible" << std::endl;
   }
   else
   {
      (*os).setf(std::ios::showpos);
      *os << std::scientific << std::showpoint << std::setprecision(14) << primalbound << std::endl;
      (*os).unsetf(std::ios::showpos);
   }

   finalDualBound = instance->convertToExternalValue(finalDualBound); // converted to original one
   SCIP_Real gap = 0.0;
   if( SCIPisEQ(scip, primalbound, finalDualBound) )
       gap = 0.0;
   else if( SCIPisZero(scip, finalDualBound)
         || SCIPisZero(scip, primalbound)
         || SCIPisInfinity(scip, REALABS(primalbound))
         || SCIPisInfinity(scip, REALABS(finalDualBound))
         || primalbound * finalDualBound < 0.0 )
      gap = SCIPinfinity(scip);
   else
      gap = REALABS((primalbound - finalDualBound)/MIN(REALABS(finalDualBound),REALABS(primalbound)));

   *os << "Dual Bound         : ";
   if( SCIPisInfinity(scip, REALABS(finalDualBound) ) )
      *os << "         -" << std::endl;
   else
   {
      (*os).setf(std::ios::showpos);
      *os <<  std::scientific << std::showpoint << std::setprecision(14)
      << finalDualBound << std::endl;
      (*os).unsetf(std::ios::showpos);
   }
   *os << "Gap                : ";
   if( SCIPisInfinity(scip, gap ) )
      *os << "  infinite" << std::endl;
   else
   {
      *os << std::fixed << std::setprecision(5) << 100.0 * gap << " %" << std::endl;
   }
}

void
ScipParaInitiator::outputProblemInfo(
      int *nNonLinearConsHdlrs
      )
{
   std::cout << "Original Problem   :" << std::endl;
   std::cout << "  Problem name     : " << SCIPgetProbName(scip) << std::endl;
   std::cout << "  Variables        : " << SCIPgetNOrigVars(scip)
         << " (" << SCIPgetNOrigBinVars(scip) << " binary, "
         << SCIPgetNOrigIntVars(scip) << " integer, "
         << SCIPgetNOrigImplVars(scip) << " implicit integer, "
         << SCIPgetNOrigContVars(scip) << " continuous)" << std::endl;
   std::cout << "  Constraints      : " << SCIPgetNOrigConss(scip) << std::endl;
   std::cout << "  Objective sense  : " <<  (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE ? "minimize" : "maximize") << std::endl;
   std::cout << "Presolved Problem  :" << std::endl;
   std::cout << "  Variables        : " << SCIPgetNVars(scip)
         << " (" << SCIPgetNBinVars(scip) << " binary, "
         << SCIPgetNIntVars(scip) << " integer, "
         << SCIPgetNImplVars(scip) << " implicit integer, "
         << SCIPgetNContVars(scip) << " continuous)" << std::endl;
   std::cout << "  Constraints      : " << SCIPgetNConss(scip) << std::endl;

   std::cout << "Constraints        : Number" << std::endl;
   for( int i = 0; i < SCIPgetNConshdlrs(scip); ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      int startnactiveconss;
      int maxnactiveconss;
      conshdlr = SCIPgetConshdlrs(scip)[i];
      startnactiveconss = SCIPconshdlrGetStartNActiveConss(conshdlr);
      maxnactiveconss = SCIPconshdlrGetMaxNActiveConss(conshdlr);
      std::cout << "  " << std::setw(17) << std::left << SCIPconshdlrGetName(conshdlr) << ": "
                <<  startnactiveconss <<  ( maxnactiveconss > startnactiveconss ? '+' : ' ') << std::endl;
      if( startnactiveconss > 0 
          && std::string(SCIPconshdlrGetName(conshdlr)) != std::string("linear") )
      {
         *nNonLinearConsHdlrs += startnactiveconss;
      }
   }
}

bool 
ScipParaInitiator::onlyLinearConsHandler(
      )
{
   for( int i = 0; i < SCIPgetNConss(scip); ++i )
   {
      SCIP_CONS** conss = SCIPgetConss(scip);
      SCIP_CONSHDLR* conshdlr = SCIPconsGetHdlr(conss[i]);
      if( std::string(SCIPconshdlrGetName(conshdlr)) != std::string("linear") )
      {
         return false;
      }
   }
   return true;
}
