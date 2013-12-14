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

/**@file    parascip.cpp
 * @brief   ParaSCIP MAIN.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cfloat>
#include "ug/paraInstance.h"
#include "ug/paraLoadCoordinator.h"
#include "ug/paraParamSet.h"
#include "ug/paraRacingRampUpParamSet.h"
#include "ug/paraInitiator.h"
#include "ug/paraNodeMpi.h"
#include "ug/paraSysTimer.h"
#include "scip/scip.h"
#include "scipParaCommMpi.h"
#include "scipParaInstance.h"
#include "scipParaDeterministicTimer.h"
#include "scipParaSolver.h"
#include "scipParaInitiator.h"

using namespace UG;
using namespace ParaSCIP;

extern void
setUserPlugins(ParaInitiator *initiator);
extern void
setUserPlugins(ParaInstance *instance);
extern void
setUserPlugins(ParaSolver *solver);

void
outputCommandLineMessages(
      char **argv
      )
{
   std::cout << std::endl;
   std::cout << "syntax: " << argv[0] << "#MPI_processes(#solvers + 1) parascip_param_file problem_file_name "
             << "[-l <logfile>] [-q] [-sl <settings>] [-s <settings>] [-sr <root_settings>] [-w <prefix_warm>] [-sth <number>] [-fsol <solution_file>] [-isol <initial solution file]" << std::endl;
   std::cout << "  -l <logfile>        : copy output into log file" << std::endl;
   std::cout << "  -q                  : suppress screen messages" << std::endl;
   std::cout << "  -sl <settings>      : load parameter settings (.set) file for LC presolving" << std::endl;
   std::cout << "  -s <settings>       : load parameter settings (.set) file for solvers" << std::endl;
   std::cout << "  -sr <root_settings> : load parameter settings (.set) file for root" << std::endl;
   std::cout << "  -w <prefix_warm>    : warm start file prefix ( prefix_warm_nodes.gz and prefix_warm_solution.txt are read )" << std::endl;
   std::cout << "  -fsol <solution file> : specify output solution file" << std::endl;
   std::cout << "  -isol <intial solution file> : specify initial solution file" << std::endl;
   std::cout << "File names need to be specified by full path form." << std::endl;
}

void
outputParaParamSet(
      ParaParamSet *paraParamSet,
      ParaInitiator *paraInitiator
      )
{
   if( !paraParamSet->getBoolParamValue(Quiet) )
   {
      std::ofstream ofsParamsOutputFile;
      std::ostringstream s;
      if( paraInitiator->getPrefixWarm() )
      {
         s << paraInitiator->getPrefixWarm();
      }
      else
      {
         s << paraParamSet->getStringParamValue(LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName();
      }
      s << ".prm";
      ofsParamsOutputFile.open(s.str().c_str());
      if( !ofsParamsOutputFile ){
         std::cout << "Cannot open ParaParams output file: file name = " << s.str() << std::endl;
         exit(1);
      }
      paraParamSet->write(&ofsParamsOutputFile);
      ofsParamsOutputFile.close();
   }
}

void
outputSolverParams(
      ParaParamSet *paraParamSet,
      ParaInitiator *paraInitiator
      )
{
   if( !paraParamSet->getBoolParamValue(Quiet) )
   {
      std::ofstream ofsSolverParamsOutputFile;
      std::ostringstream s;
      if( paraInitiator->getPrefixWarm() )
      {
         s << paraInitiator->getPrefixWarm();
      }
      else
      {
         s << paraParamSet->getStringParamValue(LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName();
      }
      s << "_solver.prm";
      ofsSolverParamsOutputFile.open(s.str().c_str());
      if( !ofsSolverParamsOutputFile ){
         std::cout << "Cannot open Solver parameters output file: file name = " << s.str() << std::endl;
         exit(1);
      }
      paraInitiator->writeSolverParameters(&ofsSolverParamsOutputFile);
      ofsSolverParamsOutputFile.close();
   }
}

/**************************************************************************************
 *                                                                                    *
 * Command line see outputCommandLineMessages()                                       *
 *                                                                                    *
 **************************************************************************************/
int
main (
      int  argc,
      char **argv
     )
{
   // mtrace();
   static const int solverOrigin = 1;

   bool racingSolversExist = false;
   ParaDeterministicTimer *detTimer = 0;

   ParaSysTimer sysTimer;
   sysTimer.start();

   // ParaComm *comm = new PARA_COMM_TYPE();
   ScipParaCommMpi *comm = new ScipParaCommMpi();
   comm->init(argc,argv);

   ParaTimer *paraTimer = new ParaTimerMpi();
   paraTimer->init(comm);

   ParaParamSet *paraParamSet = comm->createParaParamSet();

#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
   SCIP_CALL_ABORT( SCIPcreateMesshdlrPThreads(1) );
   SCIPmessageSetDefaultHandler();
#endif

   if( comm->getRank() == 0 )
   {
      if( argc < 3 )
      {
         outputCommandLineMessages(argv);
         return 1;
      }
      paraParamSet->read(comm, argv[1]);

      paraParamSet->bcast(comm, 0);
      comm->lcInit(paraParamSet);
      ParaInitiator *paraInitiator = new ScipParaInitiator(comm);
      setUserPlugins(paraInitiator);
      if( paraInitiator->init(paraParamSet, argc, argv) )
      {
         if( paraInitiator->isSolvedAtInit() )
         {
            paraInitiator->outputFinalSolverStatistics(0, paraTimer->getElapsedTime());
            return 0;
         }
      }
      ParaInstance *paraInstance = paraInitiator->getParaInstance();
      if( paraParamSet->getIntParamValue(OutputParaParams) > 0 )
      {
         outputParaParamSet(paraParamSet, paraInitiator);
         outputSolverParams(paraParamSet, paraInitiator);
      }
      paraInstance->bcast(comm, 0, paraParamSet->getIntParamDefaultValue(InstanceTransferMethod) );
      paraInitiator->sendSolverInitializationMessage();  // This messages should be received in constructor of the target Solver
      ParaLoadCoordinator *paraLc;
      if( paraParamSet->getBoolParamValue(Deterministic) )
      {
          detTimer = new ScipParaDeterministicTimer();
      }
      if( paraInitiator->isSolvedAtInit() )
      {
         paraLc = new ParaLoadCoordinator(comm, paraParamSet, paraInitiator, &racingSolversExist, paraTimer, detTimer);
         delete paraLc;
         delete paraInitiator;
         delete paraParamSet;
         delete paraTimer;
         if( detTimer ) delete detTimer;

         sysTimer.stop();
         std::cout << "[ Rank: " << comm->getRank() << " ], UTime = " << sysTimer.getUTime()
               << ", STime = " << sysTimer.getSTime() << ", RTime = " << sysTimer.getRTime() << std::endl;

         delete comm;
         return 0;
      }
      else
      {
         paraLc = new ParaLoadCoordinator(comm, paraParamSet, paraInitiator, &racingSolversExist, paraTimer, detTimer);
      }
      if( paraInitiator->isWarmStarted() )
      {
         paraLc->warmStart();
      }
      else
      {
         if( paraParamSet->getIntParamValue(RampUpPhaseProcess) == 0)
         {
            ParaNode *rootNode = new ParaNodeMpi(
                  NodeId(), NodeId(), 0, -DBL_MAX, -DBL_MAX, -DBL_MAX,
                  paraInitiator->makeRootNodeDiffSubproblem());
            paraLc->run(rootNode);
         }
         else if( paraParamSet->getIntParamValue(RampUpPhaseProcess) >= 1 &&
                  paraParamSet->getIntParamValue(RampUpPhaseProcess) <= 2 )  // racing ramp-up
         {
            ParaRacingRampUpParamSet **racingRampUpParams = new ParaRacingRampUpParamSetPtr[(comm->getSize()-1)];
            paraInitiator->generateRacingRampUpParameterSets( (comm->getSize()-1), racingRampUpParams );
            for( int i = 1; i < comm->getSize(); i++ )
            {
               PARA_COMM_CALL(
                     racingRampUpParams[i-solverOrigin]->send(comm, i)
                     );
            }
            ParaNode *rootNode = new ParaNodeMpi(
                  NodeId(), NodeId(), 0, -DBL_MAX, -DBL_MAX, -DBL_MAX,
                  paraInitiator->makeRootNodeDiffSubproblem());
            paraLc->run(rootNode, (comm->getSize()-1), racingRampUpParams );
            for( int i = 1; i < comm->getSize(); i++ )
            {
               if( racingRampUpParams[i-solverOrigin] ) delete racingRampUpParams[i-solverOrigin];
            }
            delete [] racingRampUpParams;
         }
         else
         {
            THROW_LOGICAL_ERROR2("Invalid RampUpPhaseProcess: ", paraParamSet->getIntParamValue(RampUpPhaseProcess) )
         }
      }
      delete paraLc;
      if( paraInitiator ) delete paraInitiator;
   }
   else
   {
      if( argc < 3 )
      {
         return 1;
      }
      paraParamSet->bcast(comm, 0);
      comm->solverInit(paraParamSet);
      int nNonLinearConsHdlrs = 0;
      PARA_COMM_CALL(
            comm->bcast( &nNonLinearConsHdlrs, 1, ParaINT, 0 )
      );
      ParaInstance *paraInstance = comm->createParaInstance();
      setUserPlugins(paraInstance);
      if( nNonLinearConsHdlrs > 0 )
      {
         ScipParaInstanceMpi *scipParaInstanceMpi = dynamic_cast<ScipParaInstanceMpi *>(paraInstance);
         scipParaInstanceMpi->setFileName(argv[2]);
         paraParamSet->setIntParamValue(InstanceTransferMethod,2);
      }
      paraInstance->bcast(comm, 0, paraParamSet->getIntParamDefaultValue(InstanceTransferMethod) );
      if( paraParamSet->getBoolParamValue(Deterministic) )
      {
          detTimer = new ScipParaDeterministicTimer();
      }
      ParaSolver *paraSolver = new ScipParaSolver(argc, argv, comm, paraParamSet, paraInstance, detTimer);
      setUserPlugins(paraSolver);
      if( paraParamSet->getIntParamValue(RampUpPhaseProcess) == 0 || paraSolver->isWarmStarted() )
      {
         paraSolver->run();
      }
      else if( paraParamSet->getIntParamValue(RampUpPhaseProcess) >= 1 &&
               paraParamSet->getIntParamValue(RampUpPhaseProcess) <= 2 ) // racing ramp-up
      {
         int source;
         int tag;
         comm->probe(&source, &tag);
         if( tag == TagRacingRampUpParamSet )
         {
            ParaRacingRampUpParamSet *racingRampUpParamSet = new ScipParaRacingRampUpParamSetMpi();
            PARA_COMM_CALL(
                  racingRampUpParamSet->receive(comm, 0)
                  );
            paraSolver->run( racingRampUpParamSet );
            // delete racingRampUpParamSet;   // racingRampUpParamSet is set in Solver object and deleted in the object
         }
         else
         {
            if( tag == TagTerminateRequest )
            {
               PARA_COMM_CALL(
                     comm->receive( NULL, 0, ParaBYTE, source, TagTerminateRequest )
                     );
               // when solver is deleted, solver's destructor sends termination status
            }
            else
            {
               THROW_LOGICAL_ERROR2("Invalid Tag is recicv3ed in ParaSCIP solver main: ", tag )
            }
         }
      }
      else
      {
         THROW_LOGICAL_ERROR2("Invalid RampUpPhaseProcess: ", paraParamSet->getIntParamValue(RampUpPhaseProcess) )
      }
      delete paraSolver;
      //if( paraInstance ) delete paraInstance;  /** deleted in paraSolver destructor */
   }

   delete paraParamSet;

   sysTimer.stop();
   std::cout << "[ Rank: " << comm->getRank() << " ], UTime = " << sysTimer.getUTime()
         << ", STime = " << sysTimer.getSTime() << ", RTime = " << sysTimer.getRTime() << std::endl;


#ifndef SCIP_THREADSAFE_MESSAGEHDLRS
   SCIPfreeMesshdlrPThreads();
#endif

   delete paraTimer;
   if( detTimer ) delete detTimer;

   if( racingSolversExist ) comm->abort();

   delete comm;

   return 0;

} /* END main */
