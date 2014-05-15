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

/**@file    scipParaInstancePth.cpp
 * @brief   ScipParaInstance extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


// #include "paraDef.h"
#include "ug/paraCommPth.h"
#include "scipParaInstancePth.h"
#include "scip/scipdefplugins.h"

using namespace UG;
using namespace ParaSCIP;

const static char *PRESOLVED_INSTANCE = "presolved.cip";

void
ScipParaInstance::copyScipEnvironment(
      SCIP **targetscip
      )
{
   char probname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(*targetscip != NULL);
   SCIP_Bool success = TRUE;

   /* copy all plugins and settings */
#if SCIP_VERSION == 211 && SCIP_SUBVERSION == 0
   SCIP_CALL_ABORT( SCIPcopyPlugins(scip, *targetscip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, &success) );
#else
   SCIP_CALL_ABORT( SCIPcopyPlugins(scip, *targetscip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, &success) );
#endif
   SCIP_CALL_ABORT( SCIPcopyParamSettings(scip, *targetscip) );


   /* create problem in the target SCIP */
   /* get name of the original problem and add the suffix string */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_%s", SCIPgetProbName(scip), "solver");
   SCIP_CALL_ABORT( SCIPcreateProb(*targetscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create the variable mapping hash map */
   SCIP_HASHMAP* varmap;
   SCIP_CALL_ABORT( SCIPhashmapCreate(&varmap, SCIPblkmem(*targetscip), SCIPgetNVars(scip)) );
   SCIP_HASHMAP* conssmap;
   SCIP_CALL_ABORT( SCIPhashmapCreate(&conssmap, SCIPblkmem(*targetscip), SCIPgetNConss(scip)) );   

   /* copy all variables and constraints */
   SCIP_CALL_ABORT( SCIPcopyVars(scip, *targetscip, varmap, conssmap, TRUE) );
   SCIP_CALL_ABORT( SCIPcopyConss(scip, *targetscip, varmap, conssmap, TRUE, FALSE, &success) );

   /* free hash map */
   SCIPhashmapFree(&varmap);
   SCIPhashmapFree(&conssmap);

}

int
ScipParaInstancePth::bcast(
      ParaComm *comm,
      int root,
      int method
      )
{
   DEF_PARA_COMM( commPth, comm);

#if 0  // sender side creation

   if( commPth->getRank() == root )
   {
      for( int i = 0; i < commPth->getSize(); i++ )
      {
         if( i != root )
         {
            SCIP *newScip;
            SCIP_CALL_ABORT( SCIPcreate(&newScip) );
            copyScipEnvironment(&newScip);   // this copy call should be serialized
            PARA_COMM_CALL(
               commPth->uTypeSend((void *)newScip, ParaInstanceType, i, TagParaInstance)
            );
         }
      }
   }
   else
   {
      // SCIP *received;
      PARA_COMM_CALL(
         commPth->uTypeReceive((void **)&scip, ParaInstanceType, root, TagParaInstance)
      );
      // scip = received;
   }
#else  // receiver side creation
   if( commPth->getRank() == root )
   {
      if( method == 0 )
      {
         for( int i = 0; i < commPth->getSize(); i++ )
         {
            if( i != root )
            {
               PARA_COMM_CALL(
                  commPth->uTypeSend((void *)scip, ParaInstanceType, i, TagParaInstance)
               );
            }
         }
      }
      else
      {
         for( int i = 0; i < commPth->getSize(); i++ )
         {
            if( i != root )
            {
               PARA_COMM_CALL(
                  commPth->send(NULL, 0, ParaBYTE, i, TagParaInstance)
               );
            }
         }
      }
   }
   else
   {
      if( method == 0 )
      {
         SCIP *received;

         PARA_COMM_CALL(
            commPth->uTypeReceive((void **)&received, ParaInstanceType, root, TagParaInstance)
         );

         SCIP_CALL_ABORT( SCIPcreate(&scip) );
         char probname[SCIP_MAXSTRLEN];

         assert(received != NULL);
         assert(scip != NULL);
         SCIP_Bool success = TRUE;

         commPth->lockApp();
         /* copy all plugins and settings */
#if SCIP_VERSION == 211 && SCIP_SUBVERSION == 0
         SCIP_CALL_ABORT( SCIPcopyPlugins(received, scip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
               TRUE, TRUE, TRUE, TRUE, &success) );
#else
         SCIP_CALL_ABORT( SCIPcopyPlugins(received, scip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
               TRUE, TRUE, TRUE, TRUE, FALSE, &success) );
#endif
         SCIP_CALL_ABORT( SCIPcopyParamSettings(received, scip) );

         /* create problem in the target SCIP */
         /* get name of the original problem and add the suffix string */
         (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_%s", SCIPgetProbName(received), "solver");
         SCIP_CALL_ABORT( SCIPcreateProb(scip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

         /* create the variable mapping hash map */
         SCIP_HASHMAP* varmap;
         SCIP_CALL_ABORT( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), SCIPgetNVars(received)) );
         SCIP_HASHMAP* conssmap;
         SCIP_CALL_ABORT( SCIPhashmapCreate(&conssmap, SCIPblkmem(scip), SCIPgetNConss(received)) );

         // commPth->lockApp();
         /* copy all variables and constraints */
         SCIP_CALL_ABORT( SCIPcopyVars(received, scip, varmap, conssmap, TRUE) );
         SCIP_CALL_ABORT( SCIPcopyConss(received, scip, varmap, conssmap, TRUE, FALSE, &success) );
         commPth->unlockApp();

         /* free hash map */
         SCIPhashmapFree(&varmap);
         SCIPhashmapFree(&conssmap);
      }
      else
      {
         PARA_COMM_CALL(
             commPth->receive( NULL, 0, ParaBYTE, root, TagParaInstance )
         );
      }

   }
#endif

   return 0;

}

ScipParaInstancePth::ScipParaInstancePth(
      SCIP *inScip,
      int  method
      ) : ScipParaInstance(inScip)
{
   if( method == 1)
   {
      SCIP_CALL_ABORT( SCIPwriteTransProblem(scip, PRESOLVED_INSTANCE, "cip", TRUE ) );
   }
}

/** create presolved problem instance that is solved by ParaSCIP */
void
ScipParaInstance::createProblem(
     SCIP *inScip,
     int  method,
     bool noPreprocessingInLC,  // LC preprocesing settings
     char *settingsNameLC       // LC preprocesing settings
     )
{
   if( method == 1 )
   {
      scip = inScip;
      SCIP_CALL_ABORT( SCIPreadProb(scip, PRESOLVED_INSTANCE, "cip" ) );
   }
   if( method == 0 && SCIPgetStage(inScip) == SCIP_STAGE_INIT )
   {
      if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM)
      {
         SCIP_CALL_ABORT( SCIPtransformProb(scip));
      }
      if( scip == inScip ) return;
      SCIP_Bool success = TRUE;
      char probname[SCIP_MAXSTRLEN];
      /* copy all plugins and settings */
#if SCIP_VERSION == 211 && SCIP_SUBVERSION == 0
      SCIP_CALL_ABORT( SCIPcopyPlugins(scip, inScip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
          TRUE, TRUE, TRUE, TRUE, &success) );
#else
      SCIP_CALL_ABORT( SCIPcopyPlugins(scip, inScip, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
           TRUE, TRUE, TRUE, TRUE, TRUE, &success) );
#endif
      SCIP_CALL_ABORT( SCIPcopyParamSettings(scip, inScip) );

      /* create problem in the target SCIP */
      /* get name of the original problem and add the suffix string */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_%s", SCIPgetProbName(scip), "solver_created");
      SCIP_CALL_ABORT( SCIPcreateProb(inScip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* create the variable mapping hash map */
      SCIP_HASHMAP* varmap;
      SCIP_CALL_ABORT( SCIPhashmapCreate(&varmap, SCIPblkmem(inScip), SCIPgetNVars(scip)) );
      SCIP_HASHMAP* conssmap;
      SCIP_CALL_ABORT( SCIPhashmapCreate(&conssmap, SCIPblkmem(inScip), SCIPgetNConss(scip)) );

      /* copy all variables and constraints */
      SCIP_CALL_ABORT( SCIPcopyVars(scip, inScip, varmap, conssmap, TRUE) );
      SCIP_CALL_ABORT( SCIPcopyConss(scip, inScip, varmap, conssmap, TRUE, FALSE, &success) );

      /* free hash map */
      SCIPhashmapFree(&varmap);
      SCIPhashmapFree(&conssmap);
   }
   if( method == 2 )
   {
      std::cout << "You should use instance transfer method 0 or 1!" << std::endl;
      exit(0);
   }
}

