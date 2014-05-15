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

/**@file    scipDefaultUserPlugins.cpp
 * @brief   Set SCIP default user plugins.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scipUserPlugins.h"
#include "scipParaSolver.h"
#include "scipParaInitiator.h"

using namespace UG;
using namespace ParaSCIP;

void
setUserPlugins(
   ParaInitiator *inInitiator
   )
{
   ScipParaInitiator *initiator = dynamic_cast<ScipParaInitiator *>(inInitiator);
   initiator->setUserPlugins(0);
}

void
setUserPlugins(
   ParaInstance *inInstance
   )
{
   ScipParaInstance *instance = dynamic_cast<ScipParaInstance *>(inInstance);
   instance->setUserPlugins(0);
}

void
setUserPlugins(
   ParaSolver *inSolver
   )
{
   ScipParaSolver *solver = dynamic_cast<ScipParaSolver *>(inSolver);
   solver->setUserPlugins(0);
}
