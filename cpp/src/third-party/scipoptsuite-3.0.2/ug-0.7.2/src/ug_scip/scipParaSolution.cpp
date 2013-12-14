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

/**@file    scipParaSolution.cpp
 * @brief   ParaSolution extension for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scipParaSolution.h"

using namespace ParaSCIP;

void
ScipParaSolution::write(
      ogzstream &out
      )
{
   out.write((char *)&objectiveFunctionValue, sizeof(double));
   out.write((char *)&nVars, sizeof(int));
   for(int i = 0; i < nVars; i++ )
   {
      out.write((char *)&indicesAmongSolvers[i], sizeof(int));
      out.write((char *)&values[i], sizeof(SCIP_Real));
   }
}

bool
ScipParaSolution::read(
      UG::ParaComm *comm,
      igzstream &in
      )
{
   in.read((char *)&objectiveFunctionValue, sizeof(double));
   if( in.eof() ) return false;
   in.read((char *)&nVars, sizeof(int));
   indicesAmongSolvers = new int[nVars];
   values = new SCIP_Real[nVars];
   for(int i = 0; i < nVars; i++ )
   {
      in.read((char *)&indicesAmongSolvers[i], sizeof(int));
      in.read((char *)&values[i], sizeof(SCIP_Real));
   }
   return true;
}
