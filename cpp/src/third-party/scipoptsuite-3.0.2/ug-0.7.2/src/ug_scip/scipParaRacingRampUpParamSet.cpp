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

/**@file    scipParaRacingRampUpParamSet.cpp
 * @brief   ParaRacingRampUpParamSet extension for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scipParaComm.h"
#include "scipParaRacingRampUpParamSet.h"

using namespace UG;
using namespace ParaSCIP;

/** write scipParaRacingRampUpParamSet */
void 
ScipParaRacingRampUpParamSet::write(
    ogzstream &out
    )
{
   out.write((char *)&scipRacingParamSeed, sizeof(int));
   out.write((char *)&permuteProbSeed, sizeof(int));
   out.write((char *)&generateBranchOrderSeed, sizeof(int));
   out.write((char *)&scipDiffParamSetInfo, sizeof(int));
   if( scipDiffParamSetInfo )
   {
      scipDiffParamSet->write(out);
   }
}

/** read scipParaRacingRampUpParamSet */
bool 
ScipParaRacingRampUpParamSet::read(
     ParaComm *comm,
     igzstream &in
     )
{
   in.read((char *)&scipRacingParamSeed, sizeof(int));
   if( in.eof() ) return false;
   in.read((char *)&permuteProbSeed, sizeof(int));
   in.read((char *)&generateBranchOrderSeed, sizeof(int));
   in.read((char *)&scipDiffParamSetInfo, sizeof(int));
   if( scipDiffParamSetInfo )
   {
      DEF_SCIP_PARA_COMM(scipParaComm, comm);
      scipDiffParamSet = scipParaComm->createScipDiffParamSet();
      scipDiffParamSet->read(comm, in);
   }
   return true;
}
