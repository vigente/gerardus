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

/**@file    scipParaRacingRampUpParamSet.h
 * @brief   ParaRacingRampUpParamSet extension for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_RACING_RAMP_UP_PARAM_SET_H__
#define __SCIP_PARA_RACING_RAMP_UP_PARAM_SET_H__

#include "scip/scip.h"
#include "ug/paraRacingRampUpParamSet.h"
#include "ug/paraComm.h"
#include "scipDiffParamSet.h"

namespace ParaSCIP
{

/** The racing ramp-up parameter set for SCIP solver */
class ScipParaRacingRampUpParamSet : public UG::ParaRacingRampUpParamSet
{
protected:
   int scipRacingParamSeed;             /**< seed to generate SCIP racing parameter */
   int permuteProbSeed;                 /**< seed to permute problem */
   int generateBranchOrderSeed;         /**< seed to generate branching order */
   int scipDiffParamSetInfo;            /**< 1: with scipDiffParamSet, 0: no scipDiffParamSet */
   ScipDiffParamSet *scipDiffParamSet;  /**< scip parameter set different from default values for racing ramp-up */
public:
   /** default constructor */
   ScipParaRacingRampUpParamSet(
         )
         : ParaRacingRampUpParamSet(), scipRacingParamSeed(-1),
           permuteProbSeed(0), generateBranchOrderSeed(0), scipDiffParamSetInfo(0), scipDiffParamSet(0)
   {
   }

   ScipParaRacingRampUpParamSet(
         int inTerminationCriteria,
         int inNNodesLeft,
         double inTimeLimit,
         int inScipRacingParamSeed,
         int inPermuteProbSeed,
         int inGenerateBranchOrderSeed,
         ScipDiffParamSet *inScipDiffParamSet
         )
         : ParaRacingRampUpParamSet(inTerminationCriteria, inNNodesLeft, inTimeLimit),
           scipRacingParamSeed(inScipRacingParamSeed),permuteProbSeed(inPermuteProbSeed),
           generateBranchOrderSeed(inGenerateBranchOrderSeed), scipDiffParamSetInfo(0), scipDiffParamSet(inScipDiffParamSet)
   {
      if( inScipDiffParamSet ) scipDiffParamSetInfo = 1;
   }

   /** destructor */
   virtual ~ScipParaRacingRampUpParamSet()
   {
      if( scipDiffParamSet ) delete scipDiffParamSet;
   }


   /** getter of permuteProbSeed */
   int getPermuteProbSeed(
         )
   {
      return permuteProbSeed;
   }

   /** getter of generateBranchOrderSeed */
   int getGenerateBranchOrderSeed(
         )
   {
      return generateBranchOrderSeed;
   }

   /** getter of ScipDiffParamSet */
   ScipDiffParamSet *getScipDiffParamSet(
         )
   {
      return scipDiffParamSet;
   }

   int getScipRacingParamSeed(
         )
   {
      return scipRacingParamSeed;
   }

   /** write scipParaRacingRampUpParamSet */
   void write(
         ogzstream &out
         );

   /** read scipParaRacingRampUpParamSet */
   bool read(
         UG::ParaComm *comm,
         igzstream &in
         );

   /** stringfy ScipParaRacingRampUpParamSet */
   const std::string toString(
         )
   {
      std::ostringstream s;
      s << "[ SCIP racing parameter seed; " << scipRacingParamSeed;
      s << ", Permutate problem seed: " << permuteProbSeed << ", Generate branch order seed: " << generateBranchOrderSeed << " ]" << std::endl;
      if( scipDiffParamSetInfo )
      {
         s << scipDiffParamSet->toString();
      }
      return s.str();
   }

   int getStrategy(
         )
   {
      return scipRacingParamSeed;
   }
};

}



#endif    // __SCIP_PARA_RACING_RAMP_UP_PARAM_SET_H__

