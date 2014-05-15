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

/**@file    scipParaInstancePth.h
 * @brief   ScipParaInstance extension for Pthreads communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_INSTANCE_PTH_H__
#define __SCIP_PARA_INSTANCE_PTH_H__

#include <cassert>
#include "ug/paraInstance.h"
#include "ug/paraCommPth.h"
#include "scipUserPlugins.h"
#include "scip/scip.h"

namespace ParaSCIP
{

class ScipParaInstance : public UG::ParaInstance
{
protected:
   SCIP *scip;     // this pointer should point to the scip environment of LoadCoordinator
public:
   /** constructor */
   ScipParaInstance(
         )
   {
   }

   ScipParaInstance(
		 SCIP *inScip
         ) : scip(inScip)
   {
   }
   /** destractor */
   virtual ~ScipParaInstance(
         )
   {
   }

   /** convert an internal value to external value */
   double convertToExternalValue(double internalValue)
   {
      return SCIPretransformObj(scip, internalValue);
   }

   /** convert an external value to internal value */
   double convertToInternalValue(double externalValue)
   {
      return SCIPtransformObj(scip, externalValue);
   }

   /** create presolved problem instance that is solved by ParaSCIP form scip environment in this object */
   void copyScipEnvironment(
         SCIP **scip
         );

   SCIP *getScip(
         )
   {
      return scip;
   }

   /** create presolved problem instance that is solved by ParaSCIP */
   void createProblem(
        SCIP *inScip,
        int  method,               // transferring method
        bool noPreprocessingInLC,  // LC preprocesing settings
        char *settingsNameLC       // LC preprocesing settings
        );

   /** stringfy ParaCalculationState */
   const std::string toString(
        )
   {
      return std::string("Should be written from scip environment.");
   }

   // int getOrigProbIndex(int index){ return mapToOriginalIndecies[index]; }
   const char *getProbName(){ return SCIPgetProbName(scip); }

   int getNVars(){ return SCIPgetNVars(scip); }
   SCIP_Real getVarLb(int i)
   {
      SCIP_VAR **vars = SCIPgetVars(scip);
      return SCIPvarGetLbGlobal(vars[i]);
   }
   SCIP_Real getVarUb(int i)
   {
      SCIP_VAR **vars = SCIPgetVars(scip);
      return SCIPvarGetUbGlobal(vars[i]);
   }
   int getVarType(int i)
   {
      SCIP_VAR **vars = SCIPgetVars(scip);
      return SCIPvarGetType(vars[i]);
   }

   /** set user plugins */
   void setUserPlugins(ScipUserPlugins *inUi) { /** maybe called, no need to do anything */ }

   /** include user plugins */
   void includeUserPlugins(SCIP *inScip){/** should not be called **/}

};

}

namespace ParaSCIP
{

/** ScipInstanceMpi */
class ScipParaInstancePth : public ScipParaInstance
{
public:
   /** constructor */
   ScipParaInstancePth(
         )
   {
   }

   /** constructor : only called from ScipInitiator */
   ScipParaInstancePth(
         SCIP *inScip,
         int method
         );

   /** destractor */
   ~ScipParaInstancePth(
         )
   {
   }

   /** broadcasts instance to all solvers */
   int bcast(UG::ParaComm *comm, int rank, int method);

};

typedef ScipParaInstancePth *ScipParaInstancePthPtr;

}

#endif  // __SCIP_PARA_INSTANCE_PTH_H__

