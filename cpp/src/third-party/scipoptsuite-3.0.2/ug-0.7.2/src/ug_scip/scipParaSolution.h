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

/**@file    scipParaSolution.h
 * @brief   ParaSolution extension for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_SOLUTION_H__
#define __SCIP_PARA_SOLUTION_H__

#include <sstream>
#include "ug/paraSolution.h"
#include "scip/scip.h"

namespace ParaSCIP
{

/** ScipParaSolution class */
class ScipParaSolution : public UG::ParaSolution
{
protected:
   double objectiveFunctionValue;
   int nVars;                       /**< number of variables */
   int *indicesAmongSolvers;        /**< array of variable indices ( unique index )  */
   SCIP_Real *values;               /**< array of bounds which the branchings     */
public:
   /** constructor */
   ScipParaSolution(
         )
         : objectiveFunctionValue(0.0), nVars(0), indicesAmongSolvers(0), values(0)
   {
   }

   ScipParaSolution(
         SCIP_Real      objval,
         int            inNvars,
         SCIP_VAR **    vars,
         SCIP_Real *    vals
         )
        : indicesAmongSolvers(0), values(0)
   {
      objectiveFunctionValue = objval;
      nVars = inNvars;
      indicesAmongSolvers = new int[nVars];
      values = new SCIP_Real[nVars];
      for(int i = 0; i < nVars; i++ )
      {
         indicesAmongSolvers[i] = SCIPvarGetIndex(vars[i]);
         values[i] = vals[i];
      }
   }

   ScipParaSolution(
         double inObjectiveFunctionValue,
         int inNVars,                       /**< number of variables */
         int *inIndicesAmongSolvers,        /**< array of variable indices ( probindex )  */
         SCIP_Real *inValues                /**< array of bounds which the branchings     */
         ) : objectiveFunctionValue(inObjectiveFunctionValue), nVars(inNVars)
   {
       indicesAmongSolvers = new int[inNVars];
       values = new SCIP_Real[inNVars];
       for( int i = 0; i < inNVars; i++ )
       {
          indicesAmongSolvers[i] = inIndicesAmongSolvers[i];
          values[i] = inValues[i];
       }
   }

   /** destructor */
   virtual ~ScipParaSolution(
         )
   {
      if( indicesAmongSolvers ) delete [] indicesAmongSolvers;
      if( values ) delete [] values;
   }

   /** get objective function value */
   double getObjectiveFuntionValue(){ return objectiveFunctionValue; }

   /** set objective function value */
   void setObjectiveFuntionValue(SCIP_Real val){ objectiveFunctionValue = val; }

   /** get number of variables */
   int getNVars(){ return nVars; }

   int indexAmongSolvers(int index){ return indicesAmongSolvers[index]; }

   SCIP_Real *getValues(){ return values; }
   void setValue(int i, SCIP_Real val){ assert(i < nVars && i >= 0 ); values[i] = val; }

   /** user should implement write method */
   void write(ogzstream &out);

   /** user should implement read method */
   bool read(UG::ParaComm *comm, igzstream &in);

   /** stringfy solution */
   const std::string toString(
         )
   {
      std::ostringstream s;
      s << "Obj = " << objectiveFunctionValue << std::endl;
      for( int i = 0; i < nVars; i++ )
      {
          s << "i = " << i;
          s << ", idx = " << indicesAmongSolvers[i];
          s << ", val = " << values[i] << std::endl;
      }
      return s.str();
   }

};

}

#endif // __SCIP_PARA_SOLUTION_H__

