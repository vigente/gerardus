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

/**@file    scipParaInstance.h
 * @brief   ParaInstance extenstion for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PARA_INSTANCE_H__
#define __SCIP_PARA_INSTANCE_H__

#if defined(_COMM_MPI_WORLD)
#include "ug/paraDef.h"
#include "ug/paraInstance.h"
#include "scipUserPlugins.h"
#include "scip/scip.h"

namespace ParaSCIP
{

/** ScipInstance */
class ScipParaInstance : public UG::ParaInstance
{
protected:
   int           lProbName;               /**< length of problem name */
   char          *probName;               /**< problem name */
   int           origObjSense;                /**< objective sense : SCIP_OBJSENSE_MAXIMIZE = -1,  SCIP_OBJSENSE_MINIMIZE = +1 */
   /** Do not set objScale and objOffset in SCIP */
   SCIP_Real     objScale;                /**< scalar applied to objective function; original objective value is extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real     objOffset;               /**< objective offset from bound shifting and fixing */
   int           nVars;                   /**< number of variables */
   SCIP_Real     *varLbs;                 /**< array of current lower bound of variable */
   SCIP_Real     *varUbs;                 /**< array of current upper bound of variable */
   SCIP_Real     *objCoefs;               /**< coefficient values: allocation size = nVars */
   int           *varTypes;               /**< array of variable type */
   int           lVarNames;               /**< length of varNames area */
   char          *varNames;               /**< variable names : names are concatenated */
   int           *posVarNames;            /**< positions of varNames */
   int           *mapToOriginalIndecies;  /**< array of indices to map to original problem's probindices */
                                          /**< NOTE: these indices are not transferred. Only valid for Initiator object */
   int           nConss;                  /**< number of constraints */
   int           lConsNames;              /**< length of consNames area */
   char          *consNames;              /**< constraint names : names are concatenated */
   int           *posConsNames;           /**< positions of consNames */


   /*************************
    * for linear constrains *
    * ***********************/
   int           nLinearConss;             /**< number of linear constrains */
   int           *idxLinearConsNames;      /**< array of indices to consName */
   SCIP_Real     *linearLhss;              /**< array of lhs */
   SCIP_Real     *linearRhss;              /**< array of rhs */
   int           *nLinearCoefs;            /**< array of number of coefficient values for linear constrains */
   SCIP_Real     **linearCoefs;            /**< array of non-zero coefficient values of linear constrains */
   int           **idxLinearCoefsVars;     /**< array of indices of on-zero coefficient values of linear constrains */

   /********************************************************
    * for setppc constrains (NOTE:: should be normalized ) *
    * *****************************************************/
   int           nSetppcConss;             /**< number of setppc constrains */
   int           *idxSetppcConsNames;      /**< array of indices to consName */
   int           *nIdxSetppcVars;          /**< array of numbers of indices of variables for setppc constrains */
   int           *setppcTypes;             /**< setppc Types */
   int           **idxSetppcVars;          /**< array of indices of variables for setppc constrains */

   /*********************************
    * for logicor constrains (NOTE:: should be normalized ) *
    *********************************/
   int           nLogicorConss;            /**< number of logical constrains */
   int           *idxLogicorConsNames;     /**< array of indices to consName */
   int           *nIdxLogicorVars;         /**< array of number of indices of variables for logicor constrains */
   int           **idxLogicorVars;         /**< array of indices of of variables for logicor constrains */

   /*********************************
    * for knapsack constrains (NOTE:: should be normalized ) *
    *********************************/
   int           nKnapsackConss;           /**< number of knapsack constrains */
   int           *idxKnapsackConsNames;    /**< array of indices to consName */
   SCIP_Longint  *capacities;              /**< array of capacities(rhs) */
   int           *nLKnapsackCoefs;         /**< array of number of coefficient values for knapsack constrains */
   SCIP_Longint  **knapsackCoefs;          /**< array of non-zero coefficient values of knapsack constrains */
   int           **idxKnapsackCoefsVars;   /**< array of indices of on-zero coefficient values of knapsack constrains */

   /**********************************
    * for varbound constrains (NOTE:: should be normalized ) *
    *********************************/
   int           nVarboundConss;           /**< number of varbound constrains */
   int           *idxVarboundConsNames;    /**< array of indices to consName */
   SCIP_Real     *varboundLhss;            /**< array of lhs */
   SCIP_Real     *varboundRhss;            /**< array of rhs */
   int           *idxVarboundCoefVar1s;    /**< array of indices of variable 1 */
   SCIP_Real     *varboundCoef2s;          /**< array of coefficient of variable 2 */
   int           *idxVarboundCoefVar2s;    /**< array of indices of variable 2 */

   /**********************************
    * for bounddisjunction constrains (NOTE:: should be normalized ) *
    *********************************/
   int           nVarBoundDisjunctionConss;   /**< number of bounddisjunction constrains */
   int           *idxBoundDisjunctionConsNames;    /**< array of indices to consName */
   int           *nVarsBoundDisjunction;           /** array of the number of variables */
   int           **idxVarBoundDisjunction;    /**< array of indices of variables */
   SCIP_BOUNDTYPE **boundTypesBoundDisjunction;      /**< array of bound types */
   SCIP_Real     **boundsBoundDisjunction;          /**< array of bounds */

   /*********************************
    * for SOS1 constraints (NOTE:: should be normalized ) *
    *********************************/
   int           nSos1Conss;               /**< number of SOS1 constraints */
   int           *idxSos1ConsNames;        /**< array of indices to consName */
   int           *nSos1Coefs;              /**< array of number of coefficient values for SOS1 constrains */
   SCIP_Real     **sos1Coefs;              /**< array of non-zero coefficient values of SOS1 constrains */
   int           **idxSos1CoefsVars;       /**< array of indices of on-zero coefficient values of SOS1 constrains */

   /*********************************
    * for SOS2 constraints (NOTE:: should be normalized ) *
    *********************************/
   int           nSos2Conss;               /**< number of SOS1 constraints */
   int           *idxSos2ConsNames;        /**< array of indices to consName */
   int           *nSos2Coefs;              /**< array of number of coefficient values for SOS2 constrains */
   SCIP_Real     **sos2Coefs;              /**< array of non-zero coefficient values of SOS2 constrains */
   int           **idxSos2CoefsVars;       /**< array of indices of on-zero coefficient values of SOS2 constrains */

   /**********************************
    * for aggregated constrains      *
    *********************************/
   int           nAggregatedConss;         /**< number of aggregated constrains = number of aggregated vars */
   int           lAggregatedVarNames;      /**< length of aggregatedVarNames area */
   char          *aggregatedVarNames;      /**< aggregated var names: names are concatenated */
   int           *posAggregatedVarNames;   /**< positions of aggregatedVarNames */
                                           /** no need lower bound and upper bound: always free */
   int           lAggregatedConsNames;     /**< length of aggregatedConsNames area */
   char          *aggregatedConsNames;     /**< aggregated cons names: names are concatenated */
   int           *posAggregatedConsNames;  /**< positions of aggregatedConsNames */
   SCIP_Real     *aggregatedLhsAndLhss;    /**< array of lhs and rhs (lhs = rhs) */
   int           *nAggregatedCoefs;        /**< array of number of coefficient values for aggregated constrains */
   SCIP_Real     **aggregatedCoefs;        /**< array of non-zero coefficient values of aggregated constrains */
   int           **idxAggregatedCoefsVars; /**< array of indices of on-zero coefficient values of aggregated constrains */

   ScipUserPlugins  *userPlugins;          /**< user plugins */

   /**********************************
    *  private functions             *
    **********************************/
   void allocateMemoryForOrdinaryConstraints();
   void addOrdinaryConstraintName(int c, SCIP_CONS *cons );
   void setLinearConstraint( SCIP *scip, int  c, SCIP_CONS *cons );
   void createLinearConstraintsInSCIP( SCIP *scip );
   void setSetppcConstraint( SCIP *scip, int  c, SCIP_CONS *cons );
   void createSetppcConstraintsInSCIP( SCIP *scip );
   void setLogicorConstraint( SCIP *scip, int  c, SCIP_CONS *cons );
   void createLogicorConstraintsInSCIP( SCIP *scip );
   void setKnapsackConstraint( SCIP *scip, int  c, SCIP_CONS *cons );
   void createKnapsackConstraintsInSCIP( SCIP *scip );
   void setVarboundConstraint( SCIP *scip, int  c, SCIP_CONS *cons );
   void createVarboundConstraintsInSCIP( SCIP *scip );
   void setBoundDisjunctionConstraint( SCIP *scip, int  c, SCIP_CONS *cons );
   void createBoundDisjunctionConstraintInSCIP( SCIP *scip );
   void setSos1Constraint( SCIP *scip, int  c, SCIP_CONS *cons, SCIP_CONS** consSOS1 );
   void createSos1ConstraintsInSCIP( SCIP *scip );
   void setSos2Constraint( SCIP *scip, int  c, SCIP_CONS *cons, SCIP_CONS** consSOS2 );
   void createSos2ConstraintsInSCIP( SCIP *scip );
   void getActiveVariables( SCIP *scip, SCIP_VAR **vars, SCIP_Real *scalars, int *nvars,
         SCIP_Real *constant,  SCIP_Bool transformed );
   void collectAggregatedVars(SCIP *scip, int nvars, SCIP_VAR** vars, int* nAggregatedVars,
         SCIP_VAR*** aggregatedVars, SCIP_HASHTABLE**   varAggregated );
   void setAggregatedConstraint(
         SCIP*                 scip,               /**< SCIP data structure */
         int                   c,                  /** aggregated constraint number */
         const char*           constName,          /**< constraint name */
         SCIP_VAR**            vars,               /**< array of variables */
         SCIP_Real*            vals,               /**< array of values */
         int                   nvars,              /**< number of variables */
         SCIP_Real             lhsAndrhs           /**< right hand side = left hand side */
         );
   void setAggregatedConstrains(
         SCIP*              scip,               /**< SCIP data structure */
         int                nvars,              /**< number of mutable variables in the problem */
         int                nAggregatedVars,    /**< number of aggregated variables */
         SCIP_VAR**         aggregatedVars      /**< array storing the aggregated variables */
         );
   void createAggregatedVarsAndConstrainsInSCIP( SCIP *scip );
   virtual char *getFileName() = 0;
public:
   /** constructor */
   ScipParaInstance(
         ) : lProbName(0), probName(0), origObjSense(0), objScale(0.0),
         objOffset(0.0), nVars(0),
         varLbs(0), varUbs(0), objCoefs(0), varTypes(0), lVarNames(0), varNames(0),
         posVarNames(0), mapToOriginalIndecies(0), nConss(0), lConsNames(0), consNames(0),
         posConsNames(0), nLinearConss(0), idxLinearConsNames(0), linearLhss(0),
         linearRhss(0), nLinearCoefs(0), linearCoefs(0), idxLinearCoefsVars(0),
         nSetppcConss(0), idxSetppcConsNames(0), nIdxSetppcVars(0), setppcTypes(0),
         idxSetppcVars(0), nLogicorConss(0), idxLogicorConsNames(0), nIdxLogicorVars(0), idxLogicorVars(0),
         nKnapsackConss(0), idxKnapsackConsNames(0), capacities(0), nLKnapsackCoefs(0),
         knapsackCoefs(0), idxKnapsackCoefsVars(0),
         nVarboundConss(0), idxVarboundConsNames(0), varboundLhss(0),varboundRhss(0), idxVarboundCoefVar1s(0),
         varboundCoef2s(0), idxVarboundCoefVar2s(0),
         nSos1Conss(0),idxSos1ConsNames(0), nSos1Coefs(0), sos1Coefs(0), idxSos1CoefsVars(0),
         nSos2Conss(0), idxSos2ConsNames(0), nSos2Coefs(0), sos2Coefs(0), idxSos2CoefsVars(0),
         nAggregatedConss(0), lAggregatedVarNames(0), aggregatedVarNames(0), posAggregatedVarNames(0),
         lAggregatedConsNames(0), aggregatedConsNames(0), posAggregatedConsNames(0),
         aggregatedLhsAndLhss(0), nAggregatedCoefs(0), aggregatedCoefs(0), idxAggregatedCoefsVars(0)
   {
   }

   /** constractor : only called from ScipInitiator */
   ScipParaInstance(
         SCIP *scip,
         int  method      // transferring method
         );

   /** destractor */
   virtual ~ScipParaInstance(
	        );

  /** convert an internal value to external value */
  double convertToExternalValue(double internalValue)
  {
     return  ( (internalValue + objOffset) * objScale * origObjSense );
  }

  /** convert an external value to internal value */
  double convertToInternalValue(double exteranlValue)
  {
     return  ( ( exteranlValue / ( objScale * origObjSense ) ) - objOffset );
  }

  /** create presolved problem instance that is solved by ParaSCIP */
  void createProblem(
        SCIP *scip,
        int  method,               // transferring method
        bool noPreprocessingInLC,  // LC preprocesing settings
        char *settingsNameLC       // LC preprocesing settings
        );

  /** stringfy ParaCalculationState */
  const std::string toString(
        );

  // int getOrigProbIndex(int index){ return mapToOriginalIndecies[index]; }
  const char *getProbName(){ return probName; }

  void freeMemory();

  int getNVars(){return nVars;}
  SCIP_Real getVarLb(int i){ return varLbs[i]; }
  SCIP_Real getVarUb(int i){ return varUbs[i]; }
  int getVarType(int i){ return varTypes[i]; }

  SCIP *getScip(){ THROW_LOGICAL_ERROR1("This function is only for FiberSCIP!!"); return 0; }

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

}

#include "scipParaInstanceMpi.h"
#endif

#if defined(_COMM_PTH)
#include "scipParaInstancePth.h"
#endif

#endif  // __SCIP_PARA_INSTANCE_H__
