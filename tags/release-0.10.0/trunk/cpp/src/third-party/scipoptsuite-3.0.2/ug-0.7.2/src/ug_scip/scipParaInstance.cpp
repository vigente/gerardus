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

/**@file    scipParaInstance.cpp
 * @brief   ParaInstance extenstion for SCIP solver.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cctype>

#include "ug/paraDef.h"
#include "ug/paraInstance.h"
#include "scip/scipdefplugins.h"
#include "scip/reader_lp.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/pub_misc.h"
#include "scipParaInstance.h"

using namespace ParaSCIP;

const static char *PRESOLVED_INSTANCE = "presolved.cip";

/** hash key retrieval function for variables */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the indices of both variables are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqVar)
{  /*lint --e{715}*/
   if ( key1 == key2 )
      return TRUE;
   return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVar)
{  /*lint --e{715}*/
   assert( SCIPvarGetIndex((SCIP_VAR*) key) >= 0 );
   return (unsigned int) SCIPvarGetIndex((SCIP_VAR*) key);
}

void
ScipParaInstance::allocateMemoryForOrdinaryConstraints(
      )
{
   /** for linear constraints */
   if( nLinearConss > 0 )
   {
      idxLinearConsNames = new int[nLinearConss];
      linearLhss = new SCIP_Real[nLinearConss];
      linearRhss = new SCIP_Real[nLinearConss];
      nLinearCoefs = new int[nLinearConss];
      linearCoefs = new SCIP_Real*[nLinearConss];
      idxLinearCoefsVars = new int*[nLinearConss];
   }
   /** for setppc constraints */
   if( nSetppcConss )
   {
      idxSetppcConsNames = new int[nSetppcConss];
      nIdxSetppcVars = new int[nSetppcConss];
      setppcTypes = new int[nSetppcConss];
      idxSetppcVars = new int*[nSetppcConss];
   }
   /** for logical constraints */
   if(nLogicorConss)
   {
      idxLogicorConsNames = new int[nLogicorConss];
      nIdxLogicorVars = new int[nLogicorConss];
      idxLogicorVars = new int*[nLogicorConss];
   }
   /** for knapsack constraints */
   if(nKnapsackConss)
   {
      idxKnapsackConsNames = new int[nKnapsackConss];
      capacities = new SCIP_Longint[nKnapsackConss];
      nLKnapsackCoefs = new int[nKnapsackConss];
      knapsackCoefs = new SCIP_Longint*[nKnapsackConss];
      idxKnapsackCoefsVars = new int*[nKnapsackConss];
   }
   /** for varbound constraints */
   if( nVarboundConss ){
      idxVarboundConsNames = new int[nVarboundConss];
      varboundLhss = new SCIP_Real[nVarboundConss];
      varboundRhss = new SCIP_Real[nVarboundConss];
      idxVarboundCoefVar1s = new int[nVarboundConss];
      varboundCoef2s = new SCIP_Real[nVarboundConss];
      idxVarboundCoefVar2s = new int[nVarboundConss];
   }
   /** for bounddisjunction constraints */
   if( nVarBoundDisjunctionConss ){
      idxBoundDisjunctionConsNames = new int[nVarBoundDisjunctionConss];
      nVarsBoundDisjunction = new int[nVarBoundDisjunctionConss];
      boundsBoundDisjunction = new SCIP_Real*[nVarBoundDisjunctionConss];
      boundTypesBoundDisjunction = new SCIP_BOUNDTYPE*[nVarBoundDisjunctionConss];
      idxVarBoundDisjunction = new int*[nVarBoundDisjunctionConss];
   }
   /** for SOS1 constraints */
   if( nSos1Conss ){
      idxSos1ConsNames = new int[nSos1Conss];
      nSos1Coefs = new int[nSos1Conss];
      sos1Coefs = new SCIP_Real*[nSos1Conss];
      idxSos1CoefsVars = new int*[nSos1Conss];
   }
   /** for SOS2 constraints */
   if( nSos2Conss ){
      idxSos2ConsNames = new int[nSos2Conss];
      nSos2Coefs = new int[nSos2Conss];
      sos2Coefs = new SCIP_Real*[nSos2Conss];
      idxSos2CoefsVars = new int*[nSos2Conss];
   }
}

void
ScipParaInstance::addOrdinaryConstraintName(
      int c,
      SCIP_CONS *cons
      )
{
   posConsNames[c] = lConsNames;
   (void) strcpy( &consNames[lConsNames], SCIPconsGetName(cons) );
   lConsNames += strlen( SCIPconsGetName(cons) )+ 1;
}

void
ScipParaInstance::setLinearConstraint(
      SCIP *scip,
      int  c,
      SCIP_CONS *cons
      )
{
   addOrdinaryConstraintName(c, cons);
   idxLinearConsNames[nLinearConss] = c;
   linearLhss[nLinearConss] = SCIPgetLhsLinear(scip, cons);
   linearRhss[nLinearConss] = SCIPgetRhsLinear(scip, cons);
   int nvars = SCIPgetNVarsLinear(scip, cons);
   nLinearCoefs[nLinearConss] = nvars;
   linearCoefs[nLinearConss] = new SCIP_Real[nvars];
   idxLinearCoefsVars[nLinearConss] = new int[nvars];
   SCIP_Real *coefs = SCIPgetValsLinear(scip, cons);
   SCIP_VAR  **vars = SCIPgetVarsLinear(scip, cons);
   for(int v = 0; v < nvars; ++v)
   {
      linearCoefs[nLinearConss][v] =  coefs[v];
      idxLinearCoefsVars[nLinearConss][v] = SCIPvarGetProbindex(vars[v]);
   }
   nLinearConss++;
}

void
ScipParaInstance::createLinearConstraintsInSCIP(
      SCIP *scip
      )
{
   SCIP_VAR **vars = SCIPgetVars(scip);
   for(int c = 0; c < nLinearConss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR **varsInCons = new SCIP_VAR*[nLinearCoefs[c]];
      for( int v = 0; v < nLinearCoefs[c]; v++)
      {
         varsInCons[v] = vars[idxLinearCoefsVars[c][v]];
      }
      SCIP_CALL_ABORT( SCIPcreateConsLinear(scip, &cons,
            &consNames[posConsNames[idxLinearConsNames[c]]], nLinearCoefs[c], varsInCons, linearCoefs[c],linearLhss[c], linearRhss[c],
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIP_CALL_ABORT( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      delete [] varsInCons;
   }
}

/*************************************************************/
/*** NOTE: setSetppcConstraint routine cannot work correctly */
/*************************************************************/
void
ScipParaInstance::setSetppcConstraint(
      SCIP *scip,
      int  c,
      SCIP_CONS *cons
      )
{
   addOrdinaryConstraintName(c, cons);
   idxSetppcConsNames[nSetppcConss] = c;
   int nvars = SCIPgetNVarsSetppc(scip,cons);
   nIdxSetppcVars[nSetppcConss] = nvars;
   setppcTypes[nSetppcConss] = SCIPgetTypeSetppc(scip,cons);
   idxSetppcVars[nSetppcConss] = new int[nvars];
   SCIP_VAR **vars = SCIPgetVarsSetppc(scip, cons);
   for( int v = 0; v < nvars; ++v)
   {
      idxSetppcVars[nSetppcConss][v] = SCIPvarGetProbindex(vars[v]);
   }
   nSetppcConss++;
}

void
ScipParaInstance::createSetppcConstraintsInSCIP(
      SCIP *scip
      )
{
   SCIP_VAR **vars = SCIPgetVars(scip);
   for(int c = 0; c < nSetppcConss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR **varsInCons = new SCIP_VAR*[nIdxSetppcVars[c]];
      for( int v = 0; v < nIdxSetppcVars[c]; v++)
      {
         varsInCons[v] = vars[idxSetppcVars[c][v]];
      }
      switch (setppcTypes[c])
      {
      case SCIP_SETPPCTYPE_PARTITIONING:
         SCIP_CALL_ABORT( SCIPcreateConsSetpart(scip, &cons,
               &consNames[posConsNames[idxSetppcConsNames[c]]], nIdxSetppcVars[c], varsInCons,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         break;
      case SCIP_SETPPCTYPE_PACKING:
         SCIP_CALL_ABORT( SCIPcreateConsSetpack(scip, &cons,
               &consNames[posConsNames[idxSetppcConsNames[c]]], nIdxSetppcVars[c], varsInCons,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         break;
      case SCIP_SETPPCTYPE_COVERING:
         SCIP_CALL_ABORT( SCIPcreateConsSetcover(scip, &cons,
               &consNames[posConsNames[idxSetppcConsNames[c]]], nIdxSetppcVars[c], varsInCons,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         break;
      default:
         THROW_LOGICAL_ERROR2("invalid setppc constraint type: type = ", setppcTypes[c]);
      }
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIP_CALL_ABORT( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      delete [] varsInCons;
   }
}

/**************************************************************/
/*** NOTE: setLogicorConstraint routine cannot work correctly */
/**************************************************************/
void
ScipParaInstance::setLogicorConstraint(
      SCIP *scip,
      int  c,
      SCIP_CONS *cons
      )
{
   addOrdinaryConstraintName(c, cons);
   idxLogicorConsNames[nLogicorConss] = c;
   int nvars = SCIPgetNVarsLogicor(scip, cons);
   nIdxLogicorVars[nLogicorConss] = nvars;
   idxLogicorVars[nLogicorConss] = new int[nvars];

   SCIP_VAR **vars =  SCIPgetVarsLogicor(scip, cons);
   for( int v = 0; v < nvars; ++v )
   {
      if( SCIPvarIsActive(vars[v]) )
         idxLogicorVars[nLogicorConss][v] = SCIPvarGetProbindex(vars[v]);
      else if( SCIPvarIsNegated(vars[v]) )
         idxLogicorVars[nLogicorConss][v] = -SCIPvarGetProbindex(SCIPvarGetNegationVar(vars[v])) - 1;
      else
         idxLogicorVars[nLogicorConss][v] = INT_MAX;
   }
   nLogicorConss++;
}

void
ScipParaInstance::createLogicorConstraintsInSCIP(
      SCIP *scip
      )
{
   SCIP_VAR **vars = SCIPgetVars(scip);
   for(int c = 0; c < nLogicorConss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR **varsInCons = new SCIP_VAR*[nIdxLogicorVars[c]];
      int v;
      for( v = 0; v < nIdxLogicorVars[c]; v++)
      {
         if( idxLogicorVars[c][v] == INT_MAX )
            break;
         /* negated variable */
         if( idxLogicorVars[c][v] < 0 )
         {
            SCIP_CALL_ABORT( SCIPgetNegatedVar(scip, vars[-(idxLogicorVars[c][v] + 1)], &varsInCons[v]) );
         }
         else
            varsInCons[v] = vars[idxLogicorVars[c][v]];
      }
      if( v == nIdxLogicorVars[c] )
      {
         SCIP_CALL_ABORT( SCIPcreateConsLogicor(scip, &cons,
               &consNames[posConsNames[idxLogicorConsNames[c]]], nIdxLogicorVars[c], varsInCons,
               TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
         SCIPdebug( SCIP_CALL_ABORT( SCIPprintCons(scip, cons, NULL) ) );
         SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      }
      delete [] varsInCons;
   }
}

/***************************************************************/
/*** NOTE: setKnapsackConstraint routine cannot work correctly */
/***************************************************************/
void
ScipParaInstance::setKnapsackConstraint(
      SCIP *scip,
      int  c,
      SCIP_CONS *cons
      )
{
   addOrdinaryConstraintName(c, cons);
   idxKnapsackConsNames[nKnapsackConss] = c;
   capacities[nKnapsackConss] = SCIPgetCapacityKnapsack(scip, cons);
   int nvars = SCIPgetNVarsKnapsack(scip, cons);
   nLKnapsackCoefs[nKnapsackConss] = nvars;
   knapsackCoefs[nKnapsackConss] = new SCIP_Longint[nvars];
   idxKnapsackCoefsVars[nKnapsackConss] = new int[nvars];
   SCIP_Longint *weights = SCIPgetWeightsKnapsack(scip, cons);
   SCIP_VAR **vars = SCIPgetVarsKnapsack(scip,cons);
   for(int v = 0; v < nvars; ++v )
   {
      knapsackCoefs[nKnapsackConss][v] = weights[v];
      idxKnapsackCoefsVars[nKnapsackConss][v] = SCIPvarGetProbindex(vars[v]);
   }
   nKnapsackConss++;
}


void
ScipParaInstance::createKnapsackConstraintsInSCIP(
      SCIP *scip
      )
{
   SCIP_VAR **vars = SCIPgetVars(scip);
   for(int c = 0; c < nKnapsackConss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR **varsInCons = new SCIP_VAR*[nLKnapsackCoefs[c]];
      for( int v = 0; v < nLKnapsackCoefs[c]; v++)
      {
         varsInCons[v] = vars[idxKnapsackCoefsVars[c][v]];
      }
      SCIP_CALL_ABORT( SCIPcreateConsKnapsack(scip, &cons,
            &consNames[posConsNames[idxKnapsackConsNames[c]]], nLKnapsackCoefs[c], varsInCons, knapsackCoefs[c], capacities[c],
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIP_CALL_ABORT( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      delete [] varsInCons;
   }
}

/***************************************************************/
/*** NOTE: setVarboundConstraint routine cannot work correctly */
/***************************************************************/
void
ScipParaInstance::setVarboundConstraint(
      SCIP *scip,
      int  c,
      SCIP_CONS *cons
      )
{
   addOrdinaryConstraintName(c, cons);
   idxVarboundConsNames[nVarboundConss] = c;
   varboundLhss[nVarboundConss] = SCIPgetLhsVarbound(scip,cons);
   varboundRhss[nVarboundConss] = SCIPgetRhsVarbound(scip,cons);
   idxVarboundCoefVar1s[nVarboundConss] =  SCIPvarGetProbindex(SCIPgetVarVarbound(scip,cons));
   varboundCoef2s[nVarboundConss] = SCIPgetVbdcoefVarbound(scip,cons);
   idxVarboundCoefVar2s[nVarboundConss] =  SCIPvarGetProbindex(SCIPgetVbdvarVarbound(scip,cons));
   nVarboundConss++;
}

void
ScipParaInstance::createVarboundConstraintsInSCIP(
      SCIP *scip
      )
{
   SCIP_VAR **vars = SCIPgetVars(scip);
   for(int c = 0; c < nVarboundConss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR *var1InCons = vars[idxVarboundCoefVar1s[c]];
      SCIP_VAR *var2InCons = vars[idxVarboundCoefVar2s[c]];
      SCIP_CALL_ABORT( SCIPcreateConsVarbound(scip, &cons,
            &consNames[posConsNames[idxVarboundConsNames[c]]], var1InCons, var2InCons, varboundCoef2s[c], varboundLhss[c], varboundRhss[c],
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIP_CALL_ABORT( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
   }
}

/***************************************************************/
/*** NOTE: setVarboundConstraint routine cannot work correctly */
/***************************************************************/
void
ScipParaInstance::setBoundDisjunctionConstraint(
      SCIP *scip,
      int  c,
      SCIP_CONS *cons
      )
{
   addOrdinaryConstraintName(c, cons);
   idxBoundDisjunctionConsNames[nVarBoundDisjunctionConss] = c;
   nVarsBoundDisjunction[nVarBoundDisjunctionConss] = SCIPgetNVarsBounddisjunction(scip,cons);
   idxVarBoundDisjunction[nVarBoundDisjunctionConss] = new int[nVarsBoundDisjunction[nVarBoundDisjunctionConss]];
   boundsBoundDisjunction[nVarBoundDisjunctionConss] = new SCIP_Real[nVarsBoundDisjunction[nVarBoundDisjunctionConss]];
   boundTypesBoundDisjunction[nVarBoundDisjunctionConss] = new SCIP_BOUNDTYPE[nVarsBoundDisjunction[nVarBoundDisjunctionConss]];
   SCIP_VAR **vars = SCIPgetVarsBounddisjunction ( scip, cons );
   SCIP_BOUNDTYPE *boundTypes = SCIPgetBoundtypesBounddisjunction( scip, cons );
   SCIP_Real *bounds = SCIPgetBoundsBounddisjunction( scip, cons );
   for( int v = 0; v < nVarsBoundDisjunction[nVarBoundDisjunctionConss]; ++v )
   {
      idxVarBoundDisjunction[nVarBoundDisjunctionConss][v] = SCIPvarGetProbindex(vars[v]);
      boundsBoundDisjunction[nVarBoundDisjunctionConss][v] = bounds[v];
      boundTypesBoundDisjunction[nVarBoundDisjunctionConss][v] = boundTypes[v];
   }
   nVarBoundDisjunctionConss++;
}

void
ScipParaInstance::createBoundDisjunctionConstraintInSCIP(
      SCIP *scip
      )
{
   SCIP_VAR **vars = SCIPgetVars(scip);
   for(int c = 0; c < nVarBoundDisjunctionConss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR **boundVars = new SCIP_VAR*[nVarsBoundDisjunction[c]];
      SCIP_BOUNDTYPE *types = new SCIP_BOUNDTYPE[nVarsBoundDisjunction[c]];
      for( int v = 0; v < nVarsBoundDisjunction[c]; ++v )
      {
         boundVars[v] = vars[idxVarBoundDisjunction[c][v]];
         types[v] = boundTypesBoundDisjunction[c][v];
      }

      SCIP_CALL_ABORT( SCIPcreateConsBounddisjunction(scip, &cons,
            &consNames[posConsNames[idxBoundDisjunctionConsNames[c]]], nVarsBoundDisjunction[c], boundVars,
            types, boundsBoundDisjunction[c],
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIP_CALL_ABORT( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      delete [] boundVars;
      delete [] types;
   }
}

void
ScipParaInstance::setSos1Constraint(
      SCIP *scip,
      int  c,
      SCIP_CONS *cons,
      SCIP_CONS** consSOS1
      )
{
   addOrdinaryConstraintName(c, cons);
   idxSos1ConsNames[nSos1Conss] = c;
   int nvars = SCIPgetNVarsSOS1(scip, cons);
   nSos1Coefs[nSos1Conss] = nvars;
   sos1Coefs[nSos1Conss] = new SCIP_Real[nvars];
   idxSos1CoefsVars[nSos1Conss] = new int[nvars];
   SCIP_Real *weights = SCIPgetWeightsSOS1(scip, cons);
   SCIP_VAR **vars = SCIPgetVarsSOS1(scip,cons);
   for( int v = 0; v < nvars; v++ )
   {
      sos1Coefs[nSos1Conss][v] = weights[v];
      idxSos1CoefsVars[nSos1Conss][v] = SCIPvarGetProbindex(vars[v]);
   }
   // store constraint
   consSOS1[nSos1Conss++] = cons;
}

void
ScipParaInstance::createSos1ConstraintsInSCIP(
      SCIP *scip
      )
{
   SCIP_VAR **vars = SCIPgetVars(scip);
   for(int c = 0; c < nSos1Conss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR **varsInCons = new SCIP_VAR*[nSos1Coefs[c]];
      for( int v = 0; v < nSos1Coefs[c]; v++)
      {
         varsInCons[v] = vars[idxSos1CoefsVars[c][v]];
      }
      SCIP_CALL_ABORT( SCIPcreateConsSOS1(scip, &cons,
            &consNames[posConsNames[idxSos1ConsNames[c]]], nSos1Coefs[c], varsInCons, sos1Coefs[c],
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE ) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIP_CALL_ABORT( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      delete [] varsInCons;
   }
}

void
ScipParaInstance::setSos2Constraint(
      SCIP *scip,
      int  c,
      SCIP_CONS *cons,
      SCIP_CONS** consSOS2
      )
{
   addOrdinaryConstraintName(c, cons);
   idxSos2ConsNames[nSos2Conss] = c;
   int nvars = SCIPgetNVarsSOS2(scip, cons);
   nSos2Coefs[nSos2Conss] = nvars;
   sos2Coefs[nSos2Conss] = new SCIP_Real[nvars];
   idxSos2CoefsVars[nSos2Conss] = new int[nvars];
   SCIP_Real *weights = SCIPgetWeightsSOS2(scip, cons);
   SCIP_VAR **vars = SCIPgetVarsSOS2(scip,cons);
   for( int v = 0; v < nvars; v++ )
   {
      sos2Coefs[nSos2Conss][v] = weights[v];
      idxSos2CoefsVars[nSos2Conss][v] = SCIPvarGetProbindex(vars[v]);
   }
   // store constraint
   consSOS2[nSos2Conss++] = cons;
}

void
ScipParaInstance::createSos2ConstraintsInSCIP(
      SCIP *scip
      )
{
   SCIP_VAR **vars = SCIPgetVars(scip);
   for(int c = 0; c < nSos2Conss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR **varsInCons = new SCIP_VAR*[nSos2Coefs[c]];
      for( int v = 0; v < nSos2Coefs[c]; v++)
      {
         varsInCons[v] = vars[idxSos2CoefsVars[c][v]];
      }
      SCIP_CALL_ABORT( SCIPcreateConsSOS2(scip, &cons,
            &consNames[posConsNames[idxSos2ConsNames[c]]], nSos2Coefs[c], varsInCons, sos2Coefs[c],
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE ) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIP_CALL_ABORT( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      delete [] varsInCons;
   }
}

/***********************************************************************************************
 *  getActiveVariables is copied from read_lp.c of SCIP code
 */
void
ScipParaInstance::getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( scalars != NULL );
   assert( nvars != NULL );
   assert( constant != NULL );

   if( transformed )
   {
      SCIP_CALL_ABORT( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if( requiredsize > *nvars )
      {
         *nvars = requiredsize;
         SCIP_CALL_ABORT( SCIPreallocBufferArray(scip, &vars, *nvars ) );
         SCIP_CALL_ABORT( SCIPreallocBufferArray(scip, &scalars, *nvars ) );

         SCIP_CALL_ABORT( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize, TRUE) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
         SCIP_CALL_ABORT( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );
   }
}

/***********************************************************************************************
 *  collectAggregatedVars is copied from read_lp.c of SCIP code
 */
void
ScipParaInstance::collectAggregatedVars(
   SCIP*              scip,               /**< SCIP data structure */
   int                nvars,              /**< number of mutable variables in the problem */
   SCIP_VAR**         vars,               /**< variable array */
   int*               nAggregatedVars,    /**< number of aggregated variables on output */
   SCIP_VAR***        aggregatedVars,     /**< array storing the aggregated variables on output */
   SCIP_HASHTABLE**   varAggregated       /**< hashtable for checking duplicates */
   )
{
   int j;

   /* check variables */
   for (j = 0; j < nvars; ++j)
   {
      SCIP_VARSTATUS status;
      SCIP_VAR* var;

      var = vars[j];
      status = SCIPvarGetStatus(var);

      /* collect aggregated variables in a list */
      if( status >= SCIP_VARSTATUS_AGGREGATED )
      {
         assert( status == SCIP_VARSTATUS_AGGREGATED ||
            status == SCIP_VARSTATUS_MULTAGGR ||
            status == SCIP_VARSTATUS_NEGATED );

         if ( ! SCIPhashtableExists(*varAggregated, (void*) var) )
         {
            (*aggregatedVars)[(*nAggregatedVars)++] = var;
            SCIP_CALL_ABORT( SCIPhashtableInsert(*varAggregated, (void*) var) );
         }
      }
   }
}

void
ScipParaInstance::setAggregatedConstraint(
      SCIP*                 scip,               /**< SCIP data structure */
      int                   c,                  /** aggregated constraint number */
      const char*           constName,          /**< constraint name */
      SCIP_VAR**            vars,               /**< array of variables */
      SCIP_Real*            vals,               /**< array of values */
      int                   nvars,              /**< number of variables */
      SCIP_Real             lhsAndrhs           /**< right hand side = left hand side */
      )
{
   posAggregatedVarNames[c] = lAggregatedVarNames;
   (void) strcpy( &aggregatedVarNames[lAggregatedVarNames],  SCIPvarGetName(vars[c]) );
   lAggregatedVarNames += strlen( SCIPvarGetName(vars[c]) );

   posAggregatedConsNames[c] = lAggregatedConsNames;
   (void) strcpy( &aggregatedConsNames[lAggregatedConsNames], constName );
   lAggregatedConsNames += strlen( constName )+ 1;

   aggregatedLhsAndLhss[c] = lhsAndrhs;
   nAggregatedCoefs[c] = nvars;
   int v;
   for(v = 0; v < nvars - 1; ++v)
   {
      aggregatedCoefs[nLinearConss][v] =  vals[v];
      idxLinearCoefsVars[nLinearConss][v] = SCIPvarGetProbindex(vars[v]);
   }
   aggregatedCoefs[nLinearConss][v] =  vals[v];
   idxLinearCoefsVars[nLinearConss][v] = c;
   nAggregatedConss++;
}

void
ScipParaInstance::setAggregatedConstrains(
      SCIP*              scip,               /**< SCIP data structure */
      int                nvars,              /**< number of mutable variables in the problem */
      int                nAggregatedVars,    /**< number of aggregated variables */
      SCIP_VAR**         aggregatedVars      /**< array storing the aggregated variables */
      )
{
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;
   char consname[UG::LpMaxNamelen];

   assert( scip != NULL );

   /* write aggregation constraints */
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &activevars, nvars) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &activevals, nvars) );

   /* compute lengths of aggregatedVarNames and aggregatedConsNames areas */
   lAggregatedVarNames = 0;
   lAggregatedConsNames = 0;
   for( int j = 0; j < nAggregatedVars; j++ )
   {
      lAggregatedVarNames += strlen(SCIPvarGetName(aggregatedVars[j])) + 1;
      (void) SCIPsnprintf(consname, UG::LpMaxNamelen, "aggr_%s", SCIPvarGetName(aggregatedVars[j]));
      lAggregatedConsNames += strlen(consname) + 1;
   }

   if( nAggregatedVars  )
   {
      /* allocate aggregatedVarNames and aggregatedConsNames areas */
      aggregatedVarNames = new char[lAggregatedVarNames];
      aggregatedConsNames = new char[lAggregatedConsNames];
   }

   /* set aggregated constraints */
   nAggregatedConss = 0;
   lAggregatedVarNames = 0;
   lAggregatedConsNames = 0;
   for (int j = 0; j < nAggregatedVars; ++j)
   {
      /* set up list to obtain substitution variables */
      nactivevars = 1;

      activevars[0] = aggregatedVars[j];
      activevals[0] = 1.0;
      activeconstant = 0.0;

      /* retransform given variables to active variables */
      getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, TRUE);

      activevals[nactivevars] = -1.0;
      activevars[nactivevars] = aggregatedVars[j];
      ++nactivevars;

      /* set constraint */
      (void) SCIPsnprintf(consname, UG::LpMaxNamelen, "aggr_%s", SCIPvarGetName(aggregatedVars[j]));
      setAggregatedConstraint(scip, j, consname, activevars, activevals, nactivevars, - activeconstant );
   }
   assert(nAggregatedConss == nAggregatedVars);

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);
}

void
ScipParaInstance::createAggregatedVarsAndConstrainsInSCIP(
      SCIP *scip
      )
{
   SCIP_VAR **vars = SCIPgetVars(scip);
   for(int c = 0; c < nAggregatedConss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR **varsInCons = new SCIP_VAR*[nAggregatedCoefs[c]];
      int v;
      for( v = 0; v < nAggregatedCoefs[c] - 1; v++)
      {
         varsInCons[v] = vars[idxAggregatedCoefsVars[c][v]];
      }
      SCIP_VAR* newvar;
      /* create new variable of the given name */
      SCIP_CALL_ABORT( SCIPcreateVar(scip, &newvar, &aggregatedVarNames[posAggregatedVarNames[c]], 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
            TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL_ABORT( SCIPaddVar(scip, newvar) );
      /* because the variable was added to the problem, it is captured by SCIP and we can safely release it right now
       * without making the returned *var invalid
       */
      varsInCons[v] = newvar;
      assert( SCIPvarGetProbindex(newvar) == (nVars + c) );
      SCIP_CALL_ABORT( SCIPreleaseVar(scip, &newvar) );

      SCIP_CALL_ABORT( SCIPcreateConsLinear(scip, &cons,
            &aggregatedConsNames[posAggregatedConsNames[c]], nAggregatedCoefs[c], varsInCons, aggregatedCoefs[c],aggregatedLhsAndLhss[c], aggregatedLhsAndLhss[c],
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL_ABORT( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIP_CALL_ABORT( SCIPprintCons(scip, cons, NULL) ) );
      SCIP_CALL_ABORT( SCIPreleaseCons(scip, &cons) );
      delete [] varsInCons;
   }
}

/** constractor : only called from ScipInitiator */
ScipParaInstance::ScipParaInstance(
      SCIP *scip,
      int  method
      ) : lProbName(0), probName(0), origObjSense(0), objScale(0.0), objOffset(0.0), nVars(0),
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
      nVarBoundDisjunctionConss(0), idxBoundDisjunctionConsNames(0), nVarsBoundDisjunction(0),
      idxVarBoundDisjunction(0), boundTypesBoundDisjunction(0), boundsBoundDisjunction(0),
      nSos1Conss(0),idxSos1ConsNames(0), nSos1Coefs(0), sos1Coefs(0), idxSos1CoefsVars(0),
      nSos2Conss(0), idxSos2ConsNames(0), nSos2Coefs(0), sos2Coefs(0), idxSos2CoefsVars(0),
      nAggregatedConss(0), lAggregatedVarNames(0), aggregatedVarNames(0), posAggregatedVarNames(0),
      lAggregatedConsNames(0), aggregatedConsNames(0), posAggregatedConsNames(0),
      aggregatedLhsAndLhss(0), nAggregatedCoefs(0), aggregatedCoefs(0), idxAggregatedCoefsVars(0)
{
   assert( scip != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVED ||
         SCIPgetStage(scip) == SCIP_STAGE_SOLVING ||
         SCIPgetStage(scip) == SCIP_STAGE_SOLVED );   // in case that presolving or root node solving solved the problem

   lProbName = strlen(SCIPgetProbName(scip));
   probName = new char[lProbName+1];
   strcpy(probName, SCIPgetProbName(scip));

   /** set objsen */
   origObjSense = SCIPgetObjsense(scip);       // always minimization problem
   objScale = SCIPgetTransObjscale(scip);   // Stefan is making now
   objOffset = SCIPgetTransObjoffset(scip); // Stefan is making now

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVED ) return;

   if( method == 0 || method == 2 )
   {
      nVars = SCIPgetNVars(scip);
      SCIP_VAR **vars = SCIPgetVars(scip);

      /* make varName and objCoefs and ovnm */
      posVarNames = new int[nVars];
      objCoefs = new SCIP_Real[nVars];
      // mapToOriginalIndecies = new int[nVars];

      lVarNames = 0;
      for(int v = 0; v < nVars; ++v)
      {
         posVarNames[v] = lVarNames;
         objCoefs[v] = SCIPvarGetObj(vars[v]);
         assert(SCIPvarGetProbindex(vars[v])!=-1);
         assert(SCIPvarGetProbindex(vars[v]) == v);
         lVarNames += strlen(SCIPvarGetName(vars[v])) + 1;
      }
      varNames = new char[lVarNames];
      varLbs = new SCIP_Real[nVars];
      varUbs = new SCIP_Real[nVars];
      varTypes = new int[nVars];
      for(int v = 0; v < nVars; ++v )
      {
         SCIP_VAR *var = vars[v];
         strcpy (&varNames[posVarNames[v]], SCIPvarGetName(var) );
         varLbs[SCIPvarGetProbindex(var)] = SCIPvarGetLbLocal(var); //* we should use global?
         varUbs[SCIPvarGetProbindex(var)] = SCIPvarGetUbLocal(var); //* we should use global?
         varTypes[SCIPvarGetProbindex(var)] = SCIPvarGetType(var);
      }
      if( method == 2 )
      {
         return;  // Solver reads given instance file.
      }

      /* make constraints */
      nConss = SCIPgetNConss(scip);
      SCIP_CONS **conss = SCIPgetConss(scip);

      posConsNames = new int[nConss];
      int tempSizeConsNames = 0;
      /* make consNames map */
      for( int c = 0; c < nConss; c++ )
      {
         tempSizeConsNames += strlen(SCIPconsGetName(conss[c])) + 1;
      }
      consNames = new char[tempSizeConsNames];

      SCIP_CONS** consSOS1;
      SCIP_CONS** consSOS2;
      nSos1Conss = 0;
      nSos2Conss = 0;

      /** collect SOS constraints in array for later output */
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &consSOS1, nConss) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &consSOS2, nConss) );

      SCIP_CONS *cons;
      const char* conshdlrname;
      SCIP_CONSHDLR* conshdlr;

      /** count number of each constraint */
      for (int c = 0; c < nConss; ++c)
      {
         cons = conss[c];
         assert( cons != NULL);

         /* in case the transformed is written only constraint are posted which are enabled in the current node */
         if( !SCIPconsIsEnabled(cons) )
            continue;

         conshdlr = SCIPconsGetHdlr(cons);
         assert( conshdlr != NULL );

         conshdlrname = SCIPconshdlrGetName(conshdlr);
         assert( SCIPconsIsTransformed(cons) );

         if( strcmp(conshdlrname, "linear") == 0 )
         {
            nLinearConss++;
         }
         else if( strcmp(conshdlrname, "setppc") == 0 )
         {
            nSetppcConss++;
         }
         else if ( strcmp(conshdlrname, "logicor") == 0 )
         {
            nLogicorConss++;
         }
         else if ( strcmp(conshdlrname, "knapsack") == 0 )
         {
            nKnapsackConss++;
         }
         else if ( strcmp(conshdlrname, "varbound") == 0 )
         {
            nVarboundConss++;
         }
         else if (  strcmp(conshdlrname, "bounddisjunction") == 0 )
         {
            nVarBoundDisjunctionConss++;
         }
         else if ( strcmp(conshdlrname, "SOS1") == 0 )
         {
            nSos1Conss++;
         }
         else if ( strcmp(conshdlrname, "SOS2") == 0 )
         {
            nSos2Conss++;
         }
         else
         {
            THROW_LOGICAL_ERROR3("constraint handler <", conshdlrname, "> can not print requested format");
         }
      }

      assert(nConss == (nLinearConss + nSetppcConss + nLogicorConss + nKnapsackConss + nVarboundConss + nVarBoundDisjunctionConss + nSos1Conss + nSos2Conss ) );

      allocateMemoryForOrdinaryConstraints();

      /** re-initialize counters for ordinary constraints */
      nLinearConss = 0;
      nSetppcConss = 0;
      nLogicorConss = 0;
      nKnapsackConss = 0;
      nVarboundConss = 0;
      nVarBoundDisjunctionConss = 0;
      nSos1Conss = 0;
      nSos2Conss = 0;

      /** initialize length of constraint names area */
      lConsNames = 0;
      for (int c = 0; c < nConss; ++c)
      {
         cons = conss[c];

         /* in case the transformed is written only constraint are posted which are enabled in the current node */
         if( !SCIPconsIsEnabled(cons) )
            continue;

         conshdlr = SCIPconsGetHdlr(cons);
         conshdlrname = SCIPconshdlrGetName(conshdlr);

         if( strcmp(conshdlrname, "linear") == 0 )
         {
            setLinearConstraint( scip, c, cons );
         }
         else if( strcmp(conshdlrname, "setppc") == 0 )
         {
            setSetppcConstraint( scip, c, cons );
         }
         else if ( strcmp(conshdlrname, "logicor") == 0 )
         {
            setLogicorConstraint( scip, c, cons );
         }
         else if ( strcmp(conshdlrname, "knapsack") == 0 )
         {
            setKnapsackConstraint( scip, c, cons );
         }
         else if ( strcmp(conshdlrname, "varbound") == 0 )
         {
            setVarboundConstraint( scip, c, cons );
         }
         else if ( strcmp(conshdlrname, "bounddisjunction") == 0 )
         {
            setBoundDisjunctionConstraint( scip, c, cons );
         }
         else if ( strcmp(conshdlrname, "SOS1") == 0 )
         {
            setSos1Constraint( scip, c, cons, consSOS1);
         }
         else if ( strcmp(conshdlrname, "SOS2") == 0 )
         {
            setSos2Constraint( scip, c, cons, consSOS2);
         }
      }

      assert(nConss == (nLinearConss + nSetppcConss + nLogicorConss + nKnapsackConss + nVarboundConss + nVarBoundDisjunctionConss + nSos1Conss + nSos2Conss ) );

      SCIP_VAR** aggregatedVars;
      int nAggregatedVars = 0;
      SCIP_HASHTABLE* varAggregated;

      // create hashtable for storing aggregated variables
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &aggregatedVars, nVars) );
      SCIP_CALL_ABORT( SCIPhashtableCreate(&varAggregated, SCIPblkmem(scip), 1000, hashGetKeyVar, hashKeyEqVar, hashKeyValVar,
            NULL) );

      // check for aggregated variables in SOS1 constraints and output aggregations as linear constraints
      for (int c = 0; c < nSos1Conss; ++c)
      {
         cons = consSOS1[c];
         SCIP_VAR **consvars = SCIPgetVarsSOS1(scip, cons);
         int nconsvars = SCIPgetNVarsSOS1(scip, cons);

         collectAggregatedVars(scip, nconsvars, consvars, &nAggregatedVars, &aggregatedVars, &varAggregated);
      }

      // check for aggregated variables in SOS2 constraints and output aggregations as linear constraints
      for (int c = 0; c < nSos2Conss; ++c)
      {
         cons = consSOS2[c];
         SCIP_VAR **consvars = SCIPgetVarsSOS2(scip, cons);
         int nconsvars = SCIPgetNVarsSOS2(scip, cons);

         collectAggregatedVars(scip, nconsvars, consvars, &nAggregatedVars, &aggregatedVars, &varAggregated);
      }

      // set aggregation constraints
      setAggregatedConstrains(scip, nVars, nAggregatedVars, aggregatedVars );

      // free space
      SCIPfreeBufferArray(scip, &aggregatedVars);
      SCIPhashtableFree(&varAggregated);

      // free space
      SCIPfreeBufferArray(scip, &consSOS1);
      SCIPfreeBufferArray(scip, &consSOS2);
   }
   else
   {
      SCIP_CALL_ABORT( SCIPwriteTransProblem(scip, PRESOLVED_INSTANCE, "cip", TRUE ) );
   }
}

ScipParaInstance::~ScipParaInstance(
      )
{
   if(probName) delete[] probName;
   if(varLbs)  delete[] varLbs;
   if(varUbs) delete[] varUbs;
   if(objCoefs) delete[] objCoefs;
   if(varTypes) delete[] varTypes;
   if(varNames) delete[] varNames;
   if(posVarNames) delete[] posVarNames;
   if(mapToOriginalIndecies) delete[] mapToOriginalIndecies;
   if(consNames) delete[] consNames;
   if(posConsNames) delete[] posConsNames;

   /** for linear constraints */
   if( nLinearConss > 0 )
   {
      if( idxLinearConsNames ) delete [] idxLinearConsNames;
      if( linearLhss ) delete [] linearLhss;
      if( linearRhss ) delete [] linearRhss;
      if( nLinearCoefs )
      {
         for(int i = 0; i < nLinearConss; i++ )
         {
            delete [] linearCoefs[i];
            delete [] idxLinearCoefsVars[i];
         }
         delete [] linearCoefs;
         delete [] idxLinearCoefsVars;
         delete [] nLinearCoefs;
      }
   }
   /** for setppc constraints */
   if( nSetppcConss )
   {
      if( idxSetppcConsNames ) delete [] idxSetppcConsNames;
      if( nIdxSetppcVars ) delete [] nIdxSetppcVars;
      if( setppcTypes ) delete [] setppcTypes;
      if( idxSetppcVars )
      {
         for(int i = 0; i < nSetppcConss; i++ )
         {
            delete idxSetppcVars[i];
         }
         delete [] idxSetppcVars;
      }
   }
   /** for logical constraints */
   if(nLogicorConss)
   {
      if( idxLogicorConsNames) delete [] idxLogicorConsNames;
      if( nIdxLogicorVars ) delete [] nIdxLogicorVars;
      if( idxLogicorVars )
      {
         for( int i = 0; i < nLogicorConss; i++ )
         {
            delete [] idxLogicorVars[i];
         }
         delete [] idxLogicorVars;
      }
   }
   /** for knapsack constraints */
   if(nKnapsackConss)
   {
      if( idxKnapsackConsNames ) delete [] idxKnapsackConsNames;
      if( capacities ) delete [] capacities;
      if( nLKnapsackCoefs ) delete [] nLKnapsackCoefs;
      if( knapsackCoefs )
      {
         for( int i = 0; i < nKnapsackConss; i++ )
         {
            delete [] knapsackCoefs[i];
            delete [] idxKnapsackCoefsVars[i];
         }
         delete [] knapsackCoefs;
         delete [] idxKnapsackCoefsVars;
      }
   }
   /** for varbound constraints */
   if( nVarboundConss )
   {
      if( idxVarboundConsNames ) delete [] idxVarboundConsNames;
      if( idxVarboundCoefVar1s ) delete [] idxVarboundCoefVar1s;
      if( varboundCoef2s ) delete [] varboundCoef2s;
      if( idxVarboundCoefVar2s ) delete [] idxVarboundCoefVar2s;
   }

   /** for bounddisjunction constraints */
   if( nVarBoundDisjunctionConss )
   {
      if( idxBoundDisjunctionConsNames ) delete [] idxBoundDisjunctionConsNames;
      if( nVarsBoundDisjunction ) delete[] nVarsBoundDisjunction;
      for( int c = 0; c < nVarBoundDisjunctionConss; c++ )
      {
         if( idxVarBoundDisjunction[c] ) delete[] idxVarBoundDisjunction[c];
         if( boundTypesBoundDisjunction[c] ) delete[] boundTypesBoundDisjunction[c];
         if( boundsBoundDisjunction[c] ) delete[] boundsBoundDisjunction[c];
      }
   }

   /** for SOS1 constraints */
   if( nSos1Conss )
   {
      if( nSos1Coefs ) delete [] nSos1Coefs;
      if( sos1Coefs )
      {
         for( int i = 0; i < nSos1Conss; i++ )
         {
            delete [] sos1Coefs[i];
            delete [] idxSos1CoefsVars[i];
         }
         delete [] sos1Coefs;
         delete [] idxSos1CoefsVars;
      }
   }
   /** for SOS2 constraints */
   if( nSos2Conss )
   {
      if( nSos2Coefs ) delete [] nSos2Coefs;
      if( sos2Coefs )
      {
         for( int i = 0; i < nSos1Conss; i++ )
         {
            delete [] sos2Coefs[i];
            delete [] idxSos2CoefsVars[i];
         }
         delete sos2Coefs;
         delete idxSos2CoefsVars;
      }
   }
   /** for agrregated constraints */
   if( nAggregatedConss )
   {
      for( int i = 0; i < nAggregatedConss; i++ )
      {
         delete aggregatedCoefs[i];
         delete idxAggregatedCoefsVars[i];
      }
      delete [] aggregatedCoefs;
      delete [] idxAggregatedCoefsVars;
      aggregatedCoefs = 0;
      idxAggregatedCoefsVars = 0;

      if( aggregatedVarNames ) delete [] aggregatedVarNames;
      if( posAggregatedVarNames ) delete[] posAggregatedVarNames;
      if( aggregatedConsNames ) delete [] aggregatedConsNames;
      if( posAggregatedConsNames ) delete[] posAggregatedConsNames;
      if( aggregatedLhsAndLhss ) delete[] aggregatedLhsAndLhss;
   }
}

/** create presolved problem instance that is solved by ParaSCIP */
void
ScipParaInstance::createProblem(
      SCIP *scip,
      int method,
      bool noPreprocessingInLC,
      char *settingsNameLC
      )
{
   SCIP_CALL_ABORT( SCIPcreateProb(scip, probName, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   switch ( method )
   {
   case 0 :
   {
      SCIP_CALL_ABORT( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
      for(int v = 0; v < nVars; ++v)
       {
          SCIP_VAR* newvar;

          /* create new variable of the given name */
          SCIP_CALL_ABORT( SCIPcreateVar(scip, &newvar, &varNames[posVarNames[v]], varLbs[v], varUbs[v], objCoefs[v],
                SCIP_Vartype(varTypes[v]),
                TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
          SCIP_CALL_ABORT( SCIPaddVar(scip, newvar) );
          assert( SCIPvarGetProbindex(newvar) == v );
          /* because the variable was added to the problem, it is captured by SCIP and we can safely release it right now
           * without making the returned *var invalid
           */
          SCIP_CALL_ABORT( SCIPreleaseVar(scip, &newvar) );
       }

       if( nLinearConss )
          createLinearConstraintsInSCIP( scip );
       if( nSetppcConss)
          createSetppcConstraintsInSCIP( scip );
       if( nLogicorConss )
          createLogicorConstraintsInSCIP( scip );
       if( nKnapsackConss )
          createKnapsackConstraintsInSCIP( scip );
       if( nVarboundConss )
          createVarboundConstraintsInSCIP( scip );
       if( nVarBoundDisjunctionConss )
          createBoundDisjunctionConstraintInSCIP( scip );
       if( nSos1Conss )
          createSos1ConstraintsInSCIP( scip );
       if( nSos2Conss )
          createSos2ConstraintsInSCIP( scip );
       if( nAggregatedConss )
          createAggregatedVarsAndConstrainsInSCIP( scip );
       break;
   }
   case 1 :
   {
      SCIP_CALL_ABORT( SCIPreadProb(scip, PRESOLVED_INSTANCE, "cip" ) );
      nVars = SCIPgetNVars(scip);
      SCIP_VAR **vars = SCIPgetVars(scip);
      varLbs = new SCIP_Real[nVars];
      varUbs = new SCIP_Real[nVars];
      for(int v = 0; v < nVars; ++v )
      {
         SCIP_VAR *var = vars[v];
         varLbs[SCIPvarGetProbindex(var)] = SCIPvarGetLbLocal(var); //* we should use global?
         varUbs[SCIPvarGetProbindex(var)] = SCIPvarGetUbLocal(var); //* we should use global?
      }
      break;
   }
   case 2 :
   {
      SCIP *tempScip = 0;
      SCIP_Bool success = TRUE;
      SCIP_CALL_ABORT( SCIPcreate(&tempScip) );

      /* include default SCIP plugins */
      SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(tempScip) );
      /** user include plugins */
      includeUserPlugins(tempScip);

#ifdef SCIP_THREADSAFE_MESSAGEHDLRS
      SCIPsetMessagehdlrQuiet(tempScip, TRUE);
#endif
      SCIP_CALL_ABORT( SCIPreadProb(tempScip, getFileName(), NULL ) );

      /* presolve problem */
      if( noPreprocessingInLC )
      {
         SCIP_CALL_ABORT( SCIPsetIntParam(tempScip, "presolving/maxrounds", 0));
      }
      else
      {
         if( settingsNameLC )
         {
            SCIP_CALL_ABORT( SCIPreadParams(tempScip, settingsNameLC) );
         }
         /*
         else
         {
            SCIP_CALL_ABORT( SCIPsetPresolving(tempScip, SCIP_PARAMSETTING_FAST, TRUE) );
         }
         */
      }
      // SCIP_CALL_ABORT( SCIPsetIntParam(tempScip, "constraints/quadratic/replacebinaryprod", 0));
      // SCIP_CALL_ABORT( SCIPsetBoolParam(tempScip, "constraints/nonlinear/reformulate", FALSE));
      SCIP_CALL_ABORT( SCIPpresolve(tempScip) );
      SCIP_STATUS scipStatus = SCIPgetStatus(tempScip);
      if( scipStatus == SCIP_STATUS_OPTIMAL ) return;  // the problem shold be solved in LC

      SCIP_CALL_ABORT( SCIPcopyParamSettings(tempScip, scip) );

      /* create the variable mapping hash map */
      SCIP_HASHMAP* varmap;
      SCIP_CALL_ABORT( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), SCIPgetNVars(tempScip)) );
      SCIP_HASHMAP* conssmap;
      SCIP_CALL_ABORT( SCIPhashmapCreate(&conssmap, SCIPblkmem(scip), SCIPgetNConss(tempScip)) );

      /* copy all variables and constraints */
      SCIP_CALL_ABORT( SCIPcopyVars(tempScip, scip, varmap, conssmap, TRUE) );
      SCIP_CALL_ABORT( SCIPcopyConss(tempScip, scip, varmap, conssmap, TRUE, FALSE, &success) );

      /* free hash map */
      SCIPhashmapFree(&varmap);
      SCIPhashmapFree(&conssmap);

      SCIPfree(&tempScip);

      nVars = SCIPgetNVars(scip);
      SCIP_VAR **vars = SCIPgetVars(scip);
      varLbs = new SCIP_Real[nVars];
      varUbs = new SCIP_Real[nVars];
      for(int v = 0; v < nVars; ++v )
      {
         SCIP_VAR *var = vars[v];
         varLbs[SCIPvarGetProbindex(var)] = SCIPvarGetLbLocal(var); //* we should use global?
         varUbs[SCIPvarGetProbindex(var)] = SCIPvarGetUbLocal(var); //* we should use global?
      }
      break;
   }
   default :
      THROW_LOGICAL_ERROR1("Undefined instance transfer method");
   }
}

void
ScipParaInstance::freeMemory()
{
   /** free memory for ParaInstance to save memory */
   if(varLbs)
   {
      delete[] varLbs;
      varLbs = 0;
   }
   if(varUbs)
   {
      delete[] varUbs;
      varUbs = 0;
   }
   if(objCoefs)
   {
      delete[] objCoefs;
      objCoefs = 0;
   }
   if(varTypes)
   {
      delete[] varTypes;
      varTypes = 0;
   }
   if(posVarNames)
   {
      delete[] posVarNames;
      posVarNames = 0;
   }

   if( nLinearConss )
   {
      delete [] idxLinearConsNames;
      idxLinearConsNames = 0;
      delete [] linearLhss;
      linearLhss = 0;
      delete [] linearRhss;
      linearRhss = 0;
      for(int i = 0; i < nLinearConss; i++ )
      {
         delete [] linearCoefs[i];
         delete [] idxLinearCoefsVars[i];
      }
      delete [] linearCoefs;
      delete [] idxLinearCoefsVars;
      delete [] nLinearCoefs;
      linearCoefs = 0;
      idxLinearCoefsVars = 0;
      nLinearCoefs = 0;
   }
   if( nSetppcConss )
   {
      delete [] idxSetppcConsNames;
      idxSetppcConsNames = 0;
      delete [] nIdxSetppcVars;
      nIdxSetppcVars = 0;
      delete [] setppcTypes;
      setppcTypes = 0;
      for(int i = 0; i < nSetppcConss; i++ )
      {
         delete idxSetppcVars[i];
      }
      delete [] idxSetppcVars;
      idxSetppcVars = 0;
   }
   if( nLogicorConss )
   {
      delete [] idxLogicorConsNames;
      idxLogicorConsNames = 0;
      delete [] nIdxLogicorVars;
      nIdxLogicorVars = 0;
      for( int i = 0; i < nLogicorConss; i++ )
      {
         delete [] idxLogicorVars[i];
      }
      delete [] idxLogicorVars;
      idxLogicorVars = 0;
   }
   if( nKnapsackConss )
   {
      delete [] idxKnapsackConsNames;
      idxKnapsackConsNames = 0;
      delete [] capacities;
      capacities = 0;
      delete [] nLKnapsackCoefs;
      nLKnapsackCoefs = 0;
      for( int i = 0; i < nKnapsackConss; i++ )
      {
         delete [] knapsackCoefs[i];
         delete [] idxKnapsackCoefsVars[i];
      }
      delete [] knapsackCoefs;
      delete [] idxKnapsackCoefsVars;
      knapsackCoefs = 0;
      idxKnapsackCoefsVars = 0;
   }
   if( nSos1Conss )
   {
      delete [] nSos1Coefs;
      nSos1Coefs = 0;
      for( int i = 0; i < nSos1Conss; i++ )
      {
         delete [] sos1Coefs[i];
         delete [] idxSos1CoefsVars[i];
      }
      delete [] sos1Coefs;
      delete [] idxSos1CoefsVars;
      sos1Coefs = 0;
      idxSos1CoefsVars = 0;
   }
   if( nSos2Conss )
   {
      delete [] nSos2Coefs;
      nSos2Coefs = 0;
      for( int i = 0; i < nSos1Conss; i++ )
      {
         delete [] sos2Coefs[i];
         delete [] idxSos2CoefsVars[i];
      }
      delete sos2Coefs;
      delete idxSos2CoefsVars;
      sos2Coefs = 0;
      idxSos2CoefsVars = 0;
   }
   if( nAggregatedConss )
   {
      for( int i = 0; i < nAggregatedConss; i++ )
      {
         delete aggregatedCoefs[i];
         delete idxAggregatedCoefsVars[i];
      }
      delete [] aggregatedCoefs;
      delete [] idxAggregatedCoefsVars;
      delete[] aggregatedLhsAndLhss;
      aggregatedCoefs = 0;
      idxAggregatedCoefsVars = 0;
      aggregatedLhsAndLhss = 0;
   }
}

/** stringfy ScipParaInstance: for debug */
const std::string
ScipParaInstance::toString(
      )
{
   std::ostringstream s;

   s << "lProbName = " <<  lProbName << std::endl;
   s << "probName = " << probName << std::endl;
   s << "objSense = " <<  origObjSense << std::endl;
   s << "objScale = " <<  objScale << std::endl;
   s << "objOffset = " << objOffset << std::endl;
   s << "nVars = " <<  nVars << ", lVarNames = " << lVarNames << std::endl;
   for( int i = 0; i < nVars; i++ )
   {
      s << "[" << i << "]: varNames = " << &varNames[posVarNames[i]]
        << ", varsLbs = " << varLbs[i] << ", varsUbs = " << varUbs[i]
        << ", objCoefs = " << objCoefs[i] << ", varTypes = " << varTypes[i];
      if( mapToOriginalIndecies )
      {
         s << ", mapToOriginalIndecies = " << mapToOriginalIndecies[i];
      }
      s << std::endl;
   }
   s << "nConss = " << nConss << ", lConsNames = " << lConsNames << std::endl;
   for( int i = 0; i < nConss; i ++ )
   {
      s<< "[" << i << "]: consNames = " << &consNames[posConsNames[i]] << std::endl;
   }

   /*************************
    * for linear constrains *
    * ***********************/
   s << "nLinearConss = " << nLinearConss << std::endl;
   for(int i = 0; i < nLinearConss; i++ )
   {
      s << "[" << i << "]: idxLinearConsNames = " << idxLinearConsNames[i]
        << ", linearLhss = " << linearLhss[i] << ", linearRhss = " << linearRhss[i]
        << ", nLinearCoefs = " << nLinearCoefs[i] << std::endl;
      for( int j = 0; j < nLinearCoefs[i]; j++ )
      {
         s << "    [" << j << "]: linearCoefs = " << linearCoefs[i][j]
           << ", idxLinearCoefsVars = " << idxLinearCoefsVars[i][j] << std::endl;
      }
   }

   /*************************
    * for setppc constrains *
    * ***********************/
   s << "nSetppcConss = " << nSetppcConss << std::endl;
   for(int i = 0; i < nSetppcConss; i++)
   {
      s << "[" << i << "]: idxSetppcConsNames = " << idxSetppcConsNames[i]
        << ", nIdxSetppcVars = " << nIdxSetppcVars[i] << ", setppcTypes = " << setppcTypes[i] << std::endl;
      for(int j = 0; j < nIdxSetppcVars[i]; j++ )
      {
         s << "   [" << j << "]: idxSetppcVars = " << idxSetppcVars[i][j] << std::endl;
      }
   }

   /*********************************
    * for logical constrains        *
    *********************************/
   s << "nLogicorConss = " << nLogicorConss << std::endl;
   for(int i = 0; i < nLogicorConss; i++ )
   {
      s << "[" << i << "]: idxLogicorConsNames = " << idxLogicorConsNames[i]
        << ", nIdxLogicorVars = " << nIdxLogicorVars[i] << std::endl;
      for(int j = 0; j < nIdxLogicorVars[i]; j++ )
      {
         s << "   [" << j << "]: idxLogicorVars = " << idxLogicorVars[i][j] << std::endl;
      }
   }

   /*********************************
    * for knapsack constrains       *
    *********************************/
   s << "nKnapsackConss = " << nKnapsackConss << std::endl;
   for(int i = 0; i < nKnapsackConss; i++)
   {
      s << "[" << i << "]: idxKnapsackConsNames = " << idxKnapsackConsNames[i]
        << ", capacities = " << capacities[i] << ", nLKnapsackCoefs = " << nLKnapsackCoefs[i] << std::endl;
      for(int j = 0; j < nLKnapsackCoefs[i]; j++ )
      {
         s << "   [" << j << "]: knapsackCoefs = " << knapsackCoefs[i][j] << ", idxKnapsackCoefsVars = " << idxKnapsackCoefsVars[i][j] << std::endl;
      }
   }

   /**********************************
    * for varbound constrains       *
    *********************************/
   s << "nVarboundConss = " << nVarboundConss << std::endl;
   for(int i = 0; i < nVarboundConss; i++)
   {
      s << "[" << i << "]: idxVarboundConsNames = " << idxVarboundConsNames[i]
        << ", varboundLhss = " << varboundLhss[i] << ", varboundRhss = " << varboundRhss[i]
        << ", idxVarboundCoefVar1s = " << idxVarboundCoefVar1s[i] << ", varboundCoef2s = " << varboundCoef2s[i]
        << ", idxVarboundCoefVar2s = " << idxVarboundCoefVar2s[i] << std::endl;
   }

   /*********************************
    * for SOS1 constraints          *
    *********************************/
   s << "nSos1Conss = " << nSos1Conss << std::endl;
   for( int i = 0; i < nSos1Conss; i++ )
   {
      s <<  "[" << i << "]: idxSos1ConsNames = " << idxSos1ConsNames[i]
        << ", nSos1Coefs = " << nSos1Coefs[i] << std::endl;
      for(int j = 0; j < nSos1Coefs[i]; j++ )
      {
         s << "   [" << j << "]: sos1Coefs = " << sos1Coefs[i][j] << ", idxSos1CoefsVars = " << idxSos1CoefsVars[i][j] << std::endl;
      }
   }

   /*********************************
    * for SOS2 constraints          *
    *********************************/
   s << "nSos2Conss = " << nSos2Conss << std::endl;
   for( int i = 0; i < nSos2Conss; i++ )
   {
      s <<  "[" << i << "]: idxSos2ConsNames = " << idxSos2ConsNames[i]
        << ", nSos2Coefs = " << nSos2Coefs[i] << std::endl;
      for(int j = 0; j < nSos2Coefs[i]; j++ )
      {
         s << "   [" << j << "]: sos2Coefs = " << sos2Coefs[i][j] << ", idxSos2CoefsVars = " << idxSos2CoefsVars[i][j] << std::endl;
      }
   }

   /**********************************
    * for aggregated constrains      *
    *********************************/
   s << "nAggregatedConss = " << nAggregatedConss << ", lAggregatedVarNames = " << lAggregatedVarNames
     << ", lAggregatedConsNames = " << lAggregatedConsNames << std::endl;
   for(int i = 0; i < nAggregatedConss; i++)
   {
      s <<  "[" << i << "]: aggregatedVarNames = " << aggregatedVarNames[posAggregatedVarNames[i]]
        << ", aggregatedConsNames = " << aggregatedConsNames[posAggregatedConsNames[i]]
        << ", aggregatedLhsAndLhss = " << aggregatedLhsAndLhss[i] << ", nAggregatedCoefs = " << nAggregatedCoefs[i] << std::endl;
      for(int j = 0; j < nAggregatedCoefs[i]; j++ )
      {
         s << "   [" << j << "]: aggregatedCoefs = " << aggregatedCoefs[i][j]
           << ", idxAggregatedCoefsVars = " << idxAggregatedCoefsVars[i][j] << std::endl;
      }
   }
   return s.str();
}

