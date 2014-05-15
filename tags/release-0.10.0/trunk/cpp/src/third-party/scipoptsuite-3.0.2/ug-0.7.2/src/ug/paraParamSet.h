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

/**@file    paraParamSet.h
 * @brief   Parameter set for UG framework.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_PARAM_SET_H__
#define __PARA_PARAM_SET_H__
#include <algorithm>
#include <string>
#include <iostream>
#include <map>
#include <cmath>
#include "paraComm.h"

#define OUTPUT_PARAM_VALUE_ERROR( msg1, msg2, msg3, msg4 ) \
   std::cout << "[PARAM VALUE ERROR] Param type = " << msg1 << ", Param name = " << msg2 \
     << ", Param value = " <<  msg3 <<  ": Param comment is as follows: " << std::endl \
     << msg4 << std::endl;  \
   return (-1)

namespace UG
{

/** Types of parameters */
static const int ParaParamTypeBool    = 1;     /**< bool values: true or false */
static const int ParaParamTypeInt     = 2;     /**< integer values             */
static const int ParaParamTypeLongint = 3;     /**< long integer values        */
static const int ParaParamTypeReal    = 4;     /**< real values                */
static const int ParaParamTypeChar    = 5;     /**< characters                 */
static const int ParaParamTypeString  = 6;     /**< arrays of characters       */

/** Bool parameters =====================================================*/
static const int ParaParamsFirst                    = 0;
static const int ParaParamsBoolFirst                = ParaParamsFirst;
//-------------------------------------------------------------------------
static const int Quiet                               = ParaParamsBoolFirst + 0;
static const int TagTrace                            = ParaParamsBoolFirst + 1;
static const int LogSolvingStatus                    = ParaParamsBoolFirst + 2;
static const int LogNodesTransfer                    = ParaParamsBoolFirst + 3;
static const int LogSubtreeInfo                      = ParaParamsBoolFirst + 4;
static const int OutputTabularSolvingStatus          = ParaParamsBoolFirst + 5;
static const int DeterministicTabularSolvingStatus   = ParaParamsBoolFirst + 6;
static const int RootNodeSolvabilityCheck            = ParaParamsBoolFirst + 7;
static const int UseRootNodeCuts                     = ParaParamsBoolFirst + 8;
static const int TransferLocalCuts                   = ParaParamsBoolFirst + 9;
static const int TransferConflicts                   = ParaParamsBoolFirst + 10;
static const int TransferBranchStats                 = ParaParamsBoolFirst + 11;
static const int CheckEffectOfRootNodePreprocesses   = ParaParamsBoolFirst + 12;
static const int Checkpoint                          = ParaParamsBoolFirst + 13;
static const int CollectOnce                         = ParaParamsBoolFirst + 14;
static const int ProvingRun                          = ParaParamsBoolFirst + 15;
static const int SetAllDefaultsAfterRacing           = ParaParamsBoolFirst + 16;
static const int DistributeBestPrimalSolution        = ParaParamsBoolFirst + 17;
static const int LightWeightRootNodeProcess          = ParaParamsBoolFirst + 18;
static const int RacingStatBranching                 = ParaParamsBoolFirst + 19;
static const int IterativeBreakDown                  = ParaParamsBoolFirst + 20;
static const int NoPreprocessingInLC                 = ParaParamsBoolFirst + 21;
static const int NoUpperBoundTransferInRacing        = ParaParamsBoolFirst + 22;
static const int MergeNodesAtRestart                 = ParaParamsBoolFirst + 23;
static const int NChangeIntoCollectingModeNSolvers   = ParaParamsBoolFirst + 24;
static const int Deterministic                       = ParaParamsBoolFirst + 25;
static const int EventWeightedDeterministic          = ParaParamsBoolFirst + 26;
static const int NoSolverPresolvingAtRoot            = ParaParamsBoolFirst + 27;
static const int NoSolverPresolvingAtRootDefaultSet  = ParaParamsBoolFirst + 28;
static const int NoAggressiveSeparatorInRacing       = ParaParamsBoolFirst + 29;
static const int StatisticsToStdout                  = ParaParamsBoolFirst + 30;
static const int AllBoundChangesTransfer             = ParaParamsBoolFirst + 31;
static const int NoAllBoundChangesTransferInRacing   = ParaParamsBoolFirst + 32;
static const int BreakFirstSubtree                   = ParaParamsBoolFirst + 33;
static const int InitialNodesGeneration              = ParaParamsBoolFirst + 34;
//-------------------------------------------------------------------------
static const int ParaParamsBoolLast                  = ParaParamsBoolFirst + 34;
static const int ParaParamsBoolN                     = ParaParamsBoolLast - ParaParamsBoolFirst + 1;
/** Int parameters ======================================================*/
static const int ParaParamsIntFirst                  = ParaParamsBoolLast  + 1;
//-------------------------------------------------------------------------
static const int OutputParaParams                    = ParaParamsIntFirst + 0;
static const int RampUpPhaseProcess                  = ParaParamsIntFirst + 1;
static const int NChangeIntoCollectingMode           = ParaParamsIntFirst + 2;
static const int NodeTransferMode                    = ParaParamsIntFirst + 3;
static const int AddDualBoundCons                    = ParaParamsIntFirst + 4;
static const int NotificationSynchronization         = ParaParamsIntFirst + 5;
static const int MinNumberOfCollectingModeSolvers    = ParaParamsIntFirst + 6;
static const int MaxNumberOfCollectingModeSolvers    = ParaParamsIntFirst + 7;
static const int SolverOrderInCollectingMode         = ParaParamsIntFirst + 8;
static const int RacingRampUpTerminationCriteria     = ParaParamsIntFirst + 9;
static const int StopRacingNumberOfNodesLeft         = ParaParamsIntFirst + 10;
static const int NumberOfInitialNodes                = ParaParamsIntFirst + 11;
static const int MaxNRacingParamSetSeed              = ParaParamsIntFirst + 12;
static const int TryNVariablegOrderInRacing          = ParaParamsIntFirst + 13;
static const int TryNBranchingOrderInRacing          = ParaParamsIntFirst + 14;
static const int NEvaluationSolversToStopRacing      = ParaParamsIntFirst + 15;
static const int NMaxCanditatesForCollecting         = ParaParamsIntFirst + 16;
static const int NSolverNodesStartBreaking           = ParaParamsIntFirst + 17;
static const int NStopBreaking                       = ParaParamsIntFirst + 18;
static const int NTransferLimitForBreaking           = ParaParamsIntFirst + 19;
static const int NStopSolvingMode                    = ParaParamsIntFirst + 20;
static const int NCollectOnce                        = ParaParamsIntFirst + 21;
static const int AggressivePresolveDepth             = ParaParamsIntFirst + 22;
static const int AggressivePresolveStopDepth         = ParaParamsIntFirst + 23;
static const int BigDualGapSubtreeHandling           = ParaParamsIntFirst + 24;
static const int InstanceTransferMethod              = ParaParamsIntFirst + 25;
static const int KeepNodesDepth                      = ParaParamsIntFirst + 26;
static const int UgCplexRunningMode                  = ParaParamsIntFirst + 27;
//-------------------------------------------------------------------------
static const int ParaParamsIntLast                   = ParaParamsIntFirst + 27;
static const int ParaParamsIntN                      = ParaParamsIntLast - ParaParamsIntFirst + 1;
/** Longint parameters ==================================================*/
static const int ParaParamsLongintFirst              = ParaParamsIntLast + 1;
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
static const int ParaParamsLongintLast               = ParaParamsLongintFirst - 1;  // No params -1
static const int ParaParamsLongintN                  = ParaParamsLongintLast - ParaParamsLongintFirst + 1;
/** Real parameters =======================================================*/
static const int ParaParamsRealFirst                 = ParaParamsLongintLast + 1;
//-------------------------------------------------------------------------
static const int MultiplierForCollectingMode           = ParaParamsRealFirst  + 0;
static const int MultiplierToDetermineThresholdValue   = ParaParamsRealFirst  + 1;
static const int BgapCollectingMode                    = ParaParamsRealFirst  + 2;
static const int MultiplierForBgapCollectingMode       = ParaParamsRealFirst  + 3;
static const int ABgapForSwitchingToBestSolver         = ParaParamsRealFirst  + 4;
static const int BgapStopSolvingMode                   = ParaParamsRealFirst  + 5;
static const int NotificationInterval                  = ParaParamsRealFirst  + 6;
static const int TimeLimit                             = ParaParamsRealFirst  + 7;
static const int CheckpointInterval                    = ParaParamsRealFirst  + 8;
static const int StopRacingTimeLimit                   = ParaParamsRealFirst  + 9;
static const int StopRacingTimeLimitMultiplier         = ParaParamsRealFirst  + 10;
static const int StopRacingNumberOfNodesLeftMultiplier = ParaParamsRealFirst  + 11;
static const int TimeToIncreaseCMS                     = ParaParamsRealFirst  + 12;
static const int TabularSolvingStatusInterval          = ParaParamsRealFirst  + 13;
static const int RatioToApplyLightWeightRootProcess    = ParaParamsRealFirst  + 14;
static const int MultiplierForBreakingTargetBound      = ParaParamsRealFirst  + 15;
static const int FixedVariablesRatioInMerging          = ParaParamsRealFirst  + 16;
static const int AllowableRegressionRatioInMerging     = ParaParamsRealFirst  + 17;
static const int CountingSolverRatioInRacing           = ParaParamsRealFirst  + 18;
static const int ProhibitCollectOnceMultiplier         = ParaParamsRealFirst  + 19;
static const int GeneratorRatio                        = ParaParamsRealFirst  + 20;
//-------------------------------------------------------------------------
static const int ParaParamsRealLast                  = ParaParamsRealFirst  + 20;
static const int ParaParamsRealN                     = ParaParamsRealLast - ParaParamsRealFirst + 1;
/** Char parameters ========================================================*/
static const int ParaParamsCharFirst                 = ParaParamsRealLast + 1;
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
static const int ParaParamsCharLast                  = ParaParamsCharFirst - 1;   // No params -1
static const int ParaParamsCharN                     = ParaParamsCharLast - ParaParamsCharFirst + 1;
/** String parameters ======================================================*/
static const int ParaParamsStringFirst               = ParaParamsCharLast      +1;
//-------------------------------------------------------------------------
static const int SolverSettingsForInitialPresolving  = ParaParamsStringFirst + 0;
static const int SolverSettingsAtRootNode            = ParaParamsStringFirst + 1;
static const int SolverSettingsExceptRootNode        = ParaParamsStringFirst + 2;
static const int SolverSettingsAtRacing              = ParaParamsStringFirst + 3;
static const int TempFilePath                        = ParaParamsStringFirst + 4;
static const int TagTraceFileName                    = ParaParamsStringFirst + 5;
static const int LogSolvingStatusFilePath            = ParaParamsStringFirst + 6;
static const int LogNodesTransferFilePath            = ParaParamsStringFirst + 7;
static const int SolutionFilePath                    = ParaParamsStringFirst + 8;
static const int CheckpointFilePath                  = ParaParamsStringFirst + 9;
//-------------------------------------------------------------------------
static const int ParaParamsStringLast                = ParaParamsStringFirst + 9;
static const int ParaParamsStringN                   = ParaParamsStringLast - ParaParamsStringFirst + 1;
static const int ParaParamsLast                      = ParaParamsStringLast;
static const int ParaParamsSize                      = ParaParamsLast + 1;

class ParaParam {
   const char *paramName;
   const char *comment;
public:
   ParaParam( const char *inParamName, const char *inComment ) : paramName(inParamName), comment(inComment) {}
   virtual ~ParaParam() {}
   const char *getParamName() const { return paramName;}
   const char *getComment() const {return comment; }
   virtual int getType() const = 0;
};

class ParaParamBool : public ParaParam {
   const bool defaultValue;
   bool       currentValue;
public:
   ParaParamBool( const char *name, const char *inComment, bool value ) : ParaParam(name, inComment), defaultValue(value), currentValue(value) {}
   ~ParaParamBool(){}
   int getType() const { return ParaParamTypeBool; }
   bool getDefaultValue() const { return defaultValue; }
   bool getValue() const { return currentValue; }
   void setDefaultValue() { currentValue = defaultValue; }
   void setValue(bool value){ currentValue = value; }
   bool isDefaultValue()  const { return defaultValue == currentValue; }
};

class ParaParamInt : public ParaParam {
   const int defaultValue;
   int       currentValue;
   const int minValue;
   const int maxValue;
public:
   ParaParamInt( const char *name, const char *inComment, int value, const int min, const int max )
      : ParaParam(name, inComment), defaultValue(value), currentValue(value), minValue(min), maxValue(max) {}
   ~ParaParamInt(){}
   int getType() const { return ParaParamTypeInt; }
   int getDefaultValue() const { return defaultValue; }
   int getValue() const { return currentValue; }
   void setDefaultValue() { currentValue = defaultValue; }
   void setValue(int value){ currentValue = value; }
   bool isDefaultValue() const { return defaultValue == currentValue; }
   int getMinValue() const { return minValue; }
   int getMaxValue() const { return maxValue; }
};

class ParaParamLongint : public ParaParam {
   const long long defaultValue;
   long long       currentValue;
   const long long minValue;
   const long long maxValue;
public:
   ParaParamLongint( const char *name, const char *inComment, long long value, const long long min, const long long max )
      : ParaParam(name, inComment), defaultValue(value), currentValue(value),  minValue(min), maxValue(max) {}
   ~ParaParamLongint(){}
   int getType() const { return ParaParamTypeLongint; }
   long long getDefaultValue() const { return defaultValue; }
   long long getValue() const { return currentValue; }
   void setDefaultValue() { currentValue = defaultValue; }
   void setValue(long long value){ currentValue = value; }
   bool isDefaultValue() const { return defaultValue == currentValue; }
   long long getMinValue() const { return minValue; }
   long long getMaxValue() const { return maxValue; }
};

class ParaParamReal : public ParaParam {
   const double defaultValue;
   double       currentValue;
   const double minValue;
   const double maxValue;
public:
   ParaParamReal( const char *name, const char *inComment, double value, const double min, const double max )
      : ParaParam(name, inComment), defaultValue(value), currentValue(value), minValue(min), maxValue(max) {}
   ~ParaParamReal(){}
   int getType() const { return ParaParamTypeReal; }
   double getDefaultValue() const { return defaultValue; }
   double getValue() const { return currentValue; }
   void setDefaultValue() { currentValue = defaultValue; }
   void setValue(double value){ currentValue = value; }
   bool isDefaultValue() const { return ( fabs( defaultValue - currentValue ) < 1e-20 ); }
   double getMinValue() const { return minValue; }
   double getMaxValue() const { return maxValue; }
};

class ParaParamChar : public ParaParam {
   const char defaultValue;
   char       currentValue;
   const char *allowedValues;
public:
   ParaParamChar( const char *name, const char *inComment, char value, const char *inAllowedValues)
      : ParaParam(name, inComment), defaultValue(value), currentValue(value), allowedValues(inAllowedValues){}
   ~ParaParamChar(){}
   int getType() const { return ParaParamTypeChar; }
   char getDefaultValue() const { return defaultValue; }
   char getValue() const { return currentValue; }
   void setDefaultValue() { currentValue = defaultValue; }
   void setValue(char value){ currentValue = value; }
   bool isDefaultValue() const { return defaultValue == currentValue; }
   const char *getAllowedValues() const { return allowedValues; }
};

class ParaParamString : public ParaParam {
   const char *defaultValue;
   const char *currentValue;
public:
   ParaParamString( const char *name, const char *inComment, const char *value ) : ParaParam(name, inComment), defaultValue(value), currentValue(value) {}
   ~ParaParamString() { if( currentValue != defaultValue ) delete [] currentValue; }
   int getType() const { return ParaParamTypeString; }
   const char *getDefaultValue() const { return defaultValue; }
   const char *getValue() const { return currentValue; }
   void setDefaultValue() { if( currentValue != defaultValue ) delete [] currentValue; currentValue = defaultValue; }
   void setValue(const char *value){ currentValue = value; }
   bool isDefaultValue() const { return ( std::string(defaultValue) == std::string(currentValue) ); }
};

class ParaComm;
class ParaParamSet {
protected:
   static ParaParam *paraParams[ParaParamsSize];

   int paramParaseBool( ParaParam *paraParam, char *valuestr );
   int paramParaseInt( ParaParam *paraParam, char *valuestr );
   int paramParaseLongint( ParaParam *paraParam, char *valuestr );
   int paramParaseReal( ParaParam *paraParam, char *valuestr );
   int paramParaseChar( ParaParam *paraParam, char *valuestr );
   int paramParaseString( ParaParam *paraParam, char *valuestr );
   int parameterParse(char *line, std::map<std::string, int> &mapStringToId);

public:
    ParaParamSet();
    virtual ~ParaParamSet();
    /** for bool params */
    bool getBoolParamValue(int param);
    void setBoolParamValue(int param, bool value);
    bool getBoolParamDefaultValue(int param);
    void setBoolParamDefaultValue(int param);
    bool isBoolParamDefaultValue(int param);
    /** for int params */
    int getIntParamValue(int param);
    void setIntParamValue(int param, int value);
    int getIntParamDefaultValue(int param);
    void setIntParamDefaultValue(int param);
    bool isIntParamDefaultValue(int param);
    /** for longint params */
    long long getLongintParamValue(int param);
    void setLongintParamValue(int param, long long value);
    long long getLongintParamDefaultValue(int param);
    void setLongintParamDefaultValue(int param);
    bool isLongintParamDefaultValue(int param);
    /** for real params */
    double getRealParamValue(int param);
    void setRealParamValue(int param, double value);
    double getRealParamDefaultValue(int param);
    void setRealParamDefaultValue(int param);
    bool isRealParamDefaultValue(int param);
    /** for char params */
    char getCharParamValue(int param);
    void setCharParamValue(int param, char value);
    char getCharParamDefaultValue(int param);
    void setCharParamDefaultValue(int param);
    bool isCharParamDefaultValue(int param);
    /** for string params */
    const char *getStringParamValue(int param);
    void setStringParamValue(int param, const char *value);
    const char *getStringParamDefaultValue(int param);
    void setStringParamDefaultValue(int param);
    bool isStringParamDefaultValue(int param);


    void read(ParaComm *comm, const char* filename);
    void write(std::ostream *os);

    virtual int bcast(ParaComm *comm, int root) = 0;
};

}

#endif  // __PARA_PARAM_SET_H__
