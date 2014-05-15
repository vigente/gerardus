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

/**@file    paraTagDef.h
 * @brief   Fundamental Tag definitions
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_TAG_DEF_H__
#define __PARA_TAG_DEF_H__

#define TAG_STR(tag) #tag

namespace UG
{
static const int TAG_FIRST = 0;
//------------------------------------------------------------------------------------------------
static const int TagNode                                = TAG_FIRST +  0;
static const int TagNodeReceived                        = TAG_FIRST +  1;
static const int TagDiffSubproblem                      = TAG_FIRST +  2;
static const int TagRampUp                              = TAG_FIRST +  3;
static const int TagRetryRampUp                         = TAG_FIRST +  4;
static const int TagSolution                            = TAG_FIRST +  5;
static const int TagIncumbentValue                      = TAG_FIRST +  6;
static const int TagGlobalBestDualBoundValueAtWarmStart = TAG_FIRST +  7;
static const int TagSolverState                         = TAG_FIRST +  8;
static const int TagCompletionOfCalculation             = TAG_FIRST +  9;
static const int TagAnotherNodeRequest                  = TAG_FIRST + 10;
static const int TagNoNodes                             = TAG_FIRST + 11;
static const int TagInCollectingMode                    = TAG_FIRST + 12;
static const int TagCollectAllNodes                     = TAG_FIRST + 13;
static const int TagOutCollectingMode                   = TAG_FIRST + 14;
static const int TagLCBestBoundValue                    = TAG_FIRST + 15;
static const int TagNotificationId                      = TAG_FIRST + 16;
static const int TagTerminateRequest                    = TAG_FIRST + 17;
static const int TagInterruptRequest                    = TAG_FIRST + 18;
static const int TagTerminated                          = TAG_FIRST + 19;
static const int TagRacingRampUpParamSet                = TAG_FIRST + 20;
static const int TagWinner                              = TAG_FIRST + 21;
static const int TagLightWeightRootNodeProcess          = TAG_FIRST + 22;
static const int TagBreaking                            = TAG_FIRST + 23;
static const int TagHardTimeLimit                       = TAG_FIRST + 24;
static const int TagInitialStat                         = TAG_FIRST + 25;
static const int TagAckCompletion                       = TAG_FIRST + 26;
static const int TagToken                               = TAG_FIRST + 27;
//-----------------------------------------------------------------------------------------------
static const int TAG_LAST                               = TAG_FIRST + 27;
static const int N_TAGS                                 = TAG_LAST - TAG_FIRST + 1;
}

#ifdef _COMM_MPI_WORLD
#include "paraTagDefMpi.h"
#endif
#ifdef _COMM_PTH
#include "paraTagDefPth.h"
#endif

#endif // __PARA_TAG_DEF_H__
