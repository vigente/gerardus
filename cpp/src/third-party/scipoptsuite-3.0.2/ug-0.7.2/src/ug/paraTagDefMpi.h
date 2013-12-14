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

/**@file    paraTagDefMpi.h
 * @brief   Tag definitions for MPI communication.
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_TAG_DEF_MPI_H__
#define __PARA_TAG_DEF_MPI_H__
#include "paraTagDef.h"


namespace UG
{
static const int TAG_MPI_FIRST            = TAG_LAST + 1;
//----------------------------------------------------------------
static const int TagSolverDiffParamSet1   = TAG_MPI_FIRST + 0;
static const int TagSolverDiffParamSet2   = TAG_MPI_FIRST + 1;
static const int TagDiffSubproblem1       = TAG_MPI_FIRST + 2;
static const int TagDiffSubproblem2       = TAG_MPI_FIRST + 3;
static const int TagSolution1             = TAG_MPI_FIRST + 4;
static const int TagInitialStat1          = TAG_MPI_FIRST + 5;
//---------------------------------------------------------------
static const int TAG_MPI_LAST             = TAG_MPI_FIRST + 5;
static const int N_MPI_TAGS               = TAG_MPI_LAST - TAG_FIRST + 1;
}

#endif // __PARA_TAG_DEF_MPI_H__
