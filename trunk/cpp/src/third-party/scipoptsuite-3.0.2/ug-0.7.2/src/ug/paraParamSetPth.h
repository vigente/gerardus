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

/**@file    paraParamSetPth.h
 * @brief   
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_PARAM_SET_PTH_H__
#define __PARA_PARAM_SET_PTH_H__
#include "paraParamSet.h"

namespace UG
{

class ParaParamSetPth : public ParaParamSet
{
public:
    ParaParamSetPth(){}
    int bcast(ParaComm *comm, int root)
    {
       THROW_LOGICAL_ERROR1("bcast is called in ParaParamSetPth");
    }
    ~ParaParamSetPth(){}
};

}

#endif  // __PARA_PARAM_SET_PTH_H__
