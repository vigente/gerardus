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

/**@file    paraDef.h
 * @brief   Defines for UG Framework
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_DEF_H__
#define __PARA_DEF_H__
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <cfloat>

namespace UG
{

#define MULTIPLIER_FOR_NODE_PSEUDQUEUE 10
#define DEFAULT_NUM_EPSILON          1e-9  /**< default upper bound for floating points to be considered zero */
#define MINEPSILON                   1e-20  /**< minimum value for any numerical epsilon */

#define THROW_LOGICAL_ERROR1( msg1 ) \
   { \
   std::ostringstream s; \
   s << "[LOGICAL ERROR:" <<  __FILE__ << "] func = " \
     << __func__ << ", line = " << __LINE__ << " - " \
     << ( msg1 );  \
   throw std::logic_error( s.str() ); \
   }

#define THROW_LOGICAL_ERROR2( msg1, msg2 ) \
   { \
   std::ostringstream s; \
   s << "[LOGICAL ERROR:" <<  __FILE__ << "] func = " \
     << __func__ << ", line = " << __LINE__ << " - " \
     << ( msg1 ) << ( msg2 );  \
   throw std::logic_error( s.str() ); \
   }

#define THROW_LOGICAL_ERROR3( msg1, msg2, msg3 ) \
   { \
   std::ostringstream s; \
   s << "[LOGICAL ERROR:" <<  __FILE__ << "] func = " \
     << __func__ << ", line = " << __LINE__ << " - " \
     << ( msg1 ) << ( msg2 ) << ( msg3 );  \
   throw std::logic_error( s.str() ); \
   }

#define THROW_LOGICAL_ERROR4( msg1, msg2, msg3, msg4 ) \
   { \
   std::ostringstream s; \
   s << "[LOGICAL ERROR:" <<  __FILE__ <<  "] func = " \
     << __func__ << ", line = " << __LINE__ << " - " \
     << ( msg1 ) << ( msg2 ) << ( msg3 ) << ( msg4 );  \
   throw std::logic_error( s.str() ); \
   }

#define THROW_LOGICAL_ERROR5( msg1, msg2, msg3, msg4, msg5) \
   { \
   std::ostringstream s; \
   s << "[LOGICAL ERROR:" <<  __FILE__ << "] func = " \
     << __func__ << ", line = " << __LINE__ << " - " \
     << ( msg1 ) << ( msg2 ) << ( msg3 ) << ( msg4 ) << ( msg5 );  \
   throw std::logic_error( s.str() ); \
   }

#define THROW_LOGICAL_ERROR6( msg1, msg2, msg3, msg4, msg5, msg6) \
   { \
   std::ostringstream s; \
   s << "[LOGICAL ERROR:" <<  __FILE__ << "] func = " \
     << __func__ << ", line = " << __LINE__ << " - " \
     << ( msg1 ) << ( msg2 ) << ( msg3 ) << ( msg4 ) << ( msg5 ) << ( msg6 );  \
   throw std::logic_error( s.str() ); \
   }

#define REALABS(x)        (fabs(x))
#define EPSEQ(x,y,eps)    (REALABS((x)-(y)) <= (eps))
#define EPSLT(x,y,eps)    ((x)-(y) < -(eps))
#define EPSLE(x,y,eps)    ((x)-(y) <= (eps))
#define EPSGT(x,y,eps)    ((x)-(y) > (eps))
#define EPSGE(x,y,eps)    ((x)-(y) >= -(eps))
#define EPSZ(x,eps)       (REALABS(x) <= (eps))
#define EPSP(x,eps)       ((x) > (eps))
#define EPSN(x,eps)       ((x) < -(eps))
#define EPSFLOOR(x,eps)   (floor((x)+(eps)))
#define EPSCEIL(x,eps)    (ceil((x)-(eps)))
#define EPSFRAC(x,eps)    ((x)-EPSFLOOR(x,eps))
#define EPSISINT(x,eps)   (EPSFRAC(x,eps) <= (eps))

static const int MaxStrLen = 1024;
static const int LpMaxNamelen = 1024;

static const int CompTerminatedNormally           = 0;
static const int CompTerminatedByAnotherNode      = 1;
static const int CompTerminatedByInterruptRequest = 2;
static const int CompTerminatedInRacingStage      = 3;
static const int CompInterruptedInRacingStage     = 4;
static const int CompInterruptedInMerging         = 5;

}

#endif // __PARA_DEF_H__
