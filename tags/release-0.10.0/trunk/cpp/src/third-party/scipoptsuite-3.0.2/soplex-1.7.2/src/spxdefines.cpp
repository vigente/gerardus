/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxdefines.cpp
 * @brief Debugging, floating point type and parameter definitions.
 */
#include "spxdefines.h"
#include "assert.h"
#include "spxout.h"

namespace soplex
{

const Real infinity                 = DEFAULT_INFINITY;

Real Param::s_epsilon               = DEFAULT_EPS_ZERO;
Real Param::s_epsilon_factorization = DEFAULT_EPS_FACTOR;
Real Param::s_epsilon_update        = DEFAULT_EPS_UPDATE;
int  Param::s_verbose               = 1;

bool msginconsistent(const char* name, const char* file, int line)
{
   assert(name != 0);
   assert(file != 0);
   assert(line >= 0);

   MSG_ERROR( spxout << file << "(" << line << ") "
   << "Inconsistency detected in " << name << std::endl; )

   return 0;
}

void Param::setEpsilon(Real eps)
{
   s_epsilon = eps;
}

void Param::setEpsilonFactorization(Real eps)
{
   s_epsilon_factorization = eps;
}

void Param::setEpsilonUpdate(Real eps)
{
   s_epsilon_update = eps;
}

void Param::setVerbose(int p_verbose)
{
#ifndef DISABLE_VERBOSITY
   s_verbose = p_verbose;
#endif
}

} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------

