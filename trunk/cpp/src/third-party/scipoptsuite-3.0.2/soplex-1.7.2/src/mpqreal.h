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

/**@file  mpqreal.h
 * @brief Wrapper for GMP types.
 */
#ifndef _MPQREAL_H_
#define _MPQREAL_H_

#include <math.h>
#include <iostream>

#include "spxdefines.h"



#ifdef SOPLEX_WITH_GMP
#include "gmp.h"
#include "gmpxx.h"
#endif

namespace soplex
{
/**@brief   Wrapper for GMP type mpq_class.
 * @ingroup Algebra
 *
 * We wrap mpq_class so that we can replace it by SoPlex's normal Real type if GMP is not available.
 */
#ifdef SOPLEX_WITH_GMP

/// If compiled with GMP support, MpqReal is defined as mpq_class.
typedef mpq_class MpqReal;

/// return whether MpqReal provides exact arithmetic
#define MpqRealIsExact() (true)

/// print MpqReal with limited floating point precision
std::ostream& operator<<(std::ostream& os, const MpqReal& q);

#else

/// If compiled without GMP support, MpqReal is defined as SoPlex's normal Real.
typedef Real MpqReal;

/// return whether MpqReal provides exact arithmetic
#define MpqRealIsExact() (false)

/// return maximal absolute value
MpqReal abs(const MpqReal& q);
#endif

/// cast MpqReal to Real
Real get_d(const MpqReal& q);

} // namespace soplex
#endif // _MPQREAL_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
