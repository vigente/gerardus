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

/**@file  mpqreal.cpp
 * @brief Wrapper for GMP types.
 */

#include <math.h>

#include "mpqreal.h"

namespace soplex
{
#ifdef SOPLEX_WITH_GMP

/// print MpqReal with limited floating point precision
std::ostream& operator<<(std::ostream& os, const MpqReal& q)
{
   os << mpf_class(q);
   return os;
}

/// cast MpqReal to Real
Real get_d(const MpqReal& q)
{
   return q.get_d();
}

#else

/// cast MpqReal to Real
Real get_d(const MpqReal& q)
{
   return q;
}

/// return maximal absolute value
MpqReal abs(const MpqReal& q)
{
   return fabs(q);
}

#endif
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
