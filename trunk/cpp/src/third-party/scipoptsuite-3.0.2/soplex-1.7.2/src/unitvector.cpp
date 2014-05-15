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

#include <assert.h>

#include "unitvector.h"

namespace soplex
{

bool UnitVector::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (mem() != themem)
      return MSGinconsistent("UnitVector");
   if (mem() + 1 != &themem[1])
      return MSGinconsistent("UnitVector");
   if (size() != 1)
      return MSGinconsistent("UnitVector");
   if (max() != 1)
      return MSGinconsistent("UnitVector");

   return SVector::isConsistent();
#else
   return true;
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
