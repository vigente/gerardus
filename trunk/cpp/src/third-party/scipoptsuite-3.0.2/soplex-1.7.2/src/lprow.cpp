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

//#define DEBUGGING 1

#include <stdlib.h>
#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "lprow.h"
#include "exceptions.h"

namespace soplex
{
LPRow::Type LPRow::type() const
{
   if (rhs() >= infinity)
      return GREATER_EQUAL;
   if (lhs() <= -infinity)
      return LESS_EQUAL;
   if (lhs() == rhs())
      return EQUAL;
   return RANGE;
}

void LPRow::setType( LPRow::Type p_type )
{
   switch (p_type)
   {
   case LESS_EQUAL:
      left = -infinity;
      break;
   case EQUAL:
      if (lhs() > -infinity)
         right = lhs();
      else
         left = rhs();
      break;
   case GREATER_EQUAL:
      right = infinity;
      break;
   case RANGE:
      MSG_ERROR( spxout << "ELPROW01 RANGE not supported in LPRow::setType()" 
                        << std::endl; )
      throw SPxInternalCodeException("XLPROW01 This should never happen.");
   default:
      throw SPxInternalCodeException("XLPROW02 This should never happen.");
   }
}

Real LPRow::value() const
{
   assert(type() != RANGE);

   return (rhs() < infinity) ? rhs() : lhs();
}

LPRow::LPRow(const SVector& p_rowVector, LPRow::Type p_type, Real p_value)
   : vec(p_rowVector)
{
   switch (p_type)
   {
   case LESS_EQUAL:
      left = -infinity;
      right = p_value;
      break;
   case EQUAL:
      left = p_value;
      right = p_value;
      break;
   case GREATER_EQUAL:
      left = p_value;
      right = infinity;
      break;
   default:
      throw SPxInternalCodeException("XLPROW03 This should never happen.");
   }

   assert(isConsistent());
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
