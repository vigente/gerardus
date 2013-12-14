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

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "lprowset.h"
#include "exceptions.h"

namespace soplex
{

void LPRowSet::add(DataKey& p_key, Real p_lhs, const SVector& vector, Real p_rhs)
{
   SVSet::add(p_key, vector);
   if (num() > left.dim())
   {
      left.reDim(num());
      right.reDim(num());
   }
   left [num() - 1] = p_lhs;
   right[num() - 1] = p_rhs;
}

void LPRowSet::add(const LPRowSet& p_set)
{
   int i = num();

   SVSet::add(p_set);
   if (num() > left.dim())
   {
      left.reDim(num());
      right.reDim(num());
   }

   for (int j = 0; i < num(); ++i, ++j)
   {
      left [i] = p_set.lhs(j);
      right[i] = p_set.rhs(j);
   }
}

void LPRowSet::add(DataKey nkey[], const LPRowSet& p_set)
{
   int i = num();
 
   add(p_set);

   for (int j = 0; i < num(); ++i, ++j)
      nkey[j] = key(i);
}

SVector& LPRowSet::create(DataKey& nkey, int nonzeros, Real p_lhs, Real p_rhs)
{
   if (num() + 1 > left.dim())
   {
      left.reDim(num() + 1);
      right.reDim(num() + 1);
   }
   left [num()] = p_lhs;
   right[num()] = p_rhs;

   return *SVSet::create(nkey, nonzeros);
}


/*      \SubSection{Shrinking}
 */
void LPRowSet::remove(int i)
{
   SVSet::remove(i);
   left[i] = left[num()];
   right[i] = right[num()];
   left.reDim (num());
   right.reDim(num());
}

void LPRowSet::remove(int perm[])
{
   int i;
   int j = num();
   SVSet::remove(perm);
   for (i = 0; i < j; ++i)
   {
      if (perm[i] >= 0 && perm[i] != i)
      {
         left[perm[i]] = left[i];
         right[perm[i]] = right[i];
      }
   }
   left.reDim (num());
   right.reDim(num());
}

void LPRowSet::remove(const int nums[], int n, int* perm)
{
   SVSet::remove(nums, n, perm);
   int i;
   int j = num();
   for (i = 0; i < j; ++i)
   {
      if (perm[i] >= 0 && perm[i] != i)
      {
         left[perm[i]] = left[i];
         right[perm[i]] = right[i];
      }
   }
   left.reDim (num());
   right.reDim(num());
}

void LPRowSet::clear()
{
   SVSet::clear();
   left.reDim (num());
   right.reDim(num());
}

// TK13OCT1998
void LPRowSet::setType(
   int i,
   LPRow::Type p_type)
{
   switch (p_type)
   {
   case LPRow::LESS_EQUAL:
      lhs_w(i) = -infinity;
      break;
   case LPRow::EQUAL:
      if (lhs_w(i) > -infinity)
         rhs_w(i) = lhs(i);
      else
         lhs_w(i) = rhs(i);
      break;
   case LPRow::GREATER_EQUAL:
      rhs_w(i) = infinity;
      break;
   case LPRow::RANGE :
      MSG_ERROR( spxout << "EROWST01 RANGE not supported in LPRowSet::setType()" 
                        << std::endl; )
      throw SPxInternalCodeException("XROWST01 This should never happen.");
   default:
      throw SPxInternalCodeException("XROWST02 This should never happen.");
   }
}

bool LPRowSet::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   const int ldim = left.dim();
   if (ldim != right.dim())
      return MSGinconsistent("LPRowSet");
   if (ldim != num())
      return MSGinconsistent("LPRowSet");

   return left.isConsistent() && right.isConsistent() && SVSet::isConsistent();
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
