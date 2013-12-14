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


/**@file  sorter.h
 * @brief Generic QuickSort implementation.
 */
#ifndef _SORTER_H_
#define _SORTER_H_

#include <assert.h>

namespace soplex
{
/// Generic QuickSort implementation.
/** This template function sorts an array \p t holding \p n elements of
    type T using \p compare for comparisons. Class COMPARATOR must provide
    an overloaded operator()(const T& t1,const T& t2) which returns
    - < 0, if \p t1 is to appear before \p t2,
    - = 0, if \p t1 and \p t2 can appear in any order, or
    - > 0, if \p t1 is to appear after \p t2.
*/

template < class T, class COMPARATOR >
static void sorter_qsort(T* t, int end, COMPARATOR& compare, int start = 0)
{
   assert(start >= 0);

   /* nothing to sort for less than two elements */
   if( end <= start + 1 )
      return;

   int i0, i1, j;
   Real c;

   T work, mid, tmp;

   work = t[start];
   t[start] = t[(start + end) / 2];
   t[(start + end) / 2] = work;
   
   mid = t[start];
   work = t[end - 1];
   
   for (i0 = i1 = start, j = end - 1; i1 < j;)
   {
      c = compare(mid, work);
      if (c > 0)
      {
         tmp = t[i0];
         t[i0] = work;
         i0++;
         i1++;
         work = t[i1];
         t[i1] = tmp;
      }
      else if (c < 0)
      {
         t[j] = work;
         --j;
         work = t[j];
      }
      else
      {
         i1++;
         tmp = t[i1];
         t[i1] = work;
         work = tmp;
      }
   }
   
   if (start < i0 - 1)
      sorter_qsort(t, i0, compare, start);
   if (i1 + 1 < end)
      sorter_qsort(t, end, compare, i1 + 1);
}


/**@brief  Generic implementation of Partial QuickSort.
 *
 * This template function sorts an array \p t holding \p n elements of
 * type T partially using \p compare for comparisons, i.e. ensures that
 * the \p size smallest elements are sorted to the front.
 *
 * Class COMPARATOR must provide an overloaded
 * operator()(const T& t1,const T& t2) which returns
 * - < 0, if \p t1 is to appear before \p t2,
 * - = 0, if \p t1 and \p t2 can appear in any order, or
 * - > 0, if \p t1 is to appear after \p t2.
 */
template < class T, class COMPARATOR >
static int sorter_qsortPart(
   T*                    t,                  /**< array of elements to be sorted between index start and end */
   COMPARATOR&           compare,            /**< comparator */
   int                   start,              /**< index of first element in range to be sorted */
   int                   end,                /**< index of last element in range to be sorted, plus 1 */
   int                   size,               /**< guaranteed number of sorted elements */
   int                   start2 = 0,         /**< auxiliary start index of sub range used for recursive call */
   int                   end2 = 0            /**< auxiliary end index of sub range used for recursive call */
   )
{
   assert(start >= 0);

   /* nothing to sort for less than two elements */
   if( end < start + 1 )
      return 0;
   else if( end <= start + 1 )
      return 1;

   T work;
   T pivotval;                               /* value of pivot element */
   T tmp;

   Real c;

   int pivotind;                             /* index of pivot element */
   int i0;                                   /* elements start, ..., i0-1 are strictly smaller than the pivot value */
   int i1;                                   /* elements i0, ..., i1-1 are equal to the pivot value */
   int j;                                    /* elements i1, ..., j-1 have arbitrary value; elements j, ..., end-1 are strictly greater than the pivot value */

   /* we assume that range {start, ..., start2-1} already contains the start2-start smallest elements in sorted order;
    * start2 has to lie in {start, ..., end-1} */
   if( start2 < start )
      start2 = start;

   /* if all remaining elements should be sorted, we simply call standard quicksort */
   if( start2+size >= end-1 )
   {
      sorter_qsort(t, end, compare, start2);
      return end > end2 ? end-1 : end2-1;
   }

   /* select middle element as pivot element */
   pivotind = (start2+end)/2;
   pivotval = t[pivotind];

   /* exchange first element with pivot element */
   t[pivotind] = t[start2];
   t[start2] = pivotval;

   /* the unsorted part covers {start2, ..., end-1} */
   i0 = start2;
   i1 = start2;
   j = end-1;

   /* we start at the end of the unsorted range */
   work = t[j];

   /* keep sorting as long as the part with arbitrary values is nonempty */
   while( i1 < j )
   {
      /* compare work to pivot element */
      c = compare(pivotval, work);

      /* work < mid */
      if( c > 0 )
      {
         tmp = t[i0];
         t[i0] = work;
         i0++;
         i1++;
         work = t[i1];
         t[i1] = tmp;
      }
      /* work > mid */
      else if (c < 0)
      {
         t[j] = work;
         j--;
         work = t[j];
      }
      /* work == mid */
      else
      {
         i1++;
         tmp = t[i1];
         t[i1] = work;
         work = tmp;
      }
   }

   /* if we only need to sort less than half of the "<" part, use partial sort again */
   if( 2*size <= i0 - start2 )
   {
      return sorter_qsortPart(t, compare, start, i0, size, start2, i1);
   }
   /* otherwise, and if we do not need to sort the ">" part, use standard quicksort on the "<" part */
   else if( size <= i1 - start2 )
   {
      sorter_qsort(t, i0, compare, start2);
      return i1-1;
   }
   /* otherwise we have to sort the "<" part fully (use standard quicksort) and the ">" part partially */
   else
   {
      sorter_qsort(t, i0, compare, start);
      return sorter_qsortPart(t, compare, start, end, size+start2-i1, i1+1, i1);
   }
}

} // namespace soplex
#endif // _SORTER_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
