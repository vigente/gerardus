/* $Id: setrange.c,v 1.20 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: setrange.c                                                    */
/*   Name....: Set Range Functions                                           */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Copyright (C) 2001-2012 by Thorsten Koch <koch@zib.de>
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "bool.h"
#include "mshell.h"
#include "numb.h"
#include "elem.h"
#include "tuple.h"
#include "mme.h"
#include "hash.h"
#include "stmt.h"
#include "set.h"
#include "set4.h"

#ifdef _MSC_VER
#pragma warning (disable: 4100) /* unreferenced formal parameter */
#endif

#define SET_RANGE_SID          0x53455452
#define SET_RANGE_ITER_SID     0x53455249

/* ------------------------------------------------------------------------- 
 * --- valid                 
 * -------------------------------------------------------------------------
 */
static Bool set_range_is_valid(const Set* set)
{
   return set != NULL
      && SID_ok2(set->range, SET_RANGE_SID)
      && set->head.refc > 0
      && set->head.dim == 1;
}

static Bool set_range_iter_is_valid(const SetIter* iter)
{
   return iter != NULL && SID_ok2(iter->range, SET_RANGE_ITER_SID)
      && iter->range.first >= 0
      && iter->range.last  >= 0
      && iter->range.now   >= iter->range.first;
}

/* ------------------------------------------------------------------------- 
 * --- set_new                 
 * -------------------------------------------------------------------------
 */
Set* set_range_new(int begin, int end, int step)
{
   Set* set;

   set = calloc(1, sizeof(*set));

   assert(set != NULL);

   set->head.refc    = 1;
   set->head.dim     = 1;
   set->head.members = 1 + (end - begin) / step;
   set->head.type    = SET_RANGE;

   set->range.begin  = begin;
   set->range.end    = end;
   set->range.step   = step;

   SID_set2(set->range, SET_RANGE_SID);

   assert(set_range_is_valid(set));
   
   return set;
}

/* ------------------------------------------------------------------------- 
 * --- copy
 * -------------------------------------------------------------------------
 */
static Set* set_range_copy(const Set* source)
{
   Set* set = (Set*)source;
   
   set->head.refc++;

   return set;
}

/* ------------------------------------------------------------------------- 
 * --- set_free                 
 * -------------------------------------------------------------------------
 */
static void set_range_free(Set* set)
{
   assert(set_range_is_valid(set));

   set->head.refc--;

   if (set->head.refc == 0)
   {
      SID_del2(set->range);

      free(set);
   }
}

/* ------------------------------------------------------------------------- 
 * --- lookup                 
 * -------------------------------------------------------------------------
 */
/* Return index number of element. -1 if not present
 */
static int idx_to_val(int begin, int step, int idx)
{
#if 0
   fprintf(stderr, "idx_to_val: %d %d %d = %d\n",
      begin, step, idx, begin + idx * step);
#endif
   return begin + idx * step;
}

static int val_to_idx(int begin, int step, int val)
{
#if 0
   fprintf(stderr, "val_to_idx: %d %d %d = %d\n",
      begin, step, val, (val - begin) / step);
#endif
   return (val - begin) / step;
}

static int set_range_lookup_idx(const Set* set, const Tuple* tuple, int offset)
{
   const Elem* elem;
   const Numb* numb;
   int         val;
   
   assert(set_range_is_valid(set));
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset <  tuple_get_dim(tuple));
   
   elem = tuple_get_elem(tuple, offset);

   /* If this fails, this means we have a set with mixed number
    * and integer types.
    */
   assert(elem_get_type(elem) == ELEM_NUMB);

   numb = elem_get_numb(elem);

   assert(numb_is_int(numb));

   val = numb_toint(numb);

   if (set->range.step > 0)
   {
      if (  val < set->range.begin 
         || val > set->range.end
         || ((val - set->range.begin) % set->range.step) != 0)
         return -1;
   }
   else
   {
      assert(set->range.step < 0);
      
      if (  val > set->range.begin 
         || val < set->range.end
         || ((set->range.begin - val) % set->range.step) != 0)
         return -1;
   }
   return val_to_idx(set->range.begin, set->range.step, val);
}

/* ------------------------------------------------------------------------- 
 * --- get_tuple                 
 * -------------------------------------------------------------------------
 */
static void set_range_get_tuple(
   const Set* set,
   int        idx,
   Tuple*     tuple,
   int        offset)
{
   int   val;
   Numb* numb;
      
   assert(set_range_is_valid(set));
   assert(idx >= 0);
   assert(idx <= set->head.members);
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset <  tuple_get_dim(tuple));

   val  = idx_to_val(set->range.begin, set->range.step, idx);
   numb = numb_new_integer(val);

   tuple_set_elem(tuple, offset, elem_new_numb(numb));

   numb_free(numb);
}

/* ------------------------------------------------------------------------- 
 * --- iter_init                 
 * -------------------------------------------------------------------------
 */
/* Initialise Iterator. Write into iter
 */
static SetIter* set_range_iter_init(
   const Set*   set,
   const Tuple* pattern,
   int          offset)
{
   const Elem*  elem;
   SetIter*     iter;
   
   assert(set_range_is_valid(set));
   assert(pattern == NULL || tuple_is_valid(pattern));
   assert(offset      >= 0);
   assert(pattern == NULL || offset <  tuple_get_dim(pattern));

   iter = calloc(1, sizeof(*iter));

   assert(iter != NULL);

   if (pattern == NULL)
   {
      iter->range.first = 0;
      iter->range.last  = val_to_idx(set->range.begin, set->range.step, set->range.end);
   }
   else
   {
      elem = tuple_get_elem(pattern, offset);

      switch(elem_get_type(elem))
      {
      case ELEM_NAME :
         iter->range.first = 0;
         iter->range.last  = val_to_idx(set->range.begin, set->range.step, set->range.end);
         break;
      case ELEM_NUMB :
         iter->range.first = set_range_lookup_idx(set, pattern, offset);

         if (iter->range.first >= 0)
            iter->range.last = iter->range.first;
         else
         {
            iter->range.first = 1;
            iter->range.last  = 0;
         }
         break;
      case ELEM_STRG :
         /* This should not happen. Probably a set with mixed
          * numbers and string was generated.
          */
      default :
         abort();
      }
   }
   iter->range.now = iter->range.first;

   SID_set2(iter->range, SET_RANGE_ITER_SID);

   assert(set_range_iter_is_valid(iter));

   return iter;
}

/* ------------------------------------------------------------------------- 
 * --- iter_next
 * -------------------------------------------------------------------------
 */
/* FALSE means, there is no further element
 */
static Bool set_range_iter_next(
   SetIter*   iter,
   const Set* set,
   Tuple*     tuple,
   int        offset)
{
   int   val;
   Numb* numb;

   assert(set_range_iter_is_valid(iter));
   assert(set_range_is_valid(set));
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset <  tuple_get_dim(tuple));

   if (iter->range.now > iter->range.last)
      return FALSE;

   val  = idx_to_val(set->range.begin, set->range.step, iter->range.now);
   numb = numb_new_integer(val);

   tuple_set_elem(tuple, offset, elem_new_numb(numb));

   numb_free(numb);

   iter->range.now++;

   return TRUE;
}

/* ------------------------------------------------------------------------- 
 * --- iter_exit
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void set_range_iter_exit(SetIter* iter, const Set* set)
{
   assert(set_range_iter_is_valid(iter));

   SID_del2(iter->range);
   
   free(iter);
}

/* ------------------------------------------------------------------------- 
 * --- iter_reset
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void set_range_iter_reset(SetIter* iter, const Set* set)
{
   assert(set_range_iter_is_valid(iter));
   
   iter->range.now = iter->range.first;
}

/* ------------------------------------------------------------------------- 
 * --- vtab_init
 * -------------------------------------------------------------------------
 */
void set_range_init(SetVTab* vtab)
{
   vtab[SET_RANGE].set_copy       = set_range_copy;
   vtab[SET_RANGE].set_free       = set_range_free;
   vtab[SET_RANGE].set_lookup_idx = set_range_lookup_idx;
   vtab[SET_RANGE].set_get_tuple  = set_range_get_tuple;
   vtab[SET_RANGE].iter_init      = set_range_iter_init;
   vtab[SET_RANGE].iter_next      = set_range_iter_next;
   vtab[SET_RANGE].iter_exit      = set_range_iter_exit;
   vtab[SET_RANGE].iter_reset     = set_range_iter_reset;
   vtab[SET_RANGE].set_is_valid   = set_range_is_valid;
}


