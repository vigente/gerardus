/* $Id: setlist.c,v 1.26 2012/07/29 15:09:29 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: setlist.c                                                     */
/*   Name....: Set List Functions                                            */
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
#include "ratlptypes.h"
#include "numb.h"
#include "elem.h"
#include "tuple.h"
#include "mme.h"
#include "entry.h"
#include "hash.h"
#include "list.h"
#include "stmt.h"
#include "set.h"
#include "set4.h"

#ifdef _MSC_VER
#pragma warning (disable: 4100) /* unreferenced formal parameter */
#endif

#define SET_LIST_SID          0x5345544c
#define SET_LIST_ITER_SID     0x53454c49

#define HASH_THRESHOLD        12

/* ------------------------------------------------------------------------- 
 * --- valid                 
 * -------------------------------------------------------------------------
 */
static Bool set_list_is_valid(const Set* set)
{
   int i;

   if (set == NULL
      || !SID_ok2(set->list, SET_LIST_SID)
      || set->head.refc    <= 0
      || set->head.dim     != 1
      || set->head.members <  0
      || set->list.size    <  set->head.members
      || set->list.member  == NULL)
      return FALSE;
   
   for(i = 0; i < set->head.members; i++)
      if (!elem_is_valid(set->list.member[i]))
         return FALSE;

   return TRUE;
}

static Bool set_list_iter_is_valid(const SetIter*iter)
{
   return iter != NULL
      && SID_ok2(iter->list, SET_LIST_ITER_SID)
      && iter->list.first >=  0
      && iter->list.last  >= -1
      && iter->list.now   >=  0
      && iter->list.now   >=  iter->list.first;
}

/* ------------------------------------------------------------------------- 
 * --- internal                 
 * -------------------------------------------------------------------------
 */
/* Return index number of element. -1 if not present
 */
static int lookup_elem_idx(const Set* set, const Elem* elem)
{
   int i;
   
   assert(set_list_is_valid(set));
   assert(elem_is_valid(elem));

   if (set->list.hash != NULL)
      return hash_lookup_elem_idx(set->list.hash, elem);

   for(i = 0; i < set->list.head.members; i++)
      if (!elem_cmp(elem, set->list.member[i]))
         return i;

   return -1;
}

/* ------------------------------------------------------------------------- 
 * --- new                 
 * -------------------------------------------------------------------------
 */
Set* set_list_new(int size, int flags)
{
   Set* set;

   assert(size > 0);
   
   set = calloc(1, sizeof(*set));

   assert(set != NULL);

   set->head.refc    = 1;
   set->head.dim     = 1;
   set->head.members = 0;
   set->head.type    = SET_LIST;
   set->list.size    = size;
   set->list.member  = calloc((size_t)size, sizeof(*set->list.member));

   assert(set->list.member != NULL);

   if ((flags & SET_NO_HASH) == 0 && size > HASH_THRESHOLD)
      set->list.hash = hash_new(HASH_ELEM_IDX, size);

   SID_set2(set->list, SET_LIST_SID);

   assert(set_list_is_valid(set));
   
   return set;
}

int set_list_add_elem(Set* set, const Elem* elem, SetCheckType check)
{
   int idx = -1;
   
   assert(set_list_is_valid(set));
   assert(elem_is_valid(elem));
   
   if (check != SET_CHECK_NONE && (idx = lookup_elem_idx(set, elem)) >= 0)
   {
      if (check != SET_CHECK_QUIET)
      {
         assert(check == SET_CHECK_WARN);

         if (stmt_trigger_warning(164))
         {
            fprintf(stderr, "--- Warning 164: Duplicate element ");
            elem_print(stderr, elem, TRUE);
            fprintf(stderr, " for set rejected\n");
         }
      }
   }
   else
   {
      idx = set->head.members;

      set->list.member[idx] = elem_copy(elem); 

      if (set->list.hash != NULL)
         hash_add_elem_idx(set->list.hash, set->list.member[idx], idx);
      
      set->head.members++;

      assert(set->head.members <= set->list.size);
   }
   assert(idx >= 0);
   
   return idx;
}

Set* set_list_new_from_elems(const List* list, SetCheckType check)
{
   ListElem*   le = NULL;
   Set*        set;
   int         n;

   assert(list_is_valid(list));

   n   = list_get_elems(list);

   assert(n > 0);

   set = set_list_new(n, SET_DEFAULT);

   while(n-- > 0)
      (void)set_list_add_elem(set, list_get_elem(list, &le), check);

   assert(set_list_is_valid(set));

   return set;
}

Set* set_list_new_from_tuples(const List* list, SetCheckType check)
{
   ListElem*    le = NULL;
   const Tuple* tuple;
   Set*         set;
   int          n;

   assert(list_is_valid(list));

   n   = list_get_elems(list);

   assert(n > 0);

   set = set_list_new(n, SET_DEFAULT);

   while(n-- > 0)
   {
      tuple = list_get_tuple(list, &le);

      assert(tuple_get_dim(tuple) == 1);
      
      (void)set_list_add_elem(set, tuple_get_elem(tuple, 0), check);
   }
   assert(set_list_is_valid(set));

   return set;
}

Set* set_list_new_from_entries(const List* list, SetCheckType check)
{
   ListElem*    le = NULL;
   const Tuple* tuple;
   Set*         set;
   int          n;

   assert(list_is_valid(list));

   n   = list_get_elems(list);

   assert(n > 0);

   set = set_list_new(n, SET_DEFAULT);

   while(n-- > 0)
   {
      tuple = entry_get_tuple(list_get_entry(list, &le));

      assert(tuple_get_dim(tuple) == 1);
      
      (void)set_list_add_elem(set, tuple_get_elem(tuple, 0), check);
   }
   assert(set_list_is_valid(set));

   return set;
}

/* ------------------------------------------------------------------------- 
 * --- copy
 * -------------------------------------------------------------------------
 */
static Set* set_list_copy(const Set* source)
{
   Set* set = (Set*)source;
   
   set->head.refc++;

   return set;
}

/* ------------------------------------------------------------------------- 
 * --- free                 
 * -------------------------------------------------------------------------
 */
static void set_list_free(Set* set)
{
   int i;
   
   assert(set_list_is_valid(set));

   set->head.refc--;

   if (set->head.refc == 0)
   {
      SID_del2(set->list);

      for(i = 0; i < set->head.members; i++)
         elem_free(set->list.member[i]);

      if (set->list.hash != NULL)
         hash_free(set->list.hash);
      
      free(set->list.member);
      free(set);
   }
}

/* ------------------------------------------------------------------------- 
 * --- lookup                 
 * -------------------------------------------------------------------------
 */
/* Return index number of element. -1 if not present
 */
static int set_list_lookup_idx(const Set* set, const Tuple* tuple, int offset)
{
   assert(set_list_is_valid(set));
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset <  tuple_get_dim(tuple));
   
   return lookup_elem_idx(set, tuple_get_elem(tuple, offset));
}

/* ------------------------------------------------------------------------- 
 * --- get_tuple                 
 * -------------------------------------------------------------------------
 */
static void set_list_get_tuple(
   const Set* set,
   int        idx,
   Tuple*     tuple,
   int        offset)
{
   assert(set_list_is_valid(set));
   assert(idx >= 0);
   assert(idx <= set->head.members);
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset <  tuple_get_dim(tuple));

   tuple_set_elem(tuple, offset, elem_copy(set->list.member[idx]));
}

/* ------------------------------------------------------------------------- 
 * --- iter_init                 
 * -------------------------------------------------------------------------
 */
/* Initialise Iterator. Write into iter
 */
static SetIter* set_list_iter_init(
   const Set*   set,
   const Tuple* pattern,
   int          offset)
{
   const Elem*  elem;
   SetIter*     iter;

   assert(set_list_is_valid(set));
   assert(pattern == NULL || tuple_is_valid(pattern));
   assert(offset        >= 0);
   assert(pattern == NULL || offset < tuple_get_dim(pattern));

   iter = calloc(1, sizeof(*iter));
   
   assert(iter != NULL);

   if (pattern == NULL)
   {
      iter->list.first = 0;
      iter->list.last  = set->head.members - 1;
   }
   else
   {
      elem = tuple_get_elem(pattern, offset);
   
      if (elem_get_type(elem) == ELEM_NAME)
      {
         iter->list.first = 0;
         iter->list.last  = set->head.members - 1;
      }
      else
      {
         iter->list.first = lookup_elem_idx(set, elem);

         if (iter->list.first >= 0)
            iter->list.last  = iter->list.first;
         else
         {
            iter->list.first = 1;
            iter->list.last  = 0;
         }
      }
   }
   iter->list.now = iter->list.first;

   SID_set2(iter->list, SET_LIST_ITER_SID);

   assert(set_list_iter_is_valid(iter));

   return iter;
}

/* ------------------------------------------------------------------------- 
 * --- iter_next
 * -------------------------------------------------------------------------
 */
/* FALSE means, there is no further element
 */
static Bool set_list_iter_next(
   SetIter*   iter,
   const Set* set,
   Tuple*     tuple,
   int        offset)
{
   assert(set_list_iter_is_valid(iter));
   assert(set_list_is_valid(set));
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset <  tuple_get_dim(tuple));
   
   if (iter->list.now > iter->list.last)
      return FALSE;

   tuple_set_elem(tuple, offset, elem_copy(set->list.member[iter->list.now]));

   iter->list.now++;

   return TRUE;
}

/* ------------------------------------------------------------------------- 
 * --- iter_exit
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void set_list_iter_exit(SetIter* iter, const Set* set)
{
   assert(set_list_iter_is_valid(iter));

   SID_del2(iter->list);
   
   free(iter);
}

/* ------------------------------------------------------------------------- 
 * --- iter_reset
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void set_list_iter_reset(SetIter* iter, const Set* set)
{
   assert(set_list_iter_is_valid(iter));
   
   iter->list.now = iter->list.first;
}

/* ------------------------------------------------------------------------- 
 * --- vtab_init
 * -------------------------------------------------------------------------
 */
void set_list_init(SetVTab* vtab)
{
   vtab[SET_LIST].set_copy       = set_list_copy;
   vtab[SET_LIST].set_free       = set_list_free;
   vtab[SET_LIST].set_lookup_idx = set_list_lookup_idx;
   vtab[SET_LIST].set_get_tuple  = set_list_get_tuple;
   vtab[SET_LIST].iter_init      = set_list_iter_init;
   vtab[SET_LIST].iter_next      = set_list_iter_next;
   vtab[SET_LIST].iter_exit      = set_list_iter_exit;
   vtab[SET_LIST].iter_reset     = set_list_iter_reset;
   vtab[SET_LIST].set_is_valid   = set_list_is_valid;
}

/* ------------------------------------------------------------------------- 
 * --- extras
 * -------------------------------------------------------------------------
 */
const Elem* set_list_get_elem(const Set* set, int idx)
{
   assert(set_list_is_valid(set));
   assert(idx >= 0);
   assert(idx <  set->head.members);

   return set->list.member[idx];
}





