/* $Id: setempty.c,v 1.14 2012/07/29 15:09:29 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: setempty.c                                                    */
/*   Name....: Set Empty Functions                                           */
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
#include "set.h"
#include "set4.h"

#ifdef _MSC_VER
#pragma warning (disable: 4100) /* unreferenced formal parameter */
#endif

#define SET_EMPTY_SID          0x5345544f
#define SET_EMPTY_ITER_SID     0x53454f49

/* ------------------------------------------------------------------------- 
 * --- valid                 
 * -------------------------------------------------------------------------
 */
static Bool set_empty_is_valid(const Set* set)
{
   return set != NULL
      && SID_ok2(set->empty, SET_EMPTY_SID)
      && set->head.refc > 0
      && set->head.members == 0;
}

static Bool set_empty_iter_is_valid(const SetIter* iter)
{
   return iter != NULL && SID_ok2(iter->empty, SET_EMPTY_ITER_SID);
}

/* ------------------------------------------------------------------------- 
 * --- set_new                 
 * -------------------------------------------------------------------------
 */
Set* set_empty_new(int dim)
{
   Set* set;

   set = calloc(1, sizeof(*set));

   assert(set != NULL);

   set->head.refc    = 1;
   set->head.dim     = dim;
   set->head.members = 0;
   set->head.type    = SET_EMPTY;

   SID_set2(set->empty, SET_EMPTY_SID);

   assert(set_empty_is_valid(set));
   
   return set;
}

/* ------------------------------------------------------------------------- 
 * --- copy
 * -------------------------------------------------------------------------
 */
static Set* set_empty_copy(const Set* source)
{
   Set* set = (Set*)source;
   
   set->head.refc++;

   return set;
}

/* ------------------------------------------------------------------------- 
 * --- set_free                 
 * -------------------------------------------------------------------------
 */
static void set_empty_free(Set* set)
{
   assert(set_empty_is_valid(set));

   set->head.refc--;

   if (set->head.refc == 0)
   {
      SID_del2(set->empty);

      free(set);
   }
}

/* ------------------------------------------------------------------------- 
 * --- lookup                 
 * -------------------------------------------------------------------------
 */
/* Return index number of element. -1 if not present
 */
/*ARGSUSED*/
static int set_empty_lookup_idx(const Set* set, const Tuple* tuple, int offset)
{
   assert(set_empty_is_valid(set));
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);

   return -1;
}

/* ------------------------------------------------------------------------- 
 * --- get_tuple                 
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void set_empty_get_tuple(
   const Set* set,
   int        idx,
   Tuple*     tuple,
   int        offset)
{
   assert(set_empty_is_valid(set));
   assert(idx >= 0);
   assert(idx <= set->head.members);
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset + set->head.dim <= tuple_get_dim(tuple));

   fprintf(stderr, "internal_error: set_empty_get_tuple called\n");
   abort();
}

/* ------------------------------------------------------------------------- 
 * --- iter_init                 
 * -------------------------------------------------------------------------
 */
/* Initialise Iterator. Write into iter
 */
/*ARGSUSED*/
static SetIter* iter_init(
   const Set*   set,
   const Tuple* pattern,
   int          offset)
{
   SetIter*        iter;
   
   assert(set_empty_is_valid(set));
   assert(pattern == NULL || tuple_is_valid(pattern));
   assert(offset      == 0);

   iter = calloc(1, sizeof(*iter));

   assert(iter != NULL);
   
   SID_set2(iter->empty, SET_EMPTY_ITER_SID);

   assert(set_empty_iter_is_valid(iter));

   return iter;
}

/* ------------------------------------------------------------------------- 
 * --- iter_next
 * -------------------------------------------------------------------------
 */
/* FALSE means, there is no further element
 */
/*ARGSUSED*/
static Bool iter_next(
   SetIter*   iter,
   const Set* set,
   Tuple*     tuple,
   int        offset)
{
   assert(set_empty_iter_is_valid(iter));
   assert(set_empty_is_valid(set));

   return FALSE;
}

/* ------------------------------------------------------------------------- 
 * --- iter_exit
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void iter_exit(SetIter* iter, const Set* set)
{
   assert(set_empty_iter_is_valid(iter));

   SID_del2(iter->empty);
   
   free(iter);
}

/* ------------------------------------------------------------------------- 
 * --- iter_reset
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void iter_reset(SetIter* iter, const Set* set)
{
   assert(set_empty_iter_is_valid(iter));
}

/* ------------------------------------------------------------------------- 
 * --- vtab_init
 * -------------------------------------------------------------------------
 */
void set_empty_init(SetVTab* vtab)
{
   vtab[SET_EMPTY].set_copy       = set_empty_copy;
   vtab[SET_EMPTY].set_free       = set_empty_free;
   vtab[SET_EMPTY].set_lookup_idx = set_empty_lookup_idx;
   vtab[SET_EMPTY].set_get_tuple  = set_empty_get_tuple;
   vtab[SET_EMPTY].iter_init      = iter_init;
   vtab[SET_EMPTY].iter_next      = iter_next;
   vtab[SET_EMPTY].iter_exit      = iter_exit;
   vtab[SET_EMPTY].iter_reset     = iter_reset;
   vtab[SET_EMPTY].set_is_valid   = set_empty_is_valid;
}




