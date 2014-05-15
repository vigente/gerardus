/* $Id: setpseudo.c,v 1.14 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: setpseudo.c                                                    */
/*   Name....: Set Pseudo Functions                                           */
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

#define SET_PSEUDO_SID          0x53455452
#define SET_PSEUDO_ITER_SID     0x53455249

/* ------------------------------------------------------------------------- 
 * --- valid                 
 * -------------------------------------------------------------------------
 */
static Bool set_pseudo_is_valid(const Set* set)
{
   return set != NULL
      && SID_ok2(set->pseudo, SET_PSEUDO_SID)
      && set->head.refc > 0
      && set->head.dim     == 0
      && set->head.members == 1;
}

static Bool set_pseudo_iter_is_valid(const SetIter* iter)
{
   return iter != NULL && SID_ok2(iter->pseudo, SET_PSEUDO_ITER_SID);
}

/* ------------------------------------------------------------------------- 
 * --- set_new                 
 * -------------------------------------------------------------------------
 */
Set* set_pseudo_new()
{
   Set* set;

   set = calloc(1, sizeof(*set));

   assert(set != NULL);

   set->head.refc    = 1;
   set->head.dim     = 0;
   set->head.members = 1;
   set->head.type    = SET_PSEUDO;

   SID_set2(set->pseudo, SET_PSEUDO_SID);

   assert(set_pseudo_is_valid(set));
   
   return set;
}

/* ------------------------------------------------------------------------- 
 * --- copy
 * -------------------------------------------------------------------------
 */
static Set* set_pseudo_copy(const Set* source)
{
   Set* set = (Set*)source;
   
   set->head.refc++;

   return set;
}

/* ------------------------------------------------------------------------- 
 * --- set_free                 
 * -------------------------------------------------------------------------
 */
static void set_pseudo_free(Set* set)
{
   assert(set_pseudo_is_valid(set));

   set->head.refc--;

   if (set->head.refc == 0)
   {
      SID_del2(set->pseudo);

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
static int set_pseudo_lookup_idx(const Set* set, const Tuple* tuple, int offset)
{
   assert(set_pseudo_is_valid(set));
   assert(tuple_is_valid(tuple));
   assert(offset == 0);
   assert(tuple_get_dim(tuple) == 0);

   return 0;
}

/* ------------------------------------------------------------------------- 
 * --- get_tuple                 
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void set_pseudo_get_tuple(
   const Set* set,
   int        idx,
   Tuple*     tuple,
   int        offset)
{
   assert(set_pseudo_is_valid(set));
   assert(idx == 0);
   assert(tuple_is_valid(tuple));
   assert(offset == 0);
   assert(tuple_get_dim(tuple) == 0);
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
   
   assert(set_pseudo_is_valid(set));
   assert(pattern == NULL || tuple_is_valid(pattern));
   assert(pattern == NULL || tuple_get_dim(pattern) == 0);
   assert(offset                 == 0);

   iter = calloc(1, sizeof(*iter));

   assert(iter != NULL);

   iter->pseudo.first = TRUE;
   
   SID_set2(iter->pseudo, SET_PSEUDO_ITER_SID);

   assert(set_pseudo_iter_is_valid(iter));

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
   assert(set_pseudo_iter_is_valid(iter));
   assert(set_pseudo_is_valid(set));

   if (!iter->pseudo.first)
      return FALSE;

   iter->pseudo.first = FALSE;
   
   return TRUE;
}

/* ------------------------------------------------------------------------- 
 * --- iter_exit
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void iter_exit(SetIter* iter, const Set* set)
{
   assert(set_pseudo_iter_is_valid(iter));

   SID_del2(iter->pseudo);

   free(iter);
}

/* ------------------------------------------------------------------------- 
 * --- iter_reset
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void iter_reset(SetIter* iter, const Set* set)
{
   assert(set_pseudo_iter_is_valid(iter));

   iter->pseudo.first = TRUE;
}

/* ------------------------------------------------------------------------- 
 * --- vtab_init
 * -------------------------------------------------------------------------
 */
void set_pseudo_init(SetVTab* vtab)
{
   vtab[SET_PSEUDO].set_copy       = set_pseudo_copy;
   vtab[SET_PSEUDO].set_free       = set_pseudo_free;
   vtab[SET_PSEUDO].set_lookup_idx = set_pseudo_lookup_idx;
   vtab[SET_PSEUDO].set_get_tuple  = set_pseudo_get_tuple;
   vtab[SET_PSEUDO].iter_init      = iter_init;
   vtab[SET_PSEUDO].iter_next      = iter_next;
   vtab[SET_PSEUDO].iter_exit      = iter_exit;
   vtab[SET_PSEUDO].iter_reset     = iter_reset;
   vtab[SET_PSEUDO].set_is_valid   = set_pseudo_is_valid;
}







