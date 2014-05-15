/* $Id: setprod.c,v 1.20 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: setprod.c                                                     */
/*   Name....: Set Product Functions                                         */
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

#define SET_PROD_SID          0x53455450
#define SET_PROD_ITER_SID     0x53455049

/* ------------------------------------------------------------------------- 
 * --- valid                 
 * -------------------------------------------------------------------------
 */
static Bool set_prod_is_valid(const Set* set)
{
   return set != NULL
      && SID_ok2(set->prod, SET_PROD_SID)
      && set->head.refc > 0
      && set->head.dim > 1
      && set_is_valid(set->prod.set_a)
      && set_is_valid(set->prod.set_b);
}

static Bool set_prod_iter_is_valid(const SetIter* iter)
{
   return iter != NULL
      && SID_ok2(iter->prod, SET_PROD_ITER_SID)
      && iter->prod.iter_a != NULL
      && iter->prod.iter_b != NULL;
}

/* ------------------------------------------------------------------------- 
 * --- set_new                 
 * -------------------------------------------------------------------------
 */
Set* set_prod_new(const Set* a, const Set* b)
{
   Set* set;

   assert(set_is_valid(a));
   assert(set_is_valid(b));
   assert(a->head.type != SET_PSEUDO);
   assert(b->head.type != SET_PSEUDO);
   
   if (a->head.type == SET_EMPTY || b->head.type == SET_EMPTY)
      return set_empty_new(a->head.dim + b->head.dim);
   
   set = calloc(1, sizeof(*set));

   assert(set != NULL);

   set->head.refc    = 1;
   set->head.dim     = a->head.dim + b->head.dim;
   set->head.members = a->head.members * b->head.members;
   set->head.type    = SET_PROD;

   set->prod.set_a   = set_copy(a);
   set->prod.set_b   = set_copy(b);

   SID_set2(set->prod, SET_PROD_SID);

   assert(set_prod_is_valid(set));

   return set;
}

/* ------------------------------------------------------------------------- 
 * --- copy
 * -------------------------------------------------------------------------
 */
static Set* set_prod_copy(const Set* source)
{
   Set* set = (Set*)source;
   
   set->head.refc++;

   (void)set_copy(set->prod.set_a);
   (void)set_copy(set->prod.set_b);
   
   return set;
}

/* ------------------------------------------------------------------------- 
 * --- set_free                 
 * -------------------------------------------------------------------------
 */
static void set_prod_free(Set* set)
{
   assert(set_prod_is_valid(set));

   set_free(set->prod.set_a);
   set_free(set->prod.set_b);

   set->head.refc--;

   if (set->head.refc == 0)
   {
      SID_del2(set->prod);
      
      free(set);
   }
}

/* ------------------------------------------------------------------------- 
 * --- lookup                 
 * -------------------------------------------------------------------------
 */
/* Return index number of element. -1 if not present
 */
static int set_prod_lookup_idx(const Set* set, const Tuple* tuple, int offset)
{
   int            idx_a;
   int            idx_b;
      
   assert(set_prod_is_valid(set));
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset <  tuple_get_dim(tuple));

   idx_a = set_lookup_idx(set->prod.set_a, tuple, offset);

   if (idx_a < 0)
      return -1;

   offset += set->prod.set_a->head.dim;

   idx_b = set_lookup_idx(set->prod.set_b, tuple, offset);

   if (idx_b < 0)
      return -1;

   return idx_a * set->prod.set_b->head.members + idx_b;
}

/* ------------------------------------------------------------------------- 
 * --- get_tuple                 
 * -------------------------------------------------------------------------
 */
static void set_prod_get_tuple(
   const Set* set,
   int        idx,
   Tuple*     tuple,
   int        offset)
{
   const Set* a;
   const Set* b;
   int        offset2;
   
   assert(set_prod_is_valid(set));
   assert(idx >= 0);
   assert(idx <= set->head.members);
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset + set->head.dim <= tuple_get_dim(tuple));

   a       = set->prod.set_a;
   b       = set->prod.set_b;
   offset2 = offset + a->head.dim;

   set_get_tuple_intern(a, idx / b->head.members, tuple, offset);
   set_get_tuple_intern(b, idx % b->head.members, tuple, offset2);
}

/* ------------------------------------------------------------------------- 
 * --- iter_init                 
 * -------------------------------------------------------------------------
 */
/* Initialise Iterator. Write into iter
 */
static SetIter* set_prod_iter_init(
   const Set*   set,
   const Tuple* pattern,
   int          offset)
{
   SetIter*     iter;

   assert(set_prod_is_valid(set));
   assert(pattern == NULL || tuple_is_valid(pattern));
   assert(offset     >= 0);
   assert(pattern == NULL || offset <  tuple_get_dim(pattern));

   iter = calloc(1, sizeof(*iter));

   assert(iter != NULL);

   iter->prod.elem = calloc(set->head.dim, sizeof(*iter->prod.elem));

   assert(iter->prod.elem != NULL);

   iter->prod.first = TRUE;
   
   iter->prod.iter_a = set_iter_init_intern(set->prod.set_a, pattern, offset);
   iter->prod.iter_b = set_iter_init_intern(set->prod.set_b, pattern,
      offset + set->prod.set_a->head.dim);

   SID_set2(iter->prod, SET_PROD_ITER_SID);

   assert(set_prod_iter_is_valid(iter));

   return iter;
}

/* ------------------------------------------------------------------------- 
 * --- iter_next
 * -------------------------------------------------------------------------
 */
/* This gets the fore part of the product and saves it and
 * also gets the back part
 */
static Bool get_both_parts(
   const Set* a,
   const Set* b,
   SetIter*   iter,
   SetIter*   iter_a,
   SetIter*   iter_b,
   Tuple*     tuple,
   int        offset,
   int        offset2)
{
   int i;
   
   if (!set_iter_next_intern(iter_a, a, tuple, offset))
      return FALSE;

   for(i = 0; i < a->head.dim; i++)
   {
      assert(iter->prod.elem[i] == NULL);
      
      iter->prod.elem[i] = elem_copy(tuple_get_elem(tuple, i + offset));

      assert(elem_is_valid(iter->prod.elem[i]));
   }
   if (!set_iter_next_intern(iter_b, b, tuple, offset2))
      return FALSE;

   return TRUE;
}


/* FALSE means, there is no further element
 */
/*ARGSUSED*/
static Bool set_prod_iter_next(
   SetIter*   iter,
   const Set* set,
   Tuple*     tuple,
   int        offset)
{
   Set*       a;
   Set*       b;
   SetIter*   iter_a;
   SetIter*   iter_b;
   int        offset2;
   int        i;
   
   assert(set_prod_iter_is_valid(iter));
   assert(set_prod_is_valid(set));
   assert(tuple_is_valid(tuple));
   assert(offset >= 0);
   assert(offset + set->head.dim <= tuple_get_dim(tuple));

   a       = set->prod.set_a;
   b       = set->prod.set_b;
   iter_a  = iter->prod.iter_a;
   iter_b  = iter->prod.iter_b;
   offset2 = offset + a->head.dim;

   if (iter->prod.first)
   {
      iter->prod.first = FALSE;

      return get_both_parts(a, b, iter, iter_a, iter_b, tuple, offset, offset2);
   }
   assert(!iter->prod.first);
   
   /* Get back part
    */
   if (set_iter_next_intern(iter_b, b, tuple, offset2))
   {
      /* copy fore part
       */
      for(i = 0; i < a->head.dim; i++)
         tuple_set_elem(tuple, i + offset, elem_copy(iter->prod.elem[i]));

      return TRUE;
   }

   /* No back part, so reset it
    */
   set_iter_reset_intern(iter_b, b);

   /* Clear elem cache
    */
   for(i = 0; i < set->head.dim; i++)
   {
      if (iter->prod.elem[i] != NULL)
      {
         elem_free(iter->prod.elem[i]);

         iter->prod.elem[i] = NULL;
      }
   }
   return get_both_parts(a, b, iter, iter_a, iter_b, tuple, offset, offset2);
}

/* ------------------------------------------------------------------------- 
 * --- iter_exit
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void set_prod_iter_exit(SetIter* iter, const Set* set)
{
   int i;
   
   assert(set_prod_iter_is_valid(iter));
   assert(set_prod_is_valid(set));

   SID_del2(iter->prod);

   set_iter_exit_intern(iter->prod.iter_a, set->prod.set_a);
   set_iter_exit_intern(iter->prod.iter_b, set->prod.set_b);

   for(i = 0; i < set->head.dim; i++)
      if (iter->prod.elem[i] != NULL)
         elem_free(iter->prod.elem[i]);

   free(iter->prod.elem);
   free(iter);
}

/* ------------------------------------------------------------------------- 
 * --- iter_reset
 * -------------------------------------------------------------------------
 */
/*ARGSUSED*/
static void set_prod_iter_reset(SetIter* iter, const Set* set)
{
   assert(set_prod_iter_is_valid(iter));
   assert(set_prod_is_valid(set));

   iter->prod.first = TRUE;
   
   set_iter_reset_intern(iter->prod.iter_a, set->prod.set_a);
   set_iter_reset_intern(iter->prod.iter_b, set->prod.set_b);
}

/* ------------------------------------------------------------------------- 
 * --- vtab_init
 * -------------------------------------------------------------------------
 */
void set_prod_init(SetVTab* vtab)
{
   vtab[SET_PROD].set_copy       = set_prod_copy;
   vtab[SET_PROD].set_free       = set_prod_free;
   vtab[SET_PROD].set_lookup_idx = set_prod_lookup_idx;
   vtab[SET_PROD].set_get_tuple  = set_prod_get_tuple;
   vtab[SET_PROD].iter_init      = set_prod_iter_init;
   vtab[SET_PROD].iter_next      = set_prod_iter_next;
   vtab[SET_PROD].iter_exit      = set_prod_iter_exit;
   vtab[SET_PROD].iter_reset     = set_prod_iter_reset;
   vtab[SET_PROD].set_is_valid   = set_prod_is_valid;
}

