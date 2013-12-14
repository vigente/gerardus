/* $Id: local.c,v 1.23 2012/07/29 15:09:27 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: local.c                                                       */
/*   Name....: Local Parameter Functions                                     */
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
#include "local.h"

#define LOCAL_STR_SIZE 100

typedef struct local Local;

struct local
{
   const char*  name;
   Elem*        element;
   Local*       next;
};

static Local anchor  = { "", NULL, NULL };

void local_init()
{
   assert(strlen(anchor.name) == 0);
   assert(anchor.element      == NULL);
   assert(anchor.next         == NULL);
}

void local_exit()
{
   assert(strlen(anchor.name) == 0);
   assert(anchor.element == NULL);
   
   anchor.next = NULL;
}

static void local_new(const char* name, const Elem* elem)
{
   Local* local;

   assert(name    != NULL);
   
   local = calloc(1, sizeof(*local));

   assert(local != NULL);
   
   local->name    = name;
   local->element = (elem == ELEM_NULL) ? ELEM_NULL : elem_copy(elem);
   local->next    = anchor.next;
   anchor.next    = local;
}

static void local_new_frame(void)
{
   local_new("", ELEM_NULL);
}

void local_drop_frame(void)
{
   Bool   frame = FALSE;
   Local* q     = NULL;
   Local* p;
   
   for(p = anchor.next; (p != NULL) && !frame; p = q)
   {
      q = p->next;

      if (p->element == ELEM_NULL)
         frame = TRUE;
      else
         elem_free(p->element);
      
      free(p);
   }
   anchor.next = q;
}

const Elem* local_lookup(const char* name)
{
   Local* local;

   assert(name != NULL);

   for(local = anchor.next; local != NULL; local = local->next)
      if (!strcmp(local->name, name))
         break;

   return local == NULL ? ELEM_NULL : local->element;
}

void local_install_tuple(const Tuple* pattern, const Tuple* values)
{
   const Elem* elem;
   const char* name;
   int         i;
   
   assert(tuple_is_valid(pattern));
   assert(tuple_is_valid(values));
   
   assert(tuple_get_dim(pattern) == tuple_get_dim(values));

   local_new_frame();
   
   for(i = 0; i < tuple_get_dim(pattern); i++)
   {
      elem = tuple_get_elem(pattern, i);

      if (elem_get_type(elem) == ELEM_NAME)
      {
         name = elem_get_name(elem);
         elem = tuple_get_elem(values, i);

         assert(elem_get_type(elem) != ELEM_NAME);

         local_new(name, elem);
      }
   }
}

void local_print_all(FILE* fp)
{
   const Local* local;

   for(local = anchor.next; local != NULL; local = local->next)
   {
      if (local->element == NULL)
         fprintf(fp, "New Frame");
      else
      {
         fprintf(fp, "%s = ", local->name);
         elem_print(fp, local->element, TRUE);
      }
      fprintf(fp, "\n");
   }
}

char* local_tostrall()
{
   const Local* local;
   unsigned int size = LOCAL_STR_SIZE;
   unsigned int len  = 1; /* fuer die '\0' */
   char*        str  = malloc(size);
   char*        selem;
   unsigned int selemlen;
   Bool         after_elem = FALSE;

   assert(str != NULL);

   str[0] = '\0';

   for(local = anchor.next; local != NULL; local = local->next)
   {
      /* Frame ?
       */
      if (local->element == NULL)
      {
         selem      = strdup(";");
         selemlen   = 1;
         after_elem = FALSE;
      }
      else
      {
         selem    = elem_tostr(local->element);
         selemlen = strlen(selem) + (after_elem ? 1 : 0);
      }
      if (len + selemlen >= size)
      {
         size += selemlen + LOCAL_STR_SIZE;
         str   = realloc(str, size);

         assert(str != NULL);
      }
      assert(len + selemlen < size);

      strcat(str, after_elem ? "@" : "");
      strcat(str, selem);

      free(selem);
      
      len += selemlen;
      
      after_elem = (local->element != NULL);
   }
   return str;
}


