/* $Id: strstore2.c,v 1.10 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: strstore2.c                                                   */
/*   Name....: String Storage Functions                                      */
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

#include "lint.h"
#include "mshell.h"
#include "strstore.h"
#include "mme.h"

typedef struct string_storage StrStore;

struct string_storage
{
   char*      begin;
   int        size;
   int        used;
   StrStore*  next;
};

#define MIN_STR_STORE_SIZE 65536      /* 64k */
#define MAX_STR_STORE_SIZE 1073741824 /* 1G  */

static StrStore* store_anchor = NULL;
static int       store_size   = MIN_STR_STORE_SIZE;

static void extend_storage(void)
{
   StrStore* store;

   assert(store_size > 0);
   
   /* Since we will not allocate anything again in this block, resize it to fit.
    * (not clear if this is really usefull)
    */
   if (store_anchor != NULL)
   {
      store_anchor->begin =
         realloc(store_anchor->begin, store_anchor->used * sizeof(*store_anchor->begin));
      store_anchor->size = store_anchor->used;
   }
   store = calloc(1, sizeof(*store_anchor));

   assert(store != NULL);

   store->size  = store_size;
   store->used  = 0;
   store->next  = store_anchor;
   store->begin = calloc(store_size, sizeof(*store->begin));

   assert(store->begin != NULL);

   store_anchor = store;
}

const char* str_new(const char* s)
{
   char* t;
   int   len;

   assert(store_anchor != NULL);
   assert(s            != NULL);

   len = strlen(s) + 1;

   if (len > MAX_STR_STORE_SIZE)
   {
      fprintf(stderr, "*** Error 803: String too long %d > %d\n",
         len + 1, MAX_STR_STORE_SIZE); 

      zpl_exit(EXIT_FAILURE);
   }
   if (store_anchor->size - store_anchor->used < len)
   {
      /* Double the store_size at least once each time,
       * but more often in case it is not big enough to hold a very long
       * string.
       */
      if (store_size < MAX_STR_STORE_SIZE)
      {
         do
         {
            store_size *= 2;
         }
         while(len > store_size);
      }
      extend_storage();
   }
   assert(store_anchor->size - store_anchor->used >= len);

   t = &store_anchor->begin[store_anchor->used];

   store_anchor->used += len;
   
   return strcpy(t, s);
}

void str_init(void)
{
   extend_storage();
}

void str_exit(void)
{
   StrStore* p;
   StrStore* q;

   for(p = store_anchor; p != NULL; p = q)
   {
      q = p->next;
      free(p->begin);
      free(p);
   }
   store_anchor = NULL;
}

unsigned int str_hash(const char* s)
{
#if 0
   return (unsigned int)s;
#else
   unsigned int sum;
   int          i;
   
   for(sum = 0, i = 0; s[i] != '\0'; i++)
      sum = sum * 31 + s[i];

   return sum;
#endif
}

