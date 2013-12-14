/* $Id: elem.c,v 1.33 2012/07/29 15:09:26 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: elem.c                                                        */
/*   Name....: Element Functions                                             */
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
#include "bool.h"
#include "mshell.h"
#include "numb.h"
#include "elem.h"
#include "mme.h"
#include "strstore.h"

#define ELEM_STORE_SIZE  1000
#define ELEM_SID         0x456c656d

typedef union element_value    ElemValue;
typedef struct element_storage ElemStore;

union element_value
{
   Numb*       numb;
   const char* strg;
   const char* name;
   Elem*       next;
};

struct element
{
   SID
   ElemType  type;
   ElemValue value;
};

struct element_storage
{
   Elem*       begin;
   ElemStore*  next;
};

static ElemStore* store_anchor = NULL;
static Elem*      store_free   = NULL;
static int        store_count  = 0;

static void extend_storage(void)
{
   ElemStore* store = calloc(1, sizeof(*store));
   Elem*      elem;
   int        i;
   
   assert(store != NULL);
   
   store->begin = malloc(ELEM_STORE_SIZE * sizeof(*store->begin));
   store->next  = store_anchor;
   store_anchor = store;

   for(i = 0; i < ELEM_STORE_SIZE - 1; i++)
   {
      elem             = &store->begin[i];
      elem->type       = ELEM_FREE;
      elem->value.next = &store->begin[i + 1];
      SID_set(elem, ELEM_SID);
      assert(elem_is_valid(elem));
   }
   elem             = &store->begin[i];
   elem->type       = ELEM_FREE;
   elem->value.next = store_free;
   SID_set(elem, ELEM_SID);
   assert(elem_is_valid(elem));
   
   store_free       = &store->begin[0];
   
   assert(store->begin != NULL);
   assert(store_anchor != NULL);
   assert(store_free   != NULL);
}

static Elem* new_elem(void)
{
   Elem* elem;
   
   if (store_free == NULL)
      extend_storage();

   assert(store_free != NULL);

   elem       = store_free;
   store_free = elem->value.next;

   store_count++;
   
   assert(elem->type == ELEM_FREE);
   assert(elem_is_valid(elem));
   
   return elem;
}

void elem_init()
{
}

void elem_exit()
{
   ElemStore* store;
   ElemStore* next;
   
   if (store_count != 0)
      printf("Elem store count %d\n", store_count);
   
   for(store = store_anchor; store != NULL; store = next)
   {
#if 0 /* only for debugging */
      int i;
      
      for(i = 0; i < ELEM_STORE_SIZE - 1; i++)
      {
         elem_print(stderr, &store->begin[i], TRUE);
         fprintf(stderr, "\n");
      }
#endif
      next = store->next;

      free(store->begin);
      free(store);
   }
   store_anchor = NULL;
   store_free   = NULL;
   store_count  = 0;
}

Elem* elem_new_numb(const Numb* numb)
{
   Elem* elem = new_elem();
   
   assert(elem != NULL);

   elem->type       = ELEM_NUMB;
   elem->value.numb = numb_copy(numb);
   
   return elem;
}

Elem* elem_new_strg(const char* strg)
{
   Elem* elem = new_elem();

   assert(strg != NULL);
   assert(elem != NULL);

   elem->type       = ELEM_STRG;
   elem->value.strg = strg;

   return elem;
}

Elem* elem_new_name(const char* name)
{
   Elem* elem = new_elem();

   assert(name != NULL);
   assert(elem != NULL);

   elem->type       = ELEM_NAME;
   elem->value.strg = name;

   return elem;
}

void elem_free(Elem* elem)
{
   assert(elem_is_valid(elem));

   if (elem->type == ELEM_NUMB)
      numb_free(elem->value.numb);
   
   elem->type       = ELEM_FREE;
   elem->value.next = store_free;
   store_free       = elem;
   store_count--;
}

Bool elem_is_valid(const Elem* elem)
{
   return elem != NULL && SID_ok(elem, ELEM_SID);
}

Elem* elem_copy(const Elem* source)
{
   Elem* elem = new_elem();

   assert(elem_is_valid(source));
   assert(elem_is_valid(elem));

   if (source->type != ELEM_NUMB)
      *elem = *source;
   else
   {
      elem->type       = ELEM_NUMB;
      elem->value.numb = numb_copy(source->value.numb);
   }
   return elem;
}

/* 0 wenn gleich, sonst != 0
 */
Bool elem_cmp(const Elem* elem_a, const Elem* elem_b)
{
   assert(elem_is_valid(elem_a));
   assert(elem_is_valid(elem_b));
   assert(elem_a->type != ELEM_ERR);
   assert(elem_b->type != ELEM_ERR);

   /* Auf die schnelle vorweg.
    */
   if (elem_a == elem_b)
      return FALSE;

   if (elem_a->type != elem_b->type)
   {
      fprintf(stderr,
         "*** Error 160: Comparison of elements with different types ");
      elem_print(stderr, elem_a, TRUE);
      fprintf(stderr, " / ");
      elem_print(stderr, elem_b, TRUE);
      fputc('\n', stderr);
      zpl_exit(EXIT_FAILURE);
   }
   assert(elem_a->type == elem_b->type);
   
   if (elem_a->type == ELEM_STRG)
      return strcmp(elem_a->value.strg, elem_b->value.strg) != 0;

   assert(elem_a->type == ELEM_NUMB);

   return !numb_equal(elem_a->value.numb, elem_b->value.numb);
}

ElemType elem_get_type(const Elem* elem)
{
   assert(elem_is_valid(elem));
   
   return elem->type;
}

const Numb* elem_get_numb(const Elem* elem)
{
   assert(elem_is_valid(elem));
   assert(elem->type == ELEM_NUMB);
   
   return elem->value.numb;
}

const char* elem_get_strg(const Elem* elem)
{
   assert(elem_is_valid(elem));
   assert(elem->type       == ELEM_STRG);
   assert(elem->value.strg != NULL);
   
   return elem->value.strg;
}

const char* elem_get_name(const Elem* elem)
{
   assert(elem_is_valid(elem));
   assert(elem->type       == ELEM_NAME);
   assert(elem->value.name != NULL);
   
   return elem->value.name;
}

void elem_print(FILE* fp, const Elem* elem, Bool use_quotes)
{
   assert(elem_is_valid(elem));

   switch(elem->type)
   {
   case ELEM_NUMB :
      fprintf(fp, "%.16g", numb_todbl(elem->value.numb));
      break;
   case ELEM_STRG :
      fprintf(fp, use_quotes ? "\"%s\"" : "%s", elem->value.strg);
      break;
   case ELEM_NAME :
      fprintf(fp, "%s", elem->value.name);
      break;
   case ELEM_FREE :
      fprintf(fp, "Unused Elem!");
      break;
   default :
      abort();
   }
}

unsigned int elem_hash(const Elem* elem)
{
   unsigned int hcode = 0;
   
   switch(elem->type)
   {
   case ELEM_NUMB :
      hcode = numb_hash(elem->value.numb);
      break;
   case ELEM_STRG :
      hcode = str_hash(elem->value.strg);
      break;
   case ELEM_NAME :
      hcode = str_hash(elem->value.name);
      break;
   case ELEM_FREE :
   default :
      abort();
   }
   return hcode;
}

char* elem_tostr(const Elem* elem)
{
   char* str;
   
   assert(elem_is_valid(elem));

   switch(elem->type)
   {
   case ELEM_NUMB :
      str = malloc(32);
      
      assert(str != NULL);
      
      sprintf(str, "%.16g", numb_todbl(elem->value.numb));
      break;
   case ELEM_STRG :
      str = strdup(elem->value.strg);
      break;
   case ELEM_NAME :
      str = strdup(elem->value.name);
      break;
   case ELEM_FREE :
   default :
      abort();
   }
   assert(str != NULL);

   return str;
}


