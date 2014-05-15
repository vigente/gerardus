/* $Id: heap.c,v 1.13 2012/07/29 15:09:26 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: heap.c                                                        */
/*   Name....: Heap Functions                                                */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Copyright (C) 2006-2012 by Thorsten Koch <koch@zib.de>
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
#include "set.h"
#include "symbol.h"
#include "entry.h"
#include "heap.h"

#define HEAP_SID  0x48656170

enum heap_type
{
   HEAP_ERR = 0, HEAP_ENTRY = 1
};

typedef enum heap_type HeapType;

struct heap
{
   SID
   HeapType  type;
   int       size;   /**< Maximale Anzahl von Elementen im Heap              */
   int       used;   /**< Aktuelle Anzahl von Elementen im Heap              */
   HeapData* data;
   HeapCmp   data_cmp;
};

static void heap_print(FILE* fp, const Heap* heap)
{
   int i;
   
   for(i = 0; i < heap->used; i++)
   {
      fprintf(fp, "%3d ", i);
      switch(heap->type)
      {
      case HEAP_ENTRY :
         entry_print(fp, heap->data[i].entry);
         break;
      default :
         abort();
      }
      fprintf(fp, "\n");
   }
}

Bool heap_is_valid(const Heap* heap)
{
   HeapData* data;
   int       i;
   
   if (  heap           == NULL
      || heap->type     == HEAP_ERR
      || heap->data     == NULL
      || heap->data_cmp == NULL
      || heap->size     <= 0
      || heap->used     <  0
      || heap->used     >  heap->size)
      return FALSE;

   data = heap->data;
   
   /* Heap property
    */
   for(i = 0; i < heap->used / 2; i++)
   {
      if ((*heap->data_cmp)(data[i], data[i + i]) > 0)
      {
         heap_print(stderr, heap);
         return FALSE;
      }
      if (i + i + 1 < heap->used && (*heap->data_cmp)(data[i], data[i + i + 1]) > 0)
      {
         heap_print(stderr, heap);
         return FALSE;
      }
   }
   return TRUE;
}

static Heap* heap_new(
   HeapType type,
   int      size,
   HeapCmp  data_cmp)
{
   Heap* heap = calloc(1, sizeof(*heap));

   assert(type     == HEAP_ENTRY);
   assert(size     >  0);
   assert(data_cmp != NULL);
   assert(heap     != NULL);

   heap->type     = type;
   heap->size     = size;
   heap->used     = 0;
   heap->data     = calloc((size_t)heap->size, sizeof(*heap->data));
   heap->data_cmp = data_cmp;
   
   SID_set(heap, HEAP_SID);
   assert(heap_is_valid(heap));
   
   return heap;
}

Heap* heap_new_entry(
   int     size,
   HeapCmp heap_entry_cmp)
{
   assert(size           >  0);
   assert(heap_entry_cmp != NULL);
   
   return heap_new(HEAP_ENTRY, size, heap_entry_cmp);
}

void heap_free(Heap* heap)
{
   int i;
   
   assert(heap_is_valid(heap));

   for(i = 0; i < heap->used; i++)
   {
      switch(heap->type)
      {
      case HEAP_ENTRY :
         entry_free(heap->data[i].entry);
         break;
      default :
         abort();
      }
   }
   free(heap->data);
   free(heap);
}

static void swap_entries(const Heap* heap, int i, int j)
{
   HeapData* data = heap->data;
   HeapData  t;

   t       = data[i];
   data[i] = data[j];
   data[j] = t;
}

/* Bewegt einen Eintrag weiter nach unten/hinten im Vektor
 * bis die Teilsortierung des Baumes wieder hergestellt ist.
 */
static void sift_down(
   const Heap* heap,
   int         current)
{
   HeapData* data = heap->data;
   int       child;

   /* Heap shift down
    * (Oberstes Element runter und korrigieren)
    */         
   child = current * 2;

   if (child + 1 < heap->used)
      if ((*heap->data_cmp)(data[child + 1], data[child]) < 0)
         child++;

   while(child < heap->used && (*heap->data_cmp)(data[current], data[child]) > 0)
   {
      swap_entries(heap, current, child);

      current = child;
      child  += child;
      
      if (child + 1 < heap->used)
         if ((*heap->data_cmp)(data[child + 1], data[child]) < 0)
            child++;
   }
}

/* Bewegt einen Eintrag weiter hoch/nach unten im Vektor
 * bis die Teilsortierung des Baumes wieder hergestellt ist.
 */
static void sift_up(
   const Heap* heap,
   int         current)
{
   HeapData* data   = heap->data;
   int       parent = current / 2;
   
   /* Heap sift up 
    */
   while(current > 0 && (*heap->data_cmp)(data[parent], data[current]) > 0)
   {
      swap_entries(heap, current, parent);
      current = parent;
      parent /= 2;
   }
}

/* Sortiert einen Eintrag in den Heap ein.
 */
void heap_push_entry(
   Heap*  heap,
   Entry* entry)
{
   assert(heap_is_valid(heap));
   assert(entry_is_valid(entry));
   assert(heap->used <  heap->size);

   heap->data[heap->used].entry = entry;
   
   heap->used++;

   sift_up(heap, heap->used - 1);

   assert(heap_is_valid(heap));
}

/* Holt den Eintrag mit dem kleinsten Wert aus dem Heap heraus.
 */
Entry* heap_pop_entry(
   Heap* heap)
{
   Entry* entry;
   
   assert(heap_is_valid(heap));
   assert(heap->used > 0);
   assert(heap->type == HEAP_ENTRY);
   
   /* Heap shift down
    * (Oberstes Element runter und korrigieren)
    */         
   entry = heap->data[0].entry;

   heap->data[0].entry = NULL;
   
   heap->used--;

   swap_entries(heap, 0, heap->used);
   
   sift_down(heap, 0);

   assert(heap_is_valid(heap));
      
   return entry;
}

const Entry* heap_top_entry(
   const Heap* heap)
{
   assert(heap_is_valid(heap));
   assert(heap->used > 0);
   assert(heap->type == HEAP_ENTRY);
   
   return heap->data[0].entry;
}

Bool heap_is_full(const Heap* heap)
{
   assert(heap_is_valid(heap));

   return heap->used == heap->size;
}

Bool heap_is_empty(const Heap* heap)
{
   assert(heap_is_valid(heap));

   return heap->used == 0;
}


