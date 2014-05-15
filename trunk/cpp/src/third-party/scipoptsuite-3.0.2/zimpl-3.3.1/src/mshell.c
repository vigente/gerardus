/* $Id: mshell.c,v 1.19 2012/07/29 15:09:28 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: mshell.c                                                      */
/*   Name....: Memory Allocation Shell                                       */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Copyright (C) 2003-2012 by Thorsten Koch <koch@zib.de>
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
/*
 * This is base on the routines from Jim Schimandle,
 * published in Dr. Dobbs Journal #167 09/90 p.110+.
 */
#include <sys/types.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "lint.h"

#define _MSHELL_C_

#include "mshell.h"

#if !defined(NO_MSHELL) 

#define ALIGN_SIZE   sizeof(double)
#define MEMTAG1      0xa55a
#define MEMTAG2      0xd88d
#define OLDTAG       0xb66b
#define ENDTAG       0xc77c

typedef struct memnod
{
   unsigned short tag1;
   unsigned long  size;
   struct memnod* next;
   struct memnod* prev;
   char const*    file;
   int            line;
   unsigned short tag2;
} MHDR;

#define MHDR_NIL        ((MHDR*)0)

#define HDR_SIZE        sizeof(MHDR)
#define RESERVE_SIZE    (((HDR_SIZE + (ALIGN_SIZE - 1)) \
                         / ALIGN_SIZE) * ALIGN_SIZE)
#define CLIENT_2_HDR(a) ((MHDR*)(((char*)(a)) - RESERVE_SIZE))
#define HDR_2_CLIENT(a) ((void*)(((char*)(a)) + RESERVE_SIZE))

static unsigned long mem_size = 0;
static unsigned long mem_maxi = 0;
static MHDR*         memlist  = MHDR_NIL;

#define Mem_tag_err(a, b, c, d)  mem_tag_err((a), (b), (c), (d), file, line)

static size_t mem_alloc_size(
   size_t size)
{
   size += RESERVE_SIZE + ALIGN_SIZE + ALIGN_SIZE - 1;
   
   return((size / ALIGN_SIZE) * ALIGN_SIZE);
}
   
static int* mem_endptr(
   MHDR* p)
{
   size_t offset;
   
   assert(p != NULL);
   
   offset = ((RESERVE_SIZE + p->size + ALIGN_SIZE - 1) / ALIGN_SIZE)
      * ALIGN_SIZE;
   
   return((int*)((char*)p + offset));
}  
   
static void mem_add_list(
   MHDR* p)
{
   assert((p != NULL) && (p->tag1 == MEMTAG1) && (p->tag2 == MEMTAG2));

   p->next = memlist;
   p->prev = MHDR_NIL;
   
   if (memlist != MHDR_NIL)
      memlist->prev = p;
      
   memlist = p;
   
   if (mem_size > mem_maxi)
      mem_maxi = mem_size;
}

static void mem_del_list(
   MHDR *p)
{
   assert((p != NULL) && (p->tag1 == MEMTAG1) && (p->tag2 == MEMTAG2));

   *mem_endptr(p) = ~ENDTAG;
   p->tag1        = OLDTAG;
   p->tag2        = OLDTAG;
   mem_size      -= p->size;
   
   if (p->next != MHDR_NIL)                  
      p->next->prev = p->prev;
   
   if (p->prev != MHDR_NIL)
      p->prev->next = p->next;
   else
      memlist = p->next;
}

static void mem_tag_err(
   MHDR*       p,
   int         typ,
   const char* file1,
   const int   line1,
   const char* file2,
   const int   line2)
{
   const char* const errtyp[] = { "pre1", "pre2", "post" };
   
   assert((typ >= 0) && (typ <= 2));
   
   (void)fprintf(stderr, "Memory tag error (%s) - %lx - %s(%d) at %s(%d)\n",
      errtyp[typ],
      (unsigned long)p,
      file1,
      line1,
      file2,
      line2);

   mem_display(stderr);

   abort(); 
}   

static void mem_valid(
   MHDR*       p,
   const char* file,
   const int   line)
{
   if (p->tag1 != MEMTAG1)
   {
      if (p->tag1 == OLDTAG)
         Mem_tag_err(p, 0, p->file, p->line);

      Mem_tag_err(p, 0, "Unknown", 0);
   }
   if (p->tag2 != MEMTAG2)
   {
      if (p->tag1 == OLDTAG)
         Mem_tag_err(p, 1, p->file, p->line);

      Mem_tag_err(p, 1, "Unknown", 0);
   }
   if (*mem_endptr(p) != ENDTAG)
      Mem_tag_err(p, 2, p->file, p->line);
}
      
void* mem_malloc(
   size_t       size,
   const char*  file,
   const int    line) 
{
   const char* errmsg1 =
      "mem_malloc(size=%u, file=%s, line=%d): out of memory\n";
   const char* errmsg2 =
      "mem_malloc(size=%u, file=%s, line=%d): zero size\n";

   MHDR* p;
   size_t alloc_size;

   if (size == 0)
   {
      fprintf(stderr, errmsg2, size, file, line);
      exit(EXIT_FAILURE);
   }
   alloc_size = mem_alloc_size(size);

   assert(alloc_size > 0);

   if ((p = malloc(alloc_size)) == MHDR_NIL)
   {
      fprintf(stderr, errmsg1, size, file, line);
      exit(EXIT_FAILURE);
   }
   p->tag1        = MEMTAG1;
   p->tag2        = MEMTAG2;
   p->size        = size;
   mem_size      += size;
   p->file        = file;
   p->line        = line;
   *mem_endptr(p) = ENDTAG;

   mem_add_list(p);

   return(HDR_2_CLIENT(p));
} 

void* mem_calloc(
   size_t       item,
   size_t       size,
   const char*  file,
   const int    line)
{
   const char* errmsg1 =
      "mem_calloc(item=%u, size=%u, file=%s, line=%d): out of memory\n";
   const char* errmsg2 =
      "mem_calloc(item=%u, size=%u, file=%s, line=%d): zero item/size\n";

   MHDR* p;
   
   if (item == 0 || size == 0)
   {
      fprintf(stderr, errmsg2, item, size, file, line);
      abort();
      exit(EXIT_FAILURE);
   }
   if ((p = calloc(mem_alloc_size(size * item), sizeof(char))) == MHDR_NIL)
   {
      fprintf(stderr, errmsg1, item, size, file, line);
      exit(EXIT_FAILURE);
   }
   p->tag1        = MEMTAG1;
   p->tag2        = MEMTAG2;
   p->size        = size * item;
   mem_size      += size * item;
   p->file        = file;
   p->line        = line;
   *mem_endptr(p) = ENDTAG;
   
   mem_add_list(p);
   
   return(HDR_2_CLIENT(p));
} 

void* mem_realloc(
   void*        ptr,
   size_t       size,
   const char*  file,
   const int    line)
{
   const char* errmsg1 =
      "mem_realloc(size=%u, file=%s, line=%d): out of memory\n";
   const char* errmsg2 =
      "mem_realloc(file=%s, line=%d): null pointer\n";
   const char* errmsg3 =
      "mem_realloc(size=%u, file=%s, line=%d): zero size\n";

   MHDR* p;

   if (ptr == NULL)
   {
      fprintf(stderr, errmsg2, file, line);
      exit(EXIT_FAILURE);
   }
   p = CLIENT_2_HDR(ptr);

   mem_valid(p, file, line);   
   mem_del_list(p);

   if (size == 0)
   {
      fprintf(stderr, errmsg3, size, file, line);
      exit(EXIT_FAILURE);
   }
   if ((p = realloc(p, mem_alloc_size(size))) == MHDR_NIL)
   {
      fprintf(stderr, errmsg1, size, file, line);
      exit(EXIT_FAILURE);
   }
   p->tag1        = MEMTAG1;
   p->tag2        = MEMTAG2;
   p->size        = size;
   mem_size      += size;
   p->file        = file;
   p->line        = line;
   *mem_endptr(p) = ENDTAG;
   
   mem_add_list(p);
   
   return(HDR_2_CLIENT(p));
}

char* mem_strdup(
   const char* str,
   const char* file,
   const int   line)
{
   const char* errmsg1 = "mem_strdup(file=%s, line=%d): null pointer\n";

   if (str == NULL)
   {
      (void)fprintf(stderr, errmsg1, file, line);
      exit(EXIT_FAILURE);
   }
   return(strcpy(mem_malloc((size_t)(strlen(str) + 1), file, line), str));
}
      
void mem_free(
   void*       ptr,
   const char* file,
   const int   line)
{
   const char *errmsg = "mem_free(file=%s, line=%d): null pointer\n";

   MHDR* p;

   if (ptr == NULL)
   {
      (void)fprintf(stderr, errmsg, file, line);
      abort();
   }
   p = CLIENT_2_HDR(ptr);

   mem_valid(p, file, line);   
   mem_del_list(p);

   free(p);
}

void mem_hide_x(
   void*       ptr,
   const char* file,
   const int   line)
{
   const char *errmsg = "mem_checkout(file=%s, line=%d): null pointer\n";

   MHDR* p;

   if (ptr == NULL)
   {
      (void)fprintf(stderr, errmsg, file, line);
      abort();
   }
   p = CLIENT_2_HDR(ptr);

   mem_valid(p, file, line);   
   mem_del_list(p);

   /* Size nicht verringern
    */
   mem_size += p->size;
}

size_t mem_used()
{
   return(mem_size);
}

void mem_maximum(
   FILE* fp)
{
   (void)fprintf(fp, "Maximum amount of memory used = %lu bytes\n",
      mem_maxi);
}

void mem_display(
   FILE* fp)
{
   MHDR* p;
      
   (void)fprintf(fp, "\nAddress     Size  File(Line) - total size %lu\n",
      mem_size);
   
   for(p = memlist; p != MHDR_NIL; p = p->next)
   {
      (void)fprintf(fp, "%8lx  %6lu  %s(%d) %s %s\n",
         (unsigned long)p,
         p->size,
         p->file,
         p->line,
         ((p->tag1 != MEMTAG1) || (p->tag2 != MEMTAG2)) ? "Invalid" : "",
         (*mem_endptr(p) == ENDTAG) ? "ok" : "Clobbered");
   } 
   mem_maximum(fp);
}      

void mem_check_x(
   const void* ptr,
   const char* file,
   const int   line)
{
   MHDR* p = CLIENT_2_HDR(ptr);
   
   mem_valid(p, file, line);
}      

void mem_check_all_x(
   const char* file,
   const int   line)
{
   MHDR* p;
      
   for(p = memlist; p != MHDR_NIL; p = p->next)
      mem_valid(p, file, line);
}      

#else /* NO_MSHELL */

void* mem_malloc(
   size_t       size,
   const char*  file,
   const int    line) 
{
   const char* errmsg1 =
      "mem_malloc(size=%u, file=%s, line=%d): out of memory\n";

   void* p;

   assert(size > 0);
   
   if (NULL == (p = malloc(size)))
   {
      fprintf(stderr, errmsg1, size, file, line);
      exit(EXIT_FAILURE);
   }
   return p;
}

void* mem_calloc(
   size_t       item,
   size_t       size,
   const char*  file,
   const int    line)
{
   const char* errmsg1 =
      "mem_calloc(item=%u, size=%u, file=%s, line=%d): out of memory\n";

   void* p;

   assert(item > 0);
   assert(size > 0);
   
   if (NULL == (p = calloc(item, size)))
   {
      fprintf(stderr, errmsg1, item, size, file, line);
      exit(EXIT_FAILURE);
   }
   return p;
}

void* mem_realloc(
   void*        ptr,
   size_t       size,
   const char*  file,
   const int    line)
{
   const char* errmsg1 =
      "mem_realloc(size=%u, file=%s, line=%d): out of memory\n";

   void* p;

   assert(ptr  != NULL);
   assert(size >  0);
   
   if (NULL == (p = realloc(ptr, size)))
   {
      fprintf(stderr, errmsg1, size, file, line);
      exit(EXIT_FAILURE);
   }
   return p;
}

char* mem_strdup(
   const char* str,
   const char* file,
   const int   line)
{
   const char* errmsg1 =
      "mem_strdup(size=%u, file=%s, line=%d): out of memory\n";

   char* s;
   
   assert(str != NULL);
   
   if (NULL == (s = strdup(str)))
   {
      fprintf(stderr, errmsg1, strlen(str), file, line);
      exit(EXIT_FAILURE);
   }
   return s;
}

void mem_free(
   void*       ptr,
   const char* file,
   const int   line)
{
   const char *errmsg = "mem_free(file=%s, line=%d): null pointer\n";

#ifndef NDEBUG
   if (ptr == NULL)
   {
      fprintf(stderr, errmsg, file, line);
      abort();
   }
#endif
   free(ptr);
}

#endif /* !NO_MSHELL */

/* ------------------------------------------------------------------------- */
/* Emacs Local Variables:                                                    */
/* Emacs mode:c                                                              */
/* Emacs c-basic-offset:3                                                    */
/* Emacs tab-width:8                                                         */
/* Emacs indent-tabs-mode:nil                                                */
/* Emacs End:                                                                */
/* ------------------------------------------------------------------------- */


