/* $Id: mono.c,v 1.15 2012/07/29 15:09:27 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: mono.c                                                        */
/*   Name....: Monom Functions                                               */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Copyright (C) 2007-2012 by Thorsten Koch <koch@zib.de>
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

/* #define TRACE 1 */

#include "bool.h"
#include "mshell.h"
#include "ratlptypes.h"
#include "numb.h"
#include "elem.h"
#include "tuple.h"
#include "mme.h"
#include "entry.h"
#include "mono.h"

#define MONO_SID         0x4d6f6e6f
#define MOEL_SID         0x4d6f456c

Mono* mono_new(const Numb* coeff, const Entry* entry, MFun fun)
{
   Mono* mono = calloc(1, sizeof(*mono));

   Trace("mono_new");

   assert(mono != NULL);
   assert(entry_is_valid(entry));
   assert(entry_get_type(entry) == SYM_VAR);

   mono->count       = 1;
   mono->coeff       = numb_copy(coeff);
   mono->fun         = fun;
   mono->first.entry = entry_copy(entry);
   mono->first.next  = NULL;

   SID_set(mono, MONO_SID);
   SID_set2(mono->first, MOEL_SID);
   
   assert(mono_is_valid(mono));

   return mono;
}

#ifndef NDEBUG
Bool mono_is_valid(const Mono* mono)
{
   const MonoElem* e;
   int             count = 1;
   
   if (mono == NULL
    || !SID_ok(mono, MONO_SID)
    || !SID_ok2(mono->first, MOEL_SID)
    || mono->count < 1)
      abort();

   mem_check(mono);

   assert(entry_is_valid(mono->first.entry));
   
   for(e = mono->first.next; e != NULL; e = e->next)
   {      
      count++;
      
      mem_check(e);

      if (!SID_ok(e, MOEL_SID))
         abort();
      
      assert(entry_is_valid(e->entry));
      assert(entry_get_type(e->entry) == SYM_VAR);
   }
   if (count != mono->count)
      abort();
   
   return TRUE;
}
#endif

void mono_free(Mono* mono)
{
   MonoElem* e;
   MonoElem* q;

   Trace("mono_free");

   assert(mono_is_valid(mono));
   
   for(e = mono->first.next; e != NULL; e = q)
   {
      q = e->next;
      entry_free(e->entry);
      SID_del(e);
      free(e);
   }   
   entry_free(mono->first.entry);
   numb_free(mono->coeff);
   SID_del2(mono->first);
   SID_del(mono);
   free(mono);
}

void mono_mul_entry(
   Mono*        mono,
   const Entry* entry)
{
   MonoElem* e;
   Var*      var;
   MonoElem* last;
   
   Trace("mono_add_elem");

   assert(mono_is_valid(mono));
   assert(entry_is_valid(entry));
   assert(entry_get_type(entry) == SYM_VAR);

   var = entry_get_var(entry);

   /* ??? This ensures that if the same variable is to come several times,
    * all of them come together, i.e. yxy is not allowed, yyx would be ok.
    * Is there any reason to do this?
    */
   for(e = &mono->first; e != NULL; e = e->next)
   {
      last = e;

      assert(entry_is_valid(e->entry));

      if (var == entry_get_var(e->entry))
         break;
   }
   assert(last != NULL);
       
   e = calloc(1, sizeof(*e));
   
   e->entry   = entry_copy(entry);
   e->next    = last->next;
   SID_set(e, MOEL_SID);
   last->next = e;
   mono->count++;
   
   assert(mono_is_valid(mono));
}

Mono* mono_copy(const Mono* mono)
{
   Mono*     mnew;
   MonoElem* e;
   
   assert(mono_is_valid(mono));

   mnew = mono_new(mono->coeff, mono->first.entry, mono->fun);
   
   for(e = mono->first.next; e != NULL; e = e->next)
      mono_mul_entry(mnew, e->entry);

   assert(mono_is_valid(mnew));

   return mnew;
}
     
void mono_mul_coeff(const Mono* mono, const Numb* value)
{
   Trace("mono_mul_coeff");

   assert(mono_is_valid(mono));
   assert(numb_is_valid(value));

   numb_mul(mono->coeff, value);
}

void mono_add_coeff(const Mono* mono, const Numb* value)
{
   Trace("mono_add_coeff");

   assert(mono_is_valid(mono));
   assert(numb_is_valid(value));

   numb_add(mono->coeff, value);
}

unsigned int mono_hash(const Mono* mono)
{
   size_t          hcode = 0;
   const MonoElem* e;
   
   assert(mono_is_valid(mono));   
   
   for(e = &mono->first; e != NULL; e = e->next)
      hcode += ((size_t)entry_get_var(e->entry)) >> 2;

   return DISPERSE((unsigned int)hcode);
}

/** Checks whether two monoms cosist of the same variables.
 */  
#if 1
/* We assume (I think it is true that there is only one distinct entry per var
 */
Bool mono_equal(const Mono* ma, const Mono* mb)
{
   const MonoElem* ea;
   const MonoElem* eb;
   const Entry*    entry_a;
   
   assert(mono_is_valid(ma));   
   assert(mono_is_valid(mb));   

   if (ma->count != mb->count)
      return FALSE;

   if (ma->count == 1 && (ma->first.entry != mb->first.entry))
      return FALSE;

   for(ea = &ma->first; ea != NULL; ea = ea->next)
   {
      assert(entry_is_valid(ea->entry));
      assert(entry_get_type(ea->entry) == SYM_VAR);
           
      entry_a = ea->entry;

      for(eb = &mb->first; eb != NULL; eb = eb->next)
         if (entry_a == eb->entry)
            break;

      if (eb == NULL)
         return FALSE;
      
      /* Now all variables of a kind are consecutive 
       */
      while(ea->next != NULL && ea->next->entry == entry_a) 
      {
         if (eb->next == NULL || eb->next->entry != entry_a)
            return FALSE;
               
         ea = ea->next; /*lint !e850 loop index variable is modified in body of the loop */
         eb = eb->next;               
      }
   } 
   return TRUE;
}
#else /* old */
Bool mono_equal(const Mono* ma, const Mono* mb)
{
   const MonoElem* ea;
   const MonoElem* eb;
   Var*            var_a;
   
   assert(mono_is_valid(ma));   
   assert(mono_is_valid(mb));   

   if (ma->count != mb->count)
      return FALSE;

   if (ma->count == 1 && (entry_get_var(ma->first.entry) != entry_get_var(mb->first.entry)))
      return FALSE;

   for(ea = &ma->first; ea != NULL; ea = ea->next)
   {
      assert(entry_is_valid(ea->entry));

      var_a = entry_get_var(ea->entry);
      
      for(eb = &mb->first; eb != NULL; eb = eb->next)
         if (var_a == entry_get_var(eb->entry))
            break;

      if (eb == NULL)
         return FALSE;
      
      /* Now all variables of a kind are consecutive 
       */
      while(ea->next != NULL && entry_get_var(ea->next->entry) == var_a)
      {
         if (eb->next == NULL || entry_get_var(eb->next->entry) != var_a)
            return FALSE;
               
         ea = ea->next;
         eb = eb->next;               
      }
   }
   return TRUE;
}
#endif

Mono* mono_mul(const Mono* ma, const Mono* mb)
{
   Mono*           mono;
   const MonoElem* eb;
   
   assert(mono_is_valid(ma));   
   assert(mono_is_valid(mb));   

   mono = mono_copy(ma);

   numb_mul(mono->coeff, mb->coeff);

   for(eb = &mb->first; eb != NULL; eb = eb->next)
   {
      assert(entry_is_valid(eb->entry));
      
      mono_mul_entry(mono, eb->entry);
   }
   assert(mono_is_valid(mono));

   return mono;
}

void mono_neg(Mono* mono)
{
   assert(mono_is_valid(mono));

   numb_neg(mono->coeff);
}

Bool mono_is_linear(const Mono* mono)
{
   assert(mono_is_valid(mono));

   return mono->count == 1 && mono->fun == MFUN_NONE;
}

int mono_get_degree(const Mono* mono)
{
   assert(mono_is_valid(mono));

   return mono->count;
}

const Numb* mono_get_coeff(const Mono* mono)
{
   assert(mono_is_valid(mono));

   return mono->coeff;
}

void mono_set_function(Mono* mono, MFun f)
{
   assert(mono_is_valid(mono));

   mono->fun = f;

   assert(mono_is_valid(mono));
}

MFun mono_get_function(const Mono* mono)
{
   assert(mono_is_valid(mono));

   return mono->fun;
}

Var* mono_get_var(const Mono* mono, int idx)
{
   const MonoElem* e = &mono->first;
   
   assert(mono_is_valid(mono));
   assert(mono->count > 0);
   assert(idx >= 0);
   assert(idx <= mono->count);

   if (idx > 0)
   {
      for(e = e->next; --idx > 0; e = e->next) /*lint !e441 loop variable 'e' not found in 2nd expression */
         assert(e != NULL);

      assert(e != NULL);
   }
   assert(entry_is_valid(e->entry));

   return entry_get_var(e->entry);
}

#ifndef NDEBUG
void mono_print(FILE* fp, const Mono* mono, Bool print_symbol_index)
{
   const MonoElem* e;

   assert(mono_is_valid(mono));
   
   if (numb_equal(mono->coeff, numb_one()))
      fputc('+', fp);
   else
   {
      if (numb_cmp(mono->coeff, numb_zero()) >= 0)
         fprintf(fp, "+ %g", numb_todbl(mono->coeff));
      else
         fprintf(fp, "- %g", -numb_todbl(mono->coeff));
   }
   fputc(' ', fp);
   
   for(e = &mono->first; e != NULL; e = e->next)
   {      
      entry_print(fp, e->entry);

      if (print_symbol_index)
         tuple_print(fp, entry_get_tuple(e->entry));

      if (e->next != NULL) 
         fprintf(fp, " * ");
   }
}
#endif /* !NDEBUG */
