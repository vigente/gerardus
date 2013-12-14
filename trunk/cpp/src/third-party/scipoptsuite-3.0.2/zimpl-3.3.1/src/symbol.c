/* $Id: symbol.c,v 1.36 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: symbol.c                                                      */
/*   Name....: Symbol Table Functions                                        */
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
#include "set.h"
#include "hash.h"
#include "stmt.h"
#include "symbol.h"

#define TEST_DUBLICATE   0

#define SYMBOL_SID  0x53796d62
#define SYMBOL_EXTEND_SIZE 100

struct symbol
{
   SID
   const char*  name;
   int          size;
   int          used;
   int          extend;
   SymbolType   type;
   Set*         set;
   Hash*        hash;
   Entry**      entry;
   Entry*       deflt;
   Symbol*      next;
};

static Symbol* anchor = NULL;

Symbol* symbol_new(
   const char*  name,
   SymbolType   type,
   const Set*   set,
   int          estimated_size,
   const Entry* deflt)
{
   Symbol* sym;

   assert(name           != NULL);
   assert(strlen(name)   >  0);
   assert(set            != NULL);
   assert(estimated_size >= 0);
   
   sym = calloc(1, sizeof(*sym));

   assert(sym != NULL);
   
   sym->name    = name;
   sym->type    = type;
   sym->size    = 1;
   sym->used    = 0;
   sym->extend  = SYMBOL_EXTEND_SIZE;
   sym->set     = set_copy(set);
   sym->hash    = hash_new(HASH_ENTRY, estimated_size);
   sym->entry   = calloc(1, sizeof(*sym->entry));
   sym->deflt   = (deflt != ENTRY_NULL) ? entry_copy(deflt) : ENTRY_NULL;
   sym->next    = anchor;
   anchor       = sym;

   assert(sym->entry != NULL);

   SID_set(sym, SYMBOL_SID);
   assert(symbol_is_valid(sym));

   return sym;
}

void symbol_exit(void)
{
   Symbol* q;
   Symbol* p;
   int     i;
   
   for(p = anchor; p != NULL; p = q)
   {
      assert(symbol_is_valid(p));

      SID_del(p);

      q = p->next;
      
      for(i = 0; i < p->used; i++)
         entry_free(p->entry[i]);

      free(p->entry);
      set_free(p->set);
      hash_free(p->hash);

      if (p->deflt != NULL)
         entry_free(p->deflt);

      free(p);
   }
   anchor = NULL;
}

Bool symbol_is_valid(const Symbol* sym)
{
   if (sym == NULL || !SID_ok(sym, SYMBOL_SID))
      return FALSE;

   mem_check(sym);
   mem_check(sym->entry);
      
   return TRUE;
}

Symbol* symbol_lookup(const char* name)
{
   Symbol* sym;

   assert(name != NULL);

   for(sym = anchor; sym != NULL; sym = sym->next)
      if (!strcmp(sym->name, name))
         break;

   return sym;
}

Bool symbol_has_entry(const Symbol* sym, const Tuple* tuple)
{
   assert(symbol_is_valid(sym));
   assert(tuple_is_valid(tuple));

   return hash_has_entry(sym->hash, tuple)
      || (sym->deflt != NULL && set_lookup(sym->set, tuple));
}

/* Liefert NULL wenn nicht gefunden.
 * Falls ein default zurueckgegeben wird, stimmt "tuple" nicht mit
 * entry->tuple ueberein.
 */
const Entry* symbol_lookup_entry(const Symbol* sym, const Tuple* tuple)
{
   const Entry* entry;
   
   assert(symbol_is_valid(sym));
   assert(tuple_is_valid(tuple));

   entry = hash_lookup_entry(sym->hash, tuple);

   if (NULL == entry && sym->deflt != NULL && set_lookup(sym->set, tuple))
      entry = sym->deflt;

   return entry;
}

/* Entry is eaten.
 * No check is done if entry->tuple is a member of sym->set !
 * This has to be done before.
 */
void symbol_add_entry(Symbol* sym, Entry* entry)
{
   const Tuple* tuple;
   
   assert(symbol_is_valid(sym));
   assert(entry_is_valid(entry));
   
   assert(sym->used <= sym->size);
   
   if (sym->used == sym->size)
   {
      sym->size   += sym->extend;
      sym->extend += sym->extend;
      sym->entry   = realloc(
         sym->entry, (size_t)sym->size * sizeof(*sym->entry));
      
      assert(sym->entry != NULL);
   }
   assert(sym->used < sym->size);

   tuple = entry_get_tuple(entry);

   /* There is no index set for the internal symbol.
    */
   assert(!strcmp(sym->name, SYMBOL_NAME_INTERNAL) || set_lookup(sym->set, tuple));

   if (hash_has_entry(sym->hash, tuple))
   {
      if (stmt_trigger_warning(166))
      {
         fprintf(stderr, "--- Warning 166: Duplicate element ");
         tuple_print(stderr, tuple);
         fprintf(stderr, " for symbol %s rejected\n", sym->name);
      }
      entry_free(entry);
   }
   else
   {
      /* Falls noch nicht geschehen, legen wir hier den Typ des
       * Symbols fest.
       */
      if ((sym->type == SYM_ERR) && (sym->used == 0))
         sym->type = entry_get_type(entry);

      assert(sym->type != SYM_ERR);
      
      hash_add_entry(sym->hash, entry);
      
      sym->entry[sym->used] = entry;      
      sym->used++;
   }
}

int symbol_get_dim(const Symbol* sym)
{
   assert(symbol_is_valid(sym));

   return set_get_dim(sym->set);
}

const Set* symbol_get_iset(const Symbol* sym)
{
   assert(symbol_is_valid(sym));

   return sym->set;
}

const char* symbol_get_name(const Symbol* sym)
{
   assert(symbol_is_valid(sym));

   return sym->name;
}

SymbolType symbol_get_type(const Symbol* sym)
{
   assert(symbol_is_valid(sym));

   return sym->type;
}

const Numb* symbol_get_numb(const Symbol* sym, int idx)
{
   assert(symbol_is_valid(sym));
   assert(idx >= 0);
   assert(idx <  sym->used);
   
   return entry_get_numb(sym->entry[idx]);
}

const char* symbol_get_strg(const Symbol* sym, int idx)
{
   assert(symbol_is_valid(sym));
   assert(idx >= 0);
   assert(idx <  sym->used);
   
   return entry_get_strg(sym->entry[idx]);
}

const Set* symbol_get_set(const Symbol* sym, int idx)
{
   assert(symbol_is_valid(sym));
   assert(idx >= 0);
   assert(idx <  sym->used);
   
   return entry_get_set(sym->entry[idx]);
}

Var* symbol_get_var(const Symbol* sym, int idx)
{
   assert(symbol_is_valid(sym));
   assert(idx >= 0);
   assert(idx <  sym->used);
   
   return entry_get_var(sym->entry[idx]);
}

void symbol_print(FILE* fp, const Symbol* sym)
{
   static const char* const type_name[] = { "Error", "Numb", "Strg", "Set", "Var" };
   
   int i;
   
   assert(symbol_is_valid(sym));

   fprintf(fp, "Name  : %s\n", sym->name);
   fprintf(fp, "Type  : %s\n", type_name[sym->type]);

   fprintf(fp, "Index : ");
   set_print(fp, sym->set);
   fprintf(fp, "\nEntries:\n");
   
   for(i = 0; i < sym->used; i++)
   {
      fprintf(fp, "\t%3d: ", i);
      entry_print(fp, sym->entry[i]);
      fprintf(fp, "\n");
   }
   fprintf(fp, "\n");
}

void symbol_print_all(FILE* fp)
{
   Symbol* sym;
   
   assert(fp != NULL);

   for(sym = anchor; sym != NULL; sym = sym->next)
      symbol_print(fp, sym);
}


