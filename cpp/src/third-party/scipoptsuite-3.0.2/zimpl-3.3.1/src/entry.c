/* $Id: entry.c,v 1.23 2012/07/29 15:09:26 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: entry.c                                                       */
/*   Name....: Symbol Table Entry Functions                                  */
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
#include "set.h"
#include "symbol.h"
#include "entry.h"
#include "zimpllib.h"

#define ENTRY_SID  0x456e7472

typedef union entry_value EntryValue;

union entry_value
{
   Numb*       numb;
   const char* strg;
   Set*        set;
   Var*        var;
};

struct entry
{
   SID
   int        refc;
   Tuple*     tuple;
   SymbolType type;
   EntryValue value;
};

Entry* entry_new_numb(const Tuple* tuple, const Numb* numb)
{
   Entry* entry = calloc(1, sizeof(*entry));

   assert(entry != NULL);
   assert(tuple != NULL);

   entry->refc       = 1;
   entry->tuple      = tuple_copy(tuple);
   entry->type       = SYM_NUMB;
   entry->value.numb = numb_copy(numb);

   SID_set(entry, ENTRY_SID);
   assert(entry_is_valid(entry));

   mem_check(entry);
   
   return entry;
}

Entry* entry_new_strg(const Tuple* tuple, const char* strg)
{
   Entry* entry = calloc(1, sizeof(*entry));

   assert(entry != NULL);
   assert(tuple != NULL);
   assert(strg  != NULL);
   
   entry->refc       = 1;
   entry->tuple      = tuple_copy(tuple);
   entry->type       = SYM_STRG;
   entry->value.strg = strg;

   SID_set(entry, ENTRY_SID);
   assert(entry_is_valid(entry));

   return entry;
}

Entry* entry_new_set(const Tuple* tuple, const Set* set)
{
   Entry* entry = calloc(1, sizeof(*entry));

   assert(entry != NULL);
   assert(tuple != NULL);
   assert(set   != NULL);
   
   entry->refc      = 1;
   entry->tuple     = tuple_copy(tuple);
   entry->type      = SYM_SET;
   entry->value.set = set_copy(set);

   SID_set(entry, ENTRY_SID);
   assert(entry_is_valid(entry));

   return entry;
}

Entry* entry_new_var(const Tuple* tuple, Var* var)
{
   Entry* entry = calloc(1, sizeof(*entry));

   assert(entry != NULL);
   assert(tuple != NULL);
   assert(var   != NULL);
   
   entry->refc      = 1;
   entry->tuple     = tuple_copy(tuple);
   entry->type      = SYM_VAR;
   entry->value.var = var;

   SID_set(entry, ENTRY_SID);
   assert(entry_is_valid(entry));

   return entry;
}

void entry_free(Entry* entry)
{
   assert(entry_is_valid(entry));

   entry->refc--;

   if (entry->refc == 0)
   {
      SID_del(entry);

      switch(entry->type)
      {
      case SYM_NUMB :
         numb_free(entry->value.numb);
         break;
      case SYM_STRG :
         break;
      case SYM_SET :
         set_free(entry->value.set);
         break;
      case SYM_VAR :
         break;
      default :
         abort();
      }
      tuple_free(entry->tuple);
   
      free(entry);
   }
}

Bool entry_is_valid(const Entry* entry)
{
   if (entry == NULL || !SID_ok(entry, ENTRY_SID))
      return FALSE;

   mem_check(entry);

   return TRUE;
}

Entry* entry_copy(const Entry* source)
{
   Entry* entry = (Entry*)source;
   
   assert(entry_is_valid(entry));

   entry->refc++;

   return entry;
}

Bool entry_cmp(const Entry* entry, const Tuple* tuple)
{
   assert(entry_is_valid(entry));
   assert(tuple_is_valid(tuple));
   
   return tuple_cmp(entry->tuple, tuple);
}

SymbolType entry_get_type(const Entry* entry)
{
   assert(entry_is_valid(entry));

   return entry->type;
}

const Tuple* entry_get_tuple(const Entry* entry)
{
   assert(entry_is_valid(entry));
   assert(tuple_is_valid(entry->tuple));
   
   return entry->tuple;
}

const Numb* entry_get_numb(const Entry* entry)
{
   assert(entry_is_valid(entry));
   assert(entry->type == SYM_NUMB);
   
   return entry->value.numb;
}

const char* entry_get_strg(const Entry* entry)
{
   assert(entry_is_valid(entry));
   assert(entry->type == SYM_STRG);
   
   return entry->value.strg;
}

const Set* entry_get_set(const Entry* entry)
{
   assert(entry_is_valid(entry));
   assert(entry->type == SYM_SET);
   
   return entry->value.set;
}

Var* entry_get_var(const Entry* entry)
{
   assert(entry_is_valid(entry));
   assert(entry->type == SYM_VAR);
   
   return entry->value.var;
}

void entry_print(FILE* fp, const Entry* entry)
{
   assert(entry_is_valid(entry));

   tuple_print(fp, entry->tuple);
   fprintf(fp, " -> ");
   
   switch(entry->type)
   {
   case SYM_NUMB :
      fprintf(fp, "%.16g", numb_todbl(entry->value.numb));
      break;
   case SYM_STRG :
      fprintf(fp, "\"%s\"", entry->value.strg);
      break;
   case SYM_SET :
      set_print(fp, entry->value.set);
      break;
   case SYM_VAR :
      zpl_var_print(fp, entry->value.var);
      break;
   default :
      fprintf(fp, "Entry-ERR");
      break;
   }
}

