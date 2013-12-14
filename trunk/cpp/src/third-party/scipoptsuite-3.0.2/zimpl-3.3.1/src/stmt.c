/* $Id: stmt.c,v 1.29 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: stmt.c                                                        */
/*   Name....: Statement Functions                                           */
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
#include "stkchk.h"
#include "ratlptypes.h"
#include "numb.h"
#include "elem.h"
#include "tuple.h"
#include "mme.h"
#include "set.h"
#include "symbol.h"
#include "entry.h"
#include "idxset.h"
#include "rdefpar.h"
#include "bound.h"
#include "define.h"
#include "mono.h"
#include "term.h"
#include "list.h"
#include "local.h"
#include "code.h"
#include "inst.h"
#include "stmt.h"

#define STMT_SID 0x53746d74


struct statement
{
   SID
   StmtType    type;
   const char* filename;
   int         lineno;
   const char* text;
   CodeNode*   node;
};

#define MAX_WARNINGS 1000

static int warning_count[MAX_WARNINGS];

static void activate_warnings(void)
{
   int i;

   for(i = 0; i < MAX_WARNINGS; i++)
      warning_count[i] = 0;
}

static void show_suppressed_warnings(void)
{
   int i;

   if (verbose > VERB_QUIET && verbose < VERB_CHATTER)
      for(i = 0; i < MAX_WARNINGS; i++)
         if (warning_count[i] > 1)
            fprintf(stderr, "--- Warning %3d: suppressed %d further message(s)\n",
               i, warning_count[i] - 1);
}

Bool stmt_trigger_warning(int no)
{
   Bool ret;
   
   assert(no >= 0);
   assert(no < MAX_WARNINGS);

   ret = (warning_count[no] == 0);

   warning_count[no]++;

   if (verbose >= VERB_CHATTER)
      ret = TRUE;

   if (verbose <= VERB_QUIET)
      ret = FALSE;
   
   return ret;
}
   
Stmt* stmt_new(
   StmtType    type,
   const char* filename,
   int         lineno,
   const char* text)
{
   Stmt* stmt = calloc(1, sizeof(*stmt));;

   assert(filename != NULL);
   assert(text     != NULL);
   assert(stmt     != NULL);
   assert(lineno   > 0);
   
   stmt->type     = type;
   stmt->filename = strdup(filename);
   stmt->lineno   = lineno;
   stmt->text     = strdup(text);
   stmt->node     = NULL;
   
   SID_set(stmt, STMT_SID);
   assert(stmt_is_valid(stmt));

   return stmt;
}

void stmt_free(Stmt* stmt)
{
   assert(stmt_is_valid(stmt));

   SID_del(stmt);
   
   if (stmt->node != NULL)
      code_free(stmt->node);

   free((void*)stmt->filename);
   free((void*)stmt->text);
   free(stmt);
}

Bool stmt_is_valid(const Stmt* stmt)
{
   return ((stmt != NULL)
      && SID_ok(stmt, STMT_SID)
      && (stmt->filename != NULL)
      && (stmt->lineno   >  0)
      && (stmt->text     != NULL));
}

const char* stmt_get_filename(const Stmt* stmt)
{
   assert(stmt_is_valid(stmt));

   return stmt->filename;
}

int stmt_get_lineno(const Stmt* stmt)
{
   assert(stmt_is_valid(stmt));

   return stmt->lineno;
}

const char* stmt_get_text(const Stmt* stmt)
{
   assert(stmt_is_valid(stmt));

   return stmt->text;
}

void stmt_parse(Stmt* stmt)
{
   Trace("stmt_parse");

   assert(stmt_is_valid(stmt));
   
   if (verbose >= VERB_VERBOSE)
      printf("Parsing %s %d\n", stmt->filename, stmt->lineno);

   parse_stmt(stmt);

   stmt->node = code_get_root();
}

void stmt_execute(const Stmt* stmt)
{
   int inst_count = code_get_inst_count();
   
   Trace("stmt_execute");
   
   assert(stmt_is_valid(stmt));

   if (verbose >= VERB_VERBOSE)
      printf("Executing %s %d\n", stmt->filename, stmt->lineno);

   activate_warnings();

   (void)code_prune_tree(stmt->node);

   /* ??? I don't think this can happen without a parse error first.
    */
   if (code_get_type(code_eval(stmt->node)) != CODE_VOID)
   {  /*             ^^^^^^^^^^^^^^^^^^^^ */
      fprintf(stderr, "*** Error 169: Execute must return void element\n");
      zpl_exit(EXIT_FAILURE);
   }
   show_suppressed_warnings();

   if (verbose >= VERB_CHATTER)
   {
      printf("Instructions evaluated: %u\n", code_get_inst_count() - inst_count);
      stkchk_maximum(stdout);
   }
}

void stmt_print(FILE* fp, const Stmt* stmt)
{
   static const char* const type_name[] =
   {
      "Unknown", "Set", "Param", "Var", "Min", "Max", "Cons", "Define", "Print", "SOS"
   };
   assert(stmt_is_valid(stmt));

   /* Lint weiss hier dass das assert immer erfuellt sein muss.
    * aber wir wollen es trotzdem.
    */
   assert((unsigned int)stmt->type
      < (sizeof(type_name) / sizeof(type_name[0]))); /*lint !e650 */

   fprintf(fp, "%s %04d %-7s [%s]\n",
      stmt->filename,
      stmt->lineno,
      type_name[(int)stmt->type],
      stmt->text);
}







