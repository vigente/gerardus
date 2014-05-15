/* $Id: zimpllib.c,v 1.32 2012/07/29 15:09:31 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: zimpllib.c                                                    */
/*   Name....: ZIMPL Library Interface                                       */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Copyright (C) 2005-2012 by Thorsten Koch <koch@zib.de>
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

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <setjmp.h>

#include "lint.h"
#include "bool.h"
#include "mshell.h"
#include "stkchk.h"
#include "blkmem.h"
#include "random.h"
#include "ratlptypes.h"
#include "numb.h"
#include "elem.h"
#include "tuple.h"
#include "mme.h"
#include "set.h"
#include "symbol.h"
#include "define.h"
#include "bound.h"
#include "mono.h"
#include "term.h"
#include "stmt.h"
#include "local.h"
#include "list.h"
#include "entry.h"
#include "xlpglue.h"
#include "prog.h"
#include "metaio.h"
#include "strstore.h"
#include "zimpllib.h"

extern int yydebug;
extern int yy_flex_debug;

int verbose = VERB_QUIET;

static jmp_buf zpl_read_env;
static Bool    is_longjmp_ok = FALSE;

void zpl_print_banner(FILE* fp, Bool with_license)
{
   const char* const banner = 
      "****************************************************\n" \
      "* Zuse Institute Mathematical Programming Language *\n" \
      "* Release %-5s Copyright (C)2012 by Thorsten Koch *\n" \
      "****************************************************\n";

   const char* const license = 
      "*   This is free software and you are welcome to   *\n" \
      "*     redistribute it under certain conditions     *\n" \
      "*      ZIMPL comes with ABSOLUTELY NO WARRANTY     *\n" \
      "****************************************************\n";

   if (verbose >= VERB_NORMAL)
   {
      fprintf(fp, banner, VERSION);

      if (with_license || verbose > VERB_NORMAL)
         fprintf(fp, license); 

      fputc('\n', fp);
   }
}

void zpl_exit(int retval)
{
   if (is_longjmp_ok)
      longjmp(zpl_read_env, retval);

#if defined(NDEBUG) || defined(COVERAGE)
   exit(retval);
#else
   abort(); /* to get a stack trace */
#endif
}

static Bool is_valid_identifier(const char* s)
{
   assert(s != NULL);

   /* Identifiers start with a letter or a '_'
    */
   if (!isalpha(*s) || *s == '_')
      return FALSE;

   /* Then letters, digits or '_' can follow.
    */
   while(isalnum(*++s) || *s == '_')
      ;

   return *s == '\0';
}

void zpl_var_print(FILE* fp, const Var* var)
{
   const char* name  = xlp_getvarname(prog_get_lp(), var);
   VarClass    class = xlp_getclass(prog_get_lp(), var);
   Bound*      lower = xlp_getlower(prog_get_lp(), var);
   Bound*      upper = xlp_getupper(prog_get_lp(), var);

   fprintf(fp, "\"%s\" ", name);

   switch(class)
   {
   case VAR_CON :
      fprintf(fp, "real [");
      break;
   case VAR_IMP :
      fprintf(fp, "implicit integer [");
      break;
   case VAR_INT :
      fprintf(fp, "integer [");
      break;
   default :
      abort();
   }
   bound_print(fp, lower);
   fprintf(fp, ",");
   bound_print(fp, upper);
   fprintf(fp, "] ");
         
   bound_free(upper);
   bound_free(lower);
}

void zpl_add_parameter(const char* def)
{
   const char* warning =
      "--- Warning 175: Illegal syntax for command line define \"%s\" -- ignored\n";
   Set*    set;
   Symbol* sym;
   Numb*   numb;
   Tuple*  tuple;
   Entry*  entry;
   char*   name;
   char*   value;

   assert(def != NULL);
   
   name  = strdup(def);
   value = strchr(name, '=');
   
   if (value == NULL)
   {
      fprintf(stderr, warning, def);
      free(name);
      return;
   }
   *value = '\0';
   value++;

   if (strlen(name) == 0 || strlen(value) == 0 || !is_valid_identifier(name))
   {
      if (verbose > VERB_QUIET)
         fprintf(stderr, warning, def);
      free(name);
      return;
   }
   set   = set_pseudo_new();
   sym   = symbol_new(str_new(name), SYM_ERR, set, 1, ENTRY_NULL);
   tuple = tuple_new(0);   

   if (!numb_is_number(value))
      entry = entry_new_strg(tuple, str_new(value));
   else
   {
      numb  = numb_new_ascii(value);
      entry = entry_new_numb(tuple, numb);
      numb_free(numb);
   }
   symbol_add_entry(sym, entry);
   
   tuple_free(tuple);
   set_free(set); 
   free(name); 
}

Bool zpl_read(const char* filename, Bool with_management, void* user_data)
{
   Prog*       prog = NULL;
   Set*        set;
   void*       lp  = NULL;
   Bool        ret = FALSE;

   stkchk_init();
   
   yydebug       = 0;
   yy_flex_debug = 0;

   zpl_print_banner(stdout, FALSE);

   blk_init();
   str_init();
   rand_init(13021967UL);
   numb_init(with_management);
   elem_init();
   set_init();
   mio_init();
   interns_init();
   local_init();
   
   if (0 == setjmp(zpl_read_env))
   {
      is_longjmp_ok = TRUE;
      
      set = set_pseudo_new();
      (void)symbol_new(SYMBOL_NAME_INTERNAL, SYM_VAR, set, 100, NULL);
      set_free(set);
   
      prog = prog_new();

      prog_load(prog, NULL, filename);

      if (prog_is_empty(prog))
         fprintf(stderr, "*** Error 168: No program statements to execute\n");
      else
      {
         if (verbose >= VERB_DEBUG)
            prog_print(stderr, prog);
   
         lp = xlp_alloc(filename, FALSE, user_data);

         prog_execute(prog, lp);

         ret = TRUE;
      }
   }
   is_longjmp_ok = FALSE;

   if (lp != NULL)
      xlp_free(lp);

   if (prog != NULL)
      prog_free(prog);

   local_exit();
   interns_exit();
   mio_exit();
   symbol_exit();
   define_exit();
   set_exit();
   elem_exit();
   numb_exit();
   str_exit();
   blk_exit();
   
   return ret;
}

Bool zpl_read_with_args(char** argv, int argc, Bool with_management, void* user_data)
{
   const char* options = "D:mP:sv:";

   unsigned long seed = 13021967UL;
   char**        param_table;
   int           param_count = 0;
   int           c;
   int           i;
   Prog*         prog = NULL;
   Set*          set;
   void*         lp  = NULL;
   Bool          ret = FALSE;
   char*         inppipe = NULL;
   Bool          use_startval = FALSE;
   
   stkchk_init();

   yydebug       = 0;
   yy_flex_debug = 0;
   param_table   = malloc(sizeof(*param_table));

   zpl_print_banner(stdout, FALSE);

   /* getopt might be called more than once
    */
   optind = 1;
   
   while((c = getopt(argc, argv, options)) != -1)
   {
      switch(c)
      {
      case 'D' :
         param_table =
            realloc(param_table, ((unsigned int)param_count + 1) * sizeof(*param_table));
         param_table[param_count] = strdup(optarg);

         if (verbose >= VERB_DEBUG)
            printf("Parameter %d [%s]\n", param_count, param_table[param_count]);

         param_count++;
         break;
      case 'm' :
         use_startval = TRUE;
         break;
      case 'P' :
         inppipe = strdup(optarg);
         break;
      case 's' :
         seed = (unsigned long)atol(optarg);
         break;
      case 'v' :
         verbose = atoi(optarg);
         break;
      case '?':
         fprintf(stderr, "Unknown option '%c'\n", c);
         return FALSE;
      default :
         abort();
      }
   }
   if ((argc - optind) < 1)
   {
      fprintf(stderr, "Filename missing\n");
      free(param_table);

      return FALSE;
   }

   blk_init();
   str_init();
   rand_init(seed);
   numb_init(with_management);
   elem_init();
   set_init();
   mio_init();
   interns_init();
   local_init();
   
   if (0 == setjmp( zpl_read_env))
   {
      is_longjmp_ok = TRUE;
      
      /* Make symbol to hold entries of internal variables
       */
      set = set_pseudo_new();
      (void)symbol_new(SYMBOL_NAME_INTERNAL, SYM_VAR, set, 100, NULL);
      set_free(set);
   
      /* Now store the param defines
       */
      for(i = 0; i < param_count; i++)
         zpl_add_parameter(param_table[i]);

      prog = prog_new();

      for(i = optind; i < argc; i++)
         prog_load(prog, inppipe, argv[i]);

      if (prog_is_empty(prog))
         fprintf(stderr, "*** Error 168: No program statements to execute\n");
      else
      {
         if (verbose >= VERB_DEBUG)
            prog_print(stderr, prog);
   
         lp = xlp_alloc(argv[optind], use_startval, user_data);

         prog_execute(prog, lp);

         ret = TRUE;
      }
   }
   is_longjmp_ok = FALSE;

   if (lp != NULL)
      xlp_free(lp);

   /* Now clean up. 
    */
   if (inppipe != NULL)
      free(inppipe);
   
   for(i = 0; i < param_count; i++)
      free(param_table[i]);
   free(param_table);

   if (prog != NULL)
      prog_free(prog);

   local_exit();
   interns_exit();
   mio_exit();
   symbol_exit();
   define_exit();
   set_exit();
   elem_exit();
   numb_exit();
   str_exit();
   blk_exit();

   return ret;
}
