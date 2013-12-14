/* $Id: ratlpfwrite.c,v 1.27 2012/07/29 15:09:28 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: lpfwrite.c                                                    */
/*   Name....: LP Format File Writer                                         */
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <gmp.h>

#include "lint.h"
#include "mshell.h"
#include "bool.h"
#include "gmpmisc.h"
#include "ratlptypes.h"
#include "numb.h"
#include "bound.h"
#include "mme.h"
#include "mono.h"
#include "term.h"
#include "ratlp.h"
#include "ratlpstore.h"
#include "random.h"

static void permute(int size, void** tab)
{
   int i;

   if (size < 3)
      return;

   assert(size >= 3);
   
   for(i = 0; i < size; i++)
   {
      void* t;
      int   a = rand_get_range(0, size - 1); /*lint !e426 Call to function violates semantic (1n<2n)*/
      int   b = rand_get_range(0, size - 1); /*lint !e426 Call to function violates semantic (1n<2n)*/

      assert(a >= 0);
      assert(a <  size);
      assert(b >= 0);
      assert(b <  size);
      
      t      = tab[a];
      tab[a] = tab[b];
      tab[b] = t;
   }
}

static void write_val(FILE* fp, LpFormat format, Bool force_sign, const mpq_t val)
{
   switch(format)
   {
   case LP_FORM_LPF :
   case LP_FORM_RLP :
   case LP_FORM_PIP :
      fprintf(fp, force_sign ? "%+.15g" : "%.15g", mpq_get_d(val));
      break;
   case LP_FORM_HUM :
      if (force_sign && (mpq_sgn(val) > 0)) /*lint !e634 Strong type mismatch (type 'Bool') */
         fprintf(fp, "+");

      mpq_out_str(fp, 10, val);
      break;
   default:
      abort();
   }
}

static void write_lhs(FILE* fp, LpFormat format, const Con* con, ConType type)
{
   assert(fp  != NULL);
   assert(con != NULL);
   
   switch(type)
   {
   case CON_RHS :
   case CON_LHS :
   case CON_EQUAL :
      break;
   case CON_RANGE :
      write_val(fp, format, FALSE, con->lhs);
      fprintf(fp, " <= ");
      break;
   default :
      abort();
   }
}

static void write_rhs(FILE* fp, LpFormat format, const Con* con, ConType type)
{
   assert(fp  != NULL);
   assert(con != NULL);
   
   switch(type)
   {
   case CON_RHS :
   case CON_RANGE :
      fprintf(fp, " <= ");
      write_val(fp, format, FALSE, con->rhs);
      break;
   case CON_LHS :
      fprintf(fp, " >= ");
      write_val(fp, format, FALSE, con->lhs);
      break;
   case CON_EQUAL :
      fprintf(fp, " = ");
      write_val(fp, format, FALSE, con->rhs);
      break;
   default :
      abort();
   }
   fprintf(fp, "\n");
}

static void write_row(
   FILE*      fp,
   LpFormat   format,
   const Con* con,
   char*      name,
   int        name_size)
{
   Nzo*  nzo;
   Nzo** nzotab;
   int   cnt;
   int   i;

   assert(fp        != NULL);
   assert(con       != NULL);
   assert(name      != NULL);

   /* Add 1 in case con->size == 0
    */
   nzotab = calloc((size_t)con->size + 1, sizeof(*con));

   assert(nzotab != NULL);

   for(cnt = 0, nzo = con->first; nzo != NULL; nzo = nzo->con_next)
      nzotab[cnt++] = nzo;

   assert(cnt == con->size);

   if (format == LP_FORM_RLP)
      permute(con->size, (void**)nzotab);

   cnt = 0;

   for(i = 0; i < con->size; i++)
   {
      nzo = nzotab[i];

      lps_makename(name, name_size, nzo->var->name, format == LP_FORM_HUM ? -1 : nzo->var->number);

      if (mpq_equal(nzo->value, const_one))
         fprintf(fp, " + %s", name);
      else if (mpq_equal(nzo->value, const_minus_one))
         fprintf(fp, " - %s", name);
      else
      {
         fprintf(fp, " ");
         write_val(fp, format, TRUE, nzo->value);      
         fprintf(fp, " %s", name);
      }
      if (++cnt % 6 == 0)
         fprintf(fp, "\n ");         
   }
   if (con->qme_first != NULL)
   {
      Qme* qme;

      if (cnt % 6 != 0)
         fprintf(fp, "\n ");         

      cnt = 0;

      if (format == LP_FORM_LPF || format == LP_FORM_RLP)
         fprintf(fp, " + [");

      for(qme = con->qme_first; qme != NULL; qme= qme->next)
      {
         lps_makename(name, name_size, qme->var1->name, format == LP_FORM_HUM ? -1 : qme->var1->number);

         if (mpq_equal(qme->value, const_one))
            fprintf(fp, " + %s", name);
         else if (mpq_equal(qme->value, const_minus_one))
            fprintf(fp, " - %s", name);
         else
         {
            fprintf(fp, " ");         
            write_val(fp, format, TRUE, qme->value);      
            fprintf(fp, " %s", name);
         }

         if (qme->var1 == qme->var2)
            fprintf(fp, "^2");
         else
         {
            lps_makename(name, name_size, qme->var2->name, format == LP_FORM_HUM ? -1 : qme->var2->number);

            fprintf(fp, " * %s", name);
         }
         if (++cnt % 6 == 0)
            fprintf(fp, "\n ");         
      }
      if (format == LP_FORM_LPF || format == LP_FORM_RLP)
      {
         fprintf(fp,  " ]\n");
         cnt = 0;
      }
   }
   if (con->term != NULL)
   {
      const Term* term  = con->term;
      Bool only_comment = FALSE;

      assert(term_get_degree(term) > 2 || !term_is_polynomial(term));

      if (format != LP_FORM_PIP)
      {
         if (verbose > 0)
         {
            fprintf(stderr, "--- Warning 600: File format can only handle linear and quadratic constraints\n");
            fprintf(stderr, "                 Constraint %s with degree %d ignored\n", 
               con->name, term_get_degree(term));
         }
         only_comment = TRUE;
      }
      assert(numb_equal(term_get_constant(term), numb_zero()));
   
      if (cnt % 6 != 0)
         fprintf(fp, "\n ");         
   
      cnt = 0;
   
      if (only_comment)
         fprintf(fp, "\\ ");

      for(i = 0; i < term_get_elements(term); i++)
      {
         const Mono* mono  = term_get_element(term, i);
         const Numb* coeff = mono_get_coeff(mono);
         MFun        fun   = mono_get_function(mono);
         int         k;

         if (fun == MFUN_NONE)
         {
            if (numb_equal(coeff, numb_one()))
               fprintf(fp, " +");
            else
            {
               mpq_t t;
               mpq_init(t);
               numb_get_mpq(coeff, t);
               fprintf(fp, " ");         
               write_val(fp, format, TRUE, t);      
               mpq_clear(t);
            }
            fputc(' ', fp);
         }
         else
         {
            switch(fun)
            {
            case MFUN_SQRT :
               fprintf(fp, " + sqrt(");
               break;
            case MFUN_LOG :
               fprintf(fp, " + log(");
               break;
            case MFUN_EXP :
               fprintf(fp, " + exp(");
               break;
            case MFUN_LN :
               fprintf(fp, " + ln(");
               break;
            case MFUN_SIN :
               fprintf(fp, " + sin(");
               break;
            case MFUN_COS :
               fprintf(fp, " + cos(");
               break;
            case MFUN_TAN :
               fprintf(fp, " + tan(");
               break;
            case MFUN_ABS :
               fprintf(fp, " + abs(");
               break;
            case MFUN_SGN :
               fprintf(fp, " + sgn(");
               break;
            case MFUN_POW :
               fprintf(fp, " + pow(");
               break;
            case MFUN_SGNPOW :
               fprintf(fp, " + sgnpow(");
               break;
            case MFUN_TRUE :
            case MFUN_FALSE :
            default :
               abort();
            } 
         }
          
         for(k = 0; k < mono_get_degree(mono); k++)
         {
            Var* var = mono_get_var(mono, k);
            int  j;
            
            if (k > 0)
               fprintf(fp, " * ");
            
            for(j = 1; k + j < mono_get_degree(mono); j++)
               if (var != mono_get_var(mono, k + j))
                  break;
               
            lps_makename(name, name_size, var->name, format == LP_FORM_HUM ? -1 : var->number);
               
            if (j == 1)
               fprintf(fp, "%s", name);
            else
            {
               fprintf(fp, "%s^%d", name, j);
               k += j - 1; /*lint !e850 loop index variable is modified in body of the loop */
            }
         }
         if (fun != MFUN_NONE)
         {
            if (fun == MFUN_POW || fun == MFUN_SGNPOW)
            {
               mpq_t t;
               mpq_init(t);
               numb_get_mpq(coeff, t);
               fprintf(fp, " ,");         
               write_val(fp, format, FALSE, t);      
               mpq_clear(t);
            }
            fprintf(fp, ") ");
         }
         if (++cnt % 6 == 0)
            fprintf(fp, "\n%s ", only_comment ? "\\" : "");         
      }
   }
   free(nzotab);
}

/* A specification for the LP file format can be found in the
 * ILOG CPLEX 7.0 Reference Manual page 527.
 * ILOG CPLEX 8.0 Reference Manual page 595.
 * The "Lazy Constraints" section seems to be undocumented.
 */
void lpf_write(
   const Lps*  lp,
   FILE*       fp,
   LpFormat    format,
   const char* text)
{
   const Var* var;
   Con*       con;
   Con**      contab;
   Bool  have_integer   = FALSE;
   Bool  have_separate  = FALSE;
   Bool  have_checkonly = FALSE;
   int   cnt;
   int   i;
   int   k;
   int   name_size;
   char* name;

   assert(lp != NULL);
   assert(fp != NULL);

   /* Store constraint pointers and permute them
    * (add 1 in case lp->cons == 0)
    */
   contab = calloc((size_t)lp->cons + 1, sizeof(*contab));

   assert(contab != NULL);

   k = 0;
   for(con = lp->con_root; con != NULL; con = con->next)
      contab[k++] = con;

   assert(k == lp->cons);

   if (format == LP_FORM_RLP)
      permute(lp->cons, (void**)contab);

   name_size = lps_getnamesize(lp, LP_FORM_LPF);
   name      = malloc((size_t)name_size);

   assert(name != NULL);
   
   if (text != NULL)
      fprintf(fp, "%s", text);   
      
   fprintf(fp, "\\Problem name: %s\n", lp->name);   
   fprintf(fp, "%s\n", (lp->direct == LP_MIN) ? "Minimize" : "Maximize");
   fprintf(fp, " %s: ", lp->objname == NULL ? "Objective" : lp->objname);
   
   for(var = lp->var_root, cnt = 0; var != NULL; var = var->next)
   {
      /* If cost is zero, do not include in objective function
       */
      if (mpq_equal(var->cost, const_zero))
         continue;

      lps_makename(name, name_size, var->name, format == LP_FORM_HUM ? -1 : var->number);
      
      if (mpq_equal(var->cost, const_one))
         fprintf(fp, " + %s", name);
      else if (mpq_equal(var->cost, const_minus_one))
         fprintf(fp, " - %s", name);
      else
      {
         fprintf(fp, " ");         
         write_val(fp, format, TRUE, var->cost);      
         fprintf(fp, " %s", name);
      }

      if (++cnt % 6 == 0)
         fprintf(fp, "\n ");
   }
   /* ---------------------------------------------------------------------- */

   /* First loop run for normal constraints, second one for
    * user cuts, thrid one for lazy constraints, if any.
    */
   for(i = 0; i < 3; i++)
   {      
      if (i == 0)
         fprintf(fp, "\nSubject to\n");      
      else if (i == 1)
      {
         if (!have_separate)
            continue;
         else
            fprintf(fp, "\nUser Cuts\n");
      }
      else if (i == 2)
      {
         if (!have_checkonly)
            continue;
         else
            fprintf(fp, "\nLazy Constraints\n"); 
      }

      for(k = 0; k < lp->cons; k++)
      {
         con = contab[k];

         if (con->size == 0 && con->qme_first == NULL && con->term == NULL)
            continue;

         if (i == 0 && ((con->flags & (LP_FLAG_CON_SEPAR | LP_FLAG_CON_CHECK)) != 0))
         {
            if ((con->flags & LP_FLAG_CON_SEPAR) == LP_FLAG_CON_SEPAR)
               have_separate = TRUE;
            if ((con->flags & LP_FLAG_CON_CHECK) == LP_FLAG_CON_CHECK)
               have_checkonly = TRUE;

            continue;
         }        
         if (i == 1 && (con->flags & LP_FLAG_CON_SEPAR) != LP_FLAG_CON_SEPAR)
            continue;

         if (i == 2 && (con->flags & LP_FLAG_CON_CHECK) != LP_FLAG_CON_CHECK)
            continue;
    
         if (con->type == CON_RANGE)
         {
            if (format == LP_FORM_HUM)
            {
               lps_makename(name, name_size, con->name, -1);
               fprintf(fp, " %s:\n ", name);

               write_lhs(fp, format, con, CON_RANGE);
               write_row(fp, format, con, name, name_size); 
               write_rhs(fp, format, con, CON_RANGE);
            }
            else
            {
               /* Split ranges, because LP format can't handle them.
                */
               lps_makename(name, name_size, con->name, con->number);
               fprintf(fp, " %sR:\n ", name);
   
               write_row(fp, format, con, name, name_size); /* changes name */
               write_rhs(fp, format, con, CON_RHS);
   
               lps_makename(name, name_size, con->name, con->number);
               fprintf(fp, " %sL:\n ", name);
   
               write_row(fp, format, con, name, name_size); /* changes name */
               write_rhs(fp, format, con, CON_LHS);
            }
         }
         else
         {
            lps_makename(name, name_size, con->name, format == LP_FORM_HUM ? -1 : con->number);

            if (i == 0)
               fprintf(fp, " %s:\n ", name);
            else if (i == 1)
               fprintf(fp, " %sU:\n ", name);
            else if (i == 2)
               fprintf(fp, " %sZ:\n ", name);

            if (con->ind_var != NULL)
            {
               lps_makename(name, name_size, con->ind_var->name, format == LP_FORM_HUM ? -1 : con->ind_var->number);
               fprintf(fp, "%s = %d -> ", name, con->ind_dir ? 1 : 0);
            }
            write_row(fp, format, con, name, name_size);
            write_rhs(fp, format, con, con->type);
         }
      }
   }

   /* ---------------------------------------------------------------------- */

   fprintf(fp, "Bounds\n");

   for(var = lp->var_root; var != NULL; var = var->next)
   {
      /* A variable without any entries in the matrix
       * or the objective function can be ignored.
       * (also not part of an SOS or quadratic constraint)
       */
      if (var->size == 0 && mpq_equal(var->cost, const_zero) && !var->is_used)
         continue;

      lps_makename(name, name_size, var->name, format == LP_FORM_HUM ? -1 : var->number);

      if (var->type == VAR_FIXED)
      {
         fprintf(fp, " %s = ", name);
         write_val(fp, format, FALSE, var->lower);      
         fprintf(fp, "\n");
      }
      else
      {
         /* Check if we have integer variables
          */
         if (var->vclass == VAR_INT)
            have_integer = TRUE;
         
         fprintf(fp, " ");

         if (var->type == VAR_LOWER || var->type == VAR_BOXED)
            write_val(fp, format, FALSE, var->lower);      
         else
            fprintf(fp, "-inf");
         
         fprintf(fp, " <= %s <= ", name);
         
         if (var->type == VAR_UPPER || var->type == VAR_BOXED)
         {
            write_val(fp, format, FALSE, var->upper);      
            fprintf(fp, "\n");
         }
         else
            fprintf(fp, "+inf\n");
      }
   }

   /* ---------------------------------------------------------------------- */

   if (have_integer)
   {
      fprintf(fp, "General\n");
      
      for(var = lp->var_root; var != NULL; var = var->next)
      {
         if (var->vclass != VAR_INT)
            continue;

         if (var->size == 0 && mpq_equal(var->cost, const_zero) && !var->is_used)
            continue;
         
         lps_makename(name, name_size, var->name, format == LP_FORM_HUM ? -1 : var->number);

         fprintf(fp, " %s\n", name);
      }
   }

   /* ---------------------------------------------------------------------- */

   if (lps_has_sos(lp))
   {
      const Sos* sos;
      const Sse* sse;

      fprintf(fp, "SOS\n");

      for(sos = lp->sos_root; sos != NULL; sos = sos->next)
      {
         cnt = 0;

         fprintf(fp, " %s:S%d:: ", 
            sos->name,
            sos->type == SOS_TYPE1 ? 1 : 2);
    
         for(sse = sos->first; sse != NULL; sse = sse->next)
         {
            lps_makename(name, name_size, sse->var->name, format == LP_FORM_HUM ? -1 : sse->var->number);

            fprintf(fp, " %s:", name);
            write_val(fp, format, FALSE, sse->weight);      

            if (++cnt % 6 == 0)
               fputc('\n', fp);         
         }
         if (cnt % 6 != 0)
            fputc('\n', fp);         
      }
   }
   fprintf(fp, "End\n");

   free(name);
   free(contab);
}   

/* ------------------------------------------------------------------------- */
/* Emacs Local Variables:                                                    */
/* Emacs mode:c                                                              */
/* Emacs c-basic-offset:3                                                    */
/* Emacs tab-width:8                                                         */
/* Emacs indent-tabs-mode:nil                                                */
/* Emacs End:                                                                */
/* ------------------------------------------------------------------------- */
