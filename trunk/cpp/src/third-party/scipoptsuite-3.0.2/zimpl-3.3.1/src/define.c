/* $Id: define.c,v 1.12 2012/07/29 15:09:26 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: define.c                                                      */
/*   Name....: Define Table Functions                                        */
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
#include "numb.h"
#include "elem.h"
#include "tuple.h"
#include "mme.h"
#include "define.h"

#define DEFINE_SID  0x44656669

struct define
{
   SID
   const char*  name;
   DefineType   type;
   Tuple*       param;
   CodeNode*    code;
   Define*      next;
};

#ifndef NDEBUG
static Define anchor  = { 0, "", DEF_ERR, NULL, NULL, NULL };
#else
static Define anchor  = {    "", DEF_ERR, NULL, NULL, NULL };
#endif

Define* define_new(
   const char*  name,
   DefineType   type)
{
   Define* def;

   assert(name           != NULL);
   assert(strlen(name)   >  0);
   assert(type           != DEF_ERR);
   
   def = calloc(1, sizeof(*def));

   assert(def != NULL);
   
   def->name    = name;
   def->type    = type;
   def->param   = NULL;
   def->code    = NULL;
   def->next    = anchor.next;
   anchor.next  = def;

   SID_set(def, DEFINE_SID);
   assert(define_is_valid(def));

   return def;
}

void define_set_param(
   Define*     def,
   Tuple*      param)
{
   assert(define_is_valid(def));
   assert(tuple_is_valid(param));
   
   def->param   = param;
}

void define_set_code(
   Define*     def,
   CodeNode*   code)
{
   assert(define_is_valid(def));
   assert(code != NULL);
   
   def->code = code;
}

void define_exit(void)
{
   Define* q;
   Define* p;
   
   for(p = anchor.next; p != NULL; p = q)
   {
      assert(define_is_valid(p));

      SID_del(p);

      tuple_free(p->param);
      
      q = p->next;
      
      free(p);
   }
   anchor.next = NULL;
}

Bool define_is_valid(const Define* def)
{
   if (def == NULL || !SID_ok(def, DEFINE_SID))
      return FALSE;

   mem_check(def);

   return TRUE;
}

Define* define_lookup(const char* name)
{
   Define* def;

   assert(name != NULL);

   for(def = anchor.next; def != NULL; def = def->next)
      if (!strcmp(def->name, name))
         break;

   return def;
}

const char* define_get_name(const Define* def)
{
   assert(define_is_valid(def));

   return def->name;
}

DefineType define_get_type(const Define* def)
{
   assert(define_is_valid(def));

   return def->type;
}

const Tuple* define_get_param(const Define* def)
{
   assert(define_is_valid(def));

   return def->param;
}

CodeNode* define_get_code(const Define* def)
{
   assert(define_is_valid(def));

   return def->code;
}








