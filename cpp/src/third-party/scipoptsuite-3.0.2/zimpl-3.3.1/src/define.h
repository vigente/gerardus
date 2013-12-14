/* $Id: define.h,v 1.6 2012/07/29 15:09:26 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: define.h                                                      */
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
#ifndef _DEFINE_H_
#define _DEFINE_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before define.h"
#endif
#ifndef _NUMB_H_
#error "Need to include numb.h before define.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum define_type     { DEF_ERR = 0, DEF_NUMB, DEF_STRG, DEF_BOOL, DEF_SET };

typedef enum define_type         DefineType;
typedef struct define            Define;

/*lint -sem(        define_new, nulterm(1), 1p, @p == 1) */
extern Define*      define_new(const char* name, DefineType type);
/*lint -sem(        define_set_param, custodial(2), 1p == 1 && 2p == 1) */
extern void         define_set_param(Define* def, Tuple* param);
/*lint -sem(        define_set_code, 1p == 1 && 2p == 1) */
extern void         define_set_code(Define* def, CodeNode* code);
extern void         define_exit(void);
/*lint -sem(        define_is_valid, 1p == 1) */
extern Bool         define_is_valid(const Define* def);
/*lint -sem(        define_lookup, nulterm(1), 1p, r_null) */
extern Define*      define_lookup(const char* name);
/*lint -sem(        define_get_name, 1p == 1, @p && nulterm(@)) */
extern const char*  define_get_name(const Define* def);
/*lint -sem(        define_get_type, 1p == 1) */
extern DefineType   define_get_type(const Define* def);
/*lint -sem(        define_get_param, 1p == 1) */
extern const Tuple* define_get_param(const Define* def);
/*lint -sem(        define_get_code, 1p == 1) */
extern CodeNode*    define_get_code(const Define* def);

#ifdef __cplusplus
}
#endif
#endif /* _DEFINE_H_ */
