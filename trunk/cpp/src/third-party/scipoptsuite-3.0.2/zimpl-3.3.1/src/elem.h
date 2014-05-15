/* $Id: elem.h,v 1.7 2012/07/29 15:09:26 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: elem.h                                                        */
/*   Name....: Element Functions                                             */
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
#ifndef _ELEM_H_
#define _ELEM_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before elem.h"
#endif
#ifndef _NUMB_H_
#error "Need to include numb.h before elem.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum element_type
{
   ELEM_ERR = 0, ELEM_FREE, ELEM_NUMB, ELEM_STRG, ELEM_NAME
};

typedef enum element_type        ElemType;
typedef struct element           Elem;

#define ELEM_NULL  ((Elem*)0)

extern void         elem_init(void);
extern void         elem_exit(void);
/*lint -sem(        elem_new_numb, @p == 1) */
extern Elem*        elem_new_numb(const Numb* n);
/*lint -sem(        elem_new_strg, nulterm(1), 1p, @p == 1) */
extern Elem*        elem_new_strg(const char* s);
/*lint -sem(        elem_new_name, nulterm(1), 1p, @p == 1) */
extern Elem*        elem_new_name(const char* s);
/*lint -sem(        elem_free, custodial(1), 1p == 1) */
extern void         elem_free(Elem* elem);
/*lint -sem(        elem_is_valid, 1p == 1) */
extern Bool         elem_is_valid(const Elem* elem);
/*lint -sem(        elem_copy, 1p == 1, @p == 1) */
extern Elem*        elem_copy(const Elem* elem);
/*lint -sem(        elem_cmp, 1p == 1 && 2p == 1) */
extern Bool         elem_cmp(const Elem* elem_a, const Elem* elem_b);
/*lint -sem(        elem_get_type, 1p == 1) */
extern ElemType     elem_get_type(const Elem* elem);
/*lint -sem(        elem_get_numb, 1p == 1) */
extern const Numb*  elem_get_numb(const Elem* elem);
/*lint -sem(        elem_get_strg, 1p == 1, @p && nulterm(@)) */
extern const char*  elem_get_strg(const Elem* elem);
/*lint -sem(        elem_get_name, 1p == 1, @p && nulterm(@)) */
extern const char*  elem_get_name(const Elem* elem);
/*lint -sem(        elem_print, 1p == 1 && 2p == 1) */
extern void         elem_print(FILE* fp, const Elem* elem, Bool use_quotes);
/*lint -sem(        elem_hash, 1p == 1) */
extern unsigned int elem_hash(const Elem* elem);
/*lint -sem(        elem_tostr, 1p == 1, @p && nulterm(@)) */
extern char*        elem_tostr(const Elem* elem);

#ifdef __cplusplus
}
#endif
#endif /* _ELEM_H_ */
