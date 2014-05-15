/* $Id: tuple.h,v 1.7 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: tuple.h                                                       */
/*   Name....: Tuple Functions                                               */
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
#ifndef _TUPLE_H_
#define _TUPLE_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before tuple.h"
#endif
#ifndef _ELEM_H_
#error "Need to include elem.h before tuple.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tuple Tuple;

#define TUPLE_NULL ((Tuple*)0)

/*lint -sem(        tuple_new, 1n >= 0, @p == 1) */
extern Tuple*       tuple_new(int dim);
/*lint -sem(        tuple_free, custodial(1), 1p == 1) */
extern void         tuple_free(Tuple* tuple);
/*lint -sem(        tuple_is_valid, 1p == 1) */
extern Bool         tuple_is_valid(const Tuple* tuple);
/*lint -sem(        tuple_copy, 1p == 1, @p == 1) */
extern Tuple*       tuple_copy(const Tuple* tuple);
/*lint -sem(        tuple_cmp, 1p == 1 && 2p == 1) */
extern Bool         tuple_cmp(const Tuple* tuple_a, const Tuple* tuple_b);
/*lint -sem(        tuple_get_dim, 1p == 1, @n > 0) */
extern int          tuple_get_dim(const Tuple* tuple);
/*lint -sem(        tuple_set_elem, custodial(3), 1p == 1 && 2n >= 0 && 3p == 1) */
extern void         tuple_set_elem(Tuple* tuple, int idx, Elem* elem);
/*lint -sem(        tuple_get_elem, 1p == 1 && 2n >= 0, @p == 1) */
extern const Elem*  tuple_get_elem(const Tuple* tuple, int idx);
/*lint -sem(        tuple_combine, 1p == 1 && 2p == 1, @p == 1) */
extern Tuple*       tuple_combine(const Tuple* ta, const Tuple* tb);
/*lint -sem(        tuple_print, 1p == 1 && 2p == 1) */
extern void         tuple_print(FILE* fp, const Tuple* tuple);
/*lint -sem(        tuple_hash, 1p == 1) */
extern unsigned int tuple_hash(const Tuple* tuple);
/*lint -sem(        tuple_tostr, 1p == 1, @p && nulterm(@)) */
extern char*        tuple_tostr(const Tuple* tuple);

#ifdef __cplusplus
}
#endif
#endif /* _TUPLE_H_ */

