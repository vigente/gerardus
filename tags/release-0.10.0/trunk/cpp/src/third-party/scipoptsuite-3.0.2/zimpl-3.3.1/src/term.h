/* $Id: term.h,v 1.8 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: term.h                                                        */
/*   Name....: Term Functions                                                */
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
#ifndef _TERM_H_
#define _TERM_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before term.h"
#endif

#ifndef _MONO_H_
#error "Need to include mono.h before term.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct term              Term;

/* term.c
 */
/*lint -sem(        term_new, 1n > 0, @p == 1) */
extern Term*        term_new(int size);
/*lint -sem(        term_add_elem, 1p == 1 && 2p == 1 && 3p == 1) */
extern void         term_add_elem(Term* term, const Entry* entry, const Numb* coeff, MFun mfun);
#if 0 /* ??? not used */
/*lint -sem(        term_mul_elem, 1p == 1 && 2p == 1 && 3p == 1) */
extern void         term_mul_elem(Term* term, const Entry* entry, const Numb* coeff);
#endif
/*lint -sem(        term_free, custodial(1), 1p == 1) */
extern void         term_free(Term* term);
/*lint -sem(        term_is_valid, 1p == 1) */
extern Bool         term_is_valid(const Term* term);
/*lint -sem(        term_copy, 1p == 1, @p == 1) */
extern Term*        term_copy(const Term* term);
/*lint -sem(        term_print, 1p == 1 && 2p == 1) */
extern void         term_print(FILE* fp, const Term* term, Bool print_symbol_index);
/*lint -sem(        term_append_term, 1p == 1 && 2p == 1) */
extern void         term_append_term(Term* term_a, const Term* term_b);
/*lint -sem(        term_add_term, 1p == 1 && 2p == 1, @p == 1) */
extern Term*        term_add_term(const Term* term_a, const Term* term_b);
/*lint -sem(        term_sub_term, 1p == 1 && 2p == 1, @p == 1) */
extern Term*        term_sub_term(const Term* term_a, const Term* term_b);
/*lint -sem(        term_mul_term, 1p == 1 && 2p == 1, @p == 1) */
extern Term*        term_mul_term(const Term* term_a, const Term* term_b);
/*lint -sem(        term_simplify, 1p == 1, @p == 1) */
extern Term*        term_simplify(const Term* term_org);
/*lint -sem(        term_add_constant, 1p == 1 && 2p == 1) */
extern void         term_add_constant(Term* term, const Numb* value);
/*lint -sem(        term_sub_constant, 1p == 1 && 2p == 1) */
extern void         term_sub_constant(Term* term, const Numb* value);
/*lint -sem(        term_mul_coeff, 1p == 1 && 2p == 1) */
extern void         term_mul_coeff(Term* term, const Numb* value);
/*lint -sem(        term_get_constant, 1p == 1) */
extern const Numb*  term_get_constant(const Term* term);
#if 0 /* ??? not used */
/*lint -sem(        term_negate, 1p == 1) */
extern void         term_negate(Term* term);
#endif
/*lint -sem(        term_to_objective, 1p == 1) */
extern void         term_to_objective(const Term* term);
/*lint -sem(        term_get_elements, 1p == 1, @n >= 0) */
extern int          term_get_elements(const Term* term);
/*lint -sem(        term_get_element, 1p == 1, @p == 1) */
extern Mono*        term_get_element(const Term* term, int i);
/*lint -sem(        term_get_lower_bound, 1p == 1, @p == 1) */
extern Bound*       term_get_lower_bound(const Term* term);
/*lint -sem(        term_get_upper_bound, 1p == 1, @p == 1) */
extern Bound*       term_get_upper_bound(const Term* term);
/*lint -sem(        term_is_all_integer, 1p == 1) */
extern Bool         term_is_all_integer(const Term* term);
/*lint -sem(        term_is_linear, 1p == 1) */
extern Bool         term_is_linear(const Term* term);
/*lint -sem(        term_is_polynomial, 1p == 1) */
extern Bool         term_is_polynomial(const Term* term);
/*lint -sem(        term_get_degree, 1p == 1) */
extern int          term_get_degree(const Term* term);
/*lint -sem(        term_make_conditional, 1p == 1 && 2p == 1, @p == 1) */
extern Term*        term_make_conditional(const Term* ind_term, const Term* cond_term, Bool is_true);

#ifdef __cplusplus
}
#endif
#endif /* _TERM_H_ */
