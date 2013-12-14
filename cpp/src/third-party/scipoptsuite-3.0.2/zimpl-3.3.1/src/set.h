/* $Id: set.h,v 1.8 2012/07/29 15:09:29 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: set4.c                                                        */
/*   Name....: Set Functions                                                 */
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
#ifndef _SET_H_
#define _SET_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before set.h"
#endif

#ifndef _MME_H_
#error "Need to include mme.h before set.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum set_check_type
{
   SET_CHECK_NONE, SET_CHECK_QUIET, SET_CHECK_WARN
};

typedef enum set_check_type      SetCheckType;

extern void         set_init(void);
extern void         set_exit(void);
/*lint -sem(        set_free, custodial(1), 1p == 1) */
extern void         set_free(Set* set);
/*lint -sem(        set_is_valid, 1p == 1) */
extern Bool         set_is_valid(const Set* set);
/*lint -sem(        set_copy, 1p == 1, @p == 1) */
extern Set*         set_copy(const Set* set);
/*lint -sem(        set_lookup, 1p == 1 && 2p == 1) */
extern Bool         set_lookup(const Set* set, const Tuple* tuple);
/*lint -sem(        set_iter_init, 1p == 1, @p == 1) */
extern SetIter*     set_iter_init(const Set* set, const Tuple* pattern);
/*lint -sem(        set_iter_next, 1p == 1 && 2p == 1) */
extern Tuple*       set_iter_next(SetIter* iter, const Set* set);
/*lint -sem(        set_iter_exit, 1p == 1 && 2p == 1) */
extern void         set_iter_exit(SetIter* iter, const Set* set);
/*lint -sem(        set_get_dim, 1p == 1, @n >= 0) */
extern int          set_get_dim(const Set* set);
/*lint -sem(        set_get_members, 1p == 1, @n >= 0) */
extern int          set_get_members(const Set* set);
/*lint -sem(        set_get_tuple, 1p == 1 && 2n >= 0) */
extern Tuple*       set_get_tuple(const Set* set, int idx);
/*lint -sem(        set_print, 1p == 1 && 2p == 1) */
extern void         set_print(FILE* fp, const Set* set);

/*lint -sem(        set_empty_new, 1n >= 0, @p == 1) */
extern Set*         set_empty_new(int dim);
/*lint -sem(        set_pseudo_new, @p == 1) */
extern Set*         set_pseudo_new(void);
/*lint -sem(        set_new_from_list, 1p == 1, @p == 1) */
extern Set*         set_new_from_list(const List* list, SetCheckType check);
/*lint -sem(        set_range_new, @p == 1) */
extern Set*         set_range_new(int begin, int end, int step);
/*lint -sem(        set_prod_new, 1p == 1 && 2p == 1, @p == 1) */
extern Set*         set_prod_new(const Set* a, const Set* b);
/*lint -sem(        set_union, 1p == 1 && 2p == 1, @p == 1) */
extern Set*         set_union(const Set* seta, const Set* setb);
/*lint -sem(        set_inter, 1p == 1 && 2p == 1, @p == 1) */
extern Set*         set_inter(const Set* set_a, const Set* set_b);
/*lint -sem(        set_minus, 1p == 1 && 2p == 1, @p == 1) */
extern Set*         set_minus(const Set* set_a, const Set* set_b);
/*lint -sem(        set_sdiff, 1p == 1 && 2p == 1, @p == 1) */
extern Set*         set_sdiff(const Set* set_a, const Set* set_b);
/*lint -sem(        set_proj, 1p == 1 && 2p == 1, @p == 1) */
extern Set*         set_proj(const Set* set_a, const Tuple* pattern);
/*lint -sem(        set_is_subseteq, 1p == 1 && 2p == 1) */
extern Bool         set_is_subseteq(const Set* set_a, const Set* set_b);
/*lint -sem(        set_is_subset, 1p == 1 && 2p == 1) */
extern Bool         set_is_subset(const Set* set_a, const Set* set_b);
/*lint -sem(        set_is_equal, 1p == 1 && 2p == 1) */
extern Bool         set_is_equal(const Set* set_a, const Set* set_b);
/*lint -sem(        set_subset_list, 1p == 1 && 2n > 0 && 4p == 1, @p == 1) */
extern List*        set_subsets_list(
   const Set* set, int subset_size, List* list, int* idx);

#ifdef __cplusplus
}
#endif
#endif /* _SET_H_ */
