/* $Id: list.h,v 1.6 2012/07/29 15:09:27 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: list.h                                                        */
/*   Name....: List Functions                                                */
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
#ifndef _LIST_H_
#define _LIST_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before list.h"
#endif

#ifndef _MME_H_
#error "Need to include mme.h before list.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define LIST_NULL ((List*)0)

/*lint -sem(        list_new_elem, 1p == 1, @p == 1) */
extern List*        list_new_elem(const Elem* elem);
/*lint -sem(        list_new_tuple, 1p == 1, @p == 1) */
extern List*        list_new_tuple(const Tuple* tuple);
/*lint -sem(        list_new_entry, 1p == 1, @p == 1) */
extern List*        list_new_entry(const Entry* entry);
/*lint -sem(        list_new_list, 1p == 1, @p == 1) */
extern List*        list_new_list(const List* list);
/*lint -sem(        list_free, custodial(1), 1p == 1) */
extern void         list_free(List* list);
/*lint -sem(        list_is_valid, 1p == 1) */
extern Bool         list_is_valid(const List* list);
/*lint -sem(        list_is_elemlist, 1p == 1) */
extern Bool         list_is_elemlist(const List* list);
/*lint -sem(        list_is_entrylist, 1p == 1) */
extern Bool         list_is_entrylist(const List* list);
/*lint -sem(        list_is_tuplelist, 1p == 1) */
extern Bool         list_is_tuplelist(const List* list);
/*lint -sem(        list_copy, 1p == 1, @p == 1) */
extern List*        list_copy(const List* list);
/*lint -sem(        list_add_list, 1p == 1 && 2p == 1) */
extern void         list_add_list(List* list, const List* ll);
/*lint -sem(        list_add_elem, 1p == 1 && 2p == 1) */
extern void         list_add_elem(List* list, const Elem* elem);
/*lint -sem(        list_add_tuple, 1p == 1 && 2p == 1) */
extern void         list_add_tuple(List* list, const Tuple* tuple);
/*lint -sem(        list_add_entry, 1p == 1 && 2p == 1) */
extern void         list_add_entry(List* list, const Entry* entry);
/*lint -sem(        list_insert_elem, 1p == 1 && 2p == 1) */
extern void         list_insert_elem(List* list, const Elem* elem);
/*lint -sem(        list_insert_tuple, 1p == 1 && 2p == 1) */
extern void         list_insert_tuple(List* list, const Tuple* tuple);
/*lint -sem(        list_insert_entry, 1p == 1 && 2p == 1) */
extern void         list_insert_entry(List* list, const Entry* entry);
/*lint -sem(        list_get_elems, 1p == 1, @n >= 0) */
extern int          list_get_elems(const List* list);
/*lint -sem(        list_get_elem,  1p == 1, @p == 1) */
extern const Elem*  list_get_elem(const List* list, ListElem** idxp);
/*lint -sem(        list_get_tuple, 1p == 1, @p == 1) */
extern const Tuple* list_get_tuple(const List* list, ListElem** idxp);
/*lint -sem(        list_get_entry, 1p == 1, @p == 1) */
extern const Entry* list_get_entry(const List* list, ListElem** idxp);
/*lint -sem(        list_get_list, 1p == 1, @p == 1) */
extern const List*  list_get_list(const List* list, ListElem** idxp);
/*lint -sem(        list_print, 1p == 1 && 2p == 1) */
extern void         list_print(FILE* fp, const List* list);

#ifdef __cplusplus
}
#endif
#endif /* _LIST_H_ */
