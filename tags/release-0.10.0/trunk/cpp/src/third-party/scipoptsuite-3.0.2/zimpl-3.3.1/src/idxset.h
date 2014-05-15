/* $Id: idxset.h,v 1.6 2012/07/29 15:09:27 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: idxset.h                                                      */
/*   Name....: IndexSet Functions                                            */
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
#ifndef _IDXSET_H_
#define _IDXSET_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before idxset.h"
#endif
#ifndef _TUPLE_H_
#error "Need to include tuple.h before idxset.h"
#endif
#ifndef _SET_H_
#error "Need to include set.h before idxset.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct index_set         IdxSet;

/* idxset.c
 */
/*lint -sem(        idxset_new, 1p == 1 && 2p == 1 && 3p == 1, @p == 1) */
extern IdxSet*      idxset_new(
   const Tuple* tuple, const Set* set, CodeNode* lexpr, Bool is_unrestricted);
/*lint -sem(        idxset_free, custodial(1), 1p == 1) */
extern void         idxset_free(IdxSet* idxset);
/*lint -sem(        idxset_is_valid, 1p == 1) */
extern Bool         idxset_is_valid(const IdxSet* idxset);
/*lint -sem(        idxset_copy, 1p == 1, @p == 1) */
extern IdxSet*      idxset_copy(const IdxSet* source);
/*lint -sem(        idxset_get_lexpr, 1p == 1, @p == 1) */
extern CodeNode*    idxset_get_lexpr(const IdxSet* idxset);
/*lint -sem(        idxset_get_tuple, 1p == 1, @p == 1) */
extern const Tuple* idxset_get_tuple(const IdxSet* idxset);
/*lint -sem(        idxset_get_set, 1p == 1, @p == 1) */
extern const Set*   idxset_get_set(const IdxSet* idxset);
/*lint -sem(        idxset_is_unrestricted, 1p == 1) */
extern Bool         idxset_is_unrestricted(const IdxSet* idxset);
/*lint -sem(        idxset_print, 1p == 1 && 2p == 1) */
extern void         idxset_print(FILE* fp, const IdxSet* idxset);

#ifdef __cplusplus
}
#endif
#endif /* _IDXSET_H_ */
