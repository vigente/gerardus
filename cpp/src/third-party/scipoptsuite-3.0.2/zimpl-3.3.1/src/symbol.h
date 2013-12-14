/* $Id: symbol.h,v 1.6 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: symbol.h                                                      */
/*   Name....: Symbol Table Functions                                        */
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
#ifndef _SYMBOL_H_
#define _SYMBOL_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before symbol.h"
#endif
#ifndef _TUPLE_H_
#error "Need to include tuple.h before symbol.h"
#endif
#ifndef _SET_H_
#error "Need to include set.h before symbol.h"
#endif
#ifndef _MME_H_
#error "Need to include mme.h before symbol.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*lint -sem(        symbol_new, nulterm(1), 1p && 3p == 1 && 4n >= 0, @p == 1) */
extern Symbol*      symbol_new(const char* name,
   SymbolType type, const Set* set, int estimated_size, const Entry* deflt);
extern void         symbol_exit(void);
/*lint -sem(        symbol_is_valid, 1p == 1) */
extern Bool         symbol_is_valid(const Symbol* symbol);
/*lint -sem(        symbol_lookup, nulterm(1), 1p, r_null) */
extern Symbol*      symbol_lookup(const char* name);
/*lint -sem(        symbol_has_entry, 1p == 1 && 2p == 1) */
extern Bool         symbol_has_entry(const Symbol* sym, const Tuple* tuple);
/*lint -sem(        symbol_lookup_entry, 1p == 1 && 2p == 1) */
extern const Entry* symbol_lookup_entry(const Symbol* sym, const Tuple* tuple);
/*lint -sem(        symbol_add_entry, custodial(2), 1p == 1 && 2p == 1) */
extern void         symbol_add_entry(Symbol* sym, Entry* entry);
/*lint -sem(        symbol_get_dim, 1p == 1, @n >= 0) */
extern int          symbol_get_dim(const Symbol* sym);
/*lint -sem(        symbol_get_iset, 1p == 1, @p == 1) */
extern const Set*   symbol_get_iset(const Symbol* sym);
/*lint -sem(        symbol_get_name, 1p == 1, @p && nulterm(@)) */
extern const char*  symbol_get_name(const Symbol* sym);
/*lint -sem(        symbol_get_type, 1p == 1) */
extern SymbolType   symbol_get_type(const Symbol* sym);
/*lint -sem(        symbol_get_numb, 1p == 1) */
extern const Numb*  symbol_get_numb(const Symbol* sym, int idx);
/*lint -sem(        symbol_get_strg, 1p == 1, @p && nulterm(@)) */
extern const char*  symbol_get_strg(const Symbol* sym, int idx);
/*lint -sem(        symbol_get_set, 1p == 1, @p == 1) */
extern const Set*   symbol_get_set(const Symbol* sym, int idx);
/*lint -sem(        symbol_get_var, 1p == 1, @p == 1) */
extern Var*         symbol_get_var(const Symbol* sym, int idx);
/*lint -sem(        symbol_print, 1p == 1 && 2p == 1) */
extern void         symbol_print(FILE* fp, const Symbol* sym);
/*lint -sem(        symbol_print_all, 1p == 1) */
extern void         symbol_print_all(FILE* fp);

#ifdef __cplusplus
}
#endif
#endif /* _SYMBOL_H_ */
