/* $Id: rdefpar.h,v 1.6 2012/07/29 15:09:29 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: rdefpar.c                                                     */
/*   Name....: Read Definition / Parameter                                   */
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
#ifndef _RDEFPAR_H_
#define _RDEFPAR_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct read_param        RPar;
typedef struct read_definition   RDef;

/*lint -sem(       rdef_new, nulterm(1), nulterm(2), 1p && 2p, @p == 1) */
extern RDef*       rdef_new(const char* filename, const char* pattern);
/*lint -sem(       rdef_free, custodial(1), 1p == 1) */
extern void        rdef_free(RDef* rdef);
/*lint -sem(       rdef_is_valid, 1p == 1) */
extern Bool        rdef_is_valid(const RDef* rdef);
/*lint -sem(       rdef_copy, 1p == 1, @p == 1) */
extern RDef*       rdef_copy(const RDef* rdef);
/*lint -sem(       rdef_set_param, 1p == 1 && 2p == 1) */
extern void        rdef_set_param(RDef* rdef, const RPar* rpar);
/*lint -sem(       rdef_get_filename, 1p == 1, @p && nulterm(@)) */
extern const char* rdef_get_filename(const RDef* rdef);
/*lint -sem(       rdef_get_pattern, 1p == 1, @p && nulterm(@)) */
extern const char* rdef_get_pattern(const RDef* rdef);
/*lint -sem(       rdef_get_comment, 1p == 1, @p && nulterm(@)) */
extern const char* rdef_get_comment(const RDef* rdef);
/*lint -sem(       rdef_get_match, 1p == 1) */
extern const char* rdef_get_match(const RDef* rdef);
/*lint -sem(       rdef_get_use, 1p == 1) */
extern int         rdef_get_use(const RDef* rdef);
/*lint -sem(       rdef_get_skip, 1p == 1) */
extern int         rdef_get_skip(const RDef* rdef);

/*lint -sem(       rpar_new_skip, @p == 1) */
extern RPar*       rpar_new_skip(int skip);
/*lint -sem(       rpar_new_use, @p == 1) */
extern RPar*       rpar_new_use(int use);
/*lint -sem(       rpar_new_comment, 1p && nulterm(1), @p == 1) */
extern RPar*       rpar_new_comment(const char* comment);
/*lint -sem(       rpar_new_match, 1p && nulterm(1), @p == 1) */
extern RPar*       rpar_new_match(const char* match);
/*lint -sem(       rpar_free, custodial(1), 1p == 1) */
extern void        rpar_free(RPar* rpar);
/*lint -sem(       rpar_is_valid, 1p == 1) */
extern Bool        rpar_is_valid(const RPar* rpar);
/*lint -sem(       rpar_copy, 1p == 1, @p == 1) */
extern RPar*       rpar_copy(const RPar* rpar);

#ifdef __cplusplus
}
#endif
#endif /* _RDEFPAR_H_ */

