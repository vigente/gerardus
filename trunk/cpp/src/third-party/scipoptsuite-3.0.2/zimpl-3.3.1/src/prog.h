/* $Id: prog.h,v 1.6 2012/07/29 15:09:28 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: prog.h                                                        */
/*   Name....: Program Functions                                             */
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
#ifndef _PROG_H_
#define _PROG_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before prog.h"
#endif
#ifndef _STMT_H_
#error "Need to include stmt.h before prog.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct program           Prog;

/* prog.c
 */
extern void*        prog_get_lp(void);
/*lint -sem(        prog_new, 1p) */
extern Prog*        prog_new(void);
/*lint -sem(        prog_free, custodial(1), 1p == 1) */
extern void         prog_free(Prog* prog);
/*lint -sem(        prog_is_valid, 1p == 1) */
extern Bool         prog_is_valid(const Prog* prog);
/*lint -sem(        prog_is_empty, 1p == 1) */
extern Bool         prog_is_empty(const Prog* prog);
/*lint -sem(        prog_add_stmt, custodial(2), 1p == 1 && 2p == 1) */
extern void         prog_add_stmt(Prog* prog, Stmt* stmt);
/*lint -sem(        prog_print, 1p == 1 && 2p == 1) */
extern void         prog_print(FILE* fp, const Prog* prog);
/*lint -sem(        prog_execute, 1p == 1) */
extern void         prog_execute(const Prog* prog, void* lp);
/*lint -sem(        prog_tostr, 1p == 1 && nulterm(2) && 3n > 0, nulterm(@)) */
extern char*        prog_tostr(const Prog* prog, const char* prefix, const char* title, int max_output_line_len);

/* load.c
 */
/*lint -sem(        prog_load, nulterm(3), 1p == 1) */
extern void         prog_load(Prog* prog, const char* cmd, const char* filename);

#ifdef __cplusplus
}
#endif
#endif /* _PROG_H_ */
