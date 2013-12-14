/* $Id: stmt.h,v 1.6 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: stmt.h                                                        */
/*   Name....: Statement Functions                                           */
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
#ifndef _STMT_H_
#define _STMT_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before stmt.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum statement_type
{
   STMT_ERR = 0, STMT_SET, STMT_PARAM, STMT_VAR, STMT_MIN, STMT_MAX,
   STMT_CONS, STMT_DEF, STMT_DO, STMT_SOS
};

typedef enum statement_type      StmtType;
typedef struct statement         Stmt;

/* stmt.c
 */
/*lint -sem(        stmt_new, nulterm(2), nulterm(4), 2p && 3n >= 0 && 4p, @p == 1) */
extern Stmt*        stmt_new(StmtType type, const char* filename, int lineno,
   const char* text);
/*lint -sem(        stmt_free, custodial(1), 1p == 1) */
extern void         stmt_free(Stmt* stmt);
/*lint -sem(        stmt_is_valid, 1p == 1) */
extern Bool         stmt_is_valid(const Stmt* stmt);
/*lint -sem(        stmt_get_filename, 1p == 1, @p && nulterm(@)) */
extern const char*  stmt_get_filename(const Stmt* stmt);
/*lint -sem(        stmt_get_lineno, 1p == 1, @n > 0) */
extern int          stmt_get_lineno(const Stmt* stmt);
/*lint -sem(        stmt_get_text, 1p == 1, @p && nulterm(@)) */
extern const char*  stmt_get_text(const Stmt* stmt);
/*lint -sem(        stmt_parse, 1p == 1) */
extern void         stmt_parse(Stmt* stmt);
/*lint -sem(        stmt_execute, 1p == 1) */
extern void         stmt_execute(const Stmt* stmt);
/*lint -sem(        stmt_print, 1p == 1 && 2p == 1) */
extern void         stmt_print(FILE* fp, const Stmt* stmt);
/*lint -sem(        stmt_trigger_warning, 1n >= 0) */
extern Bool         stmt_trigger_warning(int no);

/* mmlparse.y
 */
extern int          yyparse(void);

/* mmlscan.l
 */
/*lint -sem(        parse_stmt, 1p) */
extern void         parse_stmt(const Stmt* stmt);
extern const Stmt*  scan_get_stmt(void);
extern int          scan_get_column(void);

#ifdef __cplusplus
}
#endif
#endif /* _STMT_H_ */
