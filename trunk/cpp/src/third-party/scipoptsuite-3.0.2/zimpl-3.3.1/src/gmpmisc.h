/* $Id: gmpmisc.h,v 1.10 2012/07/29 15:09:26 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: gmpmisc.c                                                     */
/*   Name....: miscellenious rational arithmetic functions                   */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Copyright (C) 2003-2012 by Thorsten Koch <koch@zib.de>
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
#ifndef _GMPMISC_H_
#define _GMPMISC_H_

#ifndef __GMP_H__
#error "Need to include gmp.h before gmpmisc.h"
#endif
#ifndef _BOOL_H_
#error "Need to include bool.h before gmpmisc.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern mpq_t const_zero;
extern mpq_t const_one;
extern mpq_t const_minus_one;

/*lint -sem(gmp_str2mpq, 1p && 2p && nulterm(2)) */
extern void gmp_str2mpq(mpq_t value, const char* num);
/*lint -sem(gmp_print, 1p == 1 && 2p ) */
extern void gmp_print_mpq(FILE* fp, const mpq_t qval);
extern void gmp_init(Bool verb, Bool with_management);
extern void gmp_exit(void);

#ifdef __cplusplus
}
#endif
#endif /* _GMPMISC_H */
