/* $Id: stkchk.h,v 1.4 2012/07/29 15:09:30 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: stkchk.h                                                      */
/*   Name....: Stack usage checker                                           */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Copyright (C) 2011-2012 by Thorsten Koch <koch@zib.de>
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
#ifndef _STKCHK_H_
#define _STKCHK_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NDEBUG
extern const void* stkchk_start;
extern size_t      stkchk_maxi;

extern void   stkchk_init_x(void);
extern size_t stkchk_used_x(void);
/*lint -sem(  stkchk_maximum_x, 1p) */
extern void   stkchk_maximum_x(FILE* fp);
/*lint -sem(  stkchk_display_x, 1p) */
extern void   stkchk_display_x(FILE* fp);

#define stkchk_init()      stkchk_init_x()
#define stkchk_used()      stkchk_used_x()
#define stkchk_display(fp) stkchk_display_x(fp)
#define stkchk_maximum(fp) stkchk_maximum_x(fp)

#else 

#define stkchk_init()      /**/
#define stkchk_used()      /**/
#define stkchk_display(fp) /**/
#define stkchk_maximum(fp) /**/

#endif /* NDEBUG */

#ifdef __cplusplus
}
#endif
#endif /* _STKCHK_H_ */
