/* $Id: metaio.h,v 1.6 2012/07/29 15:09:27 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: metaio.h                                                      */
/*   Name....: Meta Input/Output                                             */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * Copyright (C) 2006-2012 by Thorsten Koch <koch@zib.de>
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
#ifndef _METAIO_H_
#define _METAIO_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before metaio.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct meta_file_ptr     MFP;

/* metaio.c
 */
/*lint -sem(        mio_add_strg_file, nulterm(1), nulterm(2), 1p && 2p) */
extern void         mio_add_strg_file(const char* name, const char* content, Bool use_copy);
extern void         mio_init(void);
extern void         mio_exit(void);
/*lint -sem(        mio_open, nulterm(2), 2p) */
/*lint -function(   fopen(1), mio_open(1)) */
/*lint -function(   fopen(r), mio_open(r)) */
extern MFP*         mio_open(const char* name, const char* ext);
/*lint -function(   fclose, mio_close) */
extern void         mio_close(MFP* mfp);
/*lint -function(   fgetc, mio_getc) */
extern int          mio_getc(const MFP* mfp);
/*lint -function(   fgets(1), mio_gets(2)) */
/*lint -function(   fgets(2), mio_gets(3)) */
/*lint -function(   fgets(3), mio_gets(1)) */
/*lint -function(   fgets(r), mio_gets(r)) */
extern char*        mio_gets(const MFP* mfp, char* buf, int len);
/*lint -sem(        mio_get_line, 1p) */
extern char*        mio_get_line(const MFP* mfp);

#ifdef __cplusplus
}
#endif
#endif /* _METAIO_H_ */
