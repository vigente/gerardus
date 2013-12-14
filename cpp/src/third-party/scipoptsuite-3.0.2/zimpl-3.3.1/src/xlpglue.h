/* $Id: xlpglue.h,v 1.27 2012/07/29 15:09:31 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: xlpglue.h                                                     */
/*   Name....: Glue between numb/term and ratlp                              */
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
#ifndef _XLPGLUE_H_
#define _XLPGLUE_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before xlpglue.h"
#endif

#ifndef _RATLPTYPES_H_
#error "Need to include ratlptypes.h before xlpglue.h"
#endif

#ifndef _MME_H_
#error "Need to include mme.h before xlpglue.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*lint -sem(    xlp_alloc, nulterm(1) && 1p, @p == 1) */
extern Lps*     xlp_alloc(const char* name, Bool need_startval, void* user_data);
/*lint -sem(    xlp_free, 1p == 1) */
extern void     xlp_free(Lps* lp);
/*lint -sem(    xlp_conname, nulterm(1), 1p == 1 && 2p && nulterm(2)) */
extern Bool     xlp_conname_exists(const Lps* lp, const char* conname);
/*lint -sem(    xlp_addcon_term, 1p == 1 && nulterm(2) && 2p && 4p == 1 && 5p == 1 && 7p == 1) */
extern Bool     xlp_addcon_term(Lps* lp, const char* name, ConType type,
   const Numb* lhs, const Numb* rhs, unsigned int flags, const Term* term);
/*lint -sem(    xlp_addvar, 1p == 1 && nulterm(2) && 2p && 4p == 1 && 5p == 1 && 6p == 1 && 7p == 1, @p == 1) */
extern Var*     xlp_addvar(Lps* lp, const char* name, VarClass usevarclass,
   const Bound* lower, const Bound* upper, const Numb* priority, const Numb* startval);
/*lint -sem(    xlp_addsos, 1p == 1 && nulterm(2) && 2p && 4p == 1 && 5p == 1, @p == 1) */
extern int      xlp_addsos_term(Lps* lp, const char* name, SosType type,
   const Numb* priority, const Term* term);
/*lint -sem(    xlp_getvarname, 1p == 1 && 2p == 1, nulterm(@) */
const char*     xlp_getvarname(const Lps* lp, const Var* var);
/*lint -sem(    xlp_getclass, 1p == 1 && 2p == 1) */
extern VarClass xlp_getclass(const Lps* lp, const Var* var);
/*lint -sem(    xlp_getlower, 1p == 1 && 2p == 1, @p == 1) */
extern Bound*   xlp_getlower(const Lps* lp, const Var* var);
/*lint -sem(    xlp_getupper, 1p == 1, @p == 1) */
extern Bound*   xlp_getupper(const Lps* lp, const Var* var);
/*lint -sem(    xlp_objname, 1p == 1 && nulterm(2) && 2p) */
extern void     xlp_objname(Lps* lp, const char* name);
/*lint -sem(    xlp_setdir, 1p == 1) */
extern void     xlp_setdir(Lps* lp, Bool minimize);
/*lint -sem(    xlp_addtocost, 1p == 1 && 2p == 1 && 3p == 1) */
extern void     xlp_addtocost(Lps* lp, Var* var, const Numb* cost);

#ifdef __cplusplus
}
#endif
#endif /* _XLPGLUE_H */








