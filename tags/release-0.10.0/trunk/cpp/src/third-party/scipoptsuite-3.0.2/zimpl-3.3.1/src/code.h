/* $Id: code.h,v 1.29 2012/07/29 15:09:26 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: code.h                                                        */
/*   Name....: Code Functions                                                */
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

#ifndef _CODE_H_
#define _CODE_H_

#ifndef _BOOL_H_
#error "Need to include bool.h before code.h"
#endif

#ifndef _RATLPTYPES_H_
#error "Need to include ratlptypes.h before code.h"
#endif

#ifndef _MME_H_
#error "Need to include mme.h before code.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*lint -sem(        code_new_inst, 1p == 1 && 2n >= 0, @p == 1) */
extern CodeNode*    code_new_inst(Inst inst, int childs, ...);
/*lint -sem(        code_new_numb, custodial(1), 1p == 1, @p == 1) */
extern CodeNode*    code_new_numb(Numb* numb);
/*lint -sem(        code_new_strg, nulterm(1), 1p, @p == 1) */
extern CodeNode*    code_new_strg(const char* strg);
/*lint -sem(        code_new_name, nulterm(1), 1p, @p == 1) */
extern CodeNode*    code_new_name(const char* name);
/*lint -sem(        code_new_size, @p == 1) */
extern CodeNode*    code_new_size(int size);
/*lint -sem(        code_new_varclass, @p == 1) */
extern CodeNode*    code_new_varclass(VarClass varclass);
/*lint -sem(        code_new_contype, @p == 1) */
extern CodeNode*    code_new_contype(ConType contype);
/*lint -sem(        code_new_bits, @p == 1) */
extern CodeNode*    code_new_bits(unsigned int bits);
/*lint -sem(        code_new_symbol, @p == 1) */
extern CodeNode*    code_new_symbol(Symbol* sym);
/*lint -sem(        code_new_define, @p == 1) */
extern CodeNode*    code_new_define(Define* def);
/*lint -sem(        code_new_bound, @p == 1) */
extern CodeNode*    code_new_bound(BoundType type);
/*lint -sem(        code_free, custodial(1), 1p == 1) */
extern void         code_free(CodeNode* node);
/*lint -sem(        code_free_value, @p == 1) */
extern void         code_free_value(CodeNode* node);
/*lint -sem(        code_is_valid, 1p == 1) */
extern Bool         code_is_valid(const CodeNode* node);
/*lint -sem(        code_get_type, 1p == 1) */
extern CodeType     code_get_type(const CodeNode* node);
/*lint -sem(        code_get_inst, 1p == 1, @p) */
extern Inst         code_get_inst(const CodeNode* node);
/*lint -sem(        code_set_root, custodial(1), 1p == 1) */
extern void         code_set_root(CodeNode* node);
/*lint -sem(        code_get_root, @p == 1) */
extern CodeNode*    code_get_root(void);
/*lint -sem(        code_set_child, custodial(3), 1p == 1 && 2n >= 0 && 3p == 1) */
extern void         code_set_child(CodeNode* node, int idx, CodeNode* child);
/*lint -sem(        code_errmsg, 1p == 1) */
extern void         code_errmsg(const CodeNode* node);
#if 0
/*lint -sem(        code_check_type, 1p == 1, @p == 1p) */
extern CodeNode*    code_check_type(CodeNode* node, CodeType expected);
#endif
/*lint -sem(        code_eval, 1p == 1, @p == 1p) */
extern CodeNode*    code_eval(CodeNode* node);
/*lint -sem(        code_prune_tree, 1p == 1) */
extern Bool         code_prune_tree(CodeNode* node);
/*lint -sem(        code_get_child, 1p == 1 && 2n >= 0, @p == 1) */
extern CodeNode*    code_get_child(const CodeNode* node, int no);
/*lint -sem(        code_get_numb, 1p == 1) */
extern const Numb*  code_get_numb(CodeNode* node);
/*lint -sem(        code_get_strg, 1p == 1, @p && nulterm(@)) */
extern const char*  code_get_strg(CodeNode* node);
/*lint -sem(        code_get_name, 1p == 1, @p && nulterm(@)) */
extern const char*  code_get_name(CodeNode* node);
extern unsigned int code_get_inst_count(void);
/*lint -sem(        code_get_tuple, 1p == 1, @p == 1) */
extern void         code_clear_inst_count(void);
extern const Tuple* code_get_tuple(CodeNode* node);
/*lint -sem(        code_get_set, 1p == 1, @p == 1) */
extern const Set*   code_get_set(CodeNode* node);
/*lint -sem(        code_get_idxset, 1p == 1, @p == 1) */
extern const IdxSet* code_get_idxset(CodeNode* node);
/*lint -sem(        code_get_entry, 1p == 1, @p == 1) */
extern const Entry* code_get_entry(CodeNode* node);
/*lint -sem(        code_get_term, 1p == 1, @p == 1) */
extern const Term*  code_get_term(CodeNode* node);
/*lint -sem(        code_get_size, 1p == 1) */
extern int          code_get_size(CodeNode* node);
/*lint -sem(        code_get_bool, 1p == 1) */
extern Bool         code_get_bool(CodeNode* node);
/*lint -sem(        code_get_list, 1p == 1, @p == 1) */
extern const List*  code_get_list(CodeNode* node);
/*lint -sem(        code_get_varclass, 1p == 1) */
extern VarClass     code_get_varclass(CodeNode* node);
/*lint -sem(        code_get_contype, 1p == 1) */
extern ConType      code_get_contype(CodeNode* node);
/*lint -sem(        code_get_rdef, 1p == 1, @p == 1) */
extern const RDef*  code_get_rdef(CodeNode* node);
/*lint -sem(        code_get_rpar, 1p == 1, @p == 1) */
extern const RPar*  code_get_rpar(CodeNode* node);
/*lint -sem(        code_get_bits, 1p == 1) */
extern unsigned int code_get_bits(CodeNode* node);
/*lint -sem(        code_get_symbol, 1p == 1) */
extern Symbol*      code_get_symbol(CodeNode* node);
/*lint -sem(        code_get_define, 1p == 1) */
extern Define*      code_get_define(CodeNode* node);
/*lint -sem(        code_get_bound, 1p == 1) */
extern const Bound* code_get_bound(CodeNode* node);
/*lint -sem(        code_value_numb, custodial(2), 1p == 1) */
extern void         code_value_numb(CodeNode* node, Numb* numb);
/*lint -sem(        code_value_strg, nulterm(2), 1p == 1 && 2p) */
extern void         code_value_strg(CodeNode* node, const char* strg);
/*lint -sem(        code_value_name, nulterm(2), 1p == 1 && 2p) */
extern void         code_value_name(CodeNode* node, const char* name);
/*lint -sem(        code_value_tuple, custodial(2), 1p == 1 && 2p == 1) */
extern void         code_value_tuple(CodeNode* node, Tuple* tuple);
/*lint -sem(        code_value_set, custodial(2), 1p == 1 && 2p == 1) */
extern void         code_value_set(CodeNode* node, Set* set);
/*lint -sem(        code_value_idxset, custodial(2), 1p == 1 && 2p == 1) */
extern void         code_value_idxset(CodeNode* node, IdxSet* idxset);
/*lint -sem(        code_value_entry, custodial(2), 1p == 1 && 2p == 1) */
extern void         code_value_entry(CodeNode* node, Entry* entry);
/*lint -sem(        code_value_term, custodial(2), 1p == 1 && 2p == 1) */
extern void         code_value_term(CodeNode* node, Term* term);
/*lint -sem(        code_value_steal_child_term, 1p == 1 && 2n >= 0, @p == 1p) */
extern Term*        code_value_steal_term(CodeNode* node, int no);
/*lint -sem(        code_value_bool, 1p == 1) */
extern void         code_value_bool(CodeNode* node, Bool bool);
/*lint -sem(        code_value_size, 1p == 1) */
extern void         code_value_size(CodeNode* node, int size);
/*lint -sem(        code_value_list, custodial(2), 1p == 1 && 2p == 1) */
extern void         code_value_list(CodeNode* node, List* list);
/*lint -sem(        code_value_varclass, 1p == 1) */
extern void         code_value_varclass(CodeNode* node, VarClass varclass);
/*lint -sem(        code_value_contype, 1p == 1) */
extern void         code_value_contype(CodeNode* node, ConType contype);
/*lint -sem(        code_value_rdef, custodial(2), 1p == 1 && 2p == 1) */
extern void         code_value_rdef(CodeNode* node, RDef* rdef);
/*lint -sem(        code_value_rpar, custodial(2), 1p == 1 && 2p == 1) */
extern void         code_value_rpar(CodeNode* node, RPar* rpar);
/*lint -sem(        code_value_bits, 1p == 1) */
extern void         code_value_bits(CodeNode* node, unsigned int bits);
/*lint -sem(        code_value_bound, custodial(2), 1p == 1 && 2p == 1) */
extern void         code_value_bound(CodeNode* node, Bound* bound);
/*lint -sem(        code_value_void, 1p == 1) */
extern void         code_value_void(CodeNode* node);
/*lint -sem(        code_copy_value, 1p == 1 && 2p == 1) */
extern void         code_copy_value(CodeNode* dst, const CodeNode* src);

/*lint -sem(        code_eval_child, 1p == 1 && 2n >= 0, @p == 1p) */
extern CodeNode*    code_eval_child(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_numb, 1p == 1 && 2n >= 0) */
extern const Numb*  code_eval_child_numb(const CodeNode* node, int no);
/*lint -sem(        char* code_eval_child_strg, 1p == 1 && 2n >= 0,
                    @p && nulterm(@)) */
extern const char*  code_eval_child_strg(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_name, 1p == 1 && 2n >= 0,
                    @p && nulterm(@)) */
extern const char*  code_eval_child_name(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_tuple, 1p == 1 && 2n >= 0, @p == 1) */
extern const Tuple* code_eval_child_tuple(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_set, 1p == 1 && 2n >= 0, @p == 1) */
extern const Set*   code_eval_child_set(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_idxset, 1p == 1 && 2n >= 0, @p == 1) */
extern const IdxSet* code_eval_child_idxset(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_entry, 1p == 1 && 2n >= 0, @p == 1) */
extern const Entry* code_eval_child_entry(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_term, 1p == 1 && 2n >= 0, @p == 1) */
extern const Term*  code_eval_child_term(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_size, 1p == 1 && 2n >= 0, @n >= 0) */
extern int          code_eval_child_size(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_bool, 1p == 1 && 2n >= 0) */
extern Bool         code_eval_child_bool(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_list, 1p == 1 && 2n >= 0, @p == 1) */
extern const List*  code_eval_child_list(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_varclass, 1p == 1 && 2n >= 0) */
extern VarClass     code_eval_child_varclass(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_contype, 1p == 1 && 2n >= 0) */
extern ConType      code_eval_child_contype(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_rdef, 1p == 1 && 2n >= 0, @p == 1) */
extern const RDef*  code_eval_child_rdef(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_rpar, 1p == 1 && 2n >= 0, @p == 1) */
extern const RPar*  code_eval_child_rpar(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_bits, 1p == 1 && 2n >= 0) */
extern unsigned int code_eval_child_bits(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_symbol, 1p == 1 && 2n >= 0) */
extern Symbol*      code_eval_child_symbol(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_define, 1p == 1 && 2n >= 0) */
extern Define*      code_eval_child_define(const CodeNode* node, int no);
/*lint -sem(        code_eval_child_bound, 1p == 1 && 2n >= 0) */
extern const Bound* code_eval_child_bound(const CodeNode* node, int no);

#ifdef __cplusplus
}
#endif
#endif /* _CODE_H_ */







