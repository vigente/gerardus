/* $Id: inst.h,v 1.52 2012/07/29 15:09:27 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   File....: inst.h                                                        */
/*   Name....: Instruction Functions                                         */
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

#ifndef _INST_H_
#define _INST_H_

#ifdef __cplusplus
extern "C" {
#endif

#define INST_NULL ((Inst)0)

/* inst.c
 */
/*lint -sem(     i_bool_and, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_and(CodeNode* self);
/*lint -sem(     i_bool_eq, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_eq(CodeNode* self);
/*lint -sem(     i_bool_exists, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_exists(CodeNode* self);
/*lint -sem(     i_bool_false, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_false(CodeNode* self);
/*lint -sem(     i_bool_ge, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_ge(CodeNode* self);
/*lint -sem(     i_bool_gt, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_gt(CodeNode* self);
/*lint -sem(     i_bool_is_elem, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_is_elem(CodeNode* self);
/*lint -sem(     i_bool_le, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_le(CodeNode* self);
/*lint -sem(     i_bool_lt, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_lt(CodeNode* self);
/*lint -sem(     i_bool_ne, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_ne(CodeNode* self);
/*lint -sem(     i_bool_not, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_not(CodeNode* self);
/*lint -sem(     i_bool_or, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_or(CodeNode* self);
/*lint -sem(     i_bool_seq, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_seq(CodeNode* self);
/*lint -sem(     i_bool_sneq, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_sneq(CodeNode* self);
/*lint -sem(     i_bool_sseq, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_sseq(CodeNode* self);
/*lint -sem(     i_bool_subs, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_subs(CodeNode* self);
/*lint -sem(     i_bool_true, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_true(CodeNode* self);
/*lint -sem(     i_bool_xor, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bool_xor(CodeNode* self);
/*lint -sem(     i_bound_new, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_bound_new(CodeNode* self);
/*lint -sem(     i_check, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_check(CodeNode* self);
/*lint -sem(     i_constraint_list, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_constraint_list(CodeNode* self);
/*lint -sem(     i_constraint, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_constraint(CodeNode* self);
/*lint -sem(     i_rangeconst, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_rangeconst(CodeNode* self);
/*lint -sem(     i_elem_list_add, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_elem_list_add(CodeNode* self);
/*lint -sem(     i_elem_list_new, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_elem_list_new(CodeNode* self);
/*lint -sem(     i_entry, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_entry(CodeNode* self);
/*lint -sem(     i_entry_list_add, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_entry_list_add(CodeNode* self);
/*lint -sem(     i_entry_list_new, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_entry_list_new(CodeNode* self);
/*lint -sem(     i_entry_list_powerset, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_entry_list_powerset(CodeNode* self);
/*lint -sem(     i_entry_list_subsets, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_entry_list_subsets(CodeNode* self);
/*lint -sem(     i_expr_abs, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_abs(CodeNode* self);
/*lint -sem(     i_expr_sgn, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_sgn(CodeNode* self);
/*lint -sem(     i_expr_add, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_add(CodeNode* self);
/*lint -sem(     i_expr_card, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_card(CodeNode* self);
/*lint -sem(     i_expr_ceil, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_ceil(CodeNode* self);
/*lint -sem(     i_expr_div, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_div(CodeNode* self);
/*lint -sem(     i_expr_exp, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_exp(CodeNode* self);
/*lint -sem(     i_expr_sqrt, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_sqrt(CodeNode* self);
/*lint -sem(     i_expr_fac, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_fac(CodeNode* self);
/*lint -sem(     i_expr_floor, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_floor(CodeNode* self);
/*lint -sem(     i_expr_if_else, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_if_else(CodeNode* self);
/*lint -sem(     i_expr_intdiv, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_intdiv(CodeNode* self);
/*lint -sem(     i_expr_length, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_length(CodeNode* self);
/*lint -sem(     i_expr_ln, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_ln(CodeNode* self);
/*lint -sem(     i_expr_log, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_log(CodeNode* self);
/*lint -sem(     i_expr_ord, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_ord(CodeNode* self);
/*lint -sem(     i_expr_prod, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_prod(CodeNode* self);
/*lint -sem(     i_expr_rand, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_rand(CodeNode* self);
/*lint -sem(     i_expr_round, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_round(CodeNode* self);
/*lint -sem(     i_expr_sum, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_sum(CodeNode* self);
/*lint -sem(     i_expr_max, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_max(CodeNode* self);
/*lint -sem(     i_expr_max2, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_max2(CodeNode* self);
/*lint -sem(     i_expr_sglmax, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_sglmax(CodeNode* self);
/*lint -sem(     i_expr_min, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_min(CodeNode* self);
/*lint -sem(     i_expr_min2, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_min2(CodeNode* self);
/*lint -sem(     i_expr_sglmin, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_sglmin(CodeNode* self);
/*lint -sem(     i_expr_mul, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_mul(CodeNode* self);
/*lint -sem(     i_expr_mod, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_mod(CodeNode* self);
/*lint -sem(     i_expr_neg, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_neg(CodeNode* self);
/*lint -sem(     i_expr_pow, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_pow(CodeNode* self);
/*lint -sem(     i_expr_sub, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_sub(CodeNode* self);
/*lint -sem(     i_expr_substr, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_expr_substr(CodeNode* self);
/*lint -sem(     i_forall, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_forall(CodeNode* self);
/*lint -sem(     i_idxset_new, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_idxset_new(CodeNode* self);
/*lint -sem(     i_idxset_pseudo_new, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_idxset_pseudo_new(CodeNode* self);
/*lint -sem(     i_list_matrix, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_list_matrix(CodeNode* self);
/*lint -sem(     i_local_deref, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_local_deref(CodeNode* self);
/*lint -sem(     i_matrix_list_new, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_matrix_list_new(CodeNode* self);
/*lint -sem(     i_matrix_list_add, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_matrix_list_add(CodeNode* self);
/*lint -sem(     i_newdef, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_newdef(CodeNode* self);
/*lint -sem(     i_newsym_para1, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_newsym_para1(CodeNode* self);
/*lint -sem(     i_newsym_para2, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_newsym_para2(CodeNode* self);
/*lint -sem(     i_newsym_set1, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_newsym_set1(CodeNode* self);
/*lint -sem(     i_newsym_set2, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_newsym_set2(CodeNode* self);
/*lint -sem(     i_newsym_var, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_newsym_var(CodeNode* self);
/*lint -sem(     i_nop, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_nop(CodeNode* self);
/*lint -sem(     i_object_max, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_object_max(CodeNode* self);
/*lint -sem(     i_object_min, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_object_min(CodeNode* self);
/*lint -sem(     i_print, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_print(CodeNode* self);
/*lint -sem(     i_set_argmax, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_argmax(CodeNode* self);
/*lint -sem(     i_set_argmin, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_argmin(CodeNode* self);
/*lint -sem(     i_set_cross, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_cross(CodeNode* self);
/*lint -sem(     i_set_empty, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_empty(CodeNode* self);
/*lint -sem(     i_set_expr, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_expr(CodeNode* self);
/*lint -sem(     i_set_idxset, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_idxset(CodeNode* self);
/*lint -sem(     i_set_indexset, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_indexset(CodeNode* self);
/*lint -sem(     i_set_inter, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_inter(CodeNode* self);
/*lint -sem(     i_set_inter2, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_inter2(CodeNode* self);
/*lint -sem(     i_set_minus, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_minus(CodeNode* self);
/*lint -sem(     i_set_new_tuple, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_new_tuple(CodeNode* self);
/*lint -sem(     i_set_new_elem, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_new_elem(CodeNode* self);
/*lint -sem(     i_set_proj, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_proj(CodeNode* self);
/*lint -sem(     i_set_pseudo, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_pseudo(CodeNode* self);
/*lint -sem(     i_set_range, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_range(CodeNode* self);
/*lint -sem(     i_set_range2, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_range2(CodeNode* self);
/*lint -sem(     i_set_sdiff, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_sdiff(CodeNode* self);
/*lint -sem(     i_set_union, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_union(CodeNode* self);
/*lint -sem(     i_set_union2, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_set_union2(CodeNode* self);
/*lint -sem(     i_sos, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_sos(CodeNode* self);
/*lint -sem(     i_soset, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_soset(CodeNode* self);
/*lint -sem(     i_subto, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_subto(CodeNode* self);
/*lint -sem(     i_symbol_deref, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_symbol_deref(CodeNode* self);
/*lint -sem(     i_define_deref, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_define_deref(CodeNode* self);
/*lint -sem(     i_term_add, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_term_add(CodeNode* self);
/*lint -sem(     i_term_coeff, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_term_coeff(CodeNode* self);
/*lint -sem(     i_term_const, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_term_const(CodeNode* self);
/*lint -sem(     i_term_expr, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_term_expr(CodeNode* self);
/*lint -sem(     i_term_mul, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_term_mul(CodeNode* self);
/*lint -sem(     i_term_quadratic, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_term_power(CodeNode* self);
/*lint -sem(     i_term_power, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_term_sub(CodeNode* self);
/*lint -sem(     i_term_sum, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_term_sum(CodeNode* self);
/*lint -sem(     i_tuple_new, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_tuple_new(CodeNode* self);
/*lint -sem(     i_tuple_empty, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_tuple_empty(CodeNode* self);
/*lint -sem(     i_tuple_list_add, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_tuple_list_add(CodeNode* self);
/*lint -sem(     i_tuple_list_new, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_tuple_list_new(CodeNode* self);


/* iread.c
 */
/*lint -sem(     i_read_new, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_read_new(CodeNode* self);
/*lint -sem(     i_read_param, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_read_param(CodeNode* self);
/*lint -sem(     i_read_comment, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_read_comment(CodeNode* self);
/*lint -sem(     i_read_match, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_read_match(CodeNode* self);
/*lint -sem(     i_read_use, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_read_use(CodeNode* self);
/*lint -sem(     i_read_skip, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_read_skip(CodeNode* self);
/*lint -sem(     i_read, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_read(CodeNode* self);

/* vinst.c
 */
/*lint -sem(     i_vabs, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vabs(CodeNode* self);
/*lint -sem(     i_vbool_and, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_and(CodeNode* self);
/*lint -sem(     i_vbool_eq, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_eq(CodeNode* self);
/*lint -sem(     i_vbool_ne, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_ne(CodeNode* self);
/*lint -sem(     i_vbool_ge, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_ge(CodeNode* self);
/*lint -sem(     i_vbool_gt, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_gt(CodeNode* self);
/*lint -sem(     i_vbool_le, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_le(CodeNode* self);
/*lint -sem(     i_vbool_lt, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_lt(CodeNode* self);
/*lint -sem(     i_vbool_not, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_not(CodeNode* self);
/*lint -sem(     i_vbool_or, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_or(CodeNode* self);
/*lint -sem(     i_vbool_xor, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vbool_xor(CodeNode* self);
/*lint -sem(     i_vexpr_fun, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vexpr_fun(CodeNode* self);
/*lint -sem(     i_vif, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vif(CodeNode* self);
/*lint -sem(     i_vif_else, 1p == 1, type(1), @p == 1p) */
extern CodeNode* i_vif_else(CodeNode* self);

#ifdef __cplusplus
}
#endif
#endif /* _INST_H_ */
