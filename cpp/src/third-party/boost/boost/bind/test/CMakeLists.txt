# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostCore)

ryppl_test_suite(
  RUN
    bind_test.cpp
    bind_dm_test.cpp
    bind_eq_test.cpp
    bind_const_test.cpp
    bind_cv_test.cpp
    bind_stateful_test.cpp
    bind_dm2_test.cpp
    bind_not_test.cpp
    bind_rel_test.cpp
    bind_function_test.cpp
    bind_lookup_problem_test.cpp
    bind_rv_sp_test.cpp
    bind_dm3_test.cpp
    bind_visit_test.cpp
    bind_placeholder_test.cpp
    bind_rvalue_test.cpp
    bind_and_or_test.cpp
    mem_fn_test.cpp
    mem_fn_void_test.cpp
    mem_fn_derived_test.cpp
    mem_fn_eq_test.cpp
    mem_fn_dm_test.cpp
    mem_fn_rv_test.cpp
    ref_fn_test.cpp
    bind_fnobj2_test.cpp
    bind_fn2_test.cpp
    bind_mf2_test.cpp
    bind_eq2_test.cpp
    mem_fn_ref_test.cpp
    bind_ref_test.cpp
    bind_eq3_test.cpp
    protect_test.cpp
    mem_fn_unary_addr_test.cpp
  COMPILE
    bind_unary_addr.cpp
  NO_SINGLE_TARGET
  )
