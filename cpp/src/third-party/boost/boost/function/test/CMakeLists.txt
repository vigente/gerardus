# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostCore)
ryppl_find_and_use_package(BoostTest)

if(MSVC)
  # <asynch-exceptions>on
endif(MSVC)

ryppl_test_suite(
  COMPILE
    function_typeof_test.cpp
    function_30.cpp
  COMPILE_FAIL
    function_test_fail1.cpp
    function_test_fail2.cpp
  RUN
    function_test.cpp
    function_n_test.cpp
    allocator_test.cpp
    stateless_test.cpp
    lambda_test.cpp
    function_arith_cxx98.cpp
    function_arith_portable.cpp
    sum_avg_cxx98.cpp
    sum_avg_portable.cpp
    mem_fun_cxx98.cpp
    mem_fun_portable.cpp
    std_bind_cxx98.cpp
    std_bind_portable.cpp
    function_ref_cxx98.cpp
    function_ref_portable.cpp
    contains_test.cpp
    contains2_test.cpp
    nothrow_swap.cpp
  LINK_LIBRARIES
    boost_test_exec_monitor
  NO_SINGLE_TARGET
  )
