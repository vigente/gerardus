# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostCore)
ryppl_find_and_use_package(BoostTest)
ryppl_find_and_use_package(BoostForeach)

if(MSVC)
  # <asynch-exceptions>on
endif(MSVC)

ryppl_test_suite(
  RUN
    stl_byval.cpp
    stl_byref.cpp
    stl_byval_r.cpp
    stl_byref_r.cpp
    array_byval.cpp
    array_byref.cpp
    array_byval_r.cpp
    array_byref_r.cpp
    cstr_byval.cpp
    cstr_byref.cpp
    cstr_byval_r.cpp
    cstr_byref_r.cpp
    pair_byval.cpp
    pair_byref.cpp
    pair_byval_r.cpp
    pair_byref_r.cpp
    user_defined.cpp
    call_once.cpp
    rvalue_const.cpp
    rvalue_nonconst.cpp
    rvalue_const_r.cpp
    rvalue_nonconst_r.cpp
    dependent_type.cpp
    misc.cpp
  COMPILE
    noncopyable.cpp
  NO_SINGLE_TARGET
  )
