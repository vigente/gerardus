# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostCore)
ryppl_find_and_use_package(BoostTest)

ryppl_test_suite(
  RUN
    format_test1.cpp
    format_test2.cpp
    format_test3.cpp
    format_test_wstring.cpp
    format_test_enum.cpp
  LINK_LIBRARIES
    boost_test_exec_monitor
  NO_SINGLE_TARGET
  )
