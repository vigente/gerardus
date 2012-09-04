# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostTest)

if(MSVC)
  # <asynch-exceptions>on
endif(MSVC)

ryppl_test_suite(
  RUN
    convenience_test.cpp
    deprecated_test.cpp
    fstream_test.cpp
    large_file_support_test.cpp
    locale_info.cpp
    operations_test.cpp
    operations_unit_test.cpp
    path_test.cpp
    path_unit_test.cpp
    ../example/simple_ls.cpp
    ../example/file_status.cpp
  LINK_LIBRARIES
    boost_prg_exec_monitor
    boost_filesystem
    boost_system
  NO_SINGLE_TARGET
  )
