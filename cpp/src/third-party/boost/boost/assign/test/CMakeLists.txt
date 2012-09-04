# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostAssign)
ryppl_find_and_use_package(BoostFunctional)
ryppl_find_and_use_package(BoostMultiIndex)
ryppl_find_and_use_package(BoostSerialization)
ryppl_find_and_use_package(BoostTest)

ryppl_test_suite(
  RUN
    basic.cpp
    std.cpp
    list_of.cpp
    ptr_list_of.cpp
    static_list_of.cpp
    tuple_list_of.cpp
    list_inserter.cpp
    ptr_list_inserter.cpp
    ptr_map_inserter.cpp
    list_of_workaround.cpp
    email_example.cpp
    my_vector_example.cpp
    multi_index_container.cpp
  LINK_LIBRARIES
    boost_unit_test_framework
  NO_SINGLE_TARGET
  )
