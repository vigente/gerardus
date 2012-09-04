# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostCore)
ryppl_find_and_use_package(BoostTest)
ryppl_find_and_use_package(BoostTimer)

# Search tests
ryppl_test_suite(
  RUN
    empty_search_test.cpp
    search_test1.cpp
    search_test2.cpp
    search_test3.cpp
  COMPILE_FAIL
    search_fail1.cpp
    search_fail2.cpp
    search_fail3.cpp
  NO_SINGLE_TARGET
  )

# Clamp tests
ryppl_test_suite(
  RUN
    clamp_test.cpp
  NO_SINGLE_TARGET
  )

# Cxx11 tests
ryppl_test_suite(
  RUN
    all_of_test.cpp
    any_of_test.cpp
    none_of_test.cpp
    one_of_test.cpp
    ordered_test.cpp
    find_if_not_test1.cpp
    copy_n_test1.cpp
    iota_test1.cpp
    is_permutation_test1.cpp
    partition_point_test1.cpp
    is_partitioned_test1.cpp
    partition_copy_test1.cpp
  NO_SINGLE_TARGET
  )

# Hex tests
ryppl_test_suite(
  RUN
    hex_test1.cpp
    hex_test2.cpp
    hex_test3.cpp
  COMPILE_FAIL
    hex_fail1.cpp
  NO_SINGLE_TARGET
  )
