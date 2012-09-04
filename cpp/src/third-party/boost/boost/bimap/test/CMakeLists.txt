# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostBimap)
ryppl_find_and_use_package(BoostSerialization)
ryppl_find_and_use_package(BoostAssign)
ryppl_find_and_use_package(BoostTest)

# XXX: some tests #include <libs/bimap/test/test_bimap.hpp>
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
set(test_bimap_hpp
  "${CMAKE_CURRENT_BINARY_DIR}/include/libs/bimap/test/test_bimap.hpp"
  )
if(NOT EXISTS "${test_bimap_hpp}")
  file(WRITE "${test_bimap_hpp}"
    "#include \"${CMAKE_CURRENT_SOURCE_DIR}/test_bimap.hpp\"\n"
    )
endif()

# tagged_test
ryppl_test_suite(test_tagged.cpp RUN)

# relation_test
ryppl_test_suite(
  RUN
    test_structured_pair.cpp
    test_mutant.cpp
    test_mutant_relation.cpp
  NO_SINGLE_TARGET
  )

# bimap_test
ryppl_test_suite(
  RUN
    # Check library user interface
    test_bimap_set_of.cpp
    test_bimap_multiset_of.cpp
    test_bimap_unordered_set_of.cpp
    test_bimap_unordered_multiset_of.cpp
    test_bimap_list_of.cpp
    test_bimap_vector_of.cpp
    # Test bimap container
    test_bimap_ordered.cpp
    test_bimap_unordered.cpp
    test_bimap_sequenced.cpp
    test_bimap_unconstrained.cpp
    test_bimap_assign.cpp
    test_bimap_property_map.cpp
    test_bimap_modify.cpp
    test_bimap_range.cpp
    test_bimap_operator_bracket.cpp
    test_bimap_lambda.cpp
    test_bimap_mutable.cpp
    test_bimap_extra.cpp
    test_bimap_convenience_header.cpp
    test_bimap_project.cpp
    test_bimap_serialization.cpp
    test_bimap_info.cpp
  LINK_LIBRARIES
    boost_serialization
  NO_SINGLE_TARGET
  )

# compile_fail_test
ryppl_test_suite(
  COMPILE_FAIL
    compile_fail/test_bimap_mutable_1.cpp
    compile_fail/test_bimap_mutable_2.cpp
    compile_fail/test_bimap_mutable_3.cpp
    compile_fail/test_bimap_info_1.cpp
    compile_fail/test_bimap_info_2.cpp
    compile_fail/test_bimap_info_3.cpp
  NO_SINGLE_TARGET
  )

# bimap_and_boost
ryppl_test_suite(
  RUN
    ../example/bimap_and_boost/property_map.cpp
    ../example/bimap_and_boost/range.cpp
    ../example/bimap_and_boost/foreach.cpp
    ../example/bimap_and_boost/lambda.cpp
    ../example/bimap_and_boost/assign.cpp
    ../example/bimap_and_boost/xpressive.cpp
    ../example/bimap_and_boost/typeof.cpp
  COMPILE
    ../example/bimap_and_boost/serialization.cpp
  LINK_LIBRARIES
    boost_serialization
  NO_SINGLE_TARGET
  )
