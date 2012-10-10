# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostStaticAssert NONE)

include(Ryppl)

# ryppl_export(
#   DEPENDS
#     BoostConfig
#     BoostCore
#   INCLUDE_DIRECTORIES
#     include
#   )

return()

boost_add_multiple_tests(
  RUN
    static_assert_test
#   static_assert_example_1
    static_assert_example_2
    static_assert_example_3
  COMPILE_FAIL
    static_assert_test_fail_1
    static_assert_test_fail_2
    static_assert_test_fail_3
    static_assert_test_fail_4
    static_assert_test_fail_5
    static_assert_test_fail_6
    static_assert_test_fail_7
    static_assert_test_fail_8
    static_assert_test_fail_9
    static_assert_test_fail_10
  )
