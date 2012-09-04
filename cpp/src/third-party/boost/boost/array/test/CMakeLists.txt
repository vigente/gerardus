# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostArray)

ryppl_test_suite(
  RUN
    array0.cpp
    array1.cpp
    array2.cpp
    array3.cpp
    array4.cpp
    array5.cpp
    array6.cpp
  NO_SINGLE_TARGET
  )
