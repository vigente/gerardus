# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostContainer)

file(GLOB sources *.cpp)
ryppl_test_suite(RUN ${sources} NO_SINGLE_TARGET)
