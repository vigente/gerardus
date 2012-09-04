# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
ryppl_find_and_use_package(Boost COMPONENTS test context NO_MODULE)
include(${Boost_USE_FILE})
include(RypplTestSuite)

ryppl_test_suite(
  RUN
    test_context.cpp
  LINK_LIBRARIES
    unit_test_framework-static
		context-static
  )
