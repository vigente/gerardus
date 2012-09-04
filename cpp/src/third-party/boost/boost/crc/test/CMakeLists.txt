# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

include(Ryppl)
include(RypplTestSuite)

ryppl_find_and_use_package(BoostCRC)
ryppl_find_and_use_package(BoostTest)
ryppl_find_and_use_package(BoostRandom)

ryppl_test_suite(
  RUN
    crc_test
  LINK_LIBRARIES
    boost_test_exec_monitor
  )
