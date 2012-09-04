# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostRandom CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostConfig)

add_library(boost_random
  src/random_device.cpp
  )
  
#ryppl_add_test_subdirectory(test)



ryppl_export(
  TARGETS
    boost_random
  DEPENDS
    BoostConfig
    BoostCore
    BoostInteger
    BoostMPL
    BoostMath
    BoostRange
    BoostStaticAssert
    BoostTypeTraits
    BoostUtility
  INCLUDE_DIRECTORIES
    include
  )
