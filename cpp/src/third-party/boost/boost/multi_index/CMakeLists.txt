# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostMultiIndex NONE)

include(Ryppl)
  
ryppl_add_test_subdirectory(test)


ryppl_export(
  DEPENDS
    BoostBind
    BoostConfig
    BoostCore
    BoostException
    BoostForeach
    BoostFunctionalHash
    BoostIterator
    BoostMPL
    BoostMath
    BoostPreprocessor
    BoostStaticAssert
    BoostTuple
    BoostTypeTraits
    BoostUtility
  RECOMMENDS
    BoostSerialization
  INCLUDE_DIRECTORIES
    include
  )
