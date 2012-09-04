# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostBimap)

include(Ryppl)

#ryppl_add_test_subdirectory(test)

ryppl_export(
  DEPENDS
    BoostConceptCheck
    BoostConfig
    BoostException
    BoostFunctionalHash
    BoostIterator
    BoostLambda
    BoostMPL
    BoostMultiIndex
    BoostPreprocessor
    BoostPropertyMap
    BoostStaticAssert
    BoostTypeTraits
    BoostUtility
  RECOMMENDS
    BoostSerialization
  INCLUDE_DIRECTORIES
    include
  )
