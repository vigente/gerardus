# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostPropertyMap NONE)

include(Ryppl)
  
#ryppl_add_test_subdirectory(test)


ryppl_export(
  DEPENDS
    BoostAny
    BoostBind
    BoostConceptCheck
    BoostConfig
    BoostConversion
    BoostCore
    BoostException
    BoostFunction
    BoostGraph
    BoostGraphParallel
    BoostIterator
    BoostMPL
    BoostMultiIndex
    BoostOptional
    BoostSmartPtr
    BoostTypeTraits
    BoostUtility
  RECOMMENDS
    BoostSerialization
  INCLUDE_DIRECTORIES
    include
  )
