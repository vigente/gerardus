# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostXpressive NONE)

include(Ryppl)

ryppl_add_test_subdirectory(test)

ryppl_export(
  DEPENDS
    BoostConfig
    BoostConversion
    BoostCore
    BoostException
    BoostFusion
    BoostInteger
    BoostIterator
    BoostMPL
    BoostNumericConversion
    BoostOptional
    BoostPreprocessor
    BoostProto
    BoostRange
    BoostSmartPtr
    BoostStaticAssert
    BoostTypeTraits
    BoostTypeof
    BoostUtility
#   BoostNumeric 
    BoostOptional 
    BoostPreprocessor 
    BoostProto 
    BoostRange 
    BoostSmartPtr 
    BoostSpirit 
    BoostTypeTraits 
    BoostTypeof 
    BoostUtility 
  INCLUDE_DIRECTORIES
    include
  )
