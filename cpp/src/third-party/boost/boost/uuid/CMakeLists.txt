# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostUUID NONE)

include(Ryppl)

ryppl_export(
  DEPENDS
    BoostConfig
    BoostException
    BoostIO
    BoostInteger
    BoostIterator
    BoostSmartPtr
    BoostStaticAssert
    BoostTypeTraits
    BoostUtility
  RECOMMENDS
    BoostRandom
    BoostSerialization
  INCLUDE_DIRECTORIES
    include
  )
