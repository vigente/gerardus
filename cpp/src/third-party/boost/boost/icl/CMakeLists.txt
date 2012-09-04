# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostICL NONE)

include(Ryppl)

ryppl_export(
  DEPENDS
    BoostConceptCheck
    BoostConfig
    BoostContainer
    BoostCore
    BoostDateTime
    BoostIterator
    BoostMPL
    BoostMove
    BoostRational
    BoostStaticAssert
    BoostTypeTraits
    BoostUtility
  INCLUDE_DIRECTORIES
    include
  )
