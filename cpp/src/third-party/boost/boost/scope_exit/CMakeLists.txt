# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostScopeExit NONE)

include(Ryppl)

ryppl_export(
  DEPENDS
    BoostConfig
    BoostCore
    BoostFunction
    BoostLocalFunction
    BoostMPL
    BoostPreprocessor
    BoostTypeTraits
    BoostTypeof
    BoostUtility
  INCLUDE_DIRECTORIES
    include
  )
