# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostException NONE)

include(Ryppl)

# add_subdirectory(test)


ryppl_export(
  DEPENDS
    BoostConfig
    BoostCore
    BoostSmartPtr
    BoostTuple
    BoostUnits
    BoostUtility
  INCLUDE_DIRECTORIES
    include
  )
