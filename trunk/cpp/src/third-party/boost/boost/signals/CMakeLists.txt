# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostSignals CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostAny)
ryppl_find_and_use_package(BoostConfig)
ryppl_find_and_use_package(BoostCore)
ryppl_find_and_use_package(BoostFunction)
ryppl_find_and_use_package(BoostIterator)
ryppl_find_and_use_package(BoostMPL)
ryppl_find_and_use_package(BoostOptional)
ryppl_find_and_use_package(BoostSmartPtr)
ryppl_find_and_use_package(BoostTypeTraits)
ryppl_find_and_use_package(BoostUtility)

include_directories(include)

if(MSVC)
  add_definitions(
    /wd4251 # class needs to have dll-interface to be used by clients
    /wd4275 # 'boost::noncopyable' used as base for dll-interface class 
    )
endif(MSVC)

add_library(boost_signals
  src/trackable.cpp
  src/connection.cpp
  src/named_slot_map.cpp
  src/signal_base.cpp
  src/slot.cpp
  )

#ryppl_add_test_subdirectory(test)

# ryppl_export(
#   TARGETS
#     boost_signals
#   DEPENDS
#     BoostAny
#     BoostConfig
#     BoostCore
#     BoostFunction
#     BoostIterator
#     BoostMPL
#     BoostOptional
#     BoostSmartPtr
#     BoostTypeTraits
#     BoostUtility
#   INCLUDE_DIRECTORIES
#     include
#   )
