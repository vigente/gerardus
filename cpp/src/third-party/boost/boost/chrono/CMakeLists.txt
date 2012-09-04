# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostChrono CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostConfig)
ryppl_find_and_use_package(BoostSystem)
ryppl_find_and_use_package(BoostRatio)

include_directories(include)

if(CMAKE_SYSTEM MATCHES "Linux")
  find_library(RT_LIBRARY rt)
endif(CMAKE_SYSTEM MATCHES "Linux") 

add_library(boost_chrono
  src/chrono.cpp
  src/thread_clock.cpp
  src/process_cpu_clocks.cpp
  )

target_link_libraries(boost_chrono
  boost_system
  ${RT_LIBRARY}
  )

# ryppl_export(
#   TARGETS
#     boost_chrono
#   DEPENDS
#     BoostAccumulators
#     BoostConfig
#     BoostCore
#     BoostException
#     BoostFormat
#     BoostInteger
#     BoostInterprocess
#     BoostMPL
#     BoostMath
#     BoostRatio
#     BoostStaticAssert
#     BoostSystem
#     BoostTypeTraits
#     BoostTypeof
#     BoostUtility
#   INCLUDE_DIRECTORIES
#     include
#   )
