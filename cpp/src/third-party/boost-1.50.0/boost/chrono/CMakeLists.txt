# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

## Changes by Ramon Casero <rcasero@gmail.com> for project Gerardus

cmake_minimum_required(VERSION 2.8.8)
project(BoostChrono CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostConfig)
ryppl_find_and_use_package(BoostSystem)
#ryppl_find_and_use_package(BoostRatio)

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

ryppl_export(
  TARGETS
    boost_chrono
  DEPENDS
#     BoostAccumulators
    BoostConfig
    BoostCore
#     BoostException
#     BoostFormat
#     BoostInteger
#     BoostInterprocess
#     BoostMPL
#     BoostMath
#     BoostRatio
#     BoostStaticAssert
    BoostSystem
#     BoostTypeTraits
#     BoostTypeof
#     BoostUtility
  INCLUDE_DIRECTORIES
    include
  )

# gerardus: Primarly, we need to generate DLLs for the boost libraries
# so that they can be linked into MEX files. However, CMake for
# Windows doesn't seem to know how to link executables to the DLLs,
# which means that we have to generate additional static libraries too
if(WIN32 AND NOT CYGWIN)
  add_library(boost_chronolib
    STATIC
    src/chrono.cpp
    src/thread_clock.cpp
    src/process_cpu_clocks.cpp
    )
  set_target_properties(boost_chronolib PROPERTIES
    OUTPUT_NAME "boost_chrono"
    )
  
  ryppl_export(
    TARGETS
    boost_chronolib
    DEPENDS
    #     BoostAccumulators
    BoostConfig
    BoostCore
    #     BoostException
    #     BoostFormat
    #     BoostInteger
    #     BoostInterprocess
    #     BoostMPL
    #     BoostMath
    #     BoostRatio
    #     BoostStaticAssert
    BoostSystem
    #     BoostTypeTraits
    #     BoostTypeof
    #     BoostUtility
    INCLUDE_DIRECTORIES
    include
    )
endif()
