# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

## Changes by Ramon Casero <rcasero@gmail.com> for project Gerardus

cmake_minimum_required(VERSION 2.8.8)
project(BoostThread C CXX) # C is required for find_package(Threads)

include(Ryppl)

ryppl_find_and_use_package(BoostChrono)
ryppl_find_and_use_package(BoostConfig)
ryppl_find_and_use_package(BoostDateTime)
ryppl_find_and_use_package(BoostSystem)
ryppl_find_package(Threads)

include_directories(
  include
  )

if(CMAKE_USE_WIN32_THREADS_INIT)
  set(sources
    src/win32/thread.cpp
    src/win32/tss_dll.cpp
    src/win32/tss_pe.cpp
    )
else()
  set(sources
    src/pthread/thread.cpp
    src/pthread/once.cpp
    )
endif()

add_library(boost_thread
  ${sources}
  )

target_link_libraries(boost_thread
  boost_chrono
  boost_date_time
  ${CMAKE_THREAD_LIBS_INIT}
  )

if(BUILD_SHARED_LIBS)
  set_property(TARGET boost_thread APPEND PROPERTY
    COMPILE_DEFINITIONS "BOOST_THREAD_BUILD_DLL=1"
    )
else()
  set_property(TARGET boost_thread APPEND PROPERTY
    COMPILE_DEFINITIONS "BOOST_THREAD_BUILD_LIB=1"
    )
  set_target_properties(boost_thread PROPERTIES
    PREFIX "lib"
    )
endif()

ryppl_export(
  TARGETS
    boost_thread
  DEPENDS
#     BoostBind
    BoostChrono
    BoostConfig
#     BoostContainer
    BoostCore
    BoostDateTime
#     BoostException
#     BoostFunction
#     BoostFunctionalHash
#     BoostIO
#     BoostInteger
#     BoostMPL
#     BoostMove
#     BoostOptional
#     BoostSmartPtr
#     BoostStaticAssert
    BoostSystem
#     BoostTypeTraits
#     BoostUtility
  INCLUDE_DIRECTORIES
    include
  )

# gerardus: Primarly, we need to generate DLLs for the boost libraries
# so that they can be linked into MEX files. However, CMake for
# Windows doesn't seem to know how to link executables to the DLLs,
# which means that we have to generate additional static libraries too
if(WIN32 AND NOT CYGWIN)
  add_library(boost_threadlib
    STATIC
    ${sources}
    )
  set_target_properties(boost_threadlib PROPERTIES
    OUTPUT_NAME "boost_thread"
    )
  
  ryppl_export(
    TARGETS
    boost_threadlib
    DEPENDS
    #     BoostBind
    BoostChrono
    BoostConfig
    #     BoostContainer
    BoostCore
    BoostDateTime
    #     BoostException
    #     BoostFunction
    #     BoostFunctionalHash
    #     BoostIO
    #     BoostInteger
    #     BoostMPL
    #     BoostMove
    #     BoostOptional
    #     BoostSmartPtr
    #     BoostStaticAssert
    BoostSystem
    #     BoostTypeTraits
    #     BoostUtility
    INCLUDE_DIRECTORIES
    include
    )
endif()
