# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostContext CXX)

include(Ryppl)

ryppl_find_and_use_package(BoostBind)
ryppl_find_and_use_package(BoostConfig)
ryppl_find_and_use_package(BoostCore)
ryppl_find_and_use_package(BoostFormat)
ryppl_find_and_use_package(BoostInteger)
ryppl_find_and_use_package(BoostMove)
ryppl_find_and_use_package(BoostPreprocessor)
ryppl_find_and_use_package(BoostProgramOptions)
ryppl_find_and_use_package(BoostSmartPtr)
ryppl_find_and_use_package(BoostSystem)
ryppl_find_and_use_package(BoostThread)
ryppl_find_and_use_package(BoostTypeTraits)
ryppl_find_and_use_package(BoostUtility)

set(abi
  aapcs
  eabi
  ms
  n32
  n64
  o32
  o64
  sysv
  )

set(binary_format
  elf
  mach-o
  pe
  )

if(WIN32)
  set(default_abi ms)
elseif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "arm")
  set(default_abi aapcs)
elseif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "mips")
  set(default_abi o32)
else()
  set(default_abi sysv)
endif()
  
if(APPLE)
  set(default_binary_format mach-o)
elseif(WIN32)
  set(default_binary_format pe)
else()
  set(default_binary_format elf)
endif()

set(BOOST_CONTEXT_ABI ${default_abi}
  CACHE STRING "ABI for Boost.Context (${abi})"
  )
set_property(CACHE BOOST_CONTEXT_ABI
  PROPERTY STRINGS "${abi}"
  )

set(BOOST_CONTEXT_BINARY_FORMAT ${default_binary_format}
  CACHE STRING "Binary format for Boost.Context (${binary_format})"
  )
set_property(CACHE BOOST_CONTEXT_BINARY_FORMAT
  PROPERTY STRINGS "${binary_format}"
  )

# TODO:
set(ASM_CONTEXT_SOURCES
  src/asm/fcontext_i386_sysv_elf_gas.S
  )

if(WIN32)
  set(CONTEXT_SOURCES
    src/stack_allocator_windows.cpp
    src/stack_utils_windows.cpp
    src/seh.cpp
    )
else()
  set(CONTEXT_SOURCES
    src/stack_allocator_posix.cpp
    src/stack_utils_posix.cpp
    )
endif()

include_directories(include)

add_library(boost_context
  ${ASM_CONTEXT_SOURCES}
  ${CONTEXT_SOURCES}
  src/fcontext.cpp
  )

# add_subdirectory(test)



ryppl_export(
  TARGETS
    boost_context
  DEPENDS
    BoostConfig
    BoostCore
    BoostInteger
    BoostUtility
  INCLUDE_DIRECTORIES
    include
  )
