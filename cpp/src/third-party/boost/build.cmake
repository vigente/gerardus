#!/usr/bin/env cmake -P

##########################################################################
# Copyright (C) 2011-2012 Daniel Pfeifer <daniel@pfeifer-mail.de>        #
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.             #
# See accompanying file LICENSE_1_0.txt or copy at                       #
#   http://www.boost.org/LICENSE_1_0.txt                                 #
##########################################################################

# cmake -DBUILDDIR=../build -DGENERATOR="Unix Makefiles" -DBUILDSTEP=configure -P build.cmake


if(NOT DEFINED BUILDDIR)
  set(BUILDDIR "${CMAKE_CURRENT_LIST_DIR}/build")
elseif(NOT IS_ABSOLUTE "${BUILDDIR}")
  set(BUILDDIR "${CMAKE_CURRENT_LIST_DIR}/${BUILDDIR}")
endif()


if(DEFINED TOOLCHAIN AND NOT DEFINED GENERATOR)
  if(TOOLCHAIN STREQUAL "vs10")
    set(GENERATOR "Visual Studio 10")
  elseif(TOOLCHAIN STREQUAL "vs9")
    set(GENERATOR "Visual Studio 9 2008")
  elseif(TOOLCHAIN STREQUAL "vs8")
    set(GENERATOR "Visual Studio 8 2005")
  elseif(TOOLCHAIN STREQUAL "vs7.1")
    set(GENERATOR "Visual Studio 7 .NET 2003")
  else()
    set(GENERATOR "Unix Makefiles")
  endif()
  message(WARNING
    "TOOLCHAIN has been set to \"${TOOLCHAIN}\".\n"
    "GENERATOR should be set instead!\n"
    "Trying \"${GENERATOR}\"..."
    )
endif(DEFINED TOOLCHAIN AND NOT DEFINED GENERATOR)


if(CMAKE_HOST_WIN32 AND NOT DEFINED GENERATOR)
  if(EXISTS "$ENV{VS100COMNTOOLS}vsvars32.bat")
    set(GENERATOR "Visual Studio 10")
  elseif(EXISTS "$ENV{VS90COMNTOOLS}vsvars32.bat")
    set(GENERATOR "Visual Studio 9 2008")
  elseif(EXISTS "$ENV{VS80COMNTOOLS}vsvars32.bat")
    set(GENERATOR "Visual Studio 8 2005")
  elseif(EXISTS "$ENV{VS71COMNTOOLS}vsvars32.bat")
    set(GENERATOR "Visual Studio 7 .NET 2003")
  endif()
endif(CMAKE_HOST_WIN32 AND NOT DEFINED GENERATOR)


if(GENERATOR MATCHES "^Visual Studio")
  set(MULTI_CONFIG TRUE)
else()
  set(MULTI_CONFIG FALSE)
endif()


if(BUILDSTEP STREQUAL "cleanconfigure")
  file(REMOVE_RECURSE "${BUILDDIR}")
endif()


if(NOT DEFINED BUILDSTEP OR BUILDSTEP MATCHES "configure$")
  if(MULTI_CONFIG)
    if(NOT EXISTS "${BUILDDIR}/CMakeCache.txt")
      file(MAKE_DIRECTORY "${BUILDDIR}")
      execute_process(COMMAND "${CMAKE_COMMAND}" "-G${GENERATOR}"
        "${toolchain_param}" "${CMAKE_CURRENT_LIST_DIR}"
        WORKING_DIRECTORY "${BUILDDIR}"
        )
    endif(NOT EXISTS "${BUILDDIR}/CMakeCache.txt")
    execute_process(COMMAND "${CMAKE_COMMAND}" "${BUILDDIR}"
      RESULT_VARIABLE result
      )
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "Configuring ${config} failed.")
    endif(NOT result EQUAL 0)
  else()
    foreach(config Debug Release)
      if(NOT EXISTS "${BUILDDIR}/${config}/CMakeCache.txt")
        file(MAKE_DIRECTORY "${BUILDDIR}/${config}")
        execute_process(COMMAND "${CMAKE_COMMAND}" "-G${GENERATOR}"
          "${toolchain_param}" -DCMAKE_BUILD_TYPE=${config} 
          "${CMAKE_CURRENT_LIST_DIR}"
          WORKING_DIRECTORY "${BUILDDIR}/${config}"
          )
      endif(NOT EXISTS "${BUILDDIR}/${config}/CMakeCache.txt")
      execute_process(COMMAND "${CMAKE_COMMAND}" "${BUILDDIR}/${config}"
        RESULT_VARIABLE result
        )
      if(NOT result EQUAL 0)
        message(FATAL_ERROR "Configuring ${config} failed.")
      endif(NOT result EQUAL 0)
    endforeach(config)
  endif()
endif(NOT DEFINED BUILDSTEP OR BUILDSTEP MATCHES "configure$")


function(build config)
  message(STATUS "Building ${ARGV1} ${config}")

  if(MULTI_CONFIG)
    set(args "${BUILDDIR}" --config ${config})
  else()
    set(args "${BUILDDIR}/${config}")
  endif()

  if(ARGV1)
    list(APPEND args --target ${ARGV1})
  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} --build ${args}
    RESULT_VARIABLE result
    )
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Building failed.")
  endif(NOT result EQUAL 0)
endfunction(build)


if(NOT DEFINED BUILDSTEP OR BUILDSTEP STREQUAL "build")
  build(Debug)
  build(Release)
endif(NOT DEFINED BUILDSTEP OR BUILDSTEP STREQUAL "build")

if(NOT DEFINED BUILDSTEP OR BUILDSTEP STREQUAL "test")
  build(Debug test)
  build(Release test)
endif(NOT DEFINED BUILDSTEP OR BUILDSTEP STREQUAL "test")

if(NOT DEFINED BUILDSTEP OR BUILDSTEP STREQUAL "documentation")
  build(Release documentation)
endif(NOT DEFINED BUILDSTEP OR BUILDSTEP STREQUAL "documentation")
