# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2011 University of Oxford
# Version: 0.3.1
# $Rev$
# $Date$
#
# University of Oxford means the Chancellor, Masters and Scholars of
# the University of Oxford, having an administrative office at
# Wellington Square, Oxford OX1 2JD, UK. 
#
# This file is part of Gerardus.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. The offer of this
# program under the terms of the License is subject to the License
# being interpreted in accordance with English Law and subject to any
# action against the University of Oxford being under the jurisdiction
# of the English Courts.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

# You need to run the following cmake command, with the appropriate
# paths, if you want to generate the project files that allow to
# develop this project from the Eclipse IDE
#
# cmake -DITK_DIR=/usr/local/lib/InsightToolkit  -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release  ./src

#############################################################################################
## Preliminary configuration, prior to defining the Gerardus project
#############################################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# for Visual Studio or other multi-configuration systems, make
# available only the Debug and Release configurations. This needs to
# go before PROJECT() or it will be ignored
SET(CMAKE_CONFIGURATION_TYPES Release Debug)
SET(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
  "Reset the configurations to Debug/Release only"
  FORCE)

# for single configuration systems (Unix-like), by default, choose the
# Release configuration
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# display type of build (e.g. Release)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# directories with ad hoc Findxxx.cmake files
SET(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/matlab)

#############################################################################################
## Matlab configuration
#############################################################################################

# hack to avoid error messages from CMake in Windows. The problem is that we must run 
# project() after we have set up the compiler version, to avoid an infinite loop in cmake.
# And we can only set up the compiler version once we know Matlab's version. Additionally,
# FindMatlab.cmake requires some variables to be set (CMAKE_FIND_LIBRARY_PREFIXES, 
# CMAKE_FIND_LIBRARY_SUFFIXES), which only happens after project(). Thus, we have a Catch22.
# This hack is simply to set the offending variables to the values they would have in Windows
if(WIN32)
   set(CMAKE_FIND_LIBRARY_PREFIXES "")
   set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
endif(WIN32)

# find Matlab
find_package(Matlab REQUIRED)
if(MATLAB_FOUND)
  message(STATUS "Matlab " ${MATLAB_VERSION} ": " ${MATLAB_ROOT})
endif(MATLAB_FOUND)

# some versions of Matlab or some libraries are not compatible with
# more advanced gcc versions
#
# the compiler version has to be set up before project(), otherwise we
# will get an infinite loop. We force the compiler version setting the
# corresponding environmental variables (equivalent to running e.g.  
# $ CC=gcc-4.4 CXX=g++-4.4 cmake ..  This is a better way than directly
# setting MAKE_C_COMPILER, CMAKE_CXX_COMPILER. For example, if we set these
# variables after project(), this creates an infinite loop.
if(NOT WIN32)
  if(BUILD_WITH_CUDA OR ${MATLAB_VERSION} MATCHES "R2012a")
    set(ENV{CC} gcc-4.4)
    set(ENV{CXX} g++-4.4)
  endif(BUILD_WITH_CUDA OR ${MATLAB_VERSION} MATCHES "R2012a")
endif(NOT WIN32)

#############################################################################################
## Start of Gerardus project
#############################################################################################

# project name
project(GERARDUS)

#############################################################################################
## Third-party libraries: GMP and MRPF
#############################################################################################

if(WIN32 AND NOT CYGWIN)

  # for Windows builds, we use the pre-compiled 64-bit libraries
  # provided by the CGAL Windows installer, that we distribute with
  # Gerardus
  set(GMP_INCLUDE_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/CGAL-4.2/auxiliary/gmp/include)
  set(MPFR_INCLUDE_DIR ${GMP_INCLUDE_DIR})
  set(GMP_LIBRARY_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/CGAL-4.2/auxiliary/gmp/lib)
  set(MPFR_LIBRARY_DIR ${GMP_LIBRARY_DIR})

else(WIN32 AND NOT CYGWIN)

  # for Unix-like builds, we have to build the GMP and MPFR libraries from the source code
  set(GMP_SOURCE_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/gmp-5.1.2)
  set(MPFR_SOURCE_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/mpfr-3.1.2)

  include(ExternalProject)

  ExternalProject_Add(GMP
    PREFIX ${GMP_SOURCE_DIR}
    STAMP_DIR ${GMP_SOURCE_DIR}/output
    SOURCE_DIR ${GMP_SOURCE_DIR}
    CONFIGURE_COMMAND ./configure --prefix=${GMP_SOURCE_DIR}/output
    BINARY_DIR ${GMP_SOURCE_DIR}
    BUILD_COMMAND make
    TEST_BEFORE_INSTALL 1
    TEST_COMMAND make check
    INSTALL_DIR ${GMP_SOURCE_DIR}/output
    INSTALL_COMMAND make install
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_TEST 1
    LOG_INSTALL 1
    )
  add_library(gmp STATIC IMPORTED GLOBAL)
  set_property(TARGET gmp PROPERTY IMPORTED_LOCATION ${GMP_SOURCE_DIR}/output/lib)
  set(GMP_FOUND TRUE 
    CACHE BOOL "GMP library is part of Gerardus, so this is always true")
  set(GMP_INCLUDE_DIR ${GMP_SOURCE_DIR}
    CACHE PATH "Path to GMP include files")
  set(GMP_LIBRARY_DIR ${GMP_SOURCE_DIR}/output/lib
    CACHE PATH "Path to GMP library")
  set(GMP_LIBRARY_DIRS ${GMP_SOURCE_DIR}/output/lib
    CACHE PATH "Path to GMP library")
  set(GMP_LIBRARIES_DIR ${GMP_SOURCE_DIR}/output/lib
    CACHE PATH "Path to GMP library")
  
  ExternalProject_Add(MPFR
    DEPENDS GMP
    PREFIX ${MPFR_SOURCE_DIR}
    STAMP_DIR ${MPFR_SOURCE_DIR}/output
    SOURCE_DIR ${MPFR_SOURCE_DIR}
    CONFIGURE_COMMAND ./configure --prefix=${MPFR_SOURCE_DIR}/output 
	--with-gmp-include=${GMP_SOURCE_DIR}/output/include --with-gmp-lib=${GMP_SOURCE_DIR}/output/lib
    BINARY_DIR ${MPFR_SOURCE_DIR}
    BUILD_COMMAND make
    TEST_BEFORE_INSTALL 1
    TEST_COMMAND make check
    INSTALL_DIR ${MPFR_SOURCE_DIR}/output
    INSTALL_COMMAND make install
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_TEST 1
    LOG_INSTALL 1
    )
  add_library(mpfr STATIC IMPORTED GLOBAL)
  set_property(TARGET mpfr PROPERTY IMPORTED_LOCATION ${MPFR_SOURCE_DIR}/output/lib)
  set(MPFR_FOUND TRUE 
    CACHE BOOL "MPFR library is part of Gerardus, so this is always true")
  set(MPFR_INCLUDE_DIR ${MPFR_SOURCE_DIR}/src
    CACHE PATH "Path to MPFR include files")
  set(MPFR_LIBRARY_DIR ${MPFR_SOURCE_DIR}/output/lib
    CACHE PATH "Path to MPFR library")
  set(MPFR_LIBRARY_DIRS ${MPFR_SOURCE_DIR}/output/lib
    CACHE PATH "Path to MPFR library")
  set(MPFR_LIBRARIES_DIR ${MPFR_SOURCE_DIR}/output/lib
    CACHE PATH "Path to MPFR library")

endif(WIN32 AND NOT CYGWIN)

#############################################################################################
## Third-party library: ITK
#############################################################################################

# find ITK
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
  MESSAGE(STATUS "ITK v" ${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH})
ENDIF(ITK_FOUND)

## flags for all targets
IF(WIN32)
  # to remove warnings from MSVC compiler asking to use fopen_s() an other unportable
  # functions exclusive to Microsoft
  ADD_DEFINITIONS("-D_CRT_SECURE_NO_DEPRECATE")
ELSE(WIN32)
  # optimise and show all warnings
  ADD_DEFINITIONS("-O2 -Wall")

  # For Mac OSX, we need to specify the system architecture.
  # If this hasn't been set, we choose 64-bit intel.
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    IF(NOT CMAKE_OSX_ARCHITECTURES)
      SET(CMAKE_OSX_ARCHITECTURES x86_64)
    ENDIF(NOT CMAKE_OSX_ARCHITECTURES)
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

ENDIF(WIN32)

#############################################################################################
## Third-party library: Boost
#############################################################################################

# Boost libraries are built as an external project. Attempts to
# integrate the build and install more closely with the Gerardus CMake
# build system have proven unsuccessful, in part because without using
# FindBoost.cmake, the linker gets confused between the local Boost
# and any system-wide installed Boost, and FindBoost.cmake runs in the
# configuration phase, before ExternalProject_Add() has the chance to
# build Boost.
#
# Hence, we use a hack:
#
#  $ cd gerardus/bin
#  $ cmake ..     ====> configure the build of Boost
#  $ make         ====> build and install Boost
#  $ cmake ..     ====> configure Gerardus and the rest of third party libraries
#  $ make         ====> build
#  $ make install ====> install

if(NOT EXISTS "${GERARDUS_SOURCE_DIR}/include/boost/version.hpp")

  message(STATUS "First run: Local Boost libraries for Gerardus need to be built")

  # build Boost libraries as an external project
  include(ExternalProject)

  if(WIN32)
    set(Boost_bootstrap "bootstrap.bat")
    set(Boost_b2 "b2")
    set(Boost_layout "versioned")
  else()
    set(Boost_bootstrap "./bootstrap.sh")
    set(Boost_b2 "./b2")
    set(Boost_layout "tagged")
  endif()
  
  ExternalProject_Add(Boost
    PREFIX ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/boost_1_53_0/
    SOURCE_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/boost_1_53_0/
    BINARY_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/boost_1_53_0/
    CONFIGURE_COMMAND ${Boost_bootstrap}
    BUILD_COMMAND ""
    TEST_BEFORE_INSTALL 0
    INSTALL_COMMAND ${Boost_b2} install --prefix=${GERARDUS_SOURCE_DIR} 
    --build-dir=./bin --layout=${Boost_layout} --build-type=complete  
    --with-chrono --with-date_time --with-filesystem --with-system --with-thread
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
    )

  return()

else()

  message(STATUS "Second run: Local Boost libraries for Gerardus found")

  # path to the Boost header files and libraries once they are built
  # and installed
  if(WIN32)
     set(Boost_INCLUDE_DIR "${GERARDUS_SOURCE_DIR}/include")
  else()
     set(Boost_INCLUDE_DIR "${GERARDUS_SOURCE_DIR}/include")
  endif()
  set(Boost_LIBRARY_DIR "${GERARDUS_SOURCE_DIR}/lib")
  set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIR}")
  set(Boost_LIBRARY_DIRS "${Boost_LIBRARY_DIR}")
  message(STATUS "Boost include path: ${Boost_INCLUDE_DIRS}")
  message(STATUS "Boost lib path:     ${Boost_LIBRARY_DIRS}")
  
  # integrate the local Boost libraries with the rest of the CMake
  # build system
  find_package(Boost 1.53.0 
    COMPONENTS chrono date_time filesystem system thread)
  
IF(Boost_FOUND)
MESSAGE("Foooooooooooound: ${Boost_LIBRARIES}")
ELSEIF()
MESSAGE("Not Foooooooooooound")
ENDIF()

  # Get Boost version number from its config.hpp file
  set(BOOST_VERSION 0)
  file(READ "${Boost_INCLUDE_DIR}/boost/version.hpp"
    _boost_VERSION_HPP_CONTENTS)
  string(REGEX REPLACE ".*#define BOOST_VERSION ([0-9]+).*" "\\1" 
    Boost_VERSION "${_boost_VERSION_HPP_CONTENTS}")
  set(Boost_VERSION ${Boost_VERSION} CACHE INTERNAL 
    "The version number for boost libraries released together with Gerardus")
  if(NOT "${Boost_VERSION}" STREQUAL "0")
    math(EXPR Boost_MAJOR_VERSION "${Boost_VERSION} / 100000")
    math(EXPR Boost_MINOR_VERSION "${Boost_VERSION} / 100 % 1000")
    math(EXPR Boost_SUBMINOR_VERSION "${Boost_VERSION} % 100")
  endif()
  
endif()

#############################################################################################
## Third-party library: CGAL
#############################################################################################

# note that CGAL depends on Boost, but it's too cumbersome to add
# an explicit dependency

# call the subproject, as CGAL is distributed with a CMAke build
ADD_SUBDIRECTORY(cpp/src/third-party/CGAL-4.2)

# after CGAL builds, we need to copy the header file compiler_config.h
# it generates to the CGAL include directory, so that it's available
# for the other sub-projects that use CGAL (e.g. itk_imfilter)
add_custom_command(
  OUTPUT ${CGAL_SOURCE_DIR}/include/CGAL/compiler_config.h
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/include/CGAL/compiler_config.h" "${CGAL_SOURCE_DIR}/include/CGAL/"
  MAIN_DEPENDENCY ${CMAKE_CURRENT_BINARY_DIR/../../include/CGAL/compiler_config.h}
  DEPENDS CGAL
  )
add_custom_target(copy_compiler_config.h DEPENDS "${CGAL_SOURCE_DIR}/include/CGAL/compiler_config.h")

#############################################################################################
## Configuration options for all Gerardus
#############################################################################################

##################################################################
## Block so that rpaths are added to both the MEX file in the build
## tree and in the installation directory. Otherwise, the rpath is
## stripped when installing, and the function will crash as it won't be
## able to find the ITK shared libraries

# auxiliary libraries built within Gerardus will not be installed in
# the system, but locally. If at some point we need to install in the
# system, the following line will have to be changed, and we'll have
# to uncomment the lines at the end of the rpath block, and use those
# for non-system installs
SET(CMAKE_INSTALL_PREFIX ${GERARDUS_SOURCE_DIR}/lib)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX})

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
#LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
#IF("${isSystemDir}" STREQUAL "-1")
#   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#ENDIF("${isSystemDir}" STREQUAL "-1")

## End rpath block
##################################################################

## flags for all targets
IF(WIN32)

  # to remove warnings from MSVC compiler asking to use fopen_s() an other unportable
  # functions exclusive to Microsoft
  ADD_DEFINITIONS("-D_CRT_SECURE_NO_DEPRECATE")

  # avoid MSVC compiler warnings like
  #    'std::copy': Function call with parameters that may be unsafe
  # given by ITK functions
  ADD_DEFINITIONS("-D_SCL_SECURE_NO_WARNINGS")

ELSE(WIN32)

  # optimise and show all warnings
  ADD_DEFINITIONS("-O2 -Wall")

ENDIF(WIN32)

#############################################################################################
## Call the subprojects 
#############################################################################################

ADD_SUBDIRECTORY(cpp/src)
ADD_SUBDIRECTORY(cpp/src/third-party/mba)
ADD_SUBDIRECTORY(matlab)
