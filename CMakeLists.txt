# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2011 University of Oxford
# Version: 0.2.2
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# for Visual Studio or other multi-configuration systems, make available only the Debug and Release configurations
# this needs to go before PROJECT() or it will be ignored
SET(CMAKE_CONFIGURATION_TYPES Release Debug)
SET(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
  "Reset the configurations to Debug/Release only"
  FORCE)

# for single configuration systems, by default, choose the Release configuration
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# directories with ad hoc Findxxx.cmake files
SET(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/matlab;${CMAKE_CURRENT_SOURCE_DIR}/cpp/src/third-party/cquammen-Clarity/CMake)

# display type of build (e.g. Release)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

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

# project name
project(GERARDUS)

# find ITK
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
  MESSAGE(STATUS "ITK v" ${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH})
ENDIF(ITK_FOUND)

# route to find Boost header files and libraries
SET(Boost_INCLUDE_DIR ${GERARDUS_SOURCE_DIR}/include 
  CACHE PATH "Path to Boost header files")
SET(Boost_LIBRARY_DIR ${GERARDUS_SOURCE_DIR}/lib
  CACHE PATH "Path to Boost libraries")
SET(Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIR}
  CACHE PATH "Path to Boost header files")
SET(Boost_LIBRARY_DIRS ${Boost_LIBRARY_DIR}
  CACHE PATH "Path to Boost libraries")
SET(Boost_LIBRARIES 
  boost_system 
  boost_filesystem 
  boost_thread
  boost_chrono
  boost_date_time
  CACHE STRING "Available Boost libraries")
SET(Boost_DEFINITIONS ""
  CACHE STRING "Definitions for Boost libraries")

# Get Boost version number from its config.hpp file
set(BOOST_VERSION 0)
file(READ "${Boost_INCLUDE_DIR}/boost/version.hpp" _boost_VERSION_HPP_CONTENTS)
string(REGEX REPLACE ".*#define BOOST_VERSION ([0-9]+).*" "\\1" 
  Boost_VERSION "${_boost_VERSION_HPP_CONTENTS}")
set(Boost_VERSION ${Boost_VERSION} CACHE INTERNAL 
  "The version number for boost libraries")
if(NOT "${Boost_VERSION}" STREQUAL "0")
  math(EXPR Boost_MAJOR_VERSION "${Boost_VERSION} / 100000")
  math(EXPR Boost_MINOR_VERSION "${Boost_VERSION} / 100 % 1000")
  math(EXPR Boost_SUBMINOR_VERSION "${Boost_VERSION} % 100")
endif()

#message(STATUS "Boost include:     ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries:   ${Boost_LIBRARIES}")
#message(STATUS "Boost definitions: ${Boost_DEFINITIONS}")
message(STATUS "Boost version: '${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}'")

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

# call the subprojects
ADD_SUBDIRECTORY(cpp/src)
ADD_SUBDIRECTORY(matlab)
