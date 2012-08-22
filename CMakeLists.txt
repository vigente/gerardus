# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2011 University of Oxford
# Version: 0.1.5
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

# find Matlab
FIND_PACKAGE(Matlab)
IF(MATLAB_FOUND)
  MESSAGE("   Matlab " ${MATLAB_VERSION} ": " ${MATLAB_ROOT})
ENDIF(MATLAB_FOUND)

# some versions of Matlab or some libraries are not compatible with
# more advanced gcc versions
#
# the compiler version has to be set up before project(), otherwise we
# will get an infinite loop
IF(NOT WIN32)
  IF(BUILD_WITH_CUDA OR ${MATLAB_VERSION} MATCHES "R2012a")
    SET(CMAKE_C_COMPILER gcc-4.4)
    SET(CMAKE_CXX_COMPILER g++-4.4)
  ENDIF(BUILD_WITH_CUDA OR ${MATLAB_VERSION} MATCHES "R2012a")
ENDIF(NOT WIN32)

# project name
PROJECT(GERARDUS)

# find ITK
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
  MESSAGE("   ITK v" ${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH})
ENDIF(ITK_FOUND)

# add Boost libraries to the include path
FIND_PACKAGE(Boost COMPONENTS filesystem system REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

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
