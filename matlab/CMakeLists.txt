# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2011 University of Oxford
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

# Note: You'll probably have to edit this file to make it compile in
# your system. Edit section "Find Matlab" with your own path to Matlab.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_BUILD_TYPE RELEASE)

# Matlab's mex compiler requires an older version of gcc than the
# default at the time of this writing, so we have to force the version
# we want
SET(CMAKE_CC_COMPILER gcc-4.3)
SET(CMAKE_CXX_COMPILER g++-4.3)

PROJECT(MATLAB_CMAKE)

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  MESSAGE(STATUS "ITK found")
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found")
ENDIF(ITK_FOUND)

# Find Matlab
SET(MATLAB_ROOT "/usr/local/matlab/R2010b/")
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
FIND_PACKAGE(Matlab REQUIRED)
IF(MATLAB_FOUND)
  MESSAGE(STATUS "Matlab found")
ELSE(MATLAB_FOUND)
  MESSAGE(FATAL_ERROR "Matlab not found")
ENDIF(MATLAB_FOUND)

# SET(MATLAB_ROOT
#   /usr/local/matlab/R2010b/
#   )
# FIND_LIBRARY(MATLAB_MEX_LIBRARY
#   mex
#   ${MATLAB_ROOT}/bin/glnxa64
#   )
# FIND_LIBRARY(MATLAB_MX_LIBRARY
#   mx
#   ${MATLAB_ROOT}/bin/glnxa64
#   )
# FIND_PATH(MATLAB_INCLUDE_DIR
#   "mex.h"
#   "${MATLAB_ROOT}/extern/include/"
#   )

INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})

# Add macros to build MEX files
INCLUDE(${CMAKE_SOURCE_DIR}/MatlabMakeMacros.cmake)

##################################################################
## Block so that rpaths are added to both the MEX file in the build
## tree and in the installation directory. Otherwise, the rpath is
## stripped when installing, and the function will crash as it won't be
## able to find the ITK shared libraries

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

## End rpath block
##################################################################

# build mex functions in the Points Toolbox
ADD_SUBDIRECTORY(PointsToolbox)
