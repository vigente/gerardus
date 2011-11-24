# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2011 University of Oxford
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

# Note: You'll probably have to edit this file to make it compile in
# your system. Edit section "Find Matlab" with your own path to Matlab.

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_BUILD_TYPE RELEASE)

# If Matlab's libraries were built with a version of gcc incompatible
# with the default in your system, you can force the correct version
# uncommenting the two lines below. However, this may cause problems
# with libraries of the wrong version trying to be used
#SET(CMAKE_CC_COMPILER gcc-4.3) 
#SET(CMAKE_CXX_COMPILER g++-4.3)

PROJECT(GERARDUS_MATLAB)

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  MESSAGE(STATUS "ITK found")
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found")
ENDIF(ITK_FOUND)

# Find Matlab
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
FIND_PACKAGE(Matlab REQUIRED)
IF(MATLAB_FOUND)
  MESSAGE(STATUS "Matlab found")
ELSE(MATLAB_FOUND)
  MESSAGE(FATAL_ERROR "Matlab not found")
ENDIF(MATLAB_FOUND)

INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})

# Add macros to build MEX files
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/MatlabMakeMacros.cmake)

# build mex functions in the toolboxes
ADD_SUBDIRECTORY(ItkToolbox)
ADD_SUBDIRECTORY(FiltersToolbox)
