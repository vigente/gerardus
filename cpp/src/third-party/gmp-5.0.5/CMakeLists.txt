# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2012 University of Oxford
# Version: 0.2.0
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

PROJECT(GMP)

# gerardus directories where GMP gets installed once it gets compiled
set(GMP_INC_DIR ${GMP_SOURCE_DIR}/output/include CACHE PATH 
  "Path to GMP include in Gerardus")
set(GMP_LIB_DIR ${GMP_SOURCE_DIR}/output/lib CACHE PATH 
  "Path to GMP library in Gerardus")

# check whether GMP has already been built and installed
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${GMP_SOURCE_DIR}/../CGAL-4.0.2/cmake/modules/)
find_package(GMP)

if(GMP_FOUND)

#  message(STATUS "GMP include found in ${GMP_INCLUDE_DIR}")
#  message(STATUS "GMP lib found in ${GMP_LIBRARIES_DIR}")

else(GMP_FOUND)

  message(STATUS "Compiling and installing GMP library to ${GMP_SOURCE_DIR}/output")

  # the GMP library is released with an autotool build system. Here,
  # instead of rewriting everything in CMake, we are just going to use
  # CMake as a wrapper of GMP's own build system
  if(UNIX OR APPLE)
    
    execute_process(COMMAND 
      mkdir -p ${GMP_SOURCE_DIR}/output)
    execute_process(COMMAND 
      ${GMP_SOURCE_DIR}/configure --prefix=${GMP_SOURCE_DIR}/output)
    execute_process(COMMAND 
      make -j3)
    execute_process(COMMAND 
      make check -j3)
    execute_process(COMMAND 
      make install)
    
  elseif(UNIX OR APPLE)
    message(SEND_ERROR "Build of GMP not implemented for Windows yet")
  endif(UNIX OR APPLE)
  
endif(GMP_FOUND)
