# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2012 University of Oxford
# Version: 0.1.0
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

################################################################
## flags for all targets
################################################################

IF(NOT WIN32)
  ADD_DEFINITIONS("-O2 -Wall")
ENDIF(NOT WIN32)

################################################################
## sparse_breakdown()
################################################################

ADD_MEX_FILE(sparse_breakdown
  sparse_breakdown.cpp)

################################################################
## mba_surface_interpolation()
################################################################

ADD_MEX_FILE(mba_surface_interpolation
  mba_surface_interpolation.cpp)

INCLUDE_DIRECTORIES(
  ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/mba/include)

LINK_DIRECTORIES(
  ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/mba/lib)

TARGET_LINK_LIBRARIES(mba_surface_interpolation
  MBA)


################################################################
## installation of targets
################################################################

IF(WIN32)
  INSTALL(TARGETS
    sparse_breakdown
    mba_surface_interpolation
    RUNTIME
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
ELSE(WIN32)
  INSTALL(TARGETS
    sparse_breakdown
    mba_surface_interpolation
    LIBRARY
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF(WIN32)
