# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2012 University of Oxford
# Version: 0.1.3
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

# MBA library

#FIND_PACKAGE(OpenGL REQUIRED)
IF(WIN32)
  ADD_LIBRARY(MBA SHARED STATIC
    src/filterRestrict2D.cpp
    src/MBA.cpp
    src/MBAdata.cpp
    src/PointAccessUtils.cpp
    src/UCBsplines.cpp
    src/UCBsplineSurface.cpp
    src/UCButils.cpp)
ELSE(WIN32)
  ADD_LIBRARY(MBA SHARED
    src/filterRestrict2D.cpp
    src/MBA.cpp
    src/MBAdata.cpp
    src/PointAccessUtils.cpp
    src/UCBsplines.cpp
    src/UCBsplineSurface.cpp
    src/UCButils.cpp)
ENDIF(WIN32)

INCLUDE_DIRECTORIES(include "${Boost_INCLUDE_DIRS}")
#TARGET_LINK_LIBRARIES(MBA glut GLU GL)

INSTALL(TARGETS MBA 
  DESTINATION 
  ${GERARDUS_SOURCE_DIR}/lib)
