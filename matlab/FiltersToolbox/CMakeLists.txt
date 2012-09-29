# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2011 University of Oxford
# Version: 0.2.1
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
## Pre-compilation for all targets
################################################################

IF(NOT WIN32)
  ADD_DEFINITIONS("-O2 -Wall")
ENDIF(NOT WIN32)

################################################################
## bwregiongrow()
################################################################

ADD_MEX_FILE(bwregiongrow bwregiongrow.cpp)
INCLUDE_DIRECTORIES(..)

################################################################
## deconvolve()
## This function has been removed for three reasons:
## 1) the Clarity library requires the FFTW library, which is
##    another dependency to install
## 2) Matlab has its own FFTW library, and I haven't found a
##    way to make both coexist in Gerardus
## 3) Matlab has very good deconvolution functions
################################################################

#LINK_DIRECTORIES(${CMAKE_INSTALL_PREFIX})
#INCLUDE_DIRECTORIES(${Clarity_BINARY_DIR})
#ADD_MEX_FILE(deconvolve Deconvolve.cpp)
#ADD_DEPENDENCIES(deconvolve Clarity)
#TARGET_LINK_LIBRARIES(deconvolve
#  ${ITK_LIBRARIES} Clarity)

################################################################
## Post-compilation for all targets
################################################################

IF(WIN32)
  INSTALL(TARGETS
    bwregiongrow
#    deconvolve
    RUNTIME
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
ELSE(WIN32)
  INSTALL(TARGETS
    bwregiongrow
#    deconvolve
    LIBRARY
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF(WIN32)
