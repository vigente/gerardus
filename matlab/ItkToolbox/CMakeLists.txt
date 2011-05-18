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

################################################################
## Post-compilation for all targets
################################################################

ADD_DEFINITIONS("-O2 -Wall")
SET(Libraries
  ITKCommon ITKIO
  )

################################################################
## itk_imfilter()
################################################################

ADD_MEX_FILE(itk_imfilter 
  NrrdImage.cpp
  BaseFilter.cpp
  ItkImFilter.cpp)
INCLUDE_DIRECTORIES(
  ../../cpp/src/third-party/itkBinaryThinningImageFilter3D/Source
  )
TARGET_LINK_LIBRARIES(itk_imfilter
  ${Libraries})

################################################################
## itk_kernel_transform
################################################################

ADD_MEX_FILE(itk_kernel_transform ItkKernelTransform.cpp)
TARGET_LINK_LIBRARIES(itk_kernel_transform
  ${Libraries})

################################################################
## Post-compilation for all targets
################################################################

INSTALL(TARGETS 
  itk_kernel_transform
  itk_imfilter
  LIBRARY DESTINATION 
  ${CMAKE_CURRENT_SOURCE_DIR})
