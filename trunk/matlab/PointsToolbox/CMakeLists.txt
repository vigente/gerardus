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

ADD_MEX_FILE(itk_kernel_transform ItkKernelTransform.cpp)
SET(Libraries
  ITKCommon ITKIO
)
TARGET_LINK_LIBRARIES(itk_kernel_transform
  ${Libraries})

ADD_DEFINITIONS("-O2 -Wall")

# where to install the MEX file
INSTALL(TARGETS itk_kernel_transform
  LIBRARY DESTINATION 
  ${CMAKE_CURRENT_SOURCE_DIR})
