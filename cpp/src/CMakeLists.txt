# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2010-2011 University of Oxford
# Version: 0.3.0
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

PROJECT(GERARDUS_CPP)

################################################################
## configuration of binary programs
################################################################

# add Boost libraries to the include and linking paths of the binary
# programs generated here
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# add TCLAP libraries to the include path
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/third-party/tclap-1.2.1/include")

# extractVoxelCoordinatesFromSegmentationMask
ADD_EXECUTABLE(extractVoxelCoordinatesFromSegmentationMask
  ExtractVoxelCoordinatesFromSegmentationMask.cxx)
TARGET_LINK_LIBRARIES(extractVoxelCoordinatesFromSegmentationMask 
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# padImageWithVoxels
ADD_EXECUTABLE(padSegmentationMaskWithVoxels
  PadSegmentationMaskWithVoxels.cxx)
TARGET_LINK_LIBRARIES(padSegmentationMaskWithVoxels
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# rotate3DImage
ADD_EXECUTABLE(rotate3DImage
  Rotate3DImage.cxx)
TARGET_LINK_LIBRARIES(rotate3DImage
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# resize3DImage
ADD_EXECUTABLE(resize3DImage
  Resize3DImage.cxx)
TARGET_LINK_LIBRARIES(resize3DImage
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# skeletonize3DSegmentation
ADD_EXECUTABLE(skeletonize3DSegmentation
  Skeletonize3DSegmentation.cxx)
INCLUDE_DIRECTORIES(third-party/itkBinaryThinningImageFilter3D/Source)
TARGET_LINK_LIBRARIES(skeletonize3DSegmentation
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# vesselness3DImage
ADD_EXECUTABLE(vesselness3DImage
  Vesselness3DImage.cxx)
TARGET_LINK_LIBRARIES(vesselness3DImage
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# rigidRegistration2D
ADD_EXECUTABLE(rigidRegistration2D
  RigidRegistration2D.cxx)
TARGET_LINK_LIBRARIES(rigidRegistration2D
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

################################################################
## installation of binary programs
################################################################

INSTALL(TARGETS
  extractVoxelCoordinatesFromSegmentationMask
  padSegmentationMaskWithVoxels
  rotate3DImage
  resize3DImage
  skeletonize3DSegmentation
  vesselness3DImage
  rigidRegistration2D
  RUNTIME
  DESTINATION ${GERARDUS_SOURCE_DIR}/programs)
