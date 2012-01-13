# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2010-2011 University of Oxford
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

# You need to run the following cmake command, with the appropriate
# paths, if you want to generate the project files that allow to
# develop this project from the Eclipse IDE
#
# cmake -DITK_DIR=/usr/local/lib/InsightToolkit  -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release  ./src

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# This project is designed to be built outside the Insight source tree
PROJECT(GerardusItkToolsProject)

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

# extractVoxelCoordinatesFromSegmentationMask
ADD_EXECUTABLE( extractVoxelCoordinatesFromSegmentationMask
  ExtractVoxelCoordinatesFromSegmentationMask.cxx )
TARGET_LINK_LIBRARIES( extractVoxelCoordinatesFromSegmentationMask 
  boost_system boost_filesystem ITKCommon ITKIO ITKNumerics )

# padImageWithVoxels
ADD_EXECUTABLE( padSegmentationMaskWithVoxels
  PadSegmentationMaskWithVoxels.cxx )
TARGET_LINK_LIBRARIES( padSegmentationMaskWithVoxels
  boost_system boost_filesystem ITKCommon ITKIO )

# rotate3DImage
ADD_EXECUTABLE( rotate3DImage
  Rotate3DImage.cxx )
TARGET_LINK_LIBRARIES( rotate3DImage
  boost_system boost_filesystem ITKCommon ITKIO )

# resize3DImage
ADD_EXECUTABLE(resize3DImage
  Resize3DImage.cxx)
TARGET_LINK_LIBRARIES(resize3DImage
  boost_system boost_filesystem ITKCommon ITKIO)

# skeletonize3DSegmentation
ADD_EXECUTABLE(skeletonize3DSegmentation
  Skeletonize3DSegmentation.cxx)
INCLUDE_DIRECTORIES(third-party/itkBinaryThinningImageFilter3D/Source)
TARGET_LINK_LIBRARIES(skeletonize3DSegmentation
  boost_system boost_filesystem ITKCommon ITKIO)

# vesselness3DImage
ADD_EXECUTABLE(vesselness3DImage
  Vesselness3DImage.cxx)
TARGET_LINK_LIBRARIES(vesselness3DImage
  boost_system boost_filesystem ITKCommon ITKIO)

# vesselness3DImage
ADD_EXECUTABLE(rigidRegistration2D
  RigidRegistration2D.cxx)
TARGET_LINK_LIBRARIES(rigidRegistration2D
  boost_system boost_filesystem ITKCommon ITKIO ITKNumerics)

# call the subprojects
ADD_SUBDIRECTORY(third-party/cquammen-Clarity)
