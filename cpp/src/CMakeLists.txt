# You need to run the following cmake command, with the appropriate
# paths, if you want to generate the project files that allow to
# develop this project from the Eclipse IDE
#
# cmake -DITK_DIR=/usr/local/lib/InsightToolkit  -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release  ./src

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_BUILD_TYPE RELEASE)

# Matlab's mex compiler requires an older version of gcc than the
# default at the time of this writing, so we have to force the version
# we want
SET(CMAKE_CC_COMPILER gcc-4.3)
SET(CMAKE_CXX_COMPILER g++-4.3)

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
boost_system ITKCommon ITKIO ITKNumerics )

# padImageWithVoxels
ADD_EXECUTABLE( padSegmentationMaskWithVoxels
PadSegmentationMaskWithVoxels.cxx )
TARGET_LINK_LIBRARIES( padSegmentationMaskWithVoxels
boost_system ITKCommon ITKIO )

# rotate3DImage
ADD_EXECUTABLE( rotate3DImage
Rotate3DImage.cxx )
TARGET_LINK_LIBRARIES( rotate3DImage
boost_system ITKCommon ITKIO )

# resize3DImage
ADD_EXECUTABLE(resize3DImage
Resize3DImage.cxx)
TARGET_LINK_LIBRARIES(resize3DImage
boost_system ITKCommon ITKIO)

# mbaSurfaceInterpolation
FIND_PACKAGE(OpenGL REQUIRED)
INCLUDE_DIRECTORIES(${OpenGL_INCLUDE_DIR})
ADD_SUBDIRECTORY(third-party/mba)
ADD_EXECUTABLE( mbaSurfaceInterpolation
MbaSurfaceInterpolation.cxx )
INCLUDE_DIRECTORIES(third-party/mba/include)
TARGET_LINK_LIBRARIES(mbaSurfaceInterpolation MBA)
