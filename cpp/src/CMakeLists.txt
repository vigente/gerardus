# You need to run the following cmake command, with the appropriate
# paths, if you want to generate the project files that allow to
# develop this project from the Eclipse IDE
#
# cmake -DITK_DIR=/usr/local/lib/InsightToolkit  -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release  ./src

# required by Eclipse CDT4
cmake_minimum_required( VERSION 2.6 )

# This project is designed to be built outside the Insight source tree
PROJECT( GerardusItkToolsProject )

# Find ITK.
FIND_PACKAGE( ITK REQUIRED )
IF( ITK_FOUND )
  INCLUDE( ${ITK_USE_FILE} )
ENDIF( ITK_FOUND )

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
ADD_EXECUTABLE( resize3DImage
Resize3DImage.cxx )
TARGET_LINK_LIBRARIES( resize3DImage
boost_system ITKCommon ITKIO )

# mbaSurfaceInterpolation
#ADD_EXECUTABLE( mbaSurfaceInterpolation
#MbaSurfaceInterpolation.cxx )
#INCLUDE_DIRECTORIES(./third-party/mba/include)
#LINK_DIRECTORIES(./third-party/mba/lib)
#TARGET_LINK_LIBRARIES(mbaSurfaceInterpolation MBA)
