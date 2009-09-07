# You need to run the following cmake command, with the appropriate
# paths, if you want to generate the project files that allow to
# develop this project from the Eclipse IDE
#
# cmake -DITK_DIR=/usr/local/lib/InsightToolkit  -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release  .

# required by Eclipse CDT4
cmake_minimum_required( VERSION 2.6 )

# This project is designed to be built outside the Insight source tree
PROJECT( GerardusItkToolsProject )

# Find ITK.
FIND_PACKAGE( ITK REQUIRED )
IF( ITK_FOUND )
  INCLUDE( ${ITK_USE_FILE} )
ENDIF( ITK_FOUND )

ADD_EXECUTABLE( extractVoxelCoordinatesFromSegmentationMask
ExtractVoxelCoordinatesFromSegmentationMask.cxx )
TARGET_LINK_LIBRARIES(  extractVoxelCoordinatesFromSegmentationMask 
ITKCommon ITKIO ITKNumerics )
