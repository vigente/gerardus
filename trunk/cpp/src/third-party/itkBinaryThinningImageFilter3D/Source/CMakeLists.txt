INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

PROJECT(BinaryThinning3D)

ENABLE_TESTING()
INCLUDE(Dart)

SET(Required_Packages
  ITK
)

FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

SET(Libraries
  ITKBasicFilters
  ITKCommon
  ITKIO
)

#IF(BUILD_TESTING)
  SET(CurrentExe "BinaryThinning3D")
  ADD_EXECUTABLE(${CurrentExe} BinaryThinning3D.cxx) 
  TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
#ENDIF(BUILD_TESTING)

ADD_TEST(BinaryThinning3DTest BinaryThinning3D ${CMAKE_SOURCE_DIR}/img_in.dcm ${CMAKE_SOURCE_DIR}/img_out.dcm)
