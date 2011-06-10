################################################################
## Pre-compilation for all targets
################################################################

ADD_DEFINITIONS("-O2 -Wall")
SET(Libraries
  ITKCommon ITKIO
  )

################################################################
## bwregiongrow()
################################################################

ADD_MEX_FILE(bwregiongrow bwregiongrow.cpp)
TARGET_LINK_LIBRARIES(bwregiongrow
  ${Libraries})

################################################################
## Post-compilation for all targets
################################################################

INSTALL(TARGETS 
  bwregiongrow
  LIBRARY 
  DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
  ARCHIVE
  DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
