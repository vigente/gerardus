################################################################
## Pre-compilation for all targets
################################################################

ADD_DEFINITIONS("-O2 -Wall")
SET(Libraries
  ITKCommon ITKIO
  )

################################################################
## itk_imfilter()
################################################################

ADD_MEX_FILE(bwregiongrow_aux bwregiongrow_aux.cpp)
TARGET_LINK_LIBRARIES(bwregiongrow_aux
  ${Libraries})

################################################################
## Post-compilation for all targets
################################################################

INSTALL(TARGETS 
  bwregiongrow_aux
  LIBRARY DESTINATION 
  ${CMAKE_CURRENT_SOURCE_DIR})
