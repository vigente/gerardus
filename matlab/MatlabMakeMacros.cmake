# This file was found in the source code of Dru F., Fillard P.,
# Vercauteren T., "An ITK Implementation of the Symmetric Log-Domain
# Diffeomorphic Demons Algorithm", Insight Journal, 2009 Jan-Jun
# http://hdl.handle.net/10380/3060

MACRO(LOAD_REQUIRED_PACKAGE Package)
  LOADPACKAGE(${Package})
  IF(NOT ${Package}_FOUND)
    MESSAGE(FATAL_ERROR "Required package ${Package} was not found.\n
    Look at Find${Package}.cmake in the CMake module directory for clues
    on what you're supposed to do to help find this package.  Good luck.\n")
  ENDIF(NOT ${Package}_FOUND)
ENDMACRO(LOAD_REQUIRED_PACKAGE)

MACRO(LOAD_OPTIONAL_PACKAGE Package)
  LOADPACKAGE(${Package} QUIET)
ENDMACRO(LOAD_OPTIONAL_PACKAGE)

MACRO(ADD_MEX_FILE Target)
  INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})
  ADD_LIBRARY(${Target} SHARED ${ARGN})
  IF(WIN32)
    TARGET_LINK_LIBRARIES(${Target} 
      ${MATLAB_MX_LIBRARY} 
      ${MATLAB_MEX_LIBRARY} 
      ${MATLAB_MAT_LIBRARY} 
      )
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES(${Target} 
      ${MATLAB_MX_LIBRARY} 
      ${MATLAB_MEX_LIBRARY} 
      ${MATLAB_MAT_LIBRARY} 
      m
      )
  ENDIF(WIN32)
  
  SET_TARGET_PROPERTIES(${Target} PROPERTIES PREFIX "")
  
  # Determine mex suffix
  IF(UNIX)
    # if this is OSX (which is UNIX) then the library suffixes depend on the architecture
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      IF(CMAKE_OSX_ARCHITECTURES MATCHES i386)
        # mac intel 32-bit
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexmaci")
      ELSEIF(CMAKE_OSX_ARCHITECTURES MATCHES x86_64)
        # mac intel 64-bit
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexmaci64")
      ELSE(CMAKE_OSX_ARCHITECTURES MATCHES i386)
        # Mac Power PC
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexmac")
      ENDIF(CMAKE_OSX_ARCHITECTURES MATCHES i386)
    ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexglx")
      ELSEIF(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexa64")
      ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
        MESSAGE(FATAL_ERROR 
          "CMAKE_SIZEOF_VOID_P (${CMAKE_SIZEOF_VOID_P}) doesn't indicate a valid platform")
      ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  ELSEIF(WIN32)
    IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
      SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexw32")
    ELSEIF(CMAKE_SIZEOF_VOID_P MATCHES "8")
      SET_TARGET_PROPERTIES(${Target} PROPERTIES SUFFIX ".mexw64")
    ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
      MESSAGE(FATAL_ERROR 
        "CMAKE_SIZEOF_VOID_P (${CMAKE_SIZEOF_VOID_P}) doesn't indicate a valid platform")
    ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
  ENDIF(UNIX)
  
  IF(MSVC)
    SET(MATLAB_FLAGS "-DMATLAB_MEX_FILE")
    SD_APPEND_TARGET_PROPERTIES(${Target} COMPILE_FLAGS ${MATLAB_FLAGS})
    SET_TARGET_PROPERTIES(${Target} PROPERTIES LINK_FLAGS "/export:mexFunction")
  ELSE(MSVC)
    IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
      SET(MATLAB_FLAGS "-fPIC" "-D_GNU_SOURCE" "-pthread"
	"-D_FILE_OFFSET_BITS=64" "-DMX_COMPAT_32")
    ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
      SET(MATLAB_FLAGS "-fPIC" "-D_GNU_SOURCE" "-pthread"
	"-D_FILE_OFFSET_BITS=64")
    ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
    SD_APPEND_TARGET_PROPERTIES(${Target} COMPILE_FLAGS ${MATLAB_FLAGS})
    
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      IF(CMAKE_OSX_ARCHITECTURES MATCHES i386)
        # mac intel 32-bit
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-L${MATLAB_ROOT}/bin/maci -Wl,-flat_namespace -undefined suppress")
      ELSEIF(CMAKE_OSX_ARCHITECTURES MATCHES x86_64)
        # mac intel 64-bit
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-L${MATLAB_ROOT}/bin/maci64 -Wl,-flat_namespace -undefined suppress")
      ELSE(CMAKE_OSX_ARCHITECTURES MATCHES i386)
        # mac powerpc?
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-L${MATLAB_SYS} -Wl,-flat_namespace -undefined suppress")
      ENDIF(CMAKE_OSX_ARCHITECTURES MATCHES i386)
    ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-Wl,-E -Wl,--no-undefined")
      ELSEIF(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${Target} PROPERTIES 
          LINK_FLAGS "-Wl,-E -Wl,--no-undefined")
      ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
        MESSAGE(FATAL_ERROR 
          "CMAKE_SIZEOF_VOID_P (${CMAKE_SIZEOF_VOID_P}) doesn't indicate a valid platform")
      ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
    ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  ENDIF(MSVC)
ENDMACRO(ADD_MEX_FILE)


MACRO(SD_APPEND_TARGET_PROPERTIES TARGET_TO_CHANGE PROP_TO_CHANGE)
  FOREACH(_newProp ${ARGN})
    GET_TARGET_PROPERTY(_oldProps ${TARGET_TO_CHANGE} ${PROP_TO_CHANGE})
    IF(_oldProps)
      IF(NOT "${_oldProps}" MATCHES "^.*${_newProp}.*$")
        SET_TARGET_PROPERTIES(${TARGET_TO_CHANGE} PROPERTIES ${PROP_TO_CHANGE} "${_newProp} ${_oldProps}")
      ENDIF(NOT "${_oldProps}" MATCHES "^.*${_newProp}.*$")
    ELSE(_oldProps)
      SET_TARGET_PROPERTIES(${TARGET_TO_CHANGE} PROPERTIES ${PROP_TO_CHANGE} ${_newProp})
    ENDIF(_oldProps)
  ENDFOREACH(_newProp ${ARGN})
ENDMACRO(SD_APPEND_TARGET_PROPERTIES TARGET_TO_CHANGE PROP_TO_CHANGE)


MACRO(SD_ADD_LINK_LIBRARIES Target)
  FOREACH (currentLib ${ARGN})
    IF (${currentLib}_LIBRARIES)
      TARGET_LINK_LIBRARIES(${Target} ${${currentLib}_LIBRARIES})
    ELSEIF (${currentLib}_LIBRARY)
      TARGET_LINK_LIBRARIES(${Target} ${${currentLib}_LIBRARY})
    ELSE (${currentLib}_LIBRARIES)
      #MESSAGE("WARNING: ${currentLib}_LIBRARY and ${currentLib}_LIBRARIES are undefined. Using ${currentLib} in linker")
      TARGET_LINK_LIBRARIES(${Target} ${currentLib})
    ENDIF (${currentLib}_LIBRARIES)
    
    IF (${currentLib}_INCLUDE_DIRS)
      INCLUDE_DIRECTORIES(${${currentLib}_INCLUDE_DIRS})
    ELSEIF (${currentLib}_INCLUDE_DIR)
      INCLUDE_DIRECTORIES(${${currentLib}_INCLUDE_DIR})
    ELSE (${currentLib}_INCLUDE_DIRS)
      #MESSAGE("WARNING: ${currentLib}_INCLUDE_DIR and ${currentLib}_INCLUDE_DIR are undefined. No specific include dir will be used for ${currentLib}")
    ENDIF (${currentLib}_INCLUDE_DIRS)
  ENDFOREACH (currentLib)
ENDMACRO(SD_ADD_LINK_LIBRARIES)


MACRO(SD_UNIT_TEST Src)
  # remove extension
  STRING(REGEX REPLACE "[.][^.]*$" "" Target ${Src})

  # parse arguments
  SET(currentPos "")
  SET(testLibs "")
  SET(testExtLibs "")
  SET(testArgs "")

  FOREACH (arg ${ARGN})
    IF (arg STREQUAL "LIBS")
      SET(currentPos "LIBS")
    ELSEIF (arg STREQUAL "EXTLIBS")
      SET(currentPos "EXTLIBS")
    ELSEIF (arg STREQUAL "ARGS")
      SET(currentPos "ARGS")
    ELSE (arg STREQUAL "LIBS")
      IF (currentPos STREQUAL "LIBS")
        SET(testLibs ${testLibs} ${arg})
      ELSEIF (currentPos STREQUAL "EXTLIBS")
        SET(testExtLibs ${testExtLibs} ${arg})
      ELSEIF (currentPos STREQUAL "ARGS")
        SET(testArgs ${testArgs} ${arg})
      ELSE (currentPos STREQUAL "ARGS")
         MESSAGE(FATAL_ERROR "Unknown argument")
      ENDIF (currentPos STREQUAL "LIBS")
    ENDIF (arg STREQUAL "LIBS")
  ENDFOREACH (arg ${ARGN})

  # setup target
  ADD_EXECUTABLE(${Target} ${Src})
  SD_ADD_LINK_LIBRARIES(${Target} ${testExtLibs})
  TARGET_LINK_LIBRARIES(${Target} ${testLibs})
  SET_TARGET_PROPERTIES(${Target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
  ADD_TEST(${Target} ${PROJECT_BINARY_DIR}/tests/${Target} ${testArgs})
ENDMACRO(SD_UNIT_TEST)


MACRO(SD_EXECUTABLE Src)
  # remove extension
  STRING(REGEX REPLACE "[.][^.]*$" "" Target ${Src})

  # parse arguments
  SET(currentPos "")
  SET(appLibs "")
  SET(appExtLibs "")

  FOREACH (arg ${ARGN})
    IF (arg STREQUAL "LIBS")
      SET(currentPos "LIBS")
    ELSEIF (arg STREQUAL "EXTLIBS")
      SET(currentPos "EXTLIBS")
    ELSE (arg STREQUAL "LIBS")
      IF (currentPos STREQUAL "LIBS")
        SET(appLibs ${appLibs} ${arg})
      ELSEIF (currentPos STREQUAL "EXTLIBS")
        SET(appExtLibs ${appExtLibs} ${arg})
      ELSE (currentPos STREQUAL "LIBS")
         MESSAGE(FATAL_ERROR "Unknown argument")
      ENDIF (currentPos STREQUAL "LIBS")
    ENDIF (arg STREQUAL "LIBS")
  ENDFOREACH (arg ${ARGN})

  # setup target
  ADD_EXECUTABLE(${Target} ${Src})
  SD_ADD_LINK_LIBRARIES(${Target} ${appExtLibs})
  TARGET_LINK_LIBRARIES(${Target} ${appLibs})
  SET_TARGET_PROPERTIES(${Target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
ENDMACRO(SD_EXECUTABLE)
