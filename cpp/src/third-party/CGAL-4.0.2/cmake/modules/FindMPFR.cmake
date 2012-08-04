# This file is provided by the CGAL project.
# Modifications by Ramon Casero  <ramon.casero@cs.ox.ac.uk> to
# integrate CGAL within project Gerardus.
#
# Version: 0.1.0
#
# - Fix bug. Variable names had not been enclosed in ${}, and where
# being ignored.

# Try to find the MPFR libraries
# MPFR_FOUND - system has MPFR lib
# MPFR_INCLUDE_DIR - the MPFR include directory
# MPFR_LIBRARIES_DIR - Directory where the MPFR libraries are located
# MPFR_LIBRARIES - the MPFR libraries
# MPFR_IN_CGAL_AUXILIARY - TRUE if the MPFR found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(CGAL_FindPackageHandleStandardArgs)
include(CGAL_GeneratorSpecificSettings)

if(MPFR_INCLUDE_DIR)
  set(MPFR_in_cache TRUE)
else()
  set(MPFR_in_cache FALSE)
endif()
if(NOT MPFR_LIBRARIES)
  set(MPFR_in_cache FALSE)
endif()

# Is it already configured?
if (MPFR_in_cache) 
   
  set(MPFR_FOUND TRUE)
  
else()  

  MESSAGE("FOO: " ${MPFR_INC_DIR})
  find_path(MPFR_INCLUDE_DIR 
    NAMES mpfr.h 
    PATHS ${ENV} ${MPFR_INC_DIR}
    ${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/include
    DOC "The directory containing the MPFR header files"
    )
  MESSAGE("Hola: " ${MPFR_INCLUDE_DIR})

  if ( MPFR_INCLUDE_DIR STREQUAL "${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/include" )
    cache_set( MPFR_IN_CGAL_AUXILIARY TRUE )
  endif()
  
  find_library(MPFR_LIBRARIES NAMES mpfr libmpfr-4 libmpfr-1
    PATHS ${ENV} ${MPFR_LIB_DIR}
    ${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/lib
    DOC "Path to the MPFR library"
    )
  
  if ( MPFR_LIBRARIES ) 
    get_filename_component(MPFR_LIBRARIES_DIR ${MPFR_LIBRARIES} PATH CACHE )
  endif()
  
  # Attempt to load a user-defined configuration for MPFR if couldn't be found
  if ( NOT MPFR_INCLUDE_DIR OR NOT MPFR_LIBRARIES_DIR )
    include( MPFRConfig OPTIONAL )
  endif()

  find_package_handle_standard_args(MPFR "DEFAULT_MSG" MPFR_LIBRARIES MPFR_INCLUDE_DIR)

endif()

