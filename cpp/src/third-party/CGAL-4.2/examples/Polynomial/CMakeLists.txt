# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Polynomial_example )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "coefficient_access.cpp" )
  create_single_source_cgal_program( "construction.cpp" )
  create_single_source_cgal_program( "degree.cpp" )
  create_single_source_cgal_program( "gcd_up_to_constant_factor.cpp" )
  create_single_source_cgal_program( "subresultants.cpp" )
  create_single_source_cgal_program( "substitute.cpp" )
  create_single_source_cgal_program( "swap_move.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

