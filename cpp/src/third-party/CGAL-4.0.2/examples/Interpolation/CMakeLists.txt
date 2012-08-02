# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Interpolation_example )

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

  create_single_source_cgal_program( "interpolation_2.cpp" )
  create_single_source_cgal_program( "linear_interpolation_2.cpp" )
  create_single_source_cgal_program( "nn_coordinates_2.cpp" )
  create_single_source_cgal_program( "nn_coordinates_3.cpp" )
  create_single_source_cgal_program( "rn_coordinates_2.cpp" )
  create_single_source_cgal_program( "sibson_interpolation_2.cpp" )
  create_single_source_cgal_program( "surface_neighbor_coordinates_3.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

