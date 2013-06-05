# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Nef_3_example )

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

  create_single_source_cgal_program( "comparison.cpp" )
  create_single_source_cgal_program( "complex_construction.cpp" )
  create_single_source_cgal_program( "exploration_SM.cpp" )
  create_single_source_cgal_program( "extended_kernel.cpp" )
  create_single_source_cgal_program( "handling_double_coordinates.cpp" )
  create_single_source_cgal_program( "interface_polyhedron.cpp" )
  create_single_source_cgal_program( "nefIO.cpp" )
  create_single_source_cgal_program( "nef_3_construction.cpp" )
  create_single_source_cgal_program( "nef_3_point_location.cpp" )
  create_single_source_cgal_program( "nef_3_simple.cpp" )
  create_single_source_cgal_program( "offIO.cpp" )
  create_single_source_cgal_program( "point_set_operations.cpp" )
  create_single_source_cgal_program( "polyline_construction.cpp" )
  create_single_source_cgal_program( "set_operations.cpp" )
  create_single_source_cgal_program( "shell_exploration.cpp" )
  create_single_source_cgal_program( "topological_operations.cpp" )
  create_single_source_cgal_program( "transformation.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

