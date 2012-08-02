# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Straight_skeleton_2_example )

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

  create_single_source_cgal_program( "Create_offset_polygons_2.cpp" )
  create_single_source_cgal_program( "Create_saop_from_polygon_with_holes_2.cpp" )
  create_single_source_cgal_program( "Create_skeleton_and_offset_polygons_2.cpp" )
  create_single_source_cgal_program( "Create_skeleton_and_offset_polygons_with_holes_2.cpp" )
  create_single_source_cgal_program( "Create_straight_skeleton_2.cpp" )
  create_single_source_cgal_program( "Create_straight_skeleton_from_polygon_with_holes_2.cpp" )
  create_single_source_cgal_program( "Low_level_API.cpp" )
  create_single_source_cgal_program( "Show_offset_polygon.cpp" )
  create_single_source_cgal_program( "Show_straight_skeleton.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

