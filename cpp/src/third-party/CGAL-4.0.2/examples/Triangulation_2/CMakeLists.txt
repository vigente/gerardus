# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Triangulation_2_example )

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

  create_single_source_cgal_program( "adding_handles.cpp" )
  create_single_source_cgal_program( "colored_face.cpp" )
  create_single_source_cgal_program( "constrained.cpp" )
  create_single_source_cgal_program( "constrained_hierarchy_plus.cpp" )
  create_single_source_cgal_program( "constrained_plus.cpp" )
  create_single_source_cgal_program( "hierarchy.cpp" )
  create_single_source_cgal_program( "info_insert_with_pair_iterator_2.cpp" )
  create_single_source_cgal_program( "info_insert_with_pair_iterator_regular_2.cpp" )
  create_single_source_cgal_program( "info_insert_with_transform_iterator_2.cpp" )
  create_single_source_cgal_program( "info_insert_with_zip_iterator_2.cpp" )
  create_single_source_cgal_program( "polygon_triangulation.cpp" )
  create_single_source_cgal_program( "print_cropped_voronoi.cpp" )
  create_single_source_cgal_program( "regular.cpp" )
  create_single_source_cgal_program( "terrain.cpp" )
  create_single_source_cgal_program( "triangulation_prog1.cpp" )
  create_single_source_cgal_program( "voronoi.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

