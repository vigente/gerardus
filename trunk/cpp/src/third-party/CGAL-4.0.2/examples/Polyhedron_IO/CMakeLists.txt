# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Polyhedron_IO_example )

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

  create_single_source_cgal_program( "iv2off.cpp" )
  create_single_source_cgal_program( "off2iv.cpp" )
  create_single_source_cgal_program( "off2off.cpp" )
  create_single_source_cgal_program( "off2stl.cpp" )
  create_single_source_cgal_program( "off2vrml.cpp" )
  create_single_source_cgal_program( "off2wav.cpp" )
  create_single_source_cgal_program( "off_bbox.cpp" )
  create_single_source_cgal_program( "off_glue.cpp" )
  create_single_source_cgal_program( "off_transform.cpp" )
  create_single_source_cgal_program( "polyhedron2vrml.cpp" )
  create_single_source_cgal_program( "polyhedron_copy.cpp" )
  create_single_source_cgal_program( "terr_trian.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

