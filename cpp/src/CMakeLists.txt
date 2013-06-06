# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2010-2011 University of Oxford
# Version: 0.2.6
# $Rev$
# $Date$
#
# University of Oxford means the Chancellor, Masters and Scholars of
# the University of Oxford, having an administrative office at
# Wellington Square, Oxford OX1 2JD, UK. 
#
# This file is part of Gerardus.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. The offer of this
# program under the terms of the License is subject to the License
# being interpreted in accordance with English Law and subject to any
# action against the University of Oxford being under the jurisdiction
# of the English Courts.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

# You need to run the following cmake command, with the appropriate
# paths, if you want to generate the project files that allow to
# develop this project from the Eclipse IDE
#
# cmake -DITK_DIR=/usr/local/lib/InsightToolkit  -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release  ./src

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# This project is designed to be built outside the Insight source tree
PROJECT(GERARDUS_CPP)

# avoid MSVC compiler warnings like
#    'std::copy': Function call with parameters that may be unsafe
# given by ITK functions
ADD_DEFINITIONS("-D_SCL_SECURE_NO_WARNINGS")

# add TCLAP libraries to the include path
INCLUDE_DIRECTORIES("${GERARDUS_SOURCE_DIR}/cpp/src/third-party/tclap-1.2.1/include")

# add Boost libraries to the include and linking paths
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})


# extractVoxelCoordinatesFromSegmentationMask
ADD_EXECUTABLE(extractVoxelCoordinatesFromSegmentationMask
  ExtractVoxelCoordinatesFromSegmentationMask.cxx )
TARGET_LINK_LIBRARIES( extractVoxelCoordinatesFromSegmentationMask 
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# padImageWithVoxels
ADD_EXECUTABLE(padSegmentationMaskWithVoxels
  PadSegmentationMaskWithVoxels.cxx)
TARGET_LINK_LIBRARIES( padSegmentationMaskWithVoxels
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# rotate3DImage
ADD_EXECUTABLE(rotate3DImage
  Rotate3DImage.cxx )
TARGET_LINK_LIBRARIES(rotate3DImage
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# resize3DImage
ADD_EXECUTABLE(resize3DImage
  Resize3DImage.cxx)
TARGET_LINK_LIBRARIES(resize3DImage
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# skeletonize3DSegmentation
ADD_EXECUTABLE(skeletonize3DSegmentation
  Skeletonize3DSegmentation.cxx)
INCLUDE_DIRECTORIES(third-party/itkBinaryThinningImageFilter3D/Source)
TARGET_LINK_LIBRARIES(skeletonize3DSegmentation
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# vesselness3DImage
ADD_EXECUTABLE(vesselness3DImage
  Vesselness3DImage.cxx)
TARGET_LINK_LIBRARIES(vesselness3DImage
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

# rigidRegistration2D
ADD_EXECUTABLE(rigidRegistration2D
  RigidRegistration2D.cxx)
TARGET_LINK_LIBRARIES(rigidRegistration2D
  boost_system boost_filesystem
  ${ITK_LIBRARIES})

if(WIN32 AND NOT CYGWIN)

  # for Windows builds, we use the pre-compiled 64-bit libraries
  # provided by the CGAL Windows installer, that we distribute with
  # Gerardus
  set(GMP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third-party/CGAL-4.2/auxiliary/gmp/include)
  set(MPFR_INCLUDE_DIR ${GMP_INCLUDE_DIR})
  set(GMP_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third-party/CGAL-4.2/auxiliary/gmp/lib)
  set(MPFR_LIBRARY_DIR ${GMP_LIBRARY_DIR})

else(WIN32 AND NOT CYGWIN)

  # for Unix-like builds, we have to build the GMP and MPFR libraries from the source code
  set(GMP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third-party/gmp-5.1.2)
  set(MPFR_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third-party/mpfr-3.1.2)

  include(ExternalProject)

  ExternalProject_Add(GMP
    PREFIX ${GMP_SOURCE_DIR}
    STAMP_DIR ${GMP_SOURCE_DIR}/output
    SOURCE_DIR ${GMP_SOURCE_DIR}
    CONFIGURE_COMMAND ./configure --prefix=${GMP_SOURCE_DIR}/output
    BINARY_DIR ${GMP_SOURCE_DIR}
    BUILD_COMMAND make
    TEST_BEFORE_INSTALL 1
    TEST_COMMAND make check
    INSTALL_DIR ${GMP_SOURCE_DIR}/output
    INSTALL_COMMAND make install
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_TEST 1
    LOG_INSTALL 1
    )
  add_library(gmp STATIC IMPORTED GLOBAL)
  set_property(TARGET gmp PROPERTY IMPORTED_LOCATION ${GMP_SOURCE_DIR}/output/lib)
  set(GMP_FOUND TRUE 
    CACHE BOOL "GMP library is part of Gerardus, so this is always true")
  set(GMP_INCLUDE_DIR ${GMP_SOURCE_DIR}
    CACHE PATH "Path to GMP include files")
  set(GMP_LIBRARY_DIR ${GMP_SOURCE_DIR}/output/lib
    CACHE PATH "Path to GMP library")
  set(GMP_LIBRARY_DIRS ${GMP_SOURCE_DIR}/output/lib
    CACHE PATH "Path to GMP library")
  set(GMP_LIBRARIES_DIR ${GMP_SOURCE_DIR}/output/lib
    CACHE PATH "Path to GMP library")
  
  ExternalProject_Add(MPFR
    DEPENDS GMP
    PREFIX ${MPFR_SOURCE_DIR}
    STAMP_DIR ${MPFR_SOURCE_DIR}/output
    SOURCE_DIR ${MPFR_SOURCE_DIR}
    CONFIGURE_COMMAND ./configure --prefix=${MPFR_SOURCE_DIR}/output --with-gmp-include=${GMP_SOURCE_DIR}/output/include --with-gmp-lib=${GMP_SOURCE_DIR}/output/lib
    BINARY_DIR ${MPFR_SOURCE_DIR}
    BUILD_COMMAND make
    TEST_BEFORE_INSTALL 1
    TEST_COMMAND make check
    INSTALL_DIR ${MPFR_SOURCE_DIR}/output
    INSTALL_COMMAND make install
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_TEST 1
    LOG_INSTALL 1
    )
  add_library(mpfr STATIC IMPORTED GLOBAL)
  set_property(TARGET mpfr PROPERTY IMPORTED_LOCATION ${MPFR_SOURCE_DIR}/output/lib)
  set(MPFR_FOUND TRUE 
    CACHE BOOL "MPFR library is part of Gerardus, so this is always true")
  set(MPFR_INCLUDE_DIR ${MPFR_SOURCE_DIR}/src
    CACHE PATH "Path to MPFR include files")
  set(MPFR_LIBRARY_DIR ${MPFR_SOURCE_DIR}/output/lib
    CACHE PATH "Path to MPFR library")
  set(MPFR_LIBRARY_DIRS ${MPFR_SOURCE_DIR}/output/lib
    CACHE PATH "Path to MPFR library")
  set(MPFR_LIBRARIES_DIR ${MPFR_SOURCE_DIR}/output/lib
    CACHE PATH "Path to MPFR library")

endif(WIN32 AND NOT CYGWIN)

# call the subprojects
ADD_SUBDIRECTORY(third-party/boost-1.50.0)

ADD_SUBDIRECTORY(third-party/mba)

ADD_SUBDIRECTORY(third-party/CGAL-4.2)

# after CGAL builds, we need to copy the header file compiler_config.h
# it generates to the CGAL include directory, so that it's available
# for the other sub-projects that use CGAL (e.g. itk_imfilter)
add_custom_command(
  OUTPUT ${CGAL_SOURCE_DIR}/include/CGAL/compiler_config.h
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/../../include/CGAL/compiler_config.h" "${CGAL_SOURCE_DIR}/include/CGAL/"
  MAIN_DEPENDENCY ${CMAKE_CURRENT_BINARY_DIR/../../include/CGAL/compiler_config.h}
  DEPENDS CGAL
  )
add_custom_target(copy_compiler_config.h DEPENDS "${CGAL_SOURCE_DIR}/include/CGAL/compiler_config.h")

################################################################
## installation of targets
################################################################

INSTALL(TARGETS
  extractVoxelCoordinatesFromSegmentationMask
  padSegmentationMaskWithVoxels
  rotate3DImage
  resize3DImage
  skeletonize3DSegmentation
  vesselness3DImage
  rigidRegistration2D
  RUNTIME
  DESTINATION ${GERARDUS_SOURCE_DIR}/programs)
