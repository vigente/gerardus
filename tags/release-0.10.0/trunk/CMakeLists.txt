# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2011-2014 University of Oxford
# Version: 0.9.1
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

#############################################################################################
## Preliminary configuration, prior to defining the Gerardus project
#############################################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# for Visual Studio or other multi-configuration systems, make
# available only the Debug and Release configurations. This needs to
# go before PROJECT() or it will be ignored
SET(CMAKE_CONFIGURATION_TYPES Release Debug)
SET(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
  "Reset the configurations to Debug/Release only"
  FORCE)

# for single configuration systems (Unix-like), by default, choose the
# Release configuration
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# display type of build (e.g. Release)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# make directory with FindMatlab.cmake module visible to cmake. This
# makes FindMatlab.cmake available to be used in the next line to find
# and configure Matlab
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/matlab)

#############################################################################################
## Matlab configuration
#############################################################################################

# hack to avoid error messages from CMake in Windows. The problem is that we must run 
# project() after we have set up the compiler version, to avoid an infinite loop in cmake.
# And we can only set up the compiler version once we know Matlab's version. Additionally,
# FindMatlab.cmake requires some variables to be set (CMAKE_FIND_LIBRARY_PREFIXES, 
# CMAKE_FIND_LIBRARY_SUFFIXES), which only happens after project(). Thus, we have a Catch22.
# This hack is simply to set the offending variables to the values they would have in Windows
if(WIN32)
   set(CMAKE_FIND_LIBRARY_PREFIXES "")
   set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
endif(WIN32)

# find Matlab
find_package(Matlab REQUIRED)
if(MATLAB_FOUND)
  message(STATUS "Matlab " ${MATLAB_VERSION} ": " ${MATLAB_ROOT})
else()
  message(FATAL_ERROR "Matlab not found")
endif(MATLAB_FOUND)

# some versions of Matlab or some libraries are not compatible with
# more advanced gcc versions
#
# the compiler version has to be set up before project(), otherwise we
# will get an infinite loop. We force the compiler version setting the
# corresponding environmental variables (equivalent to running e.g.  
# $ CC=gcc-4.4 CXX=g++-4.4 cmake ..  This is a better way than directly
# setting MAKE_C_COMPILER, CMAKE_CXX_COMPILER. For example, if we set these
# variables after project(), this creates an infinite loop.
if(NOT WIN32)

  if("${MATLAB_VERSION}" MATCHES "R2012a")
    set(GCC_VERSION_REQUIRED_BY_MATLAB "4.4")
  elseif("${MATLAB_VERSION}" MATCHES "R2012b")
    set(GCC_VERSION_REQUIRED_BY_MATLAB "4.4")
  else()
    message(FATAL_ERROR 
      "I do not know which version of gcc/g++ is needed for Matlab ${MATLAB_VERSION}.\nPlease update gerardus/CMakeLists.txt with this information.")
  endif()
#  message(STATUS "Matlab ${MATLAB_VERSION} requires gcc, g++ and gfortran version: ${GCC_VERSION_REQUIRED_BY_MATLAB}")
  message(STATUS "Matlab ${MATLAB_VERSION} requires gcc and g++ version: ${GCC_VERSION_REQUIRED_BY_MATLAB}")

################# gcc version

  # if the user has provided environmental variables with the
  # compilers, we are going to check that they have the correct
  # version
  if(DEFINED ENV{CC})
    message(STATUS "User wants to use as gcc the following binary: $ENV{CC}")
    
    # get version of the C compiler
    exec_program(
      $ENV{CC}
      ARGS                    --version
      OUTPUT_VARIABLE _compiler_output)
    string(REGEX MATCH "[0-9]\\.[0-9]\\.[0-9]"
      gcc_compiler_version ${_compiler_output})
    message(STATUS "This binary corresponds to version: ${gcc_compiler_version}")
    string(REGEX MATCH "[0-9]\\.[0-9]"
      gcc_compiler_version ${gcc_compiler_version})
    
    if(NOT gcc_compiler_version VERSION_EQUAL "${GCC_VERSION_REQUIRED_BY_MATLAB}")
      message(FATAL_ERROR "Matlab ${MATLAB_VERSION} requires gcc ${GCC_VERSION_REQUIRED_BY_MATLAB}")
      
    else()
      message(STATUS "This compiler version is valid for Matlab ${MATLAB_VERSION}")
    endif()
    
  else()
    # if the user does not provide an environmental variable with
    # the compiler, we set it up
    set(ENV{CC} "gcc-${GCC_VERSION_REQUIRED_BY_MATLAB}")
  endif()
  
################# g++ version

  if(DEFINED ENV{CXX})
    message(STATUS "User wants to use as g++ the following binary: $ENV{CXX}")
    
    # get version of the C++ compiler
    exec_program(
      $ENV{CXX}
      ARGS                    --version
      OUTPUT_VARIABLE _compiler_output)
    
    string(REGEX MATCH "[0-9]\\.[0-9]\\.[0-9]"
      gxx_compiler_version ${_compiler_output})
    message(STATUS "This binary corresponds to version: ${gxx_compiler_version}")
    string(REGEX MATCH "[0-9]\\.[0-9]"
      gxx_compiler_version ${gxx_compiler_version})
    
    if(NOT gxx_compiler_version VERSION_EQUAL "${GCC_VERSION_REQUIRED_BY_MATLAB}")
      message(FATAL_ERROR "Matlab ${MATLAB_VERSION} requires g++ ${GCC_VERSION_REQUIRED_BY_MATLAB}")
    else()
      message(STATUS "This compiler version is valid for Matlab ${MATLAB_VERSION}")
    endif()

  else()
    # if the user does not provide an environmental variable with
    # the compiler, we set it up
    set(ENV{CXX} "g++-${GCC_VERSION_REQUIRED_BY_MATLAB}")
  endif()

################# gfortran version

  # if(DEFINED ENV{F77})
  #   message(STATUS "User wants to use as gfortran the following binary: $ENV{F77}")
    
  #   # get version of the Fortran compiler
  #   exec_program(
  #     $ENV{F77}
  #     ARGS                    --version
  #     OUTPUT_VARIABLE _compiler_output)
    
  #   string(REGEX MATCH "[0-9]\\.[0-9]\\.[0-9]"
  #     gfortran_compiler_version ${_compiler_output})
  #   message(STATUS "This binary corresponds to version: ${gfortran_compiler_version}")
  #   string(REGEX MATCH "[0-9]\\.[0-9]"
  #     gfortran_compiler_version ${gfortran_compiler_version})
    
  #   if(NOT gfortran_compiler_version VERSION_EQUAL "${GCC_VERSION_REQUIRED_BY_MATLAB}")
  #     message(FATAL_ERROR "Matlab ${MATLAB_VERSION} requires gfortran ${GCC_VERSION_REQUIRED_BY_MATLAB}")
  #   else()
  #     message(STATUS "This compiler version is valid for Matlab ${MATLAB_VERSION}")
  #   endif()

  # else()
  #   # if the user does not provide an environmental variable with
  #   # the compiler, we set it up
  #   set(ENV{F77} "gfortran-${GCC_VERSION_REQUIRED_BY_MATLAB}")
  # endif()

endif(NOT WIN32)

#############################################################################################
## Start of Gerardus project
#############################################################################################

# project name
project(GERARDUS)

# CMake doesn't work for Windows when Gerardus is on a network share (e.g. cmake -E untar fails 
# without a reason, running bootstrap.bat fails too...)
if (WIN32)
  # get the two first characters of the path to Gerardus
  string(SUBSTRING "${GERARDUS_SOURCE_DIR}" 0 2 AUX)
  if (AUX STREQUAL "//")
    message(FATAL_ERROR "Gerardus cannot be built on a Windows share, e.g. \\\\my_share\\gerardus. It needs to be on C:, D:, etc. Otherwise, CMake doesn't work well.")
  endif()
endif()

#############################################################################################
## Third-party libraries: GMP and MRPF
#############################################################################################

if(WIN32 AND NOT CYGWIN)

  # for Windows builds, we use the pre-compiled 64-bit libraries
  # provided by the CGAL Windows installer, that we distribute with
  # Gerardus
  set(GMP_INCLUDE_DIR "${GERARDUS_SOURCE_DIR}/cpp/src/third-party/CGAL-4.2/auxiliary/gmp/include")
  set(MPFR_INCLUDE_DIR "${GMP_INCLUDE_DIR}")
  set(GMP_LIBRARY_DIR "${GERARDUS_SOURCE_DIR}/cpp/src/third-party/CGAL-4.2/auxiliary/gmp/lib")
  set(MPFR_LIBRARY_DIR "${GMP_LIBRARY_DIR}")

  add_library(gmp SHARED IMPORTED GLOBAL)
  set_property(TARGET gmp PROPERTY IMPORTED_LOCATION "${GMP_LIBRARY_DIR}/libgmp-10.dll")
  set_property(TARGET gmp PROPERTY IMPORTED_IMPLIB "${GMP_LIBRARY_DIR}/libgmp-10.lib")
  add_library(mpfr SHARED IMPORTED GLOBAL)
  set_property(TARGET mpfr PROPERTY IMPORTED_LOCATION "${GMP_LIBRARY_DIR}/libmpfr-4.dll")
  set_property(TARGET mpfr PROPERTY IMPORTED_IMPLIB "${GMP_LIBRARY_DIR}/libmpfr-4.lib")

  set(GMP_FOUND TRUE 
    CACHE BOOL "GMP library is part of Gerardus, so this is always true")
  set(GMP_LIBRARY_DIRS "${GMP_LIBRARY_DIR}"
    CACHE PATH "Path to GMP library")
  set(GMP_LIBRARIES_DIR "${GMP_LIBRARY_DIR}"
    CACHE PATH "Path to GMP library")
  set(MPFR_FOUND TRUE 
    CACHE BOOL "MPFR library is part of Gerardus, so this is always true")
  set(MPFR_LIBRARY_DIRS "${MPFR_LIBRARY_DIR}"
    CACHE PATH "Path to MPFR library")
  set(MPFR_LIBRARIES_DIR "${MPFR_LIBRARY_DIR}"
    CACHE PATH "Path to MPFR library")

else(WIN32 AND NOT CYGWIN) # linux builds

  # for Unix-like builds, we have to build the GMP and MPFR libraries from the source code
  set(GMP_SOURCE_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/gmp-5.1.2)
  set(MPFR_SOURCE_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/mpfr-3.1.2)

  include(ExternalProject)

  ExternalProject_Add(GMP_GERARDUS
    PREFIX ${GMP_SOURCE_DIR}
    STAMP_DIR ${GMP_SOURCE_DIR}/output
    SOURCE_DIR ${GMP_SOURCE_DIR}
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ./configure --prefix=${GMP_SOURCE_DIR}/output 
           CC=gcc-${GCC_VERSION_REQUIRED_BY_MATLAB} CXX=g++-${GCC_VERSION_REQUIRED_BY_MATLAB}
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
  add_library(gmp SHARED IMPORTED GLOBAL)
  set_property(TARGET gmp PROPERTY IMPORTED_LOCATION ${GMP_SOURCE_DIR}/output/lib/libgmp.so)
#  set(GMP_FOUND TRUE 
#    CACHE BOOL "GMP library is part of Gerardus, so this is always true")
  set(GMP_INCLUDE_DIR ${GMP_SOURCE_DIR}
    CACHE PATH "Path to GMP include files")
  set(GMP_LIBRARY_DIR ${GMP_SOURCE_DIR}/output/lib
    CACHE PATH "Path to GMP library")
  set(GMP_LIBRARY_DIRS ${GMP_SOURCE_DIR}/output/lib
    CACHE PATH "Path to GMP library")
  set(GMP_LIBRARIES_DIR ${GMP_SOURCE_DIR}/output/lib
    CACHE PATH "Path to GMP library")
  # this line needed by CGAL to avoid looking for GMP in the system
  set(GMP_LIBRARIES "${GMP_SOURCE_DIR}/output/lib/libgmp.so") 
  
  ExternalProject_Add(MPFR_GERARDUS
    DEPENDS GMP_GERARDUS
    PREFIX ${MPFR_SOURCE_DIR}
    STAMP_DIR ${MPFR_SOURCE_DIR}/output
    SOURCE_DIR ${MPFR_SOURCE_DIR}
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ./configure --prefix=${MPFR_SOURCE_DIR}/output 
           --with-gmp-include=${GMP_SOURCE_DIR}/output/include --with-gmp-lib=${GMP_SOURCE_DIR}/output/lib
           CC=gcc-${GCC_VERSION_REQUIRED_BY_MATLAB} CXX=g++-${GCC_VERSION_REQUIRED_BY_MATLAB}
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
  add_library(mpfr SHARED IMPORTED GLOBAL)
  set_property(TARGET mpfr PROPERTY IMPORTED_LOCATION ${MPFR_SOURCE_DIR}/output/lib/libmpfr.so)
#  set(MPFR_FOUND TRUE 
#    CACHE BOOL "MPFR library is part of Gerardus, so this is always true")
  set(MPFR_INCLUDE_DIR ${MPFR_SOURCE_DIR}/src
    CACHE PATH "Path to MPFR include files")
  set(MPFR_LIBRARY_DIR ${MPFR_SOURCE_DIR}/output/lib
    CACHE PATH "Path to MPFR library")
  set(MPFR_LIBRARY_DIRS ${MPFR_SOURCE_DIR}/output/lib
    CACHE PATH "Path to MPFR library")
  set(MPFR_LIBRARIES_DIR ${MPFR_SOURCE_DIR}/output/lib
    CACHE PATH "Path to MPFR library")
  # this line needed by CGAL to avoid looking for MPRF in the system
  set(MPFR_LIBRARIES "${MPFR_SOURCE_DIR}/output/lib/libmpfr.so")

endif(WIN32 AND NOT CYGWIN)

#############################################################################################
## Third-party library: IPOPT
#############################################################################################

# # currently, only linux build available in Gerardus
# if(NOT WIN32)
#   enable_language(Fortran)
  
#   set(IPOPT_SOURCE_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/Ipopt-3.11.6)
  
#   # build IPOPT
#   include(ExternalProject)
#   ExternalProject_Add(IPOPT_GERARDUS
#     PREFIX "${IPOPT_SOURCE_DIR}"
#     STAMP_DIR "${IPOPT_SOURCE_DIR}/output"
#     SOURCE_DIR "${IPOPT_SOURCE_DIR}"
#     DOWNLOAD_COMMAND ""
#     CONFIGURE_COMMAND ./configure CC=$ENV{CC} CXX=$ENV{CXX} F77=$ENV{F77} --prefix=${IPOPT_SOURCE_DIR}/output
#     BINARY_DIR "${IPOPT_SOURCE_DIR}"
#     BUILD_COMMAND make
#     TEST_BEFORE_INSTALL 1
#     # we don't run the test, because tests for any of the linear solvers
#     #  that we don't use will fail, as they won't be installed in the
#     #  system, and this will make cmake stop with an error
#     #
#     # TEST_COMMAND make test
#     TEST_COMMAND ""
#     INSTALL_DIR "${IPOPT_SOURCE_DIR}/output"
#     INSTALL_COMMAND ""
#     LOG_CONFIGURE 0
#     LOG_BUILD 0
#     LOG_TEST 0
#     LOG_INSTALL 0
#     )
# endif(NOT WIN32)

#############################################################################################
## Third-party library: SCIP
#############################################################################################

# # currently, only linux build available in Gerardus
# if(NOT WIN32)
#   set(SCIPLIB_SOURCE_DIR ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/scipoptsuite-3.0.2)
  
#   set(GERARDUSLDFLAGS_SCIPLIB "-L${GMP_SOURCE_DIR}/output/lib -L${GERARDUS_SOURCE_DIR}/cpp/src/third-party/scipoptsuite-3.0.2/zimpl-3.3.1/lib")
  
#   # build SCIP
#   include(ExternalProject)
#   ExternalProject_Add(SCIP_GERARDUS
#     PREFIX "${SCIPLIB_SOURCE_DIR}"
#     STAMP_DIR "${SCIPLIB_SOURCE_DIR}/output"
#     SOURCE_DIR "${SCIPLIB_SOURCE_DIR}"
#     DOWNLOAD_COMMAND ""
#     CONFIGURE_COMMAND ""
#     BINARY_DIR "${SCIPLIB_SOURCE_DIR}"
#     BUILD_COMMAND make CC=$ENV{CC} CXX=$ENV{CXX} LINKCC=$ENV{CC} LINKCXX=$ENV{CXX} GERARDUSLDFLAGS=${GERARDUSLDFLAGS_SCIPLIB} GERARDUSCXXFLAGS=-I${GMP_SOURCE_DIR} READLINE=false
#     TEST_BEFORE_INSTALL 1
#     TEST_COMMAND make check
#     INSTALL_DIR "${SCIPLIB_SOURCE_DIR}/output"
#     INSTALL_COMMAND ""
#     LOG_CONFIGURE 1
#     LOG_BUILD 1
#     LOG_TEST 1
#     LOG_INSTALL 1
#     )
#   add_dependencies(SCIP_GERARDUS GMP IPOPT)

#   # contrary to custom in the rest of Gerardus, we are calling the scip
#   # library "sciplib", instead of "scip". The reason is that the actual
#   # MEX Matlab function file will be called scip.mexa64, so that it has
#   # the same name as OptiToolbox gives it in Windows
#   add_library(sciplib STATIC IMPORTED GLOBAL)
#   set_property(TARGET sciplib PROPERTY IMPORTED_LOCATION ${SCIP_SOURCE_DIR}/scip-3.0.2/lib/libscip.a)
  
# endif(NOT WIN32)

#############################################################################################
## Third-party library: ITK
#############################################################################################

# find ITK
if(WIN32)
  set(ITK_DIR "C:/Program Files/ITK/lib/cmake/ITK-4.3")
endif()
find_package(ITK)
if(NOT ITK_FOUND)
  if (WIN32)
    message(FATAL_ERROR "You need to download and install ITK v4.3.1: http://downloads.sourceforge.net/sourceforge/itk/InsightApplications-4.3.0.zip")
  else(WIN32)
    # to make things easier for the user, we are going to build ITK
    
    set(ITK_TARBALL_URL "http://sourceforge.net/projects/itk/files/itk/4.3/InsightToolkit-4.3.1.tar.gz/download")
    set(ITK_TARBALL_FILENAME "${CMAKE_BINARY_DIR}/InsightToolkit-4.3.1.tar.gz")
    set(ITK_SOURCE_DIR "${GERARDUS_SOURCE_DIR}/cpp/src/third-party/InsightToolkit-4.3.1")

    if(EXISTS ${ITK_TARBALL_FILENAME})
      
      message(STATUS "ITK source code previously downloaded: ${ITK_TARBALL_FILENAME}")
      
    else(EXISTS ${ITK_TARBALL_FILENAME})

      # download the tarball with the ITK source code
      message(STATUS "Downloading ${ITK_TARBALL_URL}\nto ${ITK_TARBALL_FILENAME}")
      file(DOWNLOAD
	${ITK_TARBALL_URL}
	${ITK_TARBALL_FILENAME}
	SHOW_PROGRESS
	STATUS ITK_DOWNLOAD_STATUS)
    
      if ((ITK_DOWNLOAD_STATUS MATCHES "0;\"No error\"") 
          OR (ITK_DOWNLOAD_STATUS MATCHES "0;\"no error\""))
	message(STATUS "ITK source tarball downloaded: ${ITK_TARBALL_FILENAME}")
      else()
	message(FATAL_ERROR "${ITK_DOWNLOAD_STATUS}\nCould not download ITK source code from\n${ITK_TARBALL_URL}")
      endif()
    endif(EXISTS ${ITK_TARBALL_FILENAME})
  
    # uncompress the Boost tarball if necessary
    if(NOT EXISTS "${ITK_SOURCE_DIR}")
      execute_process(
	COMMAND ${CMAKE_COMMAND} -E tar xzf ${ITK_TARBALL_FILENAME}
	WORKING_DIRECTORY "${GERARDUS_SOURCE_DIR}/cpp/src/third-party"
	RESULT_VARIABLE ITK_UNTAR_ERROR)
      
      # it is possible that the download seems correct, but the
      # tarball is an empty file. We need to detect these cases, and
      # exit the build with an error. In those cases, the user should
      # just download the tarball, e.g. with wget from Linux, copy it
      # to the expected ${ITK_TARBALL_FILENAME}, and re-run cmake
      if(ITK_UNTAR_ERROR)
	message(FATAL_ERROR "ITK tarball cannot be uncompressed: ${ITK_TARBALL_FILENAME}")
      endif()
      
    endif(NOT EXISTS "${ITK_SOURCE_DIR}")

    # configure ITK
    message(STATUS "Configuring ITK")
    set(ITK_CMAKE_COMMAND_FLAGS -DBUILD_TESTING:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_DOCUMENTATION:BOOL=OFF -DBUILD_EXAMPLES:BOOL=OFF -DITK_USE_REVIEW=ON -DCMAKE_BUILD_TYPE=Release -DITK_LEGACY_REMOVE=ON -DCMAKE_C_COMPILER=/usr/bin/gcc-${GCC_VERSION_REQUIRED_BY_MATLAB} -DCMAKE_CXX_COMPILER=/usr/bin/g++-${GCC_VERSION_REQUIRED_BY_MATLAB})
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E make_directory "bin"
      WORKING_DIRECTORY "${ITK_SOURCE_DIR}"
      RESULT_VARIABLE ITK_MKDIR_ERROR)
    execute_process(
      COMMAND ${CMAKE_COMMAND} ${ITK_CMAKE_COMMAND_FLAGS} ..
      WORKING_DIRECTORY "${ITK_SOURCE_DIR}/bin"
      RESULT_VARIABLE ITK_CMAKE_ERROR)
    if(ITK_CMAKE_ERROR)
      message(FATAL_ERROR "Running cmake to configure the ITK source code failed")
    endif()

    # build ITK
    message(STATUS "Building ITK")
    execute_process(
      COMMAND make -j6
      WORKING_DIRECTORY "${ITK_SOURCE_DIR}/bin"
      RESULT_VARIABLE ITK_MAKE_ERROR)
    if(ITK_MAKE_ERROR)
      message(FATAL_ERROR "Running make to build the ITK source code failed")
    endif()

    # it's better if ITK gets installed in the system, but this means
    # that we need superuser privileges. We don't want the user to
    # enter his password in the script, so until we find a better
    # solution, we are just going to ask the user to go to the ITK bin
    # directory to run "sudo make install", and the exit the script
    message(STATUS "******************************************************************************")
    message(STATUS "ITK built successfully. Now you need to install it manually. For that, run")
    message(STATUS "   pushd \"${ITK_SOURCE_DIR}/bin\"")
    message(STATUS "   sudo make install; popd")
    message(STATUS "If ITK installed successfully, you can then run")
    message(STATUS "   cmake ..")
    message(STATUS "again to finish configuring Gerardus.")
    message(STATUS "******************************************************************************")
    return()

  endif(WIN32)

else(NOT ITK_FOUND)

  # ITK has been found, so we add the include directory and display its version
  include(${ITK_USE_FILE})
  message(STATUS "ITK v" ${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH})

endif(NOT ITK_FOUND)


## flags for all targets
IF(WIN32)
  # to remove warnings from MSVC compiler asking to use fopen_s() an other unportable
  # functions exclusive to Microsoft
  ADD_DEFINITIONS("-D_CRT_SECURE_NO_DEPRECATE")
ELSE(WIN32)
  # optimise and show all warnings
  ADD_DEFINITIONS("-O2 -Wall")

  # For Mac OSX, we need to specify the system architecture.
  # If this hasn't been set, we choose 64-bit intel.
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    IF(NOT CMAKE_OSX_ARCHITECTURES)
      SET(CMAKE_OSX_ARCHITECTURES x86_64)
    ENDIF(NOT CMAKE_OSX_ARCHITECTURES)
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

ENDIF(WIN32)

#############################################################################################
## Third-party library: Boost
#############################################################################################

# Boost libraries are required by the executable programs in Gerardus,
# and anything that uses CGAL. In the past, we tried to build the
# Boost libraries we needed and link to them. However, Matlab also
# uses and provides its own Boost libraries in the same directory
# where it keeps the mex libraries. This created error prone
# situations where MEX functions where linked to two different
# versions of the same Boost library.
#
# To avoid that trouble, we now avoid that duplication. Instead, we
# figure out what version of Boost Matlab is using, download the
# corresponding header files, install them in Gerardus, and link to
# the Boost libraries already distributed with Matlab.

# get the version of Boost distributed with Matlab
if(WIN32)
  set(MATLAB_BOOST_LIB_PATH "${MATLAB_ROOT}/bin/win64")
else()
  get_filename_component(MATLAB_BOOST_LIB_PATH ${MATLAB_MEX_LIBRARY} PATH)
endif()
file(GLOB AUX_BOOST_VERSION  "${MATLAB_BOOST_LIB_PATH}/*boost_system*")
get_filename_component(AUX_BOOST_VERSION ${AUX_BOOST_VERSION} NAME)
if(WIN32)
  string(REGEX MATCH "[0-9]+[._][0-9]+" AUX_BOOST_VERSION ${AUX_BOOST_VERSION}) # e.g. 1_44
#  set(AUX_BOOST_VERSION ${AUX_BOOST_VERSION}_0)                                 # e.g. 1_44_0
  string(REPLACE "_" "." AUX_BOOST_VERSION ${AUX_BOOST_VERSION})                # e.g. 1.44
else()
  string(REGEX MATCH "[0-9]+[._][0-9]+[._][0-9]+" AUX_BOOST_VERSION ${AUX_BOOST_VERSION}) # e.g. 1.44.0
                                                                                          # or 1.44
endif()
message(STATUS "Matlab is linked to Boost version ${AUX_BOOST_VERSION}")
string(REPLACE "." "_" AUX_BOOST_VERSION_ ${AUX_BOOST_VERSION})                           # e.g. 1_44_0
                                                                                          # or 1_44
# path to where the Boost header files will be found in Gerardus
if(WIN32)
  set(BOOST_SOURCE_DIR
    ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/boost_${AUX_BOOST_VERSION_}_0)
  set(Boost_INCLUDE_DIR "${GERARDUS_SOURCE_DIR}/include/boost-${AUX_BOOST_VERSION_}")
else()
  set(BOOST_SOURCE_DIR
    ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/boost_${AUX_BOOST_VERSION_})
  set(Boost_INCLUDE_DIR "${GERARDUS_SOURCE_DIR}/include")
endif()

# path to where the Boost libraries will be found in Gerardus
if(WIN32)
  set(Boost_LIBRARY_DIR "${GERARDUS_SOURCE_DIR}/lib")
else()
  set(Boost_LIBRARY_DIR "${MATLAB_BOOST_LIB_PATH}")
endif()

# we want to know whether Boost has already been installed. As
# checking for the whole project is complicated, we are going to check
# only whether file /boost/version.hpp is in its expected location. If
# it isn't, then we install the boost libraries (in linux, we just
# soft link to them; in Windows, we need to build and install them)
find_file(Boost_version_hpp version.hpp
  PATHS "${Boost_INCLUDE_DIR}/boost/"
  NO_DEFAULT_PATH)
if(NOT Boost_version_hpp)

  message(STATUS "${Boost_INCLUDE_DIR}/boost/version.hpp not found.\nPreparing to install Boost libraries")

  if(WIN32)
    set(BOOST_TARBALL_URL 
      http://downloads.sourceforge.net/project/boost/boost/${AUX_BOOST_VERSION}.0/boost_${AUX_BOOST_VERSION_}_0.tar.gz)
    set(BOOST_TARBALL_FILENAME
      ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/boost_${AUX_BOOST_VERSION_}_0.tar.gz)
  else()
    set(BOOST_TARBALL_URL 
      http://downloads.sourceforge.net/project/boost/boost/${AUX_BOOST_VERSION}/boost_${AUX_BOOST_VERSION_}.tar.gz)
    set(BOOST_TARBALL_FILENAME
      ${GERARDUS_SOURCE_DIR}/cpp/src/third-party/boost_${AUX_BOOST_VERSION_}.tar.gz)
  endif()

  # check whether we need to download the tarball with the Boost
  # source code, or whether it's already been downloaded
  if(EXISTS ${BOOST_TARBALL_FILENAME})

    message(STATUS "Boost source code previously downloaded")

  else()

    # download the tarball with the Boost source code that corresponds
    # to the Boost version in Matlab
    message(STATUS "Downloading ${BOOST_TARBALL_URL}\nto ${BOOST_TARBALL_FILENAME}")
    file(DOWNLOAD
      ${BOOST_TARBALL_URL}
      ${BOOST_TARBALL_FILENAME}
      SHOW_PROGRESS
      STATUS Boost_DOWNLOAD_STATUS)

    if ((Boost_DOWNLOAD_STATUS MATCHES "0;\"No error\"") 
        OR (Boost_DOWNLOAD_STATUS MATCHES "0;\"no error\""))
      message(STATUS "Boost source tarball downloaded")
    else()
      message(FATAL_ERROR "${Boost_DOWNLOAD_STATUS}\nCould not download Boost source code from\n${BOOST_TARBALL_URL}")
    endif()
  endif()

  # uncompress the Boost tarball if necessary
  if(NOT EXISTS "${BOOST_SOURCE_DIR}")
    execute_process(
	COMMAND ${CMAKE_COMMAND} -E tar xzf ${BOOST_TARBALL_FILENAME}
	WORKING_DIRECTORY ${GERARDUS_SOURCE_DIR}/cpp/src/third-party
	RESULT_VARIABLE BOOST_UNTAR_ERROR)

    # it is possible that the download seems correct, but the
    # tarball is an empty file. We need to detect these cases, and
    # exit the build with an error. In those cases, the user should
    # just download the tarball, e.g. with wget from Linux, copy it
    # to the expected ${BOOST_TARBALL_FILENAME}, and re-run cmake
    if(BOOST_UNTAR_ERROR)
      message(FATAL_ERROR "Boost tarball cannot be uncompressed: ${BOOST_TARBALL_FILENAME}")
    endif()

  endif()

  # In Windows, we need to build and install the libraries, because
  # Matlab only makes the .dll Boost files available, and we need the
  # .lib for linking, even for dynamic linking
  #
  # We don't use ExternalProject_Add() because we need the Boost
  # libraries to be available for FindBoost() and the rest of the
  # targets in Gerardus that link to them at configuration time
  if (WIN32)

    # for the build, we need to tell the Boost builder whether this is
    # a 32 or 64 bit system
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(Boost_address_model 64)
    else()
      set(Boost_address_model 32)
    endif()

    message(STATUS "Configuring Boost libraries")
    execute_process(
      COMMAND bootstrap.bat
      WORKING_DIRECTORY ${BOOST_SOURCE_DIR}
      RESULT_VARIABLE COMMAND_RESULT)
    if(${COMMAND_RESULT} MATCHES 0)
      message(STATUS "Boost libraries configured")
    else()
      message(FATAL_ERROR "Boost libraries configuration error:\n${COMMAND_RESULT}")
    endif()

    message(STATUS "Building and installing Boost libraries and header files")
    execute_process(
      COMMAND bjam.exe install variant=release 
        --prefix=${GERARDUS_SOURCE_DIR}
        --build-dir=./bin address-model=${Boost_address_model}
        --build-type=complete
        --with-system --with-filesystem --with-thread --with-date_time --with-serialization
      WORKING_DIRECTORY ${BOOST_SOURCE_DIR}
      RESULT_VARIABLE COMMAND_RESULT)
    if(${COMMAND_RESULT} MATCHES 0)
      message(STATUS "Boost libraries and header files built and installed")
    else()
      message(FATAL_ERROR "Boost libraries build error:\n${COMMAND_RESULT}")
    endif()

  # In Linux it's easier, because we can link directly to the .so
  #  files provided by Matlab
  else()

    # copy the header files to Gerardus' include directory
    file(COPY "${BOOST_SOURCE_DIR}/boost"
      DESTINATION "${Boost_INCLUDE_DIR}")

    # patch the Boost include files if necessary
    # version 1.44.0: https://bbs.archlinux.org/viewtopic.php?pid=1126374
    if (${AUX_BOOST_VERSION} MATCHES "1.44.0")
      message(STATUS "This version of Boost needs patching")
      file(READ "${Boost_INCLUDE_DIR}/boost/thread/xtime.hpp" AUX)
      string(REPLACE TIME_UTC TIME_UTC_ AUX "${AUX}")
      file(WRITE "${Boost_INCLUDE_DIR}/boost/thread/xtime.hpp" "${AUX}")
    endif()

    message(STATUS "Boost header files have been patched and installed")

  endif()

else()
  
  message(STATUS "${Boost_INCLUDE_DIR}/boost/version.hpp found.")
  message(STATUS "Assuming that Boost libraries and headers are installed")

endif(NOT Boost_version_hpp)

#set(BOOST_ROOT "${MATLAB_BOOST_LIB_PATH}")
set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIR}")
set(Boost_LIBRARY_DIRS "${Boost_LIBRARY_DIR}")

message(STATUS "Boost include path: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost lib path:     ${Boost_LIBRARY_DIRS}")

set(BOOST_INCLUDEDIR "${Boost_INCLUDE_DIR}")
set(BOOST_LIBRARYDIR "${Boost_LIBRARY_DIR}")
set(Boost_DEBUG                  FALSE)
set(Boost_NO_SYSTEM_PATHS        TRUE)
if (WIN32)
  set(Boost_USE_STATIC_LIBS      OFF)
else()
  set(Boost_USE_STATIC_LIBS      OFF)
endif()
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     OFF)
add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

# integrate the local Boost libraries with the rest of the CMake
# build system
if(WIN32)

  # in Windows, we can just look for the library files
  find_package(Boost ${AUX_BOOST_VERSION}
    COMPONENTS date_time filesystem serialization system thread)

else()

  # FindBoost.cmake will not find the boost libraries in linux because
  # Matlab saves them as e.g. libboost_system.so.1.44.0, so we need to
  # add them by hand

  # library suffixes get automatically added to the end of the name
  # when looking for libraries. Because we are looking for libraries
  # with names like e.g. libboost_system.so.1.44.0, we are going to
  # temporarily allow that kind of suffix
  set(AUX "${CMAKE_FIND_LIBRARY_SUFFIXES}")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so.${AUX_BOOST_VERSION}")

  find_package(Boost ${AUX_BOOST_VERSION}
    COMPONENTS date_time filesystem system thread)

  # recover the library suffixes
  set(CMAKE_FIND_LIBRARY_SUFFIXES "${AUX}")

endif()

#############################################################################################
## Third-party library: CGAL
#############################################################################################

# note that CGAL depends on Boost, but it's too cumbersome to add
# an explicit dependency

# call the subproject, as CGAL is distributed with a CMake build
add_subdirectory(cpp/src/third-party/CGAL-4.2)

# after CGAL builds, we need to copy the header file compiler_config.h
# it generates to the CGAL include directory, so that it's available
# for the other sub-projects that use CGAL (e.g. itk_imfilter)
add_custom_command(
  OUTPUT ${CGAL_SOURCE_DIR}/include/CGAL/compiler_config.h
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/include/CGAL/compiler_config.h" "${CGAL_SOURCE_DIR}/include/CGAL/"
  MAIN_DEPENDENCY ${CMAKE_CURRENT_BINARY_DIR/../../include/CGAL/compiler_config.h}
  DEPENDS CGAL
  )
add_custom_target(copy_compiler_config.h DEPENDS "${CGAL_SOURCE_DIR}/include/CGAL/compiler_config.h")

#############################################################################################
## Configuration options for all Gerardus
#############################################################################################

# make the libraries that live in the Gerardus lib directory available
# to all targets
link_directories(
  ${GERARDUS_SOURCE_DIR}/lib
  )

##################################################################
## Block so that rpaths are added to both the MEX file in the build
## tree and in the installation directory. Otherwise, the rpath is
## stripped when installing, and the function will crash as it won't be
## able to find the ITK shared libraries

# auxiliary libraries built within Gerardus will not be installed in
# the system, but locally. If at some point we need to install in the
# system, the following line will have to be changed, and we'll have
# to uncomment the lines at the end of the rpath block, and use those
# for non-system installs
SET(CMAKE_INSTALL_PREFIX ${GERARDUS_SOURCE_DIR}/lib)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX})

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
#LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
#IF("${isSystemDir}" STREQUAL "-1")
#   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#ENDIF("${isSystemDir}" STREQUAL "-1")

## End rpath block
##################################################################

## flags for all targets
IF(WIN32)

  # to remove warnings from MSVC compiler asking to use fopen_s() an other unportable
  # functions exclusive to Microsoft
  ADD_DEFINITIONS("-D_CRT_SECURE_NO_DEPRECATE")

  # avoid MSVC compiler warnings like
  #    'std::copy': Function call with parameters that may be unsafe
  # given by ITK functions
  ADD_DEFINITIONS("-D_SCL_SECURE_NO_WARNINGS")

ELSE(WIN32)

  # optimise and show all warnings
  ADD_DEFINITIONS("-O2 -Wall")

ENDIF(WIN32)

#############################################################################################
## Call the subprojects 
#############################################################################################

ADD_SUBDIRECTORY(cpp/src)
ADD_SUBDIRECTORY(cpp/src/third-party/mba)
ADD_SUBDIRECTORY(matlab)
