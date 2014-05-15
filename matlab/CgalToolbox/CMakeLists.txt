# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2012-2014 University of Oxford
# Version: 0.2.22
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

include_directories(
  "${Boost_INCLUDE_DIR}"
  "${GERARDUS_SOURCE_DIR}/cpp/src/third-party/CGAL-4.2/include"
  "${GERARDUS_SOURCE_DIR}/cpp/src/third-party/CGAL-4.2/demo/Polyhedron/include"
  "${GMP_INCLUDE_DIR}"
  "${MPFR_INCLUDE_DIR}"
  ..
  .
  )
# only for Linux/Mac. In Windows, linking to the Boost libraries 
# causes "one or more multiply defined symbols found" link errors
if(NOT WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()

################################################################
## cgal_insurftri()
################################################################

ADD_MEX_FILE(cgal_insurftri
  CgalInSurfaceTriangulation.cpp)

TARGET_LINK_LIBRARIES(cgal_insurftri
  CGAL
  CGAL_ImageIO
  ${ITK_LIBRARIES})
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_insurftri copy_compiler_config.h)

################################################################
## cgal_closest_trifacet()
################################################################

ADD_MEX_FILE(cgal_closest_trifacet
  CgalClosestTriFacet.cpp
  )

TARGET_LINK_LIBRARIES(cgal_closest_trifacet
  CGAL
  CGAL_ImageIO
  ${ITK_LIBRARIES})
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_closest_trifacet copy_compiler_config.h)

################################################################
## cgal_trifacet_area()
################################################################

ADD_MEX_FILE(cgal_trifacet_area
  CgalTriFacetArea.cpp
  )

TARGET_LINK_LIBRARIES(cgal_trifacet_area
  CGAL
  CGAL_ImageIO
  ${ITK_LIBRARIES})
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_trifacet_area copy_compiler_config.h)

################################################################
## cgal_check_self_intersect()
################################################################

add_mex_file(cgal_check_self_intersect
  CgalCheckSelfIntersect.cpp
  )

# In Windows, linking to the Boost libraries causes "one or more 
# multiply defined symbols found" link errors
if(WIN32)
  target_link_libraries(cgal_check_self_intersect
    CGAL
    CGAL_ImageIO
    ${ITK_LIBRARIES})
else()
  target_link_libraries(cgal_check_self_intersect
    ${Boost_THREAD_LIBRARY}
    CGAL
    CGAL_ImageIO
    ${ITK_LIBRARIES})
endif()
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_check_self_intersect copy_compiler_config.h)

################################################################
## cgal_alpha_shape3()
################################################################

add_mex_file(cgal_alpha_shape3
  CgalAlphaShape3.cpp
  )

target_link_libraries(cgal_alpha_shape3
  CGAL
  CGAL_ImageIO
  ${ITK_LIBRARIES})
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_alpha_shape3 copy_compiler_config.h)

################################################################
## cgal_fixed_alpha_shape3()
################################################################

add_mex_file(cgal_fixed_alpha_shape3
  CgalFixedAlphaShape3.cpp
  )

target_link_libraries(cgal_fixed_alpha_shape3
  CGAL
  CGAL_ImageIO
  ${ITK_LIBRARIES})
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_fixed_alpha_shape3 copy_compiler_config.h)

################################################################
## cgal_meshseg()
################################################################

add_mex_file(cgal_meshseg
  CgalMeshSegmentation.cpp
  )

# In Windows, linking to the Boost libraries causes "one or more 
# multiply defined symbols found" link errors
if(WIN32)
  target_link_libraries(cgal_meshseg
    CGAL
    CGAL_ImageIO
    ${ITK_LIBRARIES})
else()
  target_link_libraries(cgal_meshseg
    ${Boost_THREAD_LIBRARY}
    CGAL
    CGAL_ImageIO
    ${ITK_LIBRARIES})
endif()
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_meshseg copy_compiler_config.h)

################################################################
## cgal_surfsubdivision()
################################################################

add_mex_file(cgal_surfsubdivision
  CgalSurfaceSubdivision.cpp
  )

target_link_libraries(cgal_surfsubdivision
  CGAL
  CGAL_ImageIO
  ${ITK_LIBRARIES})
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_surfsubdivision copy_compiler_config.h)

################################################################
## cgal_tri_fillholes()
################################################################

add_mex_file(cgal_tri_fillholes
  CgalTriFillHoles.cpp
  )

target_link_libraries(cgal_tri_fillholes
  CGAL
  CGAL_ImageIO
  ${ITK_LIBRARIES})
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_tri_fillholes copy_compiler_config.h)

################################################################
## cgal_tri_simplify()
################################################################

add_mex_file(cgal_tri_simplify
  CgalTriSimplify.cpp
  )

target_link_libraries(cgal_tri_simplify
  CGAL
  CGAL_ImageIO
  ${ITK_LIBRARIES})
# we need this dependency to make sure that the CGAL library is built
# before we try to build the mex function
add_dependencies(cgal_tri_simplify copy_compiler_config.h)

################################################################
## installation of targets
################################################################

IF(WIN32)
  INSTALL(TARGETS
    cgal_insurftri
    cgal_closest_trifacet
    cgal_trifacet_area
    cgal_check_self_intersect
    cgal_alpha_shape3
    cgal_fixed_alpha_shape3
    cgal_meshseg
    cgal_surfsubdivision
    cgal_tri_fillholes
    cgal_tri_simplify
    RUNTIME
    DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}")
ELSE(WIN32)
  INSTALL(TARGETS
    cgal_insurftri
    cgal_closest_trifacet
    cgal_trifacet_area
    cgal_check_self_intersect
    cgal_alpha_shape3
    cgal_fixed_alpha_shape3
    cgal_meshseg
    cgal_surfsubdivision
    cgal_tri_fillholes
    cgal_tri_simplify
    LIBRARY
    DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}")
ENDIF(WIN32)
