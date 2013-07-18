# Author: Ramon Casero <rcasero@gmail.com>
# Copyright © 2013 University of Oxford
# Version: 0.1.0
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

add_mex_file(perform_front_propagation_2d
  mex/perform_front_propagation_2d.cpp
  mex/perform_front_propagation_2d_mex.cpp
  mex/fheap/fib.cpp)

add_mex_file(perform_front_propagation_3d
  mex/perform_front_propagation_3d.cpp
  mex/perform_front_propagation_3d_mex.cpp
  mex/fheap/fib.cpp)

add_mex_file(perform_circular_front_propagation_2d
  mex/perform_circular_front_propagation_2d.cpp 
  mex/perform_front_propagation_2d.cpp 
  mex/fheap/fib.cpp)

add_mex_file(perform_front_propagation_anisotropic
  mex/perform_front_propagation_anisotropic.cpp)

add_mex_file(fm2dAniso
  mex/anisotropic-fm-feth/fm2dAniso.cpp)

add_mex_file(skeleton
  mex/skeleton.cpp)

add_mex_file(eucdist2
  mex/eucdist2.c)

add_mex_file(perform_front_propagation_mesh
  mex/perform_front_propagation_mesh.cpp
  mex/gw/gw_core/GW_Config.cpp
  mex/gw/gw_core/GW_FaceIterator.cpp
  mex/gw/gw_core/GW_SmartCounter.cpp
  mex/gw/gw_core/GW_VertexIterator.cpp
  mex/gw/gw_core/GW_Face.cpp
  mex/gw/gw_core/GW_Mesh.cpp
  mex/gw/gw_core/GW_Vertex.cpp
  mex/gw/gw_geodesic/GW_GeodesicFace.cpp
  mex/gw/gw_geodesic/GW_GeodesicMesh.cpp
  mex/gw/gw_geodesic/GW_GeodesicPath.cpp
  mex/gw/gw_geodesic/GW_GeodesicPoint.cpp
  mex/gw/gw_geodesic/GW_TriangularInterpolation_Cubic.cpp
  mex/gw/gw_geodesic/GW_GeodesicVertex.cpp
  mex/gw/gw_geodesic/GW_TriangularInterpolation_Linear.cpp
  mex/gw/gw_geodesic/GW_TriangularInterpolation_Quadratic.cpp)

################################################################
## installation of targets
################################################################

if(WIN32)
  install(TARGETS
    perform_front_propagation_2d
    perform_front_propagation_3d
    perform_circular_front_propagation_2d
    perform_front_propagation_anisotropic
    fm2dAniso
    skeleton
    eucdist2
    perform_front_propagation_mesh
    RUNTIME
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
else(WIN32)
  install(TARGETS
    perform_front_propagation_2d
    perform_front_propagation_3d
    perform_circular_front_propagation_2d
    perform_front_propagation_anisotropic
    fm2dAniso
    skeleton
    eucdist2
    perform_front_propagation_mesh
    LIBRARY
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
endif(WIN32)
