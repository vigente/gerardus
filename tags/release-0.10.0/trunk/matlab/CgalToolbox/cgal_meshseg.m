function varargout = cgal_meshseg(varargin)
% CGAL_MESHSEG  Surface meshing of an isosurface from a segmentation or
% grayscale image.
%
% This function is a Matlab wrapper of the CGAL 3D Surface Mesh Generation
% for a grayscale input image.
%
% http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Surface_mesher/Chapter_main.html
%
% Note that even though CGAL internally associates x <-> rows, y <-> cols,
% this functions presents to the user the usual Matlab convention of 
% x <-> cols, y <-> rows.
%
% [TRI, X] = cgal_meshseg(IM, ISOVAL)
%
%   IM is a 3D array or a SCIMAT struct with a segmentation or a grayscale
%   image.
%
%   ISOVAL is a scalar value that defines the isosurface to be meshed.
%
%   Note that if you have a binary segmentation (background=0, segmented
%   voxels=1), ISOVAL=1 will give a very tight surface connecting the
%   centres of the segmented voxels on the boundary. With small values,
%   ISOVAL=0.01, the surface will be close to the centres of the adjacent
%   background voxels. With ISOVAL=0.5, the surface will be halfway between
%   the centres of the background and segmented voxels. The latter is
%   usually the desired result.
%
%   TRI is a 3-column matrix. Each row represents the indices of the three
%   vertices that form a triangle. TRI as a whole represents the closed
%   surface.
%
%   X is a 3-column matrix. Each row represents the Cartesian coordinates
%   of a vertex on the surface, indexed by TRI values.
%
% ... = cgal_meshseg(IM, ISOVAL, MINALPHA, MAXRAD, MAXD, C, MANIFOLD)
%
%   MINALPHA, MAXRAD, MAXD are scalars that implement the three meshing
%   criteria in CGAL::Surface_mesh_default_criteria_3
%
%   http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Surface_mesher_ref/Class_Surface_mesh_default_criteria_3.html#Cross_link_anchor_1519
%
%   MINALPHA is "a lower bound on the minimum angle in degrees of the
%   surface mesh facets". By default, MINALPHA = 30.0.
%
%   MAXRAD is "an upper bound on the radius of surface Delaunay balls. A
%   surface Delaunay ball is a ball circumscribing a facet, centered on the
%   surface and empty of vertices. Such a ball exists for each facet of the
%   current surface mesh. Indeed the current surface mesh is the Delaunay
%   triangulation of the current sampling restricted to the surface which
%   is just the set of facets in the three dimensional Delaunay
%   triangulation of the sampling that have a Delaunay surface ball". By
%   default, MAXRAD is 1/2 of the minimum voxel size dimension. For
%   example, if voxels have size [0.1 0.2 0.5], then by default
%   MAXRAD=0.05.
%
%   MAXD is "an upper bound on the center-center distances of the surface
%   mesh facets. The center-center distance of a surface mesh facet is the
%   distance between the facet circumcenter and the center of its surface
%   Delaunay ball". By default, MAXD is computed the same as MAXRAD.
%
%   C is a 3-vector with the coordinates of the centre of the bounding
%   sphere used by the meshing algorithm. This is an important parameter.
%   If C is close to the surface, it can produce lots of little triangles
%   in that area. If C is outside the segmentation, the computed mesh may
%   be incomplete.
%
%   MANIFOLD is a boolean flag.
%
%     MANIFOLD=false, the mesher uses CGAL::Non_manifold_tag: "When
%     instantiated with the tag Non_manifold_tag the function template
%     make_surface_mesh does not ensure that the output mesh is a manifold
%     surface. The manifold property of output mesh may nevertheless result
%     from the choice of appropriate meshing criteria".
%
%     MANIFOLD=true, the mesher uses CGAL::Manifold_tag: "When instantiated
%     with the tag Manifold_tag the function template make_surface_mesh
%     ensures that the output mesh is a manifold surface without boundary".
%
% Important!
%
% Note that this function can produce meshes with (1) stray vertices that
% belong to no triangle, (2) triangles not oriented with respect to the
% manifold, and (3) holes in the surface. These problems can be solved with
% the following functions available in Gerardus:
%
%   [tri, x] = tri_squeeze(tri, x); % remove stray vertices
%   [x, tri] = meshcheckrepair(x, tri, 'deep'); % correct orientation
%   tri = cgal_tri_fillholes(tri, x); % fill holes in the surface
%
% See also: bwmesh, tri_squeeze, meshcheckrepair, cgal_tri_fillholes.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.1.1
% $Rev$
% $Date$
%
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

error('MEX function not found')
