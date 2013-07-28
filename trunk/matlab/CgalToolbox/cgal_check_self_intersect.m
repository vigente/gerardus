function c = cgal_check_self_intersect(tri, x)
% CgalCheckSelfIntersect.cpp
%
% CGAL_CHECK_SELF_INTERSECT Check for self-intersections in a triangular
% mesh
%
% This function checks whether each triangle in a mesh intersects any other
% triangle. Finding self-intersections is useful to detect topological
% problems.
%
% C = cgal_check_self_intersect(TRI, X)
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh.
%
%   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of the
%   i-th node in the mesh.
%
%   C is a vector with one element per triangle in TRI. It gives a count of
%   the number of times TRI(I,:) causes a self-intersection in the mesh.
%
% This function uses an AABB tree component [1] to efficiently perform the
% intersection queries. However, as the CGAL documentation notes, "this
% component is not suited to the problem of finding all intersecting pairs
% of objects", so there's probably room for improvement.
%
% [1] http://www.cgal.org/Manual/latest/doc_html/cgal_manual/AABB_tree/Chapter_main.html

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.1.0
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

error('MEX file not found')
