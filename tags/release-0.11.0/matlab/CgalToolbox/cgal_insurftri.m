function varargout = cgal_insurftri(varargin)
% CGAL_INSURFACETRIANGULATION  Find whether a point is inside or outside a
% closed surface
%
%   This function evaluates whether one or more points belong inside a
%   closed surface. First, we check whether the point is on the surface
%   itself (in that case, it's considered inside). If not, a ray is
%   projected from the point and the intersections with the surface are
%   counted. An odd number means that the point is inside. This approach
%   fails if the point is not on the surface, but the ray lies on the
%   surface or crosses a vertex, because this spans many arbitrary
%   intersections. To solve this problem, a few rays are used for each
%   point, and the majority vote decides whether it's inside or outside.
%
% ISIN = cgal_insurftri(TRI, X, XI)
%
%   TRI is a 3-column matrix. Each row represents the indices of the tree
%   vertices that form a triangle. TRI as a whole represents the closed
%   surface.
%
%   X is a 3-column matrix. Each row represents the Cartesian coordinates
%   of a vertex on the surface, indexed by TRI values.
%
%   XI is a 3-column matrix. Each row represents the Carterian coordinates
%   of a point for which we want to find whether it's inside or outside the
%   closed surface. Note that if you want to test all the voxels in an
%   image, it is very slow and memory intensive to generate coordinates for
%   each voxel. In that scenario, it is much better to use the cell array
%   CI syntax shown below, and provide only values for the coordinate axes.
%
%   ISIN is a boolean vector with one element per point in XI. True means
%   that the corresponding point is inside the closed surface (or on the
%   surface), and false means that it's outside.
%
% ISIN = cgal_insurftri(TRI, X, CI)
%
%   CI is a cell array CI={XI, YI, ZI}, where XI, YI and ZI are row vectors
%   that describe a rectangular grid. For example,
%
%     CI={linspace(-.25, .25, 5), ...
%         linspace(-.25, .25, 4), ...
%         linspace(-.25, .25, 3)};
%
%   describes a sampling grid of 4 rows x 5 columns x 3 slices (note that
%   rows correspond to YI and columns to XI), of a domain 
%   [-0.25, 0.25] x [-0.25, 0.25] x [-0.25, 0.25].
%
% ISIN = cgal_insurftri(..., DIRECTIONS, TOL)
%
%   DIRECTIONS is a 3-column matrix. Each row represents a vector with a
%   ray direction. By default, 
%
%          DIRECTIONS=[ 1.0,  0.0,  0.0; ...
%                      -1.0,  1.0,  1.0; ...
%                      -1.0, -1.0, -1.0]
%
%   This default can fail with regular voxels, as rays may cross vertices.
%   A good practical alternative is to use a few random directions, e.g.
%
%          DIRECTIONS=rand(5, 3);
%
%   Warning! For the voting system to make sense, select an odd number of
%   rays.
%
%   TOL is a scalar with the distance tolerance. Points at distance <= TOL
%   are considered to be on the surface, and thus "inside". By default,
%   TOL=1e-15.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
% Version: 0.2.0
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
