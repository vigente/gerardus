function [tri, x] = delaunay_sphere(N, xi)
% DELAUNAY_SPHERE  Triangulation of a sphere or a surface topologically
% equivalent to a sphere so that it can be used as a mesh
%
% [TRI, X] = delaunay_sphere(N)
%
%   Number of samples of the sphere (in latitude and longitude).
%
%   TRI is a 3-column matrix with a triangulation of the sphere. Each
%   element in TRI is an index to a row in X. Each row represents the three
%   vertices of a triangle on the sphere.
%
%   X is a 3-column matrix with the Euclidean coordinates of the
%   triangulation vertices.
%
%   The surface can be plotted running
%
%     trisurf(tri, x(:, 1), x(:, 2), x(:, 3));
%
% ... = delaunay_sphere(N, XI)
%
%   XI is an (N+1,N+1,3)-array with the coordinates of the vertices of a
%   surface topologically equivalent to a sphere. By default, XI is the
%   concatenation of the output of sphere(N), and the surface is a unit
%   sphere.
%
%   The surface can be plotted running
%
%     surf(xi(:, :, 1), xi(:, :, 2), xi(:, :, 3))
%
% See also: scimat_tri_to_raster.

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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 2);
nargoutchk(0, 2);

if (nargin > 1 ...
        && (size(xi, 1) ~= N+1 || size(xi, 2) ~= N+1 || size(xi, 3) ~= 3))
    error('XI must be a (N+1, N+1, 3) array')
end

% sphere surface
[x, y, z] = sphere(N);

% remove repeated points
[x, idx] = unique([x(:) y(:) z(:)], 'rows');

% triangulate the sphere and extract the surface
tri = DelaunayTri(x);
tri = freeBoundary(tri);

% if an input point configuration was provided, we use it to replace the
% unit sphere points coordinates
if (nargin > 1)
    xi = reshape(xi, (N+1)*(N+1), 3);
    x = xi(idx, :);
end
