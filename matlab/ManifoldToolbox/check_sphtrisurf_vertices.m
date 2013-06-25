function ok = check_sphtrisurf_vertices(tri, uv, idx)
% CHECK_SPHTRISURF_VERTICES  Check integrity of the vertices of a spherical
% triangular mesh, to make sure facets do not overlap
%
% In a spherical triangular mesh, the angles incident to each vertex add up
% to 2*pi (360º) if the mesh has no triangle overlaps. Schematically, note
% that there are 8 (spherical) triangles incident to the central vertex
% below:
%
%    *--------*-------*
%    | \      |      /|
%    |  \     |     / |
%    |   \    |    /  |     Note: this is a just a planar representation.
%    |    \   |   /   |           In reality, we are considering spherical
%    |     \  |  /    |           triangles connecting vertices on a sphere
%    |      \ | /     |
%    |       \|/      |
%    *--------*-------*
%    |       /|\      |
%    |      / | \     |
%    |     /  |  \    |
%    |    /   |   \   |
%    |   /    |    \  |
%    |  /     |     \ |
%    | /      |      \|
%    *--------*-------*
%
% The sum of these 8 angles must add up to 2*pi (360º) if the central
% vertex is within the convex hull of its 8 neighbours. But if the vertex
% is shifted outside the convex hull of its neighbours, this produces
% triangle overlap and the sum is no longer 2*pi.
%
% OK = check_sphtrisurf_vertices(TRI, UV)
%
%   TRI is a 3-column matrix with a surface triangulation. Each row has the
%   indices of 3 vertices forming a spherical triangle.
%
%   UV is a 2-column matrix with the spherical coordinates, UV=[LAT LON] of
%   the mesh vertices.
%
% Dependencies: sss() from the SphericalTrig toolbox by Rody Oldenhuis,
% which implements formulas derived in Wertz (2001).
%
% James R. Wertz. "Mission Geometry: Orbit and Constellation Design and
% Management". Published Jointly by Microcosm Press and Kluwer Academic
% Publishers, 2001.

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

% check arguments
narginchk(2, 3);
nargoutchk(0, 1);

if (size(tri, 2) ~= 3)
    error('TRI must be 3-column matrix')
end
if (size(uv, 2) ~= 2)
    error('UV must be 2-column matrix')
end

% defaults
if (nargin < 3)
    idx = 1:size(uv, 1);
end
tol = 1e-12; % permitted error when adding the angles incident to a vertex

% simplify reading code even if we have to duplicate storage of point
% coordinates
lat = uv(:, 1);
lon = uv(:, 2);

% loop every vertex in the mesh
ok = false(length(idx), 1);
for I = 1:length(idx)
    
    % get current vertex index
    v = (idx(I));
    
    % get triangles adjacent to current vertex
    [row, ~] = find(tri == v);
    triadj = tri(row, :);
    
    % reorder vertices in the triangles, so that the current vertex is
    % always the first one
    for J = 1:size(triadj, 1)
        triadj(J, :) = [v setxor(v,triadj(J, :))];
    end
    
    % compute the lenghts of each triangle side, according to the notation
    % of the SphericalTrig toolbox
    a = distance(lat(triadj(:, 2)), lon(triadj(:, 2)), ...
        lat(triadj(:, 3)), lon(triadj(:, 3)), 'radians');
    b = distance(lat(triadj(:, 1)), lon(triadj(:, 1)), ...
        lat(triadj(:, 2)), lon(triadj(:, 2)), 'radians');
    c = distance(lat(triadj(:, 1)), lon(triadj(:, 1)), ...
        lat(triadj(:, 3)), lon(triadj(:, 3)), 'radians');
    
    % compute the spherical angle adjacent to the current vertex for each
    % triangle, and add them all. For a vertex that is within the convex
    % hull of all neighbours, the sum of angles should be 2*pi (i.e. 360º).
    % We have to allow some tolerance, because due to finite precission and
    % numerical errors the sum will not be exactly 2*pi
    ok(I) = abs(sum(sss(a, b, c)) - 2*pi) < tol;

end
