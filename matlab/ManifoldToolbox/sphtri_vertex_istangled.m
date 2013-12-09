function nTangled = sphtri_vertex_istangled(tri, latlon, vidx)
% SPHTRI_VERTEX_ISTANGLED Find tangled vertices in a spherical triangular
% mesh.
%
% NTANGLED = sphtri_vertex_istangled(TRI, LATLON)
%
%   TRI is a 3-column matrix. Each row represents the indices of the three
%   vertices that form a triangle. TRI as a whole represents the closed
%   spherical surface. Usually, TRI will come from a non-spheric closed
%   surface mesh, where the triangles are correctly oriented (i.e. normals
%   pointing outwards and positive areas). The spherical parameterization
%   of the vertices can flip triangles (normals pointing inwards and
%   negative areas).
%
%   LATLON is a 2-column matrix with the latitude and longitude coordinates
%   of the vertices' coordinates.
%
%   NTANGLED is a count vector with one element per vertex in LATLON.
%   NTANGLED(i) is the number of adjacent triangles with an inward pointing
%   normal, i.e. triangles that contribute to the tanglement. NTANGLED(i)>0
%   means that the vertex is tangled.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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

DEBUG = false;

% check arguments
narginchk(2, 3);
nargoutchk(0, 1);

% number of vertices
Nv = size(latlon, 1);

% defaults
if (nargin < 3)
    vidx = 1:Nv;
end

if (isempty(vidx))
    nTangled = [];
    return;
end

% initialize output
nTangled = zeros(size(vidx));

% compute Cartesian coordinates of the sphere points
x = zeros(size(latlon, 1), 3);
[x(:, 1), x(:, 2), x(:, 3)] = sph2cart(latlon(:, 2), latlon(:, 1), 1);

% iterate vertices to check whether they are tangled or not
for I = 1:numel(vidx)
    
    % list of triangles that are adjacent to the current vertex
    [triidx, ~] = find(ismember(tri, vidx(I)));
    
    % extract local neighbourhood from the whole mesh
    [triloc, xloc] = tri_squeeze(tri(triidx, :), x);
    
    % axis-angle representation of rotation that takes the current vertex
    % to [1, 0, 0]
    r = [cross(x(vidx(I), :), [1 0 0]) -acos(dot(x(vidx(I), :), [1 0 0]))];
    
    % 3D-rotation matrix representation of the axis-angle
    R = vrrotvec2mat(r);
    
    % rotate local neighbourhood so that the central voxel is on [1, 0, 0]
    xloc = xloc * R;
    
    % convert to latitude/longitude the centered local neighbourhood
    [lonloc, latloc] = cart2sph(xloc(:, 1), xloc(:, 2), xloc(:, 3));

    % compute signed area of each triangle in the local neighbourhood
    a = trifacet_signed_area(triloc, [lonloc latloc]);
    
    % the vertex is entangled if any of its adjacent triangles has a
    % negative area. Count the number of adjacent negative areas as a
    % measure of entanglement
    nTangled(I) = nnz(a<0);
    
    % DEBUG: plot local neighbourhood
    if (DEBUG)
        d = dmatrix_mesh(triloc);
        gplot(d, [lonloc, latloc])
    end
    
end
