function isTangled = sphtri_vertex_istangled(tri, latlon, vidx)
% SPHTRI_VERTEX_ISTANGLED Find tangled vertices in a spherical triangular
% mesh.
%
% ISTANGLED = sphtri_vertex_istangled(TRI, LATLON)
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
%   ISTANGLED is a boolean vector with one element per vertex in LATLON.
%   ISTANGLED(i)==true means that the corresponding vertex is tangled. We
%   define a vertex as tangled if any of its incident triangles has a
%   negative area (i.e. if the triangle normal points inwards).

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
    isTangled = [];
    return;
end

% initialize output
isTangled = false(size(vidx));

% iterate vertices to check whether they are tangled or not
for I = 1:numel(vidx)
    
    % list of triangles that are adjacent to the current vertex
    [triidx, ~] = find(ismember(tri, vidx(I)));
    
    % extract local neighbourhood from the whole mesh
    [triloc, latlonloc] = tri_squeeze(tri(triidx, :), latlon);
    
    % centroid of the local neighbourhood
    latlonm = meanm(latlonloc(:, 1), latlonloc(:, 1), 'radians');
    
    % center local neighbourhood around (lat=0, lon=0). This will avoid
    % wrap-around problems around lon=pi
    latlonloc(:, 1) = latlonloc(:, 1) - latlonm(1);
    latlonloc(:, 2) = latlonloc(:, 2) - latlonm(2);
    
    % compute signed area of each triangle in the local neighbourhood
    a = trifacet_signed_area(triloc, latlonloc);
    
    % the vertex is entangled if any of its adjacent triangles has a
    % negative area
    isTangled(I) = any(a<0);
    
    % DEBUG: plot local neighbourhood
    if (DEBUG)
        trisurf(triloc, latlonloc(:, 1), latlonloc(:, 2), ...
            zeros(size(latlonloc, 1), 1), uint8(a > 0));
        view(2)
    end
    
end
