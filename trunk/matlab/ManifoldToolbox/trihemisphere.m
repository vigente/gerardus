function [tri, latlon] = trihemisphere(N)
% TRIHEMISPHERE  Triangular mesh on hemisphere with roughly uniformly
% distributed points
%
% [TRI, LANLOT] = trihemisphere(N)
%
%   This function uses Anton Semechko's implementation of Reisz s-energy
%   minimisation to find 2*N points approximately uniformly distributed on
%   the sphere and their triangulation. Then, it keeps roughly the lower
%   hemisphere, and moves a bit the points at the top edge so that the top
%   edge has constant latitude. One of the vertices will be at the south
%   pole.
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh.
%
%   LANLOT is a 2-column matrix with the latitude and longitude coordinates
%   (in radians) of the mesh vertices. There will be approximately N
%   points.
%
%   The mesh can be visualised running
%
%     [tri, latlon] = trihemisphere(51);
%    [x, y, z] = sph2cart(latlon(:, 2), latlon(:, 1), 1);
%    hold off
%    trisurf(tri, x, y, z)
%    axis equal

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
narginchk(1, 1);
nargoutchk(0, 2);

% uniformly distribute particles across the surface of the unit sphere
% using Anton Semechko's implementation of Reisz s-energy minimisation
[xyz, tri] = ParticleSampleSphere('N', 2*N);

% rotation matrix that will move the first point to the north pole
A = vec2rotmat(xyz(1, :)');

% rotate all points accordingly
xyz = xyz * A;

% invert sphere so that the point will be at the south pole
xyz(:, 3) = -xyz(:, 3);

% mesh distances between points, only connected or not
d = dmatrix_mesh(tri);

% number of jumps from south pole to every other vertex
d = dijkstra(d, 1);

% number of jumps to the antipodes
nmax = max(d);

% we want to remove the vertices that closer to the north pole than to the
% south pole
idx = find(d > nmax/2);

% remove triangles with any of the rejected vertices
for I = 1:3
    ir = 1;
    while ~isempty(ir) % this loop is necessary because if more than one
                       % vertex to remove, only the last triangle will be
                       % removed
        [~, ir] = intersect(tri(:, I), idx);
        tri(ir, :) = [];
    end
end

while(1)
    
    % remove disconnected vertices
    [tri, xyz] = tri_squeeze(tri, xyz);
    
%     % DEBUG: plot mesh
%     hold off
%     trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
%     axis equal
    
    % find the vertices at the top edge
    trep = TriRep(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3));
    trep = freeBoundary(trep);
    
%     % DEBUG: plot top edge
%     hold on
%     for I = 1:size(trep, 1)
%         plot3(xyz(trep(I, :), 1), xyz(trep(I, :), 2), xyz(trep(I, :), 3), 'k*')
%         plot3(xyz(trep(I, :), 1), xyz(trep(I, :), 2), xyz(trep(I, :), 3), 'k', 'LineWidth', 3)
%     end
    
    % compute the spherical coordinates of the vertices at the top edge
    idx = unique(trep(:));
    [lon, lat, sphrad] = cart2sph(xyz(idx, 1), xyz(idx, 2), xyz(idx, 3));
    
    % move points at the top edge to a median latitude
    lat(:) = median(lat);
    
    % recompute Cartesian coordinates
    [xyz(idx, 1), xyz(idx, 2), xyz(idx, 3)] = sph2cart(lon, lat, sphrad);
    
    % compute area of the triangles
    a = cgal_trifacet_area(tri, xyz);
    
    % find very small triangles on the edge
    idx = find(a / median(a) < .33);
    
    % if no triangles have to be removed, finish
    if isempty(idx)
        break
    end
    
    % remove small triangles
    tri(idx, :) = [];
    
end

% % DEBUG: plot mesh
% hold off
% trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
% axis equal

% format output
[lon, lat] = cart2sph(xyz(:, 1), xyz(:, 2), xyz(:, 3));
latlon = [lat lon];
