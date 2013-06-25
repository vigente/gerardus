% test_check_sphtrisurf_vertices.m

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vertex surrounded by other 5 vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% uv = [lat lon] = [elevation azimuth]
uv = [
     90    0
     70    0
     70  100
     70  120
     70  180
     70 -160
    ] / 180 * pi;
lat = uv(:, 1);
lon = uv(:, 2);

% mesh
tri = [
    1 2 3
    1 3 4
    1 4 5
    1 5 6
    1 6 2
    ];

% compute xyz coordinates of the mesh points
[x, y, z] = sph2cart(lon, lat, 1);

% plot mesh
d = dmatrix_mesh([x y z], tri);
gplot3d(d, [x y z])

% check that the spherical triangles around the central point add up to
% 360º
ok = check_sphtrisurf_vertices(tri, uv, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vertex shifted creating triangle overlap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% uv = [lat lon] = [elevation azimuth]
uv = [
     60    50
     70    0
     70  100
     70  120
     70  180
     70 -160
    ] / 180 * pi;
lat = uv(:, 1);
lon = uv(:, 2);

% mesh
tri = [
    1 2 3
    1 3 4
    1 4 5
    1 5 6
    1 6 2
    ];

% compute xyz coordinates of the mesh points
[x, y, z] = sph2cart(lon, lat, 1);

% plot mesh
d = dmatrix_mesh([x y z], tri);
gplot3d(d, [x y z])

% check that the spherical triangles around the central point add up to
% 360º
ok = check_sphtrisurf_vertices(tri, uv, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bunch of random points, no overlaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% random distribution of points on the sphere
lat = (rand(40, 1) - .5) * pi;
lon = (rand(40, 1) - .5) * 2 * pi;
uv = [lat lon];

% compute xyz coordinates of the mesh points
[x, y, z] = sph2cart(lon, lat, 1);

% mesh with a delaunay triangulation
tri = DelaunayTri([x y z]);
tri = freeBoundary(tri);

% plot mesh
d = dmatrix_mesh([x y z], tri);
gplot3d(d, [x y z])

% check that no vertices produce overlap of triangles
ok = check_sphtrisurf_vertices(tri, uv);
ok'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% move a vertex a lot in the previous example, and see whether we can 
%% identify it as producing overlap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% move one of the vertices
lat(1) = lat(1) + pi/2;
lon(1) = lon(1) + pi/2;
uv = [lat lon];

% recompute xyz coordinates of the mesh points
[x, y, z] = sph2cart(lon, lat, 1);

% replot mesh
d = dmatrix_mesh([x y z], tri);
gplot3d(d, [x y z])
hold on
plot3(x(1), y(1), z(1), 'ro')

% see whether the overlap is detected
ok = check_sphtrisurf_vertices(tri, uv);
ok'

% highlight vertices affected by the overlap (we expect to see all vertices
% connected to the one we moved)
plot3(x(~ok), y(~ok), z(~ok), 'k*')

