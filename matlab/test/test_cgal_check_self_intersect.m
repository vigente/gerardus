% test_cgal_check_self_intersect.m

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
%% nice mesh without self intersections

xyz = [
    -1     1    0
     0     1    0
    -0.5   0.5  0
    -1     0    0
     0     0    0
     1     0    0
     0    -1    0
    ];

x = xyz(:, 1);
y = xyz(:, 2);
z = xyz(:, 3);

tri = [
    1 2 3
    1 4 3
    3 2 5
    3 4 5
    2 5 6
    4 5 7
    5 7 6
    ];

% interior vertices that are nearest neighbors
I = 3;
J = 5;

% plot mesh, with each face in a different colour
hold off
gplot3d(dmatrix_mesh(tri), xyz)
hold on
plot3(x, y, z, 'o')
axis equal xy
view(2)

% check triangle intersections
cross = cgal_check_self_intersect(tri, xyz);
cross'
% expected: all zeros

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% move vertex I down so that triangle 1 overlaps with triangles 4 and 6

y(I) = -1/4.5;
xyz = [x y z];

% plot mesh, with each face in a different colour
hold off
gplot3d(dmatrix_mesh(tri), xyz)
hold on
plot3(x, y, z, 'o')
axis equal xy
view(2)

% check triangle intersections
cross = cgal_check_self_intersect(tri, xyz);
cross'
% expected: 2     2     2     4     0     4     0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spherical mesh

% uniformly distribute particles across the surface of the unit sphere
% using Anton Semechko's implementation of Reisz s-energy minimisation
% [xyz, tri] = ParticleSampleSphere('N', 20);

% use always the same points for reproducibility
xyz = [
    0.1865    0.1323   -0.9735
    0.8520    0.0405   -0.5219
    0.8562    0.4596    0.2358
    0.5982   -0.1356    0.7897
    0.2115    0.5448    0.8115
    0.4822    0.7204   -0.4986
   -0.9091    0.3984   -0.1217
    0.1588    0.9741    0.1609
   -0.5982    0.1356   -0.7897
   -0.5127    0.6897    0.5113
   -0.3189    0.8227   -0.4706
   -0.8520   -0.0405    0.5219
   -0.1866   -0.1323    0.9735
    0.5127   -0.6897   -0.5113
   -0.1588   -0.9741   -0.1609
   -0.2115   -0.5448   -0.8115
   -0.4822   -0.7204    0.4986
   -0.8562   -0.4596   -0.2358
    0.3189   -0.8227    0.4706
    0.9091   -0.3984    0.1217    
];

tri = [
     7    10    11
     8    11    10
     1    11     6
     1     6     2
     8     6    11
    20     2     3
    20     3     4
     2     6     3
     8     3     6
     8    10     5
     8     5     3
     4     3     5
    13     5    10
    13     4     5
    20    14     2
     1     2    14
    20     4    19
    20    19    14
    15    14    19
    13    19     4
     7    11     9
     1     9    11
     7    12    10
    13    10    12
    15    16    14
     1    14    16
     1    16     9
    15    18    16
     9    16    18
     7     9    18
     7    18    12
    15    19    17
    15    17    18
    12    18    17
    13    17    19
    13    12    17  
     ];

% check triangle intersections
cross = cgal_check_self_intersect(tri, xyz);
cross'

% plot mesh, colouring the facets according to whether they are overlapping
% or not
subplot(1, 2, 1)
hold off
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3), cross);
axis equal

%% create some triangle overlap

% select one vertex
I = 17;

% select a nearest neighbour
J = 19;

% plot them on the mesh
hold on
plot3(xyz(I, 1), xyz(I, 2), xyz(I, 3), 'ro')
plot3(xyz(J, 1), xyz(J, 2), xyz(J, 3), 'r*')

% move the nearest neighbor in a way such that I is still within the convex
% hull, but the nearest neighbour causes overlap
lat(J) = 13/180*pi;
lon(J) = -120/180*pi;

[xyz(J, 1), xyz(J, 2), xyz(J, 3)] = sph2cart(lon(J), lat(J), 1);

% check triangle intersections
cross = cgal_check_self_intersect(tri, xyz);
cross'

% triangles that contain the vertex we moved
find(sum(tri == J, 2)>0)
find(cross)

% plot mesh, colouring the facets according to whether they are overlapping
% or not
subplot(1, 2, 2)
hold off
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3), double(cross~=0));
axis equal
hold on
plot3(xyz(I, 1), xyz(I, 2), xyz(I, 3), 'ro')
plot3(xyz(J, 1), xyz(J, 2), xyz(J, 3), 'r*')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% two 3D triangles orthogonal and intersecting, intersection=segment

xyz = [
    -1 -1 0   % point for horizontal triangle
    1 -1 0    % point for horizontal triangle
    -1 1 0    % point for horizontal triangle
    -.5 0 1   % point for vertical triangle
    .5 0 1    % point for vertical triangle
    -.5 0 -1  % point for vertical triangle
    ];

tri = [
    1 2 3     % horizontal
    4 5 6     % vertical
    ];


% check triangle intersections
cross = cgal_check_self_intersect(tri, xyz);
cross'

% plot mesh, colouring the facets according to whether they are overlapping
% or not
subplot(1, 1, 1)
hold off
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3), cross);
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% two 3D triangles orthogonal and intersecting, intersection=point

xyz = [
    -1 -1 0   % point for horizontal triangle
    1 -1 0    % point for horizontal triangle
    -1 1 0    % point for horizontal triangle
    -.5 0 2   % point for vertical triangle
    .5 0 2    % point for vertical triangle
    -.5 0 0  % point for vertical triangle
    ];

tri = [
    1 2 3     % horizontal
    4 5 6     % vertical
    ];


% check triangle intersections
cross = cgal_check_self_intersect(tri, xyz);
cross'

% plot mesh, colouring the facets according to whether they are overlapping
% or not
subplot(1, 1, 1)
hold off
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3), cross);
axis equal

