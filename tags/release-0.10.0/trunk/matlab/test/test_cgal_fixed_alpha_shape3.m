% test_cgal_fixed_alpha_shape3.m

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
%% tetrahedron

xyz = [
    0 0 0
    1 0 0
    0 1 0
    0 0 1
    ];

% compute alpha shape
tri = cgal_fixed_alpha_shape3(xyz, 1);

% plot meshes
subplot(1, 1, 1)
hold off
trisurf(tri{1}, xyz(:, 1), xyz(:, 2), xyz(:, 3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pyramid with an indentation in the base

xyz = [
    0 0 0
    1 0 0
    0 1 0
    .25 .25 0
    0 0 1
    ];

% compute alpha shape
tri = cgal_fixed_alpha_shape3(xyz, [0 0.5625    1.5469]);

% plot meshes
subplot(1, 2, 1)
hold off
trisurf(tri{2}, xyz(:, 1), xyz(:, 2), xyz(:, 3))
view(78, 40)
subplot(1, 2, 2)
hold off
trisurf(tri{3}, xyz(:, 1), xyz(:, 2), xyz(:, 3))
view(78, 40)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% crescent shape

xyz = [
    1.5 3 0
    3.5 3 0
    0 2 0
    1 2 0
    4 2 0
    5 2 0
    0 1 0
    1 1 0
    4 1 0
    5 1 0
    1 0 0
    4 0 0
    ...
    1.5 3 1
    3.5 3 1
    0 2 1
    1 2 1
    4 2 1
    5 2 1
    0 1 1
    1 1 1
    4 1 1
    5 1 1
    1 0 1
    4 0 1
    ];

% compute minimal alpha shape that creates one connected object
tri = cgal_fixed_alpha_shape3(xyz, 2.5157);
tri = tri{1};

% plot mesh
subplot(1, 1, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))

% compute convex hull
tri = cgal_fixed_alpha_shape3(xyz, Inf);
tri = tri{1};

% plot mesh
subplot(1, 1, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))

% two connected components
tri = cgal_fixed_alpha_shape3(xyz, 1.2656);
tri = tri{1};

% plot mesh
subplot(1, 1, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
