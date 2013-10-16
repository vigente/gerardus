% test_remove_vertex_from_tri.m

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
%% Remove opposite vertices from a cube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a surface mesh that is a cube
x = [
    0 0 0
    0 0 1
    0 1 0
    0 1 1
    1 0 0
    1 0 1
    1 1 0
    1 1 1
    ];

tri = [
    1     3     5
    2     1     5
    2     3     1
    4     3     2
    6     2     5
    4     2     6
    7     6     5
    3     7     5
    4     7     3
    8     6     7
    4     8     7
    4     6     8
    ];

% plot mesh
hold off
trisurf(tri, x(:,1), x(:,2), x(:,3), ones(1, size(x, 1)))
axis equal

% vertices to remove
idx = [2 7];

% plot vertices that are going to be removed
hold on
plot3(x(idx, 1), x(idx, 2), x(idx, 3), 'ro')

% remove vertices
[tri2, x2] = remove_vertex_from_tri(tri, x, idx);

% plot result
hold off
trisurf(tri2, x2(:,1), x2(:,2), x2(:,3), ones(1, size(x2, 1)))
hold on
plot3(x(idx, 1), x(idx, 2), x(idx, 3), 'ro')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove neighbour vertices from a cube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a surface mesh that is a cube
x = [
    0 0 0
    0 0 1
    0 1 0
    0 1 1
    1 0 0
    1 0 1
    1 1 0
    1 1 1
    ];

tri = [
    1     3     5
    2     1     5
    2     3     1
    4     3     2
    6     2     5
    4     2     6
    7     6     5
    3     7     5
    4     7     3
    8     6     7
    4     8     7
    4     6     8
    ];

% plot mesh
hold off
trisurf(tri, x(:,1), x(:,2), x(:,3), ones(1, size(x, 1)))
axis equal

% vertices to remove
idx = [1 2];

% plot vertices that are going to be removed
hold on
plot3(x(idx, 1), x(idx, 2), x(idx, 3), 'ro')

% remove vertices
[tri2, x2] = remove_vertex_from_tri(tri, x, idx);

% plot result
hold off
trisurf(tri2, x2(:,1), x2(:,2), x2(:,3), ones(1, size(x2, 1)))
hold on
plot3(x(idx, 1), x(idx, 2), x(idx, 3), 'ro')
axis equal
