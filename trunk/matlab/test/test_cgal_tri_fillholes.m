% test_cgal_tri_fillholes.m

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
%% Basic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empty mesh
[tri, n] = cgal_tri_fillholes([], [])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cube with one missing triangle
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
trisurf(tri, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

% make a hole by removing one of the triangles
tri2 = tri;
tri2(6, :) = [];

% plot mesh
hold off
trisurf(tri2, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

% fill hole
[tri2, n] = cgal_tri_fillholes(tri2, x);
disp(['Number of holes filled = ' num2str(n)])

% check that the hole has been filled (expected result, and empty matrix)
setdiff(sort(tri, 2), sort(tri2, 2), 'rows')

% plot mesh
hold off
trisurf(tri2, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cube with two triangles that form the upper face missing
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
trisurf(tri, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

% make a hole by removing one of the triangles
tri2 = tri;
tri2([6 12], :) = [];

% plot mesh
hold off
trisurf(tri2, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

% fill hole
[tri2, n] = cgal_tri_fillholes(tri2, x);
disp(['Number of holes filled = ' num2str(n)])

% check that the hole has been filled (expected result, and empty matrix)
setdiff(sort(tri, 2), sort(tri2, 2), 'rows')

% plot mesh
hold off
trisurf(tri2, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cube with two holes (missing two triangles)
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
trisurf(tri, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

% make a hole by removing one of the triangles
tri2 = tri;
tri2([1 6], :) = [];

% plot mesh
hold off
trisurf(tri2, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

% fill hole
[tri2, n] = cgal_tri_fillholes(tri2, x);
disp(['Number of holes filled = ' num2str(n)])

% check that the hole has been filled (expected result, and empty matrix)
setdiff(sort(tri, 2), sort(tri2, 2), 'rows')

% plot mesh
hold off
trisurf(tri2, x(:,1), x(:,2), x(:,3))
axis([0 1 0 1 0 1])

