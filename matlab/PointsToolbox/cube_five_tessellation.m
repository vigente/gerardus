function [tri, v] = cube_five_tessellation()
% CUBE_FIVE_TESSELLATION  Tessellation in 5 tetrahedra of a cube with
% perpendicular diagonals on opposite sides of the cube (non-stackable)
%
% TRI = CUBE_FIVE_TESSELLATION()
%
%   TRI is a 5-by-4 matrix, as the one produced by Matlab's function
%   DELAUNAY3. Each row of TRI contains the indices of 4 vertices of
%   the cube that form a tetrahedron.
%
%   TRI has 5 rows, i.e. the cube in tessellated by 5 tetrahedra. This
%   tessellation makes contiguous tetrahedra _not_ match when cubes are
%   stacked together. 4 of the tetraheadra have volume 1/6, while the one
%   "in the middle" has volume 2/6.
%
%   This tessellation makes opposite sides of the cube be split by
%   diagonals that are perpendicular to each other.
%
%   For more information, see e.g. 
%   http://www.iue.tuwien.ac.at/phd/wessner/node32.html.
%
%   The cube's vertices have been labelled as:
%
%            6 ------------ 8
%            / |          /|
%           /  |         / |
%          /   |        /  |
%         /    |       /   |
%        5 ----------- 7   |
%         |    |------|----|
%         |   / 2     |   / 4
%         |  /        |  /
%         | /         | /
%         |/          |/
%        1 ----------- 3
%
% [TRI, V] = CUBE_FIVE_TESSELLATION()
%
%   V is a 3-column matrix with the coordinates of the vertices of a cube
%   in [0, 1].
%
%   You can visualize the tessellation running
%
%   >> [tri, v] = cube_five_tessellation();
%   >> tetramesh( tri, v )

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009 University of Oxford
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

tri = [ ...
    1, 3, 4, 7; ...
    1, 4, 6, 7; ...
    4, 6, 7, 8; ...
    1, 5, 6, 7; ...
    1, 2, 4, 6 ];

if ( nargout > 1 )
v = [ ...
    0     0     0; ...
    0     1     0; ...
    1     0     0; ...
    1     1     0; ...
    0     0     1; ...
    0     1     1; ...
    1     0     1; ...
    1     1     1 ];
end
