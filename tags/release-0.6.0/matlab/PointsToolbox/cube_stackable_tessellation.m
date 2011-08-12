function [tri, v] = cube_stackable_tessellation()
% CUBE_STACKABLE_TESSELLATION  Tessellation in tetrahedra of a cube that
% can be stacked in a regular mesh
%
% TRI = CUBE_STACKABLE_TESSELLATION()
%
%   TRI is a numt-by-4 matrix, as the one produced by Matlab's function
%   DELAUNAY3. Each row of TRI contains the indices of 4 vertices of
%   the cube that form a tetrahedron.
%
%   TRI has 6 rows, i.e. the cube in tessellated by 6 tetrahedra. This
%   tessellation makes contiguous tetrahedra match when cubes are stacked
%   together.
%
%   A cube can also be tessellated by 5 tetrahedra, or 6 in a different
%   configuration, but then the stackability is lost. For more information,
%   see e.g. http://www.iue.tuwien.ac.at/phd/wessner/node32.html.
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
% [TRI, V] = CUBE_STACKABLE_TESSELLATION()
%
%   V is a 3-column matrix with the coordinates of the vertices of a cube
%   in [0, 1].
%
%   You can visualize the tessellation running
%
%   >> [tri, v] = cube_stackable_tessellation();
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
    1, 7, 6, 5; ...
    1, 3, 7, 2; ...
    1, 2, 7, 6; ...
    2, 3, 4, 8; ...
    2, 3, 7, 8; ...
    2, 7, 8, 6 ];

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
