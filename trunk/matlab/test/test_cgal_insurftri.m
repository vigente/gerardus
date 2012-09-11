% test_cgal_insurftri.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
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

% create a mesh for a cube
x = [...
    0 0 0; ...
    0 0 1; ...
    0 1 0; ...
    0 1 1; ...
    1 0 0; ...
    1 0 1; ...
    1 1 0; ...
    1 1 1; ...
    ]/2;

y = x(:, 2);
z = x(:, 3);
x = x(:, 1);

tri = [...
    1 2 4; ...
    1 3 4; ...
    2 4 6; ...
    4 6 8; ...
    7 8 6; ...
    7 6 5; ...
    7 5 3; ...
    1 3 5; ...
    3 4 7; ...
    7 4 8; ...
    1 2 5; ...
    2 5 6; ...
    ];

% plot cube
trisurf(tri, x, y, z, 1:12, 'EdgeColor', 'none')
axis xy equal

%% test individual arbitrary points inside and outside the cube

% with the default directions, the last test point fails, because as it's
% in the middle of the cube, the rays go through vertices/edges
cgal_insurftri(tri, [x y z], ...
    [0 0 0; ... % on vertex, a ray will be on adjacent edge
    .5 .5 .5; ... % on vertex, a ray won't touch anything else from surface
    1 0 0; ... % outside, a ray won't touch surface
    -1 0 0; ... % outside, a ray on edge
    .25 .3 .13; ... % inside, off to one side
    .25 .25 .25; ... % in the middle of the cube
    ])

% choose four random ray directions
directions = rand(5, 3);
cgal_insurftri(tri, [x y z], ...
    [0 0 0; ... % on vertex, a ray will be on adjacent edge
    .5 .5 .5; ... % on vertex, a ray won't touch anything else from surface
    1 0 0; ... % outside, a ray won't touch surface
    -1 0 0; ... % outside, a ray on edge
    .25 .3 .13; ... % inside, off to one side
    .25 .25 .25; ... % in the middle of the cube
    .1 .25 0; ... % on a facet, away from an edge
    ], ...
    directions)

%% test a rectangular volume

% choose four random ray directions
directions = rand(5, 3);
cgal_insurftri(tri, [x y z], ...
    {linspace(-.25, .25, 5), linspace(-.25, .25, 4), linspace(-.25, .25, 3)}, ...
    directions)
