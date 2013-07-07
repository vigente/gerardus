% test_dmatrix_sphmesh.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.1.1
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

%% spherical triangle with vertices in the north pole, equator at lon=0, 
%% and equator at lon=90º
uv = [
    90 0
    0  0
    0  90
    ] / 180 * pi;
tri = [1 2 3];

% compute distance matrix
d = dmatrix_sphmesh(tri, uv);

% ans =
% 
%          0    1.5708    1.5708
%     1.5708         0    1.5708
%     1.5708    1.5708         0

%% add another spherical triangle, this one with a common edge and touching
%% the south pole
uv = [
    90 0
    0  0
    0  90
    -90 0
    ] / 180 * pi;
tri = [
    1 2 3
    2 3 4
    ];

% plot
[x, y, z] = sph2cart(uv(:, 2), uv(:, 1), 1);
hold off
trisurf(tri, x, y, z)
axis equal

% compute distance matrix
d = dmatrix_sphmesh(tri, uv);
full(d)

% ans =
% 
%          0    1.5708    1.5708         0
%     1.5708         0    1.5708    1.5708
%     1.5708    1.5708         0    1.5708
%          0    1.5708    1.5708         0
         