% test_tri_qcqp_smacof_nofold_2d_pip.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.0.1
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

%% Small toy mesh

% 5 points in the perimeter, 2 free vertices, possible untangled solution
x = [
    0 3
    1 0
    3 1
    4 .5
    3.5 4
    1.5 2
    3 2
    ];

tri = [
    6 1 2
    6 2 3
    7 3 4
    4 5 7
    5 1 7
    7 1 6
    7 6 3
    ];

% maximum area of all the triangles in the mesh
amax = max(trifacet_signed_area(tri, x));

% plot mesh
d = dmatrix_mesh(tri);
hold off
gplot(d, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

%% typical case: all vertices are free vertices (i.e. all vertices are 
%% unknowns)

ymin = [-1 -0.5];
ymax = [5 5.5e7];
amin = 0.25;
amax = 2 * amax;

% constraints and bounds in PIP format
[con, bnd] = tri_qcqp_smacof_nofold_2d_pip(tri, ymin, ymax, amin, amax);

% view constraints
char(con)

% view bounds
char(bnd)
