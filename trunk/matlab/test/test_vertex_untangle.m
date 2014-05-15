% test_vertex_untangle.m

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create data neighbourhoods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Untangled vertex, to check that function will not try to untangle it

x = [
    0 0
    1 0
    1 1
    0 1
    0.25 0.25 % central vertex
    ];

tri = [
    1 2 5
    2 3 5
    3 4 5
    4 1 5
    ];


%% Very basic square neighbourhood

x = [
    0 0
    1 0
    1 1
    0 1
    0.75 1.25 % central vertex
    ];

tri = [
    1 2 5
    2 3 5
    3 4 5
    4 1 5
    ];

%% More complicated; one large triangle with a small rectangular bit 
%% attached. Solution exists

x = [
    0.0086    0.0662
    0.0122    0.0653
    0.0065    0.0342
    0.0181    0.0263
    0.0024   -0.0578
    -0.0305   -0.1080
    ];

tri = [
    4     3     1
    4     1     2
    3     4     6
    2     5     4
    4     5     6
    ];

%% Knife shape. Solution exists

x = [
    -0.0113   -0.0679
    -0.0093   -0.0465
    -0.0046   -0.0083
    -0.0006   -0.0103
    -0.0144   -0.0174
    0.0111    0.0251
    0.0247    0.1172
    ];

tri = [
    6     4     3
    3     7     6
    5     7     3
    4     2     3
    3     1     5
    2     1     3
    ];

%% Knife shape, folded over to create fixed point overlap: No solution exists

x = [
    -0.0113   -0.0679
    -0.0093   -0.0465
    -0.0046   -0.0083
    -0.0006   -0.0103
    -0.0144   -0.0174
    0.0111    0.0251
    0.0247    -0.02
    ];

tri = [
    6     4     3
    3     7     6
    5     7     3
    4     2     3
    3     1     5
    2     1     3
    ];

%% Polygon with no existing solution

x = [
    0 0
    6 0
    6 1
    3.5 1
    3.5 5
    6 5
    6 6
    0 6
    0 5
    2.5 5
    2.5 1
    0 1
    3 3
    ];

tri = [
    1 2 13
    2 3 13
    3 4 13
    4 5 13
    5 6 13
    6 7 13
    7 8 13
    8 9 13
    9 10 13
    10 11 13
    11 12 13
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Untangling code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% local distance/adjacency matrix
d = dmatrix_mesh(tri);

% index of the central vertex and neighbours
idx0 = find(sum(d ~= 0) == size(d, 1)-1);
idxn = [1:idx0-1 idx0+1:size(d, 1)];

% plot tangled neighbourhood
subplot(1, 2, 1)
hold off
gplot(d, x)
hold on
plot(x(idx0, 1), x(idx0, 2), 'ro')
gplot(d(idxn, idxn), x(idxn, :), 'r')

% untangle vertex
[x(idx0, :), exitflag] = vertex_untangle(tri, x);
exitflag

% plot untangled neighbourhood
subplot(1, 2, 2)
hold off
gplot(d, x)
hold on
plot(x(idx0, 1), x(idx0, 2), 'ro')
gplot(d(idxn, idxn), x(idxn, :), 'r')

if (exitflag < 1)
    error('Assertion error: A valid solution expected, but the algorithm did not converge')
end
