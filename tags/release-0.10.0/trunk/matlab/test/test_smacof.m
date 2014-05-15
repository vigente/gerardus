% test_smacof.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.0.2
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

%% Toy example. Full distance matrix.

% because the distance matrix is full and Euclidean distances are
% congruent, the solution can be the same as obtained with classical MDS,
% (almost identical to the ground truth). But this is not always the case,
% because depending on Y0, the algorithm can fall into a local minimum. I
% think that this is because X is quite small. Testing random
% configurations X with more points, this problem does not happen, and the
% global maximum is found.

x = [
    0 3
    1 0
    3 1
    4 .5
    3.5 4
    1.5 2
    3 2
    ];

N = size(x, 1);

% full Euclidean distance matrix
d = dmatrix(x');

% random starting configuration
y0 = rand(size(x)) * 4;

% plot points
subplot(2, 1, 1)
hold off
gplot(d, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

% SMACOF algorithm parameters
opts.MaxIter = 500;
opts.Epsilon = 1e-6;
opts.Display = 'iter';
opts.TolFun = 1e-6;

% solve MDS problem with SMACOF
[y, stopCondition, sigma, t] = smacof(d, y0, [], opts);

% rigid registration with true solution
[~, y] = procrustes(x, y, 'Scaling', false);

% plot solution
gplot(d, y, 'r')
plot(y0(:, 1), y0(:, 2), '.k')
plot(y(:, 1), y(:, 2), 'xr')

stopCondition

% plot stress evolution
subplot(2, 1, 2)
plot(t, sigma)
xlabel('Time (sec)')
ylabel('\sigma')

%% Toy example. Sparse distance matrix.
% because of the sparsity of the distance matrix, the solution will usually
% have fold-overs

x = [
    0 3
    1 0
    3 1
    4 .5
    3.5 4
    1.5 2
    3 2
    ];

N = size(x, 1);

tri = [
    6 1 2
    6 2 3
    7 3 4
    4 5 7
    5 1 7
    7 1 6
    7 6 3
    ];

% sparse Euclidean distance matrix
d = dmatrix_mesh(tri, x);

% plot points
subplot(2, 1, 1)
hold off
gplot(d, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

% random starting configuration
y0 = rand(7, 2) * 4;

% SMACOF algorithm parameters
opts.MaxIter = 500;
opts.Epsilon = 1e-3;
opts.Display = 'iter';
opts.TolFun = 1e-6;

% solve MDS problem with SMACOF
[y, stopCondition, sigma, t] = smacof(d, y0, [], opts);

% rigid registration with true solution
[~, y] = procrustes(x, y, 'Scaling', false);

% plot solution
gplot(d, y, 'r')
plot(y0(:, 1), y0(:, 2), '.k')
plot(y(:, 1), y(:, 2), 'xr')

stopCondition

% plot stress evolution
subplot(2, 1, 2)
plot(t, sigma)
xlabel('Time (sec)')
ylabel('\sigma')

