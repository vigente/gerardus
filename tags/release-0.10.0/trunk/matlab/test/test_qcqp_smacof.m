% test_qcqp_smacof.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.0.3
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

%% Toy example. Small triangular mesh.
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
y0 = rand(N, 2) * 4 - 2;

%% Solve toy example using QCQP without constraints. This is equivalent to 
%% using smacof() to solve the problem.

% as the distance matrix is sparse, it is expected that the algorithm will
% minimise the stress to almost 0, but fold-overs can still appear

% SMACOF algorithm parameters
smacof_opts.MaxIter = 100;
smacof_opts.Epsilon = 1e-3;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
scip_opts.display = 'off';
scip_opts.warnings = 'off';
scip_opts.maxiter = 1500;
scip_opts.maxnodes = 10000;
scip_opts.maxtime = 100;
%scip_opts.gamsfile = 'rcasero-problem.gams';

% solve MDS problem with QCQP-SMACOF without constraints
lb = [-4 * ones(N, 1), -4 * ones(N, 1)];
ub = [2 * ones(N, 1), 2 * ones(N, 1)];
[y, stopCondition, sigma, t] ...
    = qcqp_smacof(d, y0, lb, ub, [], [], [], [], [], ...
    smacof_opts, scip_opts);

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

%% Solve same toy example, now with constraints that prevent fold-overs 
%% (i.e. triangles with negative area)

% this will run slowly because we start from a random configuration

% boundaries: maximum area of all the triangles in the mesh
amax = max(trifacet_signed_area(tri, x));

% boundaries: box within which coordinates must be located (we center the
% box around zero, even though the toy example is not)
ymin = [-1 -1] - 3;
ymax = [5 5] - 3;

% matrices and vectors for QCQP (assuming all vertices are free)
[lb, ub, w, A, rl, ru, qc] ...
    =  tri_qcqp_smacof_nofold_2d(tri, ymin, ymax, 0.25, 1.5 * amax);

% SMACOF algorithm parameters
smacof_opts.MaxIter = 100;
smacof_opts.Epsilon = 1e-2;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
scip_opts.display = 'off';
scip_opts.warnings = 'off';
scip_opts.maxiter = 1500;
scip_opts.maxnodes = 10000;
scip_opts.maxtime = 5000;
%opts.gamsfile = 'rcasero-problem.gams';

% solve MDS problem with QCQP-SMACOF with constraints
[y, stopCondition, sigma, t] ...
    = qcqp_smacof(d, y0, lb, ub, w, A, rl, ru, qc, ...
    smacof_opts, scip_opts);

% rigid registration with true solution
[~, y] = procrustes(x, y, 'Scaling', false);

% plot points
subplot(2, 1, 1)
hold off
gplot(d, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

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

%% Toy example where the optimal solution produces triangle overlap, and we
%% use constraints to avoid the fold-over

% toy example without fold-overs
x = [
    0 0
    3 0
    3 1
    0 1
    1 .6
    2 .5
    ];

tri = [
    2 6 1
    2 3 6
    6 3 5
    3 4 5
    5 4 1
    6 5 1
    ];

d = dmatrix_mesh(tri, x);

% compute areas to make sure that we have fold-overs
a = trifacet_signed_area(tri, x);

% plot toy example
subplot(1, 1, 1)
hold off
gplot(d, x)
axis equal

% switch two vertices to create the fold-overs
aux = x(5, :);
x(5, :) = x(6, :);
x(6, :) = aux;
d = dmatrix_mesh(tri, x);

% plot toy example
subplot(2, 1, 1)
hold off
gplot(d, x)
hold on
plot(x(:, 1), x(:, 2), 'o')
axis equal

% compute areas to make sure that we have fold-overs
a = trifacet_signed_area(tri, x);

% boundaries: maximum area of all the triangles in the mesh
amax = max(a);

% boundaries: box within which coordinates must be located
ymin = [-1 -1];
ymax = [4 3];
amin = .25;
amax = 1.5 * amax;

% matrices and vectors for QCQP (assuming all vertices are free)
[lb, ub, w, A, rl, ru, qc] ...
    = tri_qcqp_smacof_nofold_2d(tri, ymin, ymax, amin, amax);

% SMACOF algorithm parameters
smacof_opts.MaxIter = 100;
smacof_opts.Epsilon = 1e-3;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
scip_opts.display = 'off';
scip_opts.warnings = 'off';
scip_opts.maxiter = 1500;
scip_opts.maxnodes = 10000;
scip_opts.maxtime = 100;

% solve MDS problem with QCQP-SMACOF with constraints
% we use the mesh with fold-up as starting point
[y, stopCondition, sigma, t] ...
    = qcqp_smacof(d, x, lb, ub, w, A, rl, ru, qc, ...
    smacof_opts, scip_opts);

% plot solution
subplot(2, 1, 2)
hold off
gplot(d, y)
hold on
plot(y(:, 1), y(:, 2), 'o')
axis equal

stopCondition

% plot stress evolution
subplot(2, 1, 2)
plot(t, sigma)
xlabel('Time (sec)')
ylabel('\sigma')
