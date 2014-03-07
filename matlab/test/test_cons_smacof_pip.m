% test_cons_smacof_pip.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.2.4
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PLANE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Toy example, small triangulated rectangular grid. This is a basic sanity
%% check. Here, we use the ground truth as the initial estimate, and we use
%% the full matrix of true distances. The algorithm should realize it's
%% already at an optimum (the global optimum, in fact), and that all
%% constraints are true
x = [
    0 0
    1 0
    2 0 
    3 0
    0 0.5
    1 0.5
    2 0.5
    3 0.5
    0 1
    1 1
    2 1 
    3 1
    0 1.5
    1 1.5
    2 1.5
    3 1.5
    ];

N = size(x, 1);
tri = delaunay(x);
w = dmatrix_mesh(tri);

% areas of the triangles
trifacet_signed_area(tri, x)

% plot triangulation
hold off
gplot(w, x)

% compute full matrix of true distances
dtot = dmatrix(x');

% boundaries: maximum area of all the triangles in the mesh
amin = 0.1;
amax = 0.4;

% boundaries: box within which coordinates must be located
ymin = [-1 -1];
ymax = [4 2.5];

% constraints and boundaries for QCQP (assuming all vertices are free)
[con, bnd] = tri_qcqp_smacof_nofold_2d_pip(tri, ymin, ymax, amin, amax);

% SMACOF algorithm parameters
smacof_opts.MaxIter = 100;
smacof_opts.Epsilon = 1e-2;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
scip_opts.limits_time = 100;
%scip_opts.limits_gap = 8/100;
%scip_opts.limits_solutions = 5;
scip_opts.display_verblevel = 4;

% solve MDS problem with constrained SMACOF
[y, stopCondition, sigma, t] ...
    = cons_smacof_pip(dtot, x, [], bnd, [], con, ...
    smacof_opts, scip_opts);

% check that this is a valid solution
if (all(trifacet_signed_area(tri, y) >= 0))
    disp('Valid solution')
else
    error('Invalid solution')
end

% rigid registration with true solution
[~, y] = procrustes(x, y, 'Scaling', false);

% plot points
subplot(2, 1, 1)
hold off
gplot(w, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

% plot solution
gplot(w, y, 'r')
plot(y(:, 1), y(:, 2), 'xr')

stopCondition

% plot stress evolution
subplot(2, 1, 2)
plot(t, sigma)
xlabel('Time (sec)')
ylabel('\sigma')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same toy example, but now we start from a random distribution of points

rng(0); % this initialization produces a very good solution
rng(1); % this initialization produces a poor solution with fold-overs,
        % that at the same time has all areas positive

% boundaries: maximum area of all the triangles in the mesh
amin = 0.1;
amax = 0.4;

% constraints and boundaries for QCQP (assuming all vertices are free)
[con, bnd] = tri_qcqp_smacof_nofold_2d_pip(tri, ymin, ymax, amin, amax);

% random initialization
y0 = rand(N, 2);
y0(:, 1) = y0(:, 1) * 3;
y0(:, 2) = y0(:, 2) * 1.5;

% SMACOF algorithm parameters
smacof_opts.MaxIter = 100;
smacof_opts.Epsilon = 1e-2;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
scip_opts.limits_time = 100;
scip_opts.limits_solutions = 1;
scip_opts.display_verblevel = 0;

% solve MDS problem with constrained SMACOF
[y, stopCondition, sigma, t] ...
    = cons_smacof_pip(dtot, y0, [], bnd, [], con, ...
    smacof_opts, scip_opts);

% check that this is a valid solution
if (all(trifacet_signed_area(tri, y) >= 0))
    disp('Valid solution')
else
    error('Invalid solution')
end

% rigid registration with true solution
[~, y] = procrustes(x, y, 'Scaling', false);

% plot points
subplot(2, 1, 1)
hold off
gplot(w, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

% plot solution
gplot(w, y, 'r')
plot(y(:, 1), y(:, 2), 'xr')

stopCondition

% plot stress evolution
subplot(2, 1, 2)
plot(t, sigma)
xlabel('Time (sec)')
ylabel('stress')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same as before, but now we move some vertices and recompute the distance
%% matrix. Our goal is to create an MDS solution that creates fold-overs.
%% Then, we can test whether the constraints protect the constrained SMACOF
%% solution from those fold-overs

rng(0);

% create fold-overs
x(10, :) = [2.25 0.7];
x(11, :) = [2.75 1.4];

% recompute distance matrix
dtot = dmatrix(x');

% SMACOF algorithm parameters
smacof_opts.MaxIter = 100;
smacof_opts.Epsilon = 1e-2;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
scip_opts.limits_time = 100;
%scip_opts.limits_gap = 8/100;
scip_opts.limits_solutions = 1;
scip_opts.display_verblevel = 0;

% solve MDS problem with constrained SMACOF
tic
[y, stopCondition, sigma, t] ...
    = cons_smacof_pip(dtot, x, isFree, bnd, [], con, ...
    smacof_opts, scip_opts);
toc

% plot points
subplot(2, 1, 1)
hold off
gplot(w, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

% plot solution
subplot(2, 1, 2)
gplot(w, y, 'r')
hold on
plot(y(:, 1), y(:, 2), 'xr')

stopCondition

% plot stress evolution
subplot(2, 1, 2)
hold off
plot(t, sigma)
xlabel('Time (sec)')
ylabel('stress')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same toy example, using fixed points

x = [
    0 0
    1 0
    2 0 
    3 0
    0 0.5
    1 0.5
    2 0.5
    3 0.5
    0 1
    1 1
    2 1 
    3 1
    0 1.5
    1 1.5
    2 1.5
    3 1.5
    ];

N = size(x, 1);
tri = delaunay(x);
w = dmatrix_mesh(tri);

% areas of the triangles
trifacet_signed_area(tri, x)

% select free vertices
isFree = true(size(x, 1), 1);
isFree([1 4 7 14 15]) = false;

% plot triangulation
hold off
gplot(w, x)
hold on
plot(x(~isFree, 1), x(~isFree, 2), 'ro')

% compute full matrix of true distances
dtot = dmatrix(x');

% boundaries: maximum area of all the triangles in the mesh
amin = 0.1;
amax = 0.4;

% boundaries: box within which coordinates must be located
ymin = [-1 -1];
ymax = [4 2.5];

% initialize
y0 = x;
y0(isFree, :) = y0(isFree, :) + rand(nnz(isFree), 2);

% constraints and boundaries for QCQP (assuming all vertices are free)
[con, bnd] = tri_qcqp_smacof_nofold_2d_pip(tri, ymin, ymax, amin, amax, ...
    isFree, y0);

% SMACOF algorithm parameters
smacof_opts.MaxIter = 100;
smacof_opts.Epsilon = 1e-3;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
scip_opts.limits_time = 100;
%scip_opts.limits_gap = 8/100;
scip_opts.limits_solutions = 1;
scip_opts.display_verblevel = 0;

% solve MDS problem with constrained SMACOF
[y, stopCondition, sigma, t] ...
    = cons_smacof_pip(dtot, y0, isFree, bnd, [], con, ...
    smacof_opts, scip_opts);

% check that this is a valid solution
if (all(trifacet_signed_area(tri, y) >= 0))
    disp('Valid solution')
else
    error('Invalid solution')
end

% rigid registration with true solution
[~, y] = procrustes(x, y, 'Scaling', false);

% plot points
subplot(2, 1, 1)
hold off
gplot(w, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

% plot solution
gplot(w, y, 'r')
plot(y(:, 1), y(:, 2), 'xr')

stopCondition

% plot stress evolution
subplot(2, 1, 2)
hold off
plot(t, sigma)
xlabel('Time (sec)')
ylabel('stress')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SPHERE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Toy example for the sphere. Sanity check: The initialization is the 
%% optimal solution

% uniform sampling of the sphere
[x, tri] = ParticleSampleSphere('N', 20);
N = size(x, 1);

% plot mesh
hold off
trisurf(tri, x(:, 1), x(:, 2), x(:, 3))
axis equal

% compute volume of each tetrahedron
vol = zeros(size(tri, 1), 1);
for I = 1:size(tri, 1)
    vol(I) = det(x(tri(I, :), :))/6;
end
if (all(vol >= 0))
    disp('Valid solution')
else
    error('Invalid solution')
end
vmin = .05;
vmax = .15;
R = 1;

% compute bounds and constraints for the spherical problem
[con, bnd] = tri_ccqp_smacof_nofold_sph_pip(tri, R, vmin, vmax);

% matrix of distances between the points on the sphere
dtot = dmatrix(x');

% SMACOF algorithm parameters
clear smacof_opts
smacof_opts.MaxIter = 1000;
smacof_opts.Epsilon = 1e-2;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
clear scip_opts
% scip_opts.limits_time = 100;
scip_opts.limits_solutions = 1;
scip_opts.display_verblevel = 4;

% solve MDS problem with constrained SMACOF
[y, stopCondition, sigma, t] ...
    = cons_smacof_pip(dtot, x, [], bnd, [], con, ...
    smacof_opts, scip_opts);

sum(y.^2, 2)'

% check that this is a valid solution
vol = zeros(size(tri, 1), 1);
for I = 1:size(tri, 1)
    vol(I) = det(y(tri(I, :), :))/6;
end
if (all(vol >= 0))
    disp('Valid solution')
else
    error('Invalid solution')
end

% plot points
subplot(2, 1, 1)
hold off
trisurf(tri, y(:, 1), y(:, 2), y(:, 3))
axis equal

stopCondition

% plot stress evolution
subplot(2, 1, 2)
plot(t, sigma)
xlabel('Time (sec)')
ylabel('stress')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same sphere toy example. Start with the optimal solution, except for
%% some vertices that have been perturbed to produce fold-overs. We don't
%% know which vertices are the ones with overlap, so we consider all of
%% them to be free

rng(0)

% create some overlaps
y0 = x;
y0(end-4:end, :) = y0(end-4:end, :) + 2*rand(5, 3);
y0(end-4:end, :) = y0(end-4:end, :) ./ repmat(sqrt(sum(y0(end-4:end, :).^2, 2)), 1, 3);

% compute volume of each tetrahedron
vol0 = zeros(size(tri, 1), 1);
for I = 1:size(tri, 1)
    vol0(I) = det(y0(tri(I, :), :))/6;
end
if (all(vol0 >= 0))
    disp('Valid solution')
else
    disp('Invalid solution')
end

% plot initialization
subplot(2, 1, 1)
trisurf(tri, y0(:, 1), y0(:, 2), y0(:, 3))
axis equal

% SMACOF algorithm parameters
clear smacof_opts
smacof_opts.MaxIter = 50;
smacof_opts.Epsilon = 1e-2;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
clear scip_opts
% scip_opts.limits_time = 100;
scip_opts.limits_solutions = 1;
scip_opts.display_verblevel = 0;

% solve MDS problem with constrained SMACOF
[y, stopCondition, sigma, t] ...
    = cons_smacof_pip(dtot, y0, [], bnd, [], con, ...
    smacof_opts, scip_opts);

sum(y.^2, 2)'

% check that this is a valid solution
vol = zeros(size(tri, 1), 1);
for I = 1:size(tri, 1)
    vol(I) = det(y(tri(I, :), :))/6;
end
if (all(vol >= 0))
    disp('Valid solution')
else
    error('Invalid solution')
end

% plot points
subplot(2, 1, 1)
hold off
trisurf(tri, y(:, 1), y(:, 2), y(:, 3))
axis equal

% plot points
subplot(2, 1, 2)
hold off
trisurf(tri, y0(:, 1), y0(:, 2), y0(:, 3))
axis equal

stopCondition

% plot stress evolution
subplot(2, 1, 2)
plot(t, sigma)
xlabel('Time (sec)')
ylabel('stress')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same as previous one, but now we assume that we know which vertices are
%% fixed

% free vertices
isFree = false(size(x, 1), 1);
isFree(end-4:end) = true;

% recompute bounds and constraints for the spherical problem
[con, bnd] = tri_ccqp_smacof_nofold_sph_pip(tri, R, vmin, vmax, isFree, x);

% SMACOF algorithm parameters
clear smacof_opts
smacof_opts.MaxIter = 100;
smacof_opts.Epsilon = 1e-2;
smacof_opts.Display = 'iter';
smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
clear scip_opts
% scip_opts.limits_time = 100;
scip_opts.limits_solutions = 1;
scip_opts.display_verblevel = 4;

% solve MDS problem with constrained SMACOF
[y, stopCondition, sigma, t] ...
    = cons_smacof_pip(dtot, y0, isFree, bnd, [], con, ...
    smacof_opts, scip_opts);

sum(y.^2, 2)'

% check that this is a valid solution
vol = zeros(size(tri, 1), 1);
for I = 1:size(tri, 1)
    vol(I) = det(y(tri(I, :), :))/6;
end
if (all(vol >= 0))
    disp('Valid solution')
else
    error('Invalid solution')
end

% plot points
subplot(2, 1, 1)
hold off
trisurf(tri, y(:, 1), y(:, 2), y(:, 3))
axis equal

% plot points
subplot(2, 1, 2)
hold off
trisurf(tri, y0(:, 1), y0(:, 2), y0(:, 3))
axis equal

stopCondition

% plot stress evolution
subplot(2, 1, 2)
plot(t, sigma)
xlabel('Time (sec)')
ylabel('stress')
