% test_tri_sphparam.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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

%% create toy example

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

%% Classical MDS parametrization

% algorithm parameters
clear sphparam_opts
sphparam_opts.sphrad = R;
sphparam_opts.TopologyCheck = true;
sphparam_opts.Display = 'iter';

% matrix of true chord-length distances between all pairs of points on the
% sphere
dtot = dmatrix(x');

% we need arc length geodesic distances for the input of tri_sphparam()
darclen = chord2arclen(dtot, sphparam_opts.sphrad);

% compute MDS parametrization
[y, t] = tri_sphparam(tri, x, 'cmdscale', darclen, [], sphparam_opts);

% register solution to ground truch
err = procrustes(x, y, 'Scaling', false);

%% Unconstrained SMACOF parametrization

% algorithm parameters
clear sphparam_opts
sphparam_opts.sphrad = R;
sphparam_opts.TopologyCheck = true;
sphparam_opts.Display = 'iter';

% SMACOF parameters
clear smacof_opts
smacof_opts.MaxIter = 10;

% with SMACOF, we can use a sparse distance matrix
d = dmatrix_mesh(tri, x);

% we need arc length geodesic distances for the input of tri_sphparam()
d = chord2arclen(d, sphparam_opts.sphrad);

% random perturbation of the optimal solution as initialization
y0 = x + 0.5 * rand(size(y));
y0 = y0 ./ repmat(sum(y0.^2, 2), 1, 3);

% compute SMACOF parametrization
[y, yIsValid, stopCondition, sigma, t] ...
    = tri_sphparam(tri, x, 'smacof', d, y0, sphparam_opts, smacof_opts);

% plot parametrization result
subplot(2, 1, 1)
hold off
trisurf(tri, y(:, 1), y(:, 2), y(:, 3))
axis equal
subplot(2, 1, 2)
hold off
plot(t, sigma)

%% Constrained SMACOF parametrization, untangling of individual components

% algorithm parameters
clear sphparam_opts
sphparam_opts.sphrad = R;
sphparam_opts.TopologyCheck = true;
sphparam_opts.Display = 'iter';

% SMACOF algorithm parameters
clear smacof_opts
smacof_opts.MaxIter = 15;
% smacof_opts.Epsilon = 1e-2;
smacof_opts.Display = 'iter';
% smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
clear scip_opts
% scip_opts.limits_time = 100;
scip_opts.limits_solutions = 1;
scip_opts.display_verblevel = 0;

% with SMACOF, we can use a sparse distance matrix
d = dmatrix_mesh(tri, x);

% we need arc length geodesic distances for the input of tri_sphparam()
d = chord2arclen(d, sphparam_opts.sphrad);

% random perturbation of the optimal solution as initialization
rng(0)
y0 = x;
y0(end-4:end, :) = y0(end-4:end, :) + 2 * rand(5, 3);
y0 = y0 ./ repmat(sqrt(sum(y0.^2, 2)), 1, 3);

% solve MDS problem with constrained SMACOF
[y, yIsValid, stopCondition, sigma, t] ...
    = tri_sphparam(tri, x, 'consmacof-local', d, y0, ...
    sphparam_opts, smacof_opts, scip_opts);

% plot parametrization result
subplot(2, 1, 1)
hold off
trisurf(tri, y(:, 1), y(:, 2), y(:, 3))
axis equal
subplot(2, 1, 2)
hold off
plot(t{1}, sigma{1})

%% Constrained SMACOF parametrization, untangling whole mesh in one go

% algorithm parameters
clear sphparam_opts
sphparam_opts.sphrad = R;
sphparam_opts.TopologyCheck = true;
sphparam_opts.Display = 'iter';

% SMACOF algorithm parameters
clear smacof_opts
smacof_opts.MaxIter = 15;
% smacof_opts.Epsilon = 1e-2;
smacof_opts.Display = 'iter';
% smacof_opts.TolFun = 1e-6;

% SCIP algorithm parameters
clear scip_opts
% scip_opts.limits_time = 100;
scip_opts.limits_solutions = 1;
scip_opts.display_verblevel = 0;

% with SMACOF, we can use a sparse distance matrix
d = dmatrix_mesh(tri, x);

% we need arc length geodesic distances for the input of tri_sphparam()
d = chord2arclen(d, sphparam_opts.sphrad);

% random perturbation of the optimal solution as initialization
rng(0)
y0 = x;
y0(end-4:end, :) = y0(end-4:end, :) + 2 * rand(5, 3);
y0 = y0 ./ repmat(sqrt(sum(y0.^2, 2)), 1, 3);

% solve MDS problem with constrained SMACOF
[y, yIsValid, stopCondition, sigma, t] ...
    = tri_sphparam(tri, x, 'consmacof-global', d, y0, ...
    sphparam_opts, smacof_opts, scip_opts);

% plot parametrization result
subplot(2, 1, 1)
hold off
trisurf(tri, y(:, 1), y(:, 2), y(:, 3))
axis equal
subplot(2, 1, 2)
hold off
plot(t, sigma)
