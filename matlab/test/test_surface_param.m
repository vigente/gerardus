% test_surface_param.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Open surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('data/valve-annula-points-009.mat');

% plot points
hold off
subplot(2, 2, 1)
plot3(x(:, 1), x(:, 2), x(:, 3), '*')
axis equal

% compute simple XY parametrization
uv0 = surface_param(x, 'xy');
d0 = dmatrix(uv', uv', 'euclidean');

% plot parametrization
subplot(2, 2, 2)
hold off
plot(uv0(:, 1), uv0(:, 2), '.')
axis equal
title('XY')


% compute PCA parametrization
uv = surface_param(x, 'pca');

% plot parametrization, aligned with the XY case
[~, uv] = procrustes(uv0, uv, 'Scaling', false);
subplot(2, 2, 3)
hold off
plot(uv(:, 1), uv(:, 2), '.')
axis equal
title('PCA')


% compute Isomap parametrization
param.type = 'isomap';
param.neigh = 'epsilon';
param.size = .0020;
param.d = dmatrix(x', x', 'euclidean');
uv = surface_param(x, param);

% plot parametrization
[~, uv] = procrustes(uv0, uv, 'Scaling', false);
subplot(2, 2, 4)
hold off
plot(uv(:, 1), uv(:, 2), '.')
axis equal
title('Isomap')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Closed surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% spherical Isomap with an alpha-shape distance matrix

% load discrete point set
scimat = scinrrd_load('data/thick-slice-points-seg.mat');
idx = find(scimat.data);
[r, c, s] = ind2sub(size(scimat.data), idx);
x = scinrrd_index2world([r c s], scimat.axis);

% plot point set
subplot(2, 2, 1)
hold off
plot3(x(:, 1), x(:, 2), x(:, 3), '.')
axis equal

% mesh the points using an alpha-shape
[~, s] = alphavol(x, 0.002);
tri = s.bnd;
param.d = dmatrix_mesh(x, tri);

% plot mesh
subplot(2, 2, 2)
hold off
trisurf(tri, x(:, 1), x(:, 2), x(:, 3));
axis equal

% compute spherical Isomap parametrization, no need the constrain the
% distance matrix more
param.type = 'sphisomap';
param.neigh = 'epsilon';
param.size = Inf;
param.init = 'random';
[uv, out] = surface_param(x, param);
lat = uv(:, 1);
lon = uv(:, 2);

% plot parametrization
subplot(2, 2, 3)
hold off
[xsph, ysph, zsph] = sph2cart(lon, lat, 1);
[~, xyzsph] = procrustes(x, [xsph ysph zsph], 'Scaling', false);
trisurf(tri, xsph, ysph, zsph);
axis equal

% plot error
subplot(2, 2, 4)
hold off
plot(out.err)
ylabel('||D-D_{param}||_{Frob}')
xlabel('Iteration (each point movement)')
title('Isometry error')
