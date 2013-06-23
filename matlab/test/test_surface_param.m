% test_surface_param.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.2.2
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
d0 = dmatrix(uv0', uv0', 'euclidean');

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
%% Closed surface from scattered point set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% plot mesh
subplot(2, 2, 2)
hold off
trisurf(tri, x(:, 1), x(:, 2), x(:, 3));
axis equal

%% spherical Isomap

% compute spherical Isomap parametrization, no need the constrain the
% distance matrix more
param.d = dmatrix_mesh(x, tri);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Closed surface from segmentation mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load segmentation
scimat = scinrrd_load('data/008-lvhull-downsampled-4.mha');
scimat = scinrrd_load('data/008-rvhull-downsampled-4.mha');

% compute surface mesh from segmentation
opt = .004;
method = 'simplify';

tic
[x, tri] = v2s(single(scimat.data), 1, opt, method);
tri = tri(:, 1:3);
toc

% plot mesh
subplot(2, 1, 1)
hold off
plotmesh(x, tri)
axis equal

% normalize the mesh so that it's more spherical (CALD does not need this,
% but this way the comparison is fairer)
x = pca_normalize(x);


%% CALD (Control Area and Length Distortions)

% % first using the SPHARM-MAT toolbox GUI
% faces = tri;
% vertices = x;
% save('/tmp/bar_obj.mat', 'faces', 'vertices')
% SPHARM_MAT
% load('/tmp/bar_CALD_smo.mat')
% 
% % plot parametrization
% subplot(2, 2, 3)
% hold off
% [~, xyzsph] = procrustes(x, sph_verts, 'Scaling', false);
% trisurf(tri, xyzsph(:, 1), xyzsph(:, 2), xyzsph(:, 3));
% title('CALD')
% axis equal

% compute CALD parametrization using our surface_param() function
param.type = 'cald';
param.tri = tri;
[uv, out] = surface_param(x, param);
latCALD = uv(:, 1);
lonCALD = uv(:, 2);

% plot parametrization
subplot(2, 2, 3)
hold off
lat = uv(:, 1);
lon = uv(:, 2);
[xsph, ysph, zsph] = sph2cart(lon, lat, 1);
[~, xyzsph] = procrustes(x, [xsph, ysph, zsph], 'Scaling', false);
trisurf(tri, xyzsph(:, 1), xyzsph(:, 2), xyzsph(:, 3));
title('CALD')
axis equal

%% spherical Isomap

% compute spherical Isomap parametrization
param.d = dmatrix_mesh(x, tri);
param.type = 'sphisomap';
param.neigh = 'epsilon';
param.size = Inf;
param.init = 'random';
param.maxiter = 50;
[uv, out] = surface_param(x, param);
latIso = uv(:, 1);
lonIso = uv(:, 2);

% plot parametrization
subplot(2, 2, 4)
hold off
[xsph, ysph, zsph] = sph2cart(lonIso, latIso, 1);
[~, xyzsph] = procrustes(x, [xsph ysph zsph], 'Scaling', false);
trisurf(tri, xsph, ysph, zsph);
title('Spherical Isomap')
axis equal

% plot error
subplot(1, 1, 1)
hold off
plot(out.err)
ylabel('||D-D_{param}||_{Frob}')
xlabel('Iteration (each point movement)')
title('Isometry error')

% plot distances to optimize vs spherical distances
d = dijkstra(sparse(param.d), 1:length(param.d));
subplot(2, 1, 2)
hold off
plot(d(:), out.dsph(:), '.')
hold on
plot([0 max(d(:))], [0 max(d(:))], 'r', 'LineWidth', 2)
title('Spherical Isomap')
xlabel('Distance on the mesh')
ylabel('Distance on the spherical parametrization')
axis equal

%% compare spherical Isomap and CALD

% compute great circle distances for CALD
dsphCALD = zeros(size(out.dsph));
for I = 1:length(dsphCALD)
    dsphCALD(:, I) = distance(latCALD(I), lonCALD(I), latCALD, lonCALD, ...
        [out.sphrad 0], 'radians');
end

% plot distances to optimize vs spherical distances for CALD
subplot(2, 1, 1)
hold off
plot(d(:), dsphCALD(:), '.')
hold on
plot([0 max(d(:))], [0 max(d(:))], 'r', 'LineWidth', 2)
title('CALD')
xlabel('Distance on the mesh')
ylabel('Distance on the spherical parametrization')
axis equal


% plot boxplots with the distance ratios
N = numel(d);
subplot(1, 1, 1)
hold off
boxplot([d(:)-dsphCALD(:) ; d(:)-out.dsph(:)], ...
    [ones(N, 1); 2*ones(N, 1)], 'labels', {'CALD', 'Sph Isomap'}, ...
    'notch', 'on')
title('Target distance - distance on spherical parametrization')
