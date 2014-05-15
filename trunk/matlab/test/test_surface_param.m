% test_surface_param.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.6.0
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
%% Open surface for 4 cardiac valves
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
%% Open planar triangular mesh to compare isomap and MDS map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% quadrangular mesh split into asymmetric triangles

% create rectangle with triangular mesh
[tri, x] = surface_tridomain('rect', 'step', 1/10, [0 0], [1 1/5*3]);
x(:, end+1) = 0;

% plot original mesh and parametrization
con = dmatrix_mesh(tri);
hold off
gplot(con, x(:, 1:2))

%% Isomap on a triangular mesh
% note that we need to compute the local distances on the mesh, and then
% apply Dijkstra. If we don't provide the distance matrix, by default
% surface_param() will compute distances between all pairs of vertices
clear param
param.type = 'isomap';
param.d = dmatrix_mesh(tri, x);
param.d = dijkstra(param.d, 1:length(param.d));
uv = surface_param(x, param);

% rigid registration to overlap the solution with the original data
[~, uv] = procrustes(x, uv, 'Scaling',false);

% add parametrization to the plot
hold on
gplot(con, uv, 'r')


%% Open classic MDSmap on the same triangular mesh
clear param
param.type = 'cmdsmap';
param.options.constraint_map = 0.4 * ones(size(x, 1), 1); % trigger warning
param.options.end_points = [1 2 3]; % trigger warning
param.options.nb_iter_max = 4; % trigger warning
param.tri = tri;
uv = surface_param(x, param);

% rigid registration to overlap the solution with the original data
[~, uv] = procrustes(x, uv, 'Scaling',false);

% add parametrization to the plot
hold on
gplot(con, uv, 'g')

%% check that MDSmap with Dijkstra is the same as Isomap
clear param
param.type = 'cmdsmap';
param.options.constraint_map = 0.4 * ones(1, size(x, 1)); % trigger warning
param.options.end_points = [1 2 3]; % trigger warning
param.options.nb_iter_max = 4; % trigger warning
[~, param.d] = dmatrix_mesh(tri, x, 'dijkstra');
uv = surface_param(x, param);

% rigid registration to overlap the solution with the original data
[~, uv] = procrustes(x, uv, 'Scaling',false);

% add parametrization to the plot
% it is expected that this perfectly overlaps the Isomap result computed
% above
hold on
gplot(con, uv, 'k')

%% lmdscale (MDSmap with local neighbourhood)
clear param
param.type = 'lmdscale';
param.tri = tri;
param.options.constraint_map = 0.4 * ones(size(x, 1), 1);
param.options2.MaxIter = 100;
param.options2.MaxInc = .01;
[uv, out] = surface_param(x, param);
out.stopCondition

% rigid registration to overlap the solution with the original data
[~, uv] = procrustes(x, uv, 'Scaling',false);

% add parametrization to the plot
% it is expected that this perfectly overlaps the Isomap result computed
% above
hold on
gplot(con, uv, 'c')

% plot error
figure
plot(out.err.maxinc)
hold on
plot(out.err.medinc)

% plot stress
hold off
plot(out.err.stress1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Closed surface from scattered point set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create a closed surface triangular mesh of a left ventricle

% seed the random generator so that we always get the same points
rng(0)

% load low resolution segmentation of LV, with detail removed
scimat = scinrrd_load('data/008-lv-resampled-0_31.mha');

% keep perimeter only
scimat.data = bwperim(scimat.data);

% coordinates of perimeter points
idx = find(scimat.data);
[r, c, s] = ind2sub(size(scimat.data), idx);
x = scinrrd_index2world([r c s], scimat.axis);

% subsample boundary points
x = x(randi(size(x, 1), round(size(x, 1)/5), 1), :);

% order a bit the vertices, so that nearby vertices are close
x = sortrows(x);

% plot point set
subplot(2, 2, 1)
hold off
plot3(x(:, 1), x(:, 2), x(:, 3), '.')
axis equal

% mesh the points using an alpha-shape
[~, s] = alphavol(x, scimat.axis(1).spacing * 40);
tri = s.bnd;

% check that we don't have non-manifold vertices
idx = tri_find_nonmanifold_vertex(tri, x, scimat.axis);
if (nnz(idx))
    error('Assertion fail: mesh has non-manifold vertices')
end

% remove points not connected to any triangle
[tri, x] = tri_squeeze(tri, x);

% plot mesh
subplot(2, 2, 2)
hold off
trisurf(tri, x(:, 1), x(:, 2), x(:, 3));
axis equal

% PCA-normalization of mesh
x = pca_normalize(x);

% plot mesh
subplot(2, 2, 3)
hold off
trisurf(tri, x(:, 1), x(:, 2), x(:, 3));
axis equal

% median radius of the points
sphrad = median(sqrt(sum(x.^2, 2)))

%% sphproj: simple projection on a sphere around the centroid
clear param
param.type = 'sphproj';
[latlon, out] = surface_param(x, param);
out.medrad

% rigid registration to overlap the solution with the original data
[xx, yy, zz] = sph2cart(latlon(:, 2), latlon(:, 1), sphrad);
[~, x1] = procrustes(x, [xx, yy, zz], 'Scaling',false);

% plot result
subplot(1, 2, 1)
hold off
trisurf(tri, x(:, 1), x(:, 2), x(:, 3));
axis equal
subplot(1, 2, 2)
hold off
trisurf(tri, x1(:, 1), x1(:, 2), x1(:, 3))
title('sphproj')
axis equal

%% CALD

% % using the SPHARM-MAT toolbox GUI
% faces = tri;
% vertices = x;
% save('/tmp/bar_obj.mat', 'faces', 'vertices')
% SPHARM_MAT
% load('/tmp/bar_CALD_smo.mat')

clear param
param.type = 'cald';
param.tri = tri;
param.options.MeshGridSize = 50;
param.options.MaxSPHARMDegree = 6;
param.options.Tolerance = 2;
param.options.Smoothing = 2;
param.options.Iteration = 100;
param.options.LocalIteration = 10;
latlon = surface_param(x, param);

% rigid registration to overlap the solution with the original data
[xx, yy, zz] = sph2cart(latlon(:, 2), latlon(:, 1), sphrad);
[~, x1] = procrustes(x, [xx, yy, zz], 'Scaling',false);

% plot parametrization
subplot(1, 2, 2)
hold off
trisurf(tri, x1(:, 1), x1(:, 2), x1(:, 3))
title('CALD')
axis equal

%% smdscale

% compute spherical Isomap parametrization, no need the constrain the
% distance matrix more
clear param
param.type = 'smdscale';
param.tri = tri;
param.dmethod = 'fastmarching';
% param.options.constraint_map = .2 * sphrad * ones(size(x, 1), 1);
param.sphrad = sqrt(sum(cgal_trifacet_area(tri, x))/4/pi);
param.init = 'random';
param.options2.MaxIter = 15;
param.options2.MaxAlpha = .01 / 180 * pi;
[uv, out] = surface_param(x, param);
out.stopCondition
lat = uv(:, 1);
lon = uv(:, 2);

% plot parametrization
subplot(1, 2, 2)
hold off
[xx, yy, zz] = sph2cart(lon, lat, param.sphrad);
[~, x1] = procrustes(x, [xx, yy, zz], 'Scaling', false);
trisurf(param.tri, x1(:, 1), x1(:, 2), x1(:, 3));
title('smdscale')
axis equal

% check for self-intersections
bad = cgal_check_self_intersect(param.tri, x1);
nnz(bad>0)

% plot error
figure
subplot(2, 2, 1)
hold off
plot(out.err.rawstress)
title('Raw stress')
subplot(2, 2, 2)
hold off
plot(out.err.stress1)
title('Stress1')
subplot(2, 1, 2)
hold off
plot(out.err.maxalpha)
hold on
plot(out.err.medalpha)
title('Alpha')

% TODO:
% %% compare spherical Isomap and CALD
% 
% % compute great circle distances for CALD
% dsphCALD = zeros(size(out.dsph));
% for I = 1:length(dsphCALD)
%     dsphCALD(:, I) = distance(latCALD(I), lonCALD(I), latCALD, lonCALD, ...
%         [out.sphrad 0], 'radians');
% end
% 
% % plot distances to optimize vs spherical distances for CALD
% subplot(2, 1, 1)
% hold off
% plot(d(:), dsphCALD(:), '.')
% hold on
% plot([0 max(d(:))], [0 max(d(:))], 'r', 'LineWidth', 2)
% title('CALD')
% xlabel('Distance on the mesh')
% ylabel('Distance on the spherical parametrization')
% axis equal
% 
% 
% % plot boxplots with the distance ratios
% N = numel(d);
% subplot(1, 1, 1)
% hold off
% boxplot([d(:)-dsphCALD(:) ; d(:)-out.dsph(:)], ...
%     [ones(N, 1); 2*ones(N, 1)], 'labels', {'CALD', 'Sph Isomap'}, ...
%     'notch', 'on')
% title('Target distance - distance on spherical parametrization')
