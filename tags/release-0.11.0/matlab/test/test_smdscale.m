% test_smdscale.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012-2013 University of Oxford
% Version: 0.3.0
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
%% Uniform sampling of points on a sphere, with a triangular mesh, using
%% true geodesic distances, random initialization

% uniformly distribute particles across the surface of the unit sphere
% using Anton Semechko's implementation of Reisz s-energy minimisation
[xyz, tri] = ParticleSampleSphere('N', 50);

% plot mesh
hold off
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
axis equal

% compute true lat, lon coordinates of each point
[lon, lat] = cart2sph(xyz(:, 1), xyz(:, 2), xyz(:, 3));

% compute full matrix of distances between all pairs of points
lat2 = repmat(lat, 1, size(lat, 1));
lon2 = repmat(lon, 1, size(lon, 1));
d = distance(lat2, lon2, lat2', lon2', 'radians');

% compute spherical MDS using the full matrix on a random initialization
tic
opt.MaxIter = 50;
opt.MaxAlpha = 1/180*pi; % stop if points don't move more than 1 degree
sphrad = 1;
[lat2, lon2, stopCondition, err, dout] = ...
    smdscale(d, sphrad, [], [], opt);
toc

% plot errors
cla
subplot(2, 2, 1)
plot(err.rawstress)
title('Raw stress')
subplot(2, 2, 2)
plot(err.stress1)
title('Stress-1')
subplot(2, 1, 2)
hold off
plot(err.maxalpha/pi*180)
hold on
plot(err.medalpha/pi*180, '--')
ylabel('degrees')
title('Alpha')

% compute the Euclidean coordinates of the projected points
[xsph, ysph, zsph] = sph2cart(lon2, lat2, sphrad);

% use a rigid Procrustes to find a rotation that aligns the sphere with the
% LV points
[~, xyzsph] = procrustes(xyz, [xsph ysph zsph], 'Scaling', false);

% plot the aligned points
hold off
subplot(1, 2, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
axis equal
subplot(1, 2, 2)
trisurf(tri, xyzsph(:, 1), xyzsph(:, 2), xyzsph(:, 3))
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uniform sampling of points on a sphere, with a triangular mesh, using
%% Fast Marching distances, random initialization

% compute Fast Marching full matrix of distances between all pairs of
% points
for I = 1:size(lat, 1)
    d(:, I) = perform_fast_marching_mesh(xyz, tri, I);
end

% compute spherical MDS using the full matrix on a random initialization
tic
opt.MaxIter = 50;
opt.MaxAlpha = 1/180*pi; % stop if points don't move more than 1 degree
sphrad = 1;
[lat2, lon2, stopCondition, err, dout] = ...
    smdscale(d, sphrad, [], [], opt);
toc

% plot errors
cla
subplot(2, 2, 1)
plot(err.rawstress)
title('Raw stress')
subplot(2, 2, 2)
plot(err.stress1)
title('Stress-1')
subplot(2, 1, 2)
hold off
plot(err.maxalpha/pi*180)
hold on
plot(err.medalpha/pi*180, '--')
ylabel('degrees')
title('Alpha')

% compute the Euclidean coordinates of the projected points
[xsph, ysph, zsph] = sph2cart(lon2, lat2, sphrad);

% use a rigid Procrustes to find a rotation that aligns the sphere with the
% LV points
[~, xyzsph] = procrustes(xyz, [xsph ysph zsph], 'Scaling', false);

% plot the aligned points
hold off
subplot(1, 2, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
axis equal
subplot(1, 2, 2)
trisurf(tri, xyzsph(:, 1), xyzsph(:, 2), xyzsph(:, 3))
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uniform sampling of points on a sphere, with a triangular mesh, using
%% Fast Marching distances, local neighbourhood, random initialization

% compute Fast Marching full matrix of distances between all pairs of
% points
for I = 1:size(lat, 1)
    d(:, I) = perform_fast_marching_mesh(xyz, tri, I);
end

% create local neighbourhoods by removing connections in distance matrix
% larger than 45 degrees
d(d > 135/180*pi) = 0;
d = sparse(d);

% compute spherical MDS using the full matrix on a random initialization
tic
opt.MaxIter = 100;
opt.MaxAlpha = 1/180*pi; % stop if points don't move more than 1 degree
sphrad = 1;
[lat2, lon2, stopCondition, err, dout] = ...
    smdscale(d, sphrad, [], [], opt);
toc

% check whether the result is a valid parametrization
ok = check_sphtrisurf_vertices(tri, [lat2 lon2]);
all(ok)

% plot errors
cla
subplot(2, 2, 1)
plot(err.rawstress)
title('Raw stress')
subplot(2, 2, 2)
plot(err.stress1)
title('Stress-1')
subplot(2, 1, 2)
hold off
plot(err.maxalpha/pi*180)
hold on
plot(err.medalpha/pi*180, '--')
ylabel('degrees')
title('Alpha')

% compute the Euclidean coordinates of the projected points
[xsph, ysph, zsph] = sph2cart(lon2, lat2, sphrad);

% use a rigid Procrustes to find a rotation that aligns the sphere with the
% LV points
[~, xyzsph] = procrustes(xyz, [xsph ysph zsph], 'Scaling', false);

% plot the aligned points
hold off
subplot(1, 2, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
axis equal
subplot(1, 2, 2)
trisurf(tri, xyzsph(:, 1), xyzsph(:, 2), xyzsph(:, 3))
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Left Ventricle mesh without PCA normalization

clear

% load segmentation mask of LV hull
im = scinrrd_load('data/008-lvhull-downsampled-4.mha');

% compute mesh from segmentation
opt = .004;
method = 'simplify';

tic
[xyz, tri] = v2s(single(im.data), 1, opt, method);
toc

% extract boundary mesh
tri = TriRep(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3));
tri = freeBoundary(tri);

% remove unused nodes
[tri, xyz] = tri_squeeze(tri, xyz);

% keep a note of the direct directions for later
con = dmatrix_mesh(tri)~=0;

% plot mesh
hold off
subplot(1, 1, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
axis equal

% initial guess for the sphere embedding
[lat, lon] = proj_on_sphere(xyz);

% compute Fast Marching full matrix of distances between all pairs of
% points
d = zeros(size(lat, 1));
for I = 1:size(lat, 1)
    d(:, I) = perform_fast_marching_mesh(xyz, tri, I);
end

% embbed the point set on the sphere
tic
opt.MaxIter = 50;
opt.MaxAlpha = 1/180*pi; % stop if points don't move more than 1 degree
[lat2, lon2, stopCondition, err, dout, sphrad] = ...
    smdscale(sparse(d), [], lat, lon, opt);
toc

% check whether the result is a valid parametrization
ok = check_sphtrisurf_vertices(tri, [lat2 lon2]);
all(ok)

% plot errors
cla
subplot(2, 2, 1)
plot(err.rawstress)
title('Raw stress')
subplot(2, 2, 2)
plot(err.stress1)
title('Stress-1')
subplot(2, 1, 2)
hold off
plot(err.maxalpha/pi*180)
hold on
plot(err.medalpha/pi*180, '--')
ylabel('degrees')
title('Alpha')

% compute the Euclidean coordinates of the projected points
[xsph, ysph, zsph] = sph2cart(lon2, lat2, sphrad);

% use a rigid Procrustes to find a rotation that aligns the sphere with the
% LV points
[~, xyzsph] = procrustes(xyz, [xsph ysph zsph], 'Scaling', false);

% plot the aligned points
hold off
subplot(1, 2, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
axis equal
subplot(1, 2, 2)
trisurf(tri, xyzsph(:, 1), xyzsph(:, 2), xyzsph(:, 3))
axis equal

% plot the normalised distance matrix error as a boxplot
subplot(1, 2, 1)
hold off
boxplot(abs(d(con)-dout(con))./d(con))

% plot the normalised distance matrix error as distance/distance scattered
% plot
subplot(1, 2, 2)
hold off
plot(d(con), (d(con)-dout(con)./d(con)), '.')
hold on
plot([0 30], [0 30], 'r', 'LineWidth', 2)
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Left Ventricle mesh with PCA normalization

clear

% load segmentation mask of LV hull
im = scinrrd_load('data/008-lvhull-downsampled-4.mha');

% compute mesh from segmentation
opt = .004;
method = 'simplify';

tic
[xyz, tri] = v2s(single(im.data), 1, opt, method);
toc

% extract boundary mesh
tri = TriRep(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3));
tri = freeBoundary(tri);

% remove unused nodes
[tri, xyz] = tri_squeeze(tri, xyz);

% keep a note of the direct directions for later
con = dmatrix_mesh(tri)~=0;

% plot mesh
hold off
subplot(1, 1, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
axis equal

% PCA normalization of the mesh
xyz = pca_normalize(xyz);

% initial guess for the sphere embedding
[lat, lon, sphrad] = proj_on_sphere(xyz);

% compute Fast Marching full matrix of distances between all pairs of
% points
d = zeros(size(lat, 1));
for I = 1:size(lat, 1)
    d(:, I) = perform_fast_marching_mesh(xyz, tri, I);
end

% embbed the point set on the sphere
tic
opt.MaxIter = 50;
opt.MaxAlpha = 1/180*pi; % stop if points don't move more than 1 degree
[lat2, lon2, stopCondition, err, dout] = ...
    smdscale(sparse(d), sphrad, lat, lon, opt);
toc

% check whether the result is a valid parametrization
ok = check_sphtrisurf_vertices(tri, [lat2 lon2]);
all(ok)

% plot errors
cla
subplot(2, 2, 1)
plot(err.rawstress)
title('Raw stress')
subplot(2, 2, 2)
plot(err.stress1)
title('Stress-1')
subplot(2, 1, 2)
hold off
plot(err.maxalpha/pi*180)
hold on
plot(err.medalpha/pi*180, '--')
ylabel('degrees')
title('Alpha')

% compute the Euclidean coordinates of the projected points
[xsph, ysph, zsph] = sph2cart(lon2, lat2, sphrad);

% use a rigid Procrustes to find a rotation that aligns the sphere with the
% LV points
[~, xyzsph] = procrustes(xyz, [xsph ysph zsph], 'Scaling', false);

% plot the aligned points
hold off
subplot(1, 2, 1)
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
axis equal
subplot(1, 2, 2)
trisurf(tri, xyzsph(:, 1), xyzsph(:, 2), xyzsph(:, 3))
axis equal

% plot the normalised distance matrix error as a boxplot
subplot(1, 2, 1)
hold off
boxplot(abs(d(con)-dout(con))./d(con))

% plot the normalised distance matrix error as distance/distance scattered
% plot
subplot(1, 2, 2)
hold off
plot(d(con), (d(con)-dout(con)./d(con)), '.')
hold on
plot([-1 2], [-1 2], 'r', 'LineWidth', 2)
axis equal
