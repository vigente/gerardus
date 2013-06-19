% test_smdscale.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012-2013 University of Oxford
% Version: 0.2.3
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
%% We have a sparse distance matrix that defines a local neighbourhood on a 
%% scattered set of points
%%
%% We want to project it onto the sphere

% load discrete point set with an associated local neighbourhood
load('data/thick-slice-points-xyz-d.mat')

% plot neighbourhood
hold off
subplot(1, 1, 1)
gplot3d(d, xyz)
axis equal

% initial guess for the sphere embedding
[lat, lon] = proj_on_sphere(xyz);

% embbed the point set on the sphere using the initial guess
tic
opt.maxiter = 100;
[lat, lon, err, stopCondition, dsph, sphrad] = ...
    smdscale(sparse(d), [], lat, lon, opt);
toc

% plot error
hold off
plot(err)
ylabel('||D-D_{param}||_{Frob}')
xlabel('Iteration (each point movement)')
title('Isometry error')

% compute the Euclidean coordinates of the projected points
[xsph, ysph, zsph] = sph2cart(lon, lat, sphrad);

% use a rigid Procrustes to find a rotation that aligns the sphere with the
% LV points
[~, xyzsph] = procrustes(xyz, [xsph ysph zsph], 'Scaling', false);

% plot the aligned points
hold off
subplot(1, 1, 1)
gplot3d(d, xyz)
axis equal
hold on
gplot3d(d, xyzsph, 'r')
axis equal

% plot the normalised distance matrix error
hold off
idx = d>0;
boxplot( abs(d(idx)-dsph(idx))./d(idx) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In this case we have no good initial guess, so we use a random one

% random distribution of points on the sphere
N = length(d);
lat = (rand(1, N)-.5) * pi;
lon = rand(1, N) * pi * 2;
[x, y, z] = sph2cart(lon, lat, sphrad);
xyz = cat(1, x, y, z)';
clear x y z

% plot initial guess
hold off
subplot(1, 1, 1)
gplot3d(d, xyz)
axis equal

% embbed the point set on the sphere using the random points. We use the
% Dijkstra'd full distance matrix, because local optimisation reduces the
% error but doesn't produce a visually good result
tic
opt.maxiter = 100;
[lat, lon, err0, stopCondition, dsph] = ...
    smdscale(dijkstra(sparse(d), 1:N), sphrad, lat, lon, opt);
toc

% convert sphere coordinates to Euclidean coordinates
[x, y, z] = sph2cart(lon, lat, sphrad);
xyz = [x y z];
clear x y z

% plot result
hold on
gplot3d(d, xyz, 'r')
axis equal

% plot the error functions
hold off
plot(err0, '--')
hold on
plot(err)
title('Distance matrix approximation error')
xlabel('No. of optimisation steps (moving one point counts as 1 step)')
ylabel('Frobenius norm of distance matrix error')
legend('Random initialisation', 'Sphere projection intialisation')

% now we run MDS again, but this time starting from the previous result,
% and using the sparse distance matrix, to fine tune the result
tic
opt.maxiter = 100;
[lat, lon, err1, stopCondition, dsph] = ...
    smdscale(d, sphrad, lat, lon, opt);
toc

% convert sphere coordinates to Euclidean coordinates
[x, y, z] = sph2cart(lon, lat, sphrad);
xyz = cat(1, x, y, z)';
clear x y z

% plot result
hold off
gplot3d(d, xyz, 'r')
axis equal

% plot the error functions
hold off
plot(err0, '--')
hold on
plot(err)
plot(err1, ':')
title('Distance matrix approximation error')
xlabel('No. of optimisation steps (moving one point counts as 1 step)')
ylabel('Frobenius norm of distance matrix error')
legend('Random initialisation', 'Sphere projection intialisation', ...
    'Random init + Local neighbourhood')
