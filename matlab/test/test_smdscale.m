% test_smdscale.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012-2013 University of Oxford
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

% embbed the point set on the sphere
tic
opt.MaxAlpha = 1/180*pi; % stop if points don't move more than 1 degree
[lat, lon, stopCondition, err, dout, sphrad] = ...
    smdscale(sparse(d), [], lat, lon, opt);
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

% plot the normalised distance matrix error as a boxplot
hold off
idx = d>0;
boxplot(abs(d(idx)-dout(idx))./d(idx))

% plot the normalised distance matrix error as distance/distance scattered
% plot
hold off
plot(d(idx), (d(idx)-dout(idx)./d(idx)), '.')
hold on
plot([0 3.5e-3], [0 0])
axis([0 3.5e-3 -2.5 .5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% In this case we have no good initial guess, so we use a random one

% random distribution of points on the sphere
N = length(d);
lat = (rand(N, 1)-.5) * pi;
lon = rand(N, 1) * pi * 2;
[x, y, z] = sph2cart(lon, lat, sphrad);
xyz = [x(:), y(:), z(:)];
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
opt.MaxAlpha = 1/180*pi; % stop if points don't move more than 1 degree
[lat, lon, stopCondition, err0, dout] = ...
    smdscale(dijkstra(sparse(d), 1:N), sphrad, lat, lon, opt);
toc

% convert sphere coordinates to Euclidean coordinates
[x, y, z] = sph2cart(lon, lat, sphrad);
xyz = [x(:), y(:), z(:)];
clear x y z

% plot result
hold on
gplot3d(d, xyz, 'r')
axis equal

% plot the normalised distance matrix error as a boxplot
hold off
idx = d>0;
boxplot(abs(d(idx)-dout(idx))./d(idx))

% plot the normalised distance matrix error as distance/distance scattered
% plot
hold off
plot(d(idx), (d(idx)-dout(idx)./d(idx)), '.')
hold on
plot([0 3.5e-3], [0 0])
axis([0 3.5e-3 -2.5 .5])

% now we run MDS again, but this time starting from the previous result,
% and using the sparse distance matrix, to fine tune the result
tic
[lat, lon, stopCondition, err1, dout] = ...
    smdscale(d, sphrad, lat, lon, opt);
toc

% convert sphere coordinates to Euclidean coordinates
[x, y, z] = sph2cart(lon, lat, sphrad);
xyz = [x(:), y(:), z(:)];
clear x y z

% plot result
hold off
gplot3d(d, xyz, 'r')
axis equal

% plot the normalised distance matrix error as a boxplot
hold off
idx = d>0;
boxplot(abs(d(idx)-dout(idx))./d(idx))

% plot the normalised distance matrix error as distance/distance scattered
% plot
hold off
plot(d(idx), (d(idx)-dout(idx)./d(idx)), '.')
hold on
plot([0 3.5e-3], [0 0])
axis([0 3.5e-3 -2.5 .5])
