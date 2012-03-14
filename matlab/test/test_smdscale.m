% test_smdscale.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
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

% number of points
N = 100;

% elevation or latitude
lat = (rand(1, N)-.5) * pi;

% longitude
lon = rand(1, N) * pi * 2;

% constrain all points to be on the surface of a sphere of radius 1
r = ones(1, N);

% convert angles from radians to degrees
lat = lat / pi * 180;
lon = lon / pi * 180;

% compute true geodesic distance between each pair of points, in degrees of
% arc
d = gdmatrix(lat, lon);

% convert latitude/longitude to Cartesian coordinates
[x, y, z] = sph2cart(lon/180*pi, lat/180*pi, r);
x = cat(1, x, y, z);
clear y z

% plot points
close all
plot3(x(1, :), x(2, :), x(3, :), 'x')
axis equal

% save ground truth coordinates for validation
lat0 = lat;
lon0 = lon;
x0 = x;

% % add noise to the points
% lat = lat + randn(size(lat)) * 100;
% lon = lon + randn(size(lon)) * 100;
% 
% % convert to Cartesian coordinates
% [x, y, z] = sph2cart(lon/180*pi, lat/180*pi, r);
% x = cat(1, x, y, z);
% clear y z

% compute MDS on the surface of the matrix
opt.fronorm = 1e-3;
% opt.frorel = .001;
[lat, lon, err, stopCondition, x] = smdscale(d, [], opt);

% plot initialization of the algorithm
hold on
plot3(x(1, :), x(2, :), x(3, :), 'r.')
for I = 1:N
    plot3([x(1, I) x0(1, I)], [x(2, I) x0(2, I)], [x(3, I) x0(3, I)], 'r')
end
axis equal

% convert solution's latitude/longitude to Cartesian coordinates
[x, y, z] = sph2cart(lon/180*pi, lat/180*pi, r);
x = cat(1, x, y, z);
clear y z

% align points with procrustes
[~, x] = procrustes(x0', x');
x = x';

% plot solution
hold on
plot3(x(1, :), x(2, :), x(3, :), 'bo')
axis equal

% plot error
figure
plot(err)
xlabel('iteration')
ylabel('Frobenius norm of distance matrix error')
