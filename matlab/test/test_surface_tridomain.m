% test_surface_tridomain.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Planar rectangular domain giving number of points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% without extending the domain size, odd number of points
[tri, uv] = surface_tridomain('rect', 'num', [5 3], [-1 2], [1.5 3]);

% plot mesh
hold off
trisurf(tri, uv(:, 1), uv(:, 2), zeros(size(uv, 1), 1))
view(2)

% extending the domain size. The spacing should not change, so there will
% be more points
[tri2, uv2] = surface_tridomain('rect', 'num', [5 3], [-1 2], [1.5 3], [2 3]);

% plot mesh
hold off
trisurf(tri2, uv2(:, 1), uv2(:, 2), zeros(size(uv2, 1), 1))
hold on
trisurf(tri, uv(:, 1), uv(:, 2), zeros(size(uv, 1), 1), 1)
view(2)


% without extending the domain size, even number of points
[tri, uv] = surface_tridomain('rect', 'num', [6 4], [-1 2], [1.5 3]);

% plot mesh
hold off
trisurf(tri, uv(:, 1), uv(:, 2), zeros(size(uv, 1), 1))
view(2)

% extending the domain size. The spacing should not change, so there will
% be more points
[tri2, uv2] = surface_tridomain('rect', 'num', [6 4], [-1 2], [1.5 3], [2 3]);

% plot mesh
hold off
trisurf(tri2, uv2(:, 1), uv2(:, 2), zeros(size(uv2, 1), 1))
hold on
trisurf(tri, uv(:, 1), uv(:, 2), zeros(size(uv, 1), 1), 1)
view(2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Planar rectangular domain giving step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% without extending the domain size, and with the Y step size not fitting
% exactly into the provided space
[tri, uv] = surface_tridomain('rect', 'step', [.5 .3], [-1 2], [1.5 3]);

% plot mesh
hold off
trisurf(tri, uv(:, 1), uv(:, 2), zeros(size(uv, 1), 1))
view(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spherical domain with uniform angle sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% without extending the domain size
[tri, uv] = surface_tridomain('sphang', 'num', [4 4]);

% Cartesian coordinates to the points on the sphere
[x, y, z] = sph2cart(uv(:, 2), uv(:, 1), 1);

% plot mesh
trisurf(tri, x, y, z)


%% finer sampling

% without extending the domain size
[tri, uv] = surface_tridomain('sphang', 'num', 100);

% Cartesian coordinates to the points on the sphere
[x, y, z] = sph2cart(uv(:, 2), uv(:, 1), 1);

% plot mesh
trisurf(tri, x, y, z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spherical domain with uniform point distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xyz,tri]=meshunitsphere(1);

% plot mesh
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))


% without extending the domain size
[tri, uv] = surface_tridomain('sphang', 'num', [4 4]);

% Cartesian coordinates to the points on the sphere
[x, y, z] = sph2cart(uv(:, 2), uv(:, 1), 1);

% plot mesh
trisurf(tri, x, y, z)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Third party libraries for uniform triangular mesh on unit sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DistMesh Toolbox by Per-Olof Persson

% it doesn't allow an easy selection of the number of vertices. The minimum
% number of points seems to be 132 for h0=0.36. Larger values of h0 make
% the algorithm oscillate
%
% http://persson.berkeley.edu/distmesh/
fd=@(p) dsphere(p,0,0,0,1);
[xyz,tri]=distmeshsurface(fd,@huniform,0.36,1.1*[-1,-1,-1;1,1,1]);

%% dist toolbox by Anton Semechko

% http://www.mathworks.co.uk/matlabcentral/fileexchange/37004-uniform-sampling-of-a-sphere

% Uniformly distribute 162 particles across the surface of the unit sphere 
[xyz,tri,~,Ue]=ParticleSampleSphere('N',20); % this operation takes ~8 sec on my machine (6GB RAM & 2.8 GHz processor)

% plot mesh
hold off
trisurf(tri, xyz(:, 1), xyz(:, 2), xyz(:, 3))
axis equal

[lon, lat] = cart2sph(xyz(:, 1), xyz(:, 2), xyz(:, 3));

% compute spherical length of edge connections
d = dmatrix_sphmesh(tri, [lat lon]);

% check deviation from equilength
boxplot(d(d~=0))
