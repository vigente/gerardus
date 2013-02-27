% test_surface_interpolation.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.2.0
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
%% PCA parameterisation for all interpolation methods on a plane

load('test/data/valve-annula-points-009.mat');

% plot points
plot3(x(:, 1), x(:, 2), x(:, 3), '*')
axis equal

% PCA parameterisation
param.type = 'pca';

% compute interpolating surface for all interpolation methods, and plot the
% results on the same figure so that obvious errors stand out
hold off
for type = {'tps', 'tsi', 'gridfit', 'mba', 'mbae'}
    
    % compute interpolating surface
    interp.type = type{1};
    interp.res = [1 1] * 0.01/50;
    interp.klim = 1;
    xi = surface_interpolation(x, param, interp);
    
    % plot result
    surf(xi(:, :, 1), xi(:, :, 2), xi(:, :, 3))
    axis equal
    hold on
    
end

% now plot separately to spot problems with individual methods
hold off
for type = {'tps', 'tsi', 'gridfit', 'mba', 'mbae'}
    
    % compute interpolating surface
    interp.type = type{1};
    interp.res = [1 1] * 0.01/50;
    interp.klim = 2;
    xi = surface_interpolation(x, param, interp);
    
    % plot result
    surf(xi(:, :, 1), xi(:, :, 2), xi(:, :, 3))
    title(['Method: ' interp.type])
    axis equal
    pause
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Isomap parameterisation for all interpolation methods on a plane

% Isomap parameterisation
param.type = 'isomap';
param.size = 20;

% compute interpolating surface for all interpolation methods, and plot the
% results on the same figure so that obvious errors stand out
hold off
for type = {'tps', 'tsi', 'gridfit', 'mba', 'mbae'}
    
    % compute interpolating surface
    interp.type = type{1};
    interp.res = [1 1] * 0.01/50;
    interp.klim = 1;
    xi = surface_interpolation(x, param, interp);
    
    % plot result
    surf(xi(:, :, 1), xi(:, :, 2), xi(:, :, 3))
    hold on
    plot3(x(:, 1), x(:, 2), x(:, 3), '*')
    axis equal
    
end

% plot separately to spot problems with individual methods
for type = {'tps', 'tsi', 'gridfit', 'mba', 'mbae'}
    
    % compute interpolating surface
    interp.type = type{1};
    interp.res = [1 1] * 0.01/50;
    interp.klim = 1;
    xi = surface_interpolation(x, param, interp);
    
    % plot result
    hold off
    surf(xi(:, :, 1), xi(:, :, 2), xi(:, :, 3))
    title(['Method: ' interp.type])
    hold on
    plot3(x(:, 1), x(:, 2), x(:, 3), '*')
    axis equal
    
    pause
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We have a scattered set of points and we want to compute an 
%% interpolation surface topologically equivalent to a sphere

% load discrete point set with an associated local neighbourhood
scimat = scinrrd_load('test/data/thick-slice-points-seg.mat');

% compute a distance matrix that makes sense for data from thick-slices
[param.d, xyz] = scimat_dmatrix_thickslice(scimat, 3);

% plot neighbourhood
hold off
gplot3d(param.d, xyz)
axis equal

% interpolate surface
param.type = 'sphisomap'; % Isomap parameterisation on a sphere
param.maxiter = 50; % stopping criterion
param.dtype = 'none';
param.dtype = 'full';
param.dtype = 'sparse';
param.dtype = 'full+sparse';
interp.type = 'sphspline';
interp.res = 2*pi/20; % arc increment (in radians) when sampling the sphere

[xi, em, x, lon, lat] = ...
    surface_interpolation(xyz, param, interp);

% plot interpolation result
hold off
surf(xi(:, :, 1), xi(:, :, 2), xi(:, :, 3))
hold on
plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), '*')
axis equal
