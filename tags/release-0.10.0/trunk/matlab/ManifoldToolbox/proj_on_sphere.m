function [lat, lon, sphrad] = proj_on_sphere(x)
% PROJ_ON_SPHERE  Project a point configuration onto a sphere
%
% [LAT, LON, SPHRAD] = proj_on_sphere(X)
%
%   X is a 3-row matrix. Each colum has the (x,y,z)-coordinates of a point.
%
%   LAT, LON are vectors with the latitude and longitude of the points in X
%   projected onto a sphere. The centre of the sphere is the centroid of X.
%
%   SPHRAD is the median of the distance of the points in X to the
%   centroid. It can be used as the radius of the sphere that best contains
%   the projection of X in terms of distances between the the points.

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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 1);
nargoutchk(0, 3);

if (size(x, 2) ~= 3)
    error('X must be a matrix with 3-columns')
end

% project 3D Cartesian coordinates onto
% the surface of a sphere lon, lat given in radians, centered around 0
[lon, lat, r] = cart2sph(...
    x(:, 1) - mean(x(:, 1)), ...
    x(:, 2) - mean(x(:, 2)), ...
    x(:, 3) - mean(x(:, 3)));

% we use the median radius of the points in the configuration as the radius
% of the sphere that best can contains their projections
sphrad = median(r);
