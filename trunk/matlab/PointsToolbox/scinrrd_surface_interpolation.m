function [nrrd, uv, x] = scinrrd_surface_interpolation(nrrd, x, param, interp)
% SCINRRD_SURFACE_INTERPOLATION  Interpolate a surface and create a
% segmentation mask from a scattered set of points
%
% NRRD = scinrrd_surface_interpolation(NRRD0, X)
%
%   NRRD0 is the SCI NRRD struct that contains the image.
%
%   X is a 3-column matrix. Each row has the coordinates of a point that
%   belongs to the surface we want to interpolate.
%
%   NRRD is a SCI NRRD struct with a segmentation of the surface that
%   interpolates the points in X.
%
% [NRRD, UV, X2] = scinrrd_surface_interpolation(NRRD0, X, PARAM, INTERP)
%
%   UV is a 2-column matrix with the parameterisation of X, (U, V)->X. In
%   planar parameterisations, UV has units of meters. In spherical
%   parameterisations, UV=[LON, LAT], in units of radians.
%
%   X2: Some methods (e.g. 'mbae') add extra points to X before computing
%   the interpolants. X2 contains the points actually used for
%   interpolation.
%
%   PARAM is a struct that describes the method used to parametrise the
%   surface and the set of points X, and its parameters. For details, see
%   help surface_interpolation.m.
%
%   INTERP is a struct that describes the interpolation method, and its
%   parameters. For options, see help surface_interpolation.m.
%
% See also: surface_interpolation.m.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2013 University of Oxford
% Version: 0.7.0
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

% check arguments
narginchk(2, 4);
nargoutchk(0, 3);

% defaults
if (nargin < 3)
    param = [];
end
if (nargin < 4)
    interp = [];
end

% get voxel size
interp.res = [nrrd.axis.spacing];

%% compute interpolating surface
[y, uv] = surface_interpolation(x, param, interp);

%% map interpolated surface points to voxels

% dimensions of the volume with coordinates of surface points
[R, C, S] = size(y);

% convert real world coordinates to indices
idx = round(scinrrd_world2index(reshape(y, R*C, S), nrrd.axis));

% remove points outside the volume
badidx = isnan(sum(idx, 2));
idx = idx(~badidx, :);

% reset segmentation mask
nrrd.data = zeros([nrrd.axis.size], 'uint8');

% construct new segmentation mask
nrrd.data(sub2ind(size(nrrd.data), ...
    idx(:, 1), idx(:, 2), floor(idx(:, 3)))) = 1;
nrrd.data(sub2ind(size(nrrd.data), ...
    idx(:, 1), idx(:, 2), ceil(idx(:, 3)))) = 1;
