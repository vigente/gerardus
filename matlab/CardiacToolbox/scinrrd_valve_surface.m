function nrrd = scinrrd_valve_surface(nrrd, x, HOR)
% SCINRRD_VALVE_SURFACE  Inter atrio-ventricular surface
%
% NRRD = SCINRRD_VALVE_SURFACE(NRRD0, X)
%
%   NRRD0 is the SCI NRRD struct that contains the cardiac Magnetic
%   Resonance Image (MRI).
%
%   X is a 3-row matrix. Each column has the coordinates of a point from
%   one of the valve annula.
%
%   NRRD is a SCI NRRD struct with a segmentation of the surface that
%   interpolates the valve annula points.
%
% NRRD = SCINRRD_VALVE_SURFACE(NRRD0, X, HOR)
%
%   HOR is a flag. If HOR=true, the valve points are rotated to make the
%   first 2 eigenvectors of their principal components horizontal. By
%   default, HOR=false.
%
%   Note on SCI NRRD: Software applications developed at the University of
%   Utah Scientific Computing and Imaging (SCI) Institute, e.g. Seg3D,
%   internally use NRRD volumes to store medical data.
%
%   When label volumes (segmentation masks) are saved to a Matlab file
%   (.mat), they use a struct called "scirunnrrd" to store all the NRRD
%   information:
%
%   >>  scirunnrrd
%
%   scirunnrrd = 
%
%          data: [4-D uint8]
%          axis: [4x1 struct]
%      property: []

% Author: Ramon Casero.
% Copyright © 2010 University of Oxford
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
error( nargchk( 2, 3, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% defaults
if (nargin < 3 || isempty(HOR))
    HOR=false;
end

% extract data volume and convert to boolean to save space
mask = boolean(nrrd.data);

% we want the valve plane to be as horizontal as possible before doing the
% interpolation
if HOR
    m = mean(x, 2);
    x = x - m(:, ones(1, size(x, 2)));
    eigv = pts_pca(x);
    x = eigv' * x;
else
    m = zeros(3, 1);
    eigv = eye(3);
end

% the warp will be defined as from the xy-plane to the surface defined by
% the valves
s = x(1:2, :);
t = x(3, :);

% find the rotation on the XY axes
vx = eigv' * [1 0 0]';
vy = eigv' * [0 1 0]';

% get the projection of the rotated vectors on the z-axis. From this, we
% can compute how much we need to enlarge the pre-rotated domain so that
% the rotated domain covers all the region of interest
K = 1/min(cosd(asind(vx(3))), cosd(asind(vy(3))));

%% create grid of points that covers the axial plane

% get image vertices
vmin = [nrrd.axis.min];
vmax = [nrrd.axis.max];

delta = (vmax - vmin)/2;
vmin = vmin - delta;
vmax = vmax + delta;

% get voxel size
res = [nrrd.axis.spacing];

% generate 3D grid of coordinates: note the inversion of coordinates,
% necessary so that xg will change with columns, and yg with rows
[gy, gx] = ndgrid( vmin(1):res(1):vmax(1), ...
    vmin(2):res(2):vmax(2));

% rotate grid
gxy = [gx(:)'; gy(:)'; m(3)+zeros(1, numel(gx))];
gxy = gxy - m(:, ones(1, size(gxy, 2)));
gxy = eigv' * gxy;

% compute interpolating surface
y = pts_tps_map(s', t', gxy(1:2, :)');

% transform the interpolated surface to the original coordinate system
gxy(3, :) = y;
gxy = eigv * gxy;
gxy = gxy + m(:, ones(1, size(gxy, 2)));

% convert real world coordinates to indices
idx = round(scinrrd_world2index(gxy', nrrd.axis));

% remove points outside the volume
badidx = isnan(sum(idx, 2));
idx = idx(~badidx, :);

% construct new segmentation mask. Note that zz has type double.
mask = mask * 0;
mask( sub2ind( size( mask ), ...
    idx(:, 1 ), idx(:, 2 ), floor(idx(:, 3)) ) ) = 1;
mask( sub2ind( size( mask ), ...
    idx(:, 1 ), idx(:, 2 ), ceil(idx(:, 3)) ) ) = 1;
nrrd.data = mask;
