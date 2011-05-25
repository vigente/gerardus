function [im, gx, gy, gz] = scinrrd_intersect_plane(nrrd, m, v)
% SCINRRD_INTERSECT_PLANE  Intersection of a plane with an image volume
%
% [IM, GX, GY, GZ] = SCINRRD_INTERSECT_PLANE(NRRD, M, V)
%
%   IM is an image that displays the intersection of the plane with the
%   image volume in SCI NRRD format. Voxels that fall outside the image
%   volume are returned as NaN.
%
%   GX, GY, GZ are matrices of the same size as IM, and contain the
%   Cartesian coordinates of the voxels in IM. You can visualize the
%   resulting plane using
%
%     >> surf(gx, gy, gz, im, 'EdgeColor', 'none')
%
%   The plane is uniquely defined in 3D space using a point and a vector:
%
%     * M is a 3-vector with the coordinates of a point contained in the
%       plane. It is assumed that M is within the image boundaries
%   
%     * V is a 3-vector that represents a normalized vector orthogonal to
%       the plane
%   
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

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2011 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
error(nargchk(3, 3, nargin, 'struct'));
error(nargoutchk(0, 4, nargout, 'struct'));

% remove dummy dimension of data if necessary
nrrd = scinrrd_squeeze(nrrd, true);

% the longest distance within the image volume is the length of the largest
% diagonal. We compute it in voxel units
lmax = ceil(sqrt(sum(([nrrd.axis.size] - 1).^2)));

% create horizontal grid for the plane sampling points. We are going to have one
% sampling point per voxel. The gird is centered on 0
[gc, gr] = meshgrid(-lmax:lmax, -lmax:lmax);
gs = 0 * gc;

% compute a rotation matrix to map the Z-axis to v. Note that v is given in
% x, y, z coordinates, but we want to work with r, c, s indices, i.e. y, x,
% z coordinates
v = v([2 1 3]);
rotmat = vec2rotmat(v(:));

% rotate the grid sampling points accordingly
grcs = rotmat * [gr(:) gc(:) gs(:)]';

% compute index coordinates of real world coordinates of the plane point
idxm = scinrrd_world2index(m, nrrd.axis);

% center rotated grid on the plane point m
grcs(1, :) = grcs(1, :) + idxm(1);
grcs(2, :) = grcs(2, :) + idxm(2);
grcs(3, :) = grcs(3, :) + idxm(3);

% round coordinates, so that we are sampling at voxel centers and don't
% need to interpolate
grcs = round(grcs);

% sampling points that are outside the image domain
idxout = (grcs(1, :) < 1) | (grcs(1, :) > nrrd.axis(1).size) ...
    | (grcs(2, :) < 1) | (grcs(2, :) > nrrd.axis(2).size) ...
    | (grcs(3, :) < 1) | (grcs(3, :) > nrrd.axis(3).size);

%sampling points that are inside
idxin = ~idxout;

% rcs indices => linear indices
idx = sub2ind(size(nrrd.data), grcs(1, idxin), grcs(2, idxin), ...
    grcs(3, idxin));

% sample image volume with the rotated and translated plane
im = nan(size(grcs, 2), 1);
im(idxin) = nrrd.data(idx);

% compute real world coordinates for the sampling points
gxyz = scinrrd_index2world(grcs', nrrd.axis)';

% reshape the sampled points to get again a grid distribution
im = reshape(im, size(gc));
gx = reshape(gxyz(1, :), size(gr));
gy = reshape(gxyz(2, :), size(gc));
gz = reshape(gxyz(3, :), size(gs));

% find columns where all elements are NaNs
idxout = all(isnan(im), 1);

% remove those columns
im = im(:, ~idxout);
gx = gx(:, ~idxout);
gy = gy(:, ~idxout);
gz = gz(:, ~idxout);

% find rows where all elements are NaNs
idxout = all(isnan(im), 2);

% remove those columns
im = im(~idxout, :);
gx = gx(~idxout, :);
gy = gy(~idxout, :);
gz = gz(~idxout, :);
