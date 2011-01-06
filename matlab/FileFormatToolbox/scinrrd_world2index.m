function idx = scinrrd_world2index(x, ax)
% SCINRRD_WORLD2INDEX  Convert real world coordinates to data volume
% indices for NRRD volumes created by SCI applications (e.g. Seg3D)
% 
%   Function SCINRRD_WORLD2INDEX() maps between the real world coordinates
%   of points within the NRRD data volume and the indices of the 
%   image volume used to store the voxel intensity values.
%
%      [x, y, z] -> [r, c, s]
%
%   This function assumes that input world coordinates are given in 
%   (x, y, z)-order. However, output indices are in 
%   (row, column, slice)-order, corresponding to (y, x, z).
%
%   For points that are not within the data volume, the returned
%   indices are "NaN".
%
%   Note also that the indices are not rounded, to allow for sub-pixel
%   accuracy. If integer indices are required, then just use round(idx).
%
% IDX = SCINRRD_WORLD2INDEX(X, AXIS)
%
%   X is a 3-column matrix where each row contains the real world
%   (x,y,z)-coordinates of a point.
%
%   IDX has the same size as X, and the voxel indices in 
%   (row, column, slice)-order, that corresponds to (y, x, z)-order.
%
%   AXIS is the 4x1 struct array nrrd.axis from an SCI NRRD struct. It
%   contains the following fields:
%
%   >> nrrd.axis
%
%   4x1 struct array with fields:
%       size
%       spacing
%       min
%       max
%       center
%       label
%       unit
%
% Example:
%
% >> idx = scinrrd_world2index([.01, .011, .02], scirunnrrd.axis)
%
% idx =
%
%     55   189   780
%
%   Note on SCI NRRD: Software applications developed at the University of
%   Utah Scientific Computing and Imaging (SCI) Institute, e.g. Seg3D,
%   internally use NRRD volumes to store medical data.
%
%   When data or label volumes are saved to a Matlab file (.mat), they use
%   a struct called "scirunnrrd" to store all the NRRD information:
%
%   >>  scirunnrrd
%
%   scirunnrrd = 
%
%          data: [4-D uint8]
%          axis: [4x1 struct]
%      property: []
%
% See also: scinrrd_index2world.
    
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009-2011 University of Oxford
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
error(nargchk(2, 2, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

if (size(x, 2) ~= 3)
    error('X must be a 3-column matrix, so that each row has the 3D coordinates of a point')
end

% init output
idx = zeros(size(x));

% extract parameters
xmin = [ax.min];
dx = [ax.spacing];
n = [ax.size];
% remove dummy dimension, if present
if (length( xmin ) == 4)
    xmin = xmin(2:end);
    dx = dx(2:end);
end

% number of dimensions (we expect D=3, but in case this gets more general)
D = length(dx);

% (x, y, z) => (y, x, z)
x = x(:, [2 1 3]);

% convert real world coordinates to indices
for I = 1:D
    idx(:, I) = (x(:, I) - xmin(I)) / dx(I) + .5;
end


% find which coordinates are outside the volume
for I = 1:D
    idx( idx( :, I ) < 0.5 | idx( :, I ) > n(I)+0.5, I ) = NaN;
end
