function x = scimat_index2world(idx, ax, CHOP)
% scimat_index2world  Convert image indices to real world coordinates for
% the SCIMAT image struct that we use in Gerardus.
%
%   Function scimat_index2world() converts the coordinates of a voxel given
%   as index coordinates [row, column, slice] into its real world
%   coordinates [x, y, z].
%
%      [r, c, s] -> [x, y, z]
%
%   This agrees with Matlab's convention that images are expected to be
%   (r, c, s) <-> (y, x, z), but point coordinates are given in the
%   (x, y, z)-order.
%
%   For points that are not within the data volume, the returned
%   coordinates are "NaN".
%
% X = scimat_index2world(IDX, AXIS)
%
%   IDX has the same size as X, and the voxel indices in 
%   (row, column, slice)-order, that corresponds to (y, x, z)-order.
%
%   AXIS is the scimat.axis field from a SCIMAT struct (see "help
%   scimat_load" for details).
%
%   X is a 3-column matrix where each row contains the real world
%   (x, y, z)-coordinates of a point.
%
% IDX = scimat_index2world(..., CHOP)
%
%   CHOP is a flag to convert points outside the image volume to NaNs. By
%   default, CHOP=true.
%
%
% Example:
%
% >> x = scimat_index2world([55 189 780], scimat.axis)
%
% x =
%
%     .0100, .0110, .0200
%
% See also: scimat_world2index, scimat_load, scimat_im2scimat.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009-2014 University of Oxford
% Version: 0.3.0
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
narginchk(2, 3);
nargoutchk(0, 1);

% defaults
if (nargin < 3 || isempty(CHOP))
    CHOP = true;
end

if (size(idx, 2) ~= 3)
    error('IDX must be a 3-column matrix, so that each row has the 3 indices of a voxel')
end

% init output
x = zeros(size(idx));

% extract parameters
xmin = [ax.min];
dx = [ax.spacing];
n = [ax.size];
% remove dummy dimension, if present
if (length(xmin) == 4)
    xmin = xmin(2:end);
    dx = dx(2:end);
end

% number of dimensions (we expect D=3, but in case this gets more general)
D = length(dx);

% find coordinates that are outside the volume
if CHOP
    for I = 1:D
        idx(idx(:, I) < 0.5 | idx(:, I) > n(I)+0.5, I) = NaN;
    end
end

% convert indices to real world coordinates
for I = 1:D
    x(:, I) = (idx(:, I) - 1) * dx(I) + xmin(I) + dx(I)/2;
end

% (y, x, z) => (x, y, z)
x = x(:, [2 1 3]);
