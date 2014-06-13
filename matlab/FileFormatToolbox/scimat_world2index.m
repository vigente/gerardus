function idx = scimat_world2index(x, scimat, CHOP)
% SCIMAT_WORLD2INDEX  Convert real world coordinates to image indices for
% the SCIMAT image struct that we use in Gerardus.
% 
%   Function scimat_world2index() converts the coordinates of a voxel given
%   as real world coordinates [x, y, z] into index coordinates 
%   [row, column, slice].
%
%      [x, y, z] -> [r, c, s]
%
%   This agrees with Matlab's convention that images are expected to be
%   (r, c, s) <-> (y, x, z), but point coordinates are given in the
%   (x, y, z)-order.
%
%   For points that are not within the data volume, the returned
%   indices are "NaN".
%
%   Note also that the indices are not rounded, to allow for sub-pixel
%   accuracy. If integer indices are required, then just use round(idx).
%
% IDX = scimat_world2index(X, SCIMAT)
%
%   X is a 3-column matrix where each row contains the real world
%   (x,y,z)-coordinates of a point.
%
%   IDX has the same size as X, and the voxel indices in 
%   (row, column, slice)-order, that corresponds to (y, x, z)-order.
%
%   SCIMAT is a struct with the image space metadata, i.e. spacing, offset
%   and orientation (see "help scimat" for details). SCIMAT.data (the fild
%   that contains the image itself) is not used by the function, and thus
%   can be present or absent. Note that Matlab will pass SCIMAT.data by
%   reference, so passing the whole image does not require more memory or
%   slow the function down.
%
% IDX = scimat_world2index(..., CHOP)
%
%   CHOP is a flag to convert points outside the image volume to NaNs. By
%   default, CHOP=true.
%
%
% Example:
%
% >> idx = scimat_world2index([.01, .011, .02], scimat)
%
% idx =
%
%     55   189   780
%
% See also: scimat_index2world, scimat_load, scimat_im2scimat.
    
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009-2014 University of Oxford
% Version: 0.4.0
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

if (size(x, 2) ~= 3)
    error('X must be a 3-column matrix, so that each row has the 3D coordinates of a point')
end

% init output
idx = zeros(size(x));

% extract parameters
xmin = [scimat.axis.min];
dx = [scimat.axis.spacing];
n = [scimat.axis.size];
orig = xmin + dx/2;

% remove dummy dimension, if present
if (length(xmin) == 4)
    xmin = xmin(2:end);
    dx = dx(2:end);
    orig = orig(2:end);
end

% number of dimensions (we expect D=3, but in case this gets more general)
D = length(dx);
if (D ~= 3)
    error('Input SCIMAT volume expected to have 3 dimensions.')
end

% (x, y, z) => (y, x, z)
x = x(:, [2 1 3]);

% convert real world coordinates to indices. The reason for the loop is
% that x may have more than 1 point to convert, and this is probably more
% efficient than having to replicate orig and dx to fit the size of x
for I = 1:D
    idx(:, I) = (x(:, I) - orig(I)) / dx(I) + 1;
end

% find which coordinates are outside the volume
if CHOP
    for I = 1:D
        idx(idx(:, I) < 0.5 | idx(:, I) > n(I)+0.5, I) = NaN;
    end
end
