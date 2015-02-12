function x = scimat_index2world(idx, scimat, CHOP)
% SCIMAT_INDEX2WORLD  Convert image indices to real world coordinates for
% the SCIMAT image struct that we use in Gerardus.
%
%   Function SCIMAT_INDEX2WORLD() converts the coordinates of a voxel given
%   as index coordinates [row, column, slice, frame] into its real world
%   coordinates [x, y, z, t].
%
%      [r, c, s, f] -> [x, y, z, t]
%
%   This agrees with Matlab's convention that images are expected to be 
%   (r, c) <-> (y, x), but point coordinates are given in the (x, y)-order.
%
%   This function can also be applied to images that are not 4D, and in
%   that case, index and real world coordinates will have the same number
%   of elements as dimensions the image has.
%
%   The relation between indices IDX and real world coordinates X is
%
%     X = s.*(IDX-1)*R + t
%
%   where s is the voxel size, R the rotation matrix, and t the
%   image offset.
%
%   For points that are not within the data volume, the returned
%   coordinates are "NaN".
%
% X = SCIMAT_INDEX2WORLD(IDX, SCIMAT)
%
%   IDX has one element per image dimension. If the image is 4D, IDX is
%   given in (row, column, slice, frame)-order.
%
%   SCIMAT is a struct with the image and axis metadata, i.e. spacing,
%   offset and orientation (see "help scimat" for details). SCIMAT.data
%   (the field that contains the image itself) is not used by the function,
%   and thus can be present or absent. Note that Matlab will pass
%   SCIMAT.data by reference, so passing the whole image does not require
%   more memory and does not slow the function down.
%
%   X is a 3-column matrix where each row contains the real world
%   (x, y, z, t)-coordinates of a point.
%
% IDX = SCIMAT_INDEX2WORLD(..., CHOP)
%
%   CHOP is a flag to convert points outside the image volume to NaNs. By
%   default, CHOP=true.
%
%
% Example:
%
% >> x = scimat_index2world([55 189 780], scimat)
%
% x =
%
%     .0100, .0110, .0200
%
% See also: scimat, scimat_world2index, scimat_load, scimat_im2scimat.

% Authors: Ramon Casero <rcasero@gmail.com>, 
% Benjamin Villard <b.016434@gmail.com>,
% Christopher Kelly  <christopher.kelly28@googlemail.com>
% Copyright © 2009-2015 University of Oxford
% Version: 0.5.0
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
if (~isfield(scimat, 'rotmat'))
    scimat.rotmat = [];
end

% extract parameters
xmin = [scimat.axis.min];
dx = [scimat.axis.spacing];
n = [scimat.axis.size];
orig = xmin + dx/2;
R = scimat.rotmat;

% number of dimensions
D = length(scimat.axis);

% find coordinates that are outside the volume
if CHOP
    for I = 1:D
        idx(idx(:, I) < 0.5 | idx(:, I) > n(I)+0.5, I) = NaN;
    end
end

%% convert indices to real world coordinates

% x = (i-1) * dx
x = (idx - 1) .* repmat(dx, size(idx, 1), 1);

% (y, x) => (x, y)
x = x(:, [2 1 3:end]);

% apply rotation only to the spatial coordinates
if (~isempty(R))
    x(:, 1:size(R, 1)) = x(:, 1:size(R, 1)) * R;
end

% apply offset
x = x + repmat(orig([2 1 3:end]), size(x, 1), 1);
