function idx = bwbox(bw, m)
% BWBOX  Find tight box around segmentation (and possibly add a margin)
%
% IDX = bwbox(BW)
%
%   BW is a segmentation volume of any number of dimensions. Voxels with 
%   value 0 belong to the background, and anything else, to the object.
%
%   IDX is a matrix where each row contains the first and last indices of
%   foreground voxels along the corresponding dimension. All indices
%   together form a hyperdimensional box that tighly contains the object.
%
% IDX = bwbox(BW, M)
%
%   M is a scalar. The box will be expanded by M voxels in each dimension.
%
% >> bw =
% 
%      0     0     0     0     0
%      0     0     0     0     0
%      0     1     1     1     1
%      0     1     1     1     1
%      0     1     1     1     1
%      0     0     0     0     0
%      0     0     0     0     0
%
% >> bwbox(bw)
% 
% ans =
% 
%      3     5
%      2     5
%
% >> bwbox(bw, 1)
% 
% ans =
% 
%      2     6
%      1     5
%
% See also: regionprops(..., 'BoundingBox')

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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 2);
nargoutchk(0, 1);

% defaults
if (nargin < 2 || isempty(m))
    m = 0;
end

% compute boundary
stats = regionprops(bw, 'BoundingBox');

% first index of bounding box
idx = stats.BoundingBox(1:ndims(bw)) + 0.5;

% size of bounding box
sz = stats.BoundingBox(ndims(bw)+1:end);

% note that regionprops() output swaps the 1st and 2nd dimensions to
% present x, y, z results, instead of row, col, slice
if (ndims(bw) == 2)
    idx = idx([2 1]);
    sz = sz([2 1]);
elseif (ndims(bw) > 2)
    idx = idx([2 1 ndims(bw):end]);
    sz = sz([2 1 ndims(bw):end]);
end

% format output
idx = [idx; idx+sz-1]';

% extend the boundaries with a margin
idx(:, 1) = max(1, idx(:, 1)-m);
idx(:, 2) = min(size(bw)', idx(:, 2)+m);

