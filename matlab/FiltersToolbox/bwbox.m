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
%      0     0     0     0     0
%      0     0     0     0     0
%      0     0     0     0     0
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
% >> bwbox(bw)
% 
% ans =
% 
%      2     6
%      1     5
%
% See also: regionprops(..., 'BoundingBox')

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.1.0
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

% initialise output
idx = zeros(ndims(bw), 2);

for I = 1:ndims(bw)

    % project all voxels onto the I-th dimension, so that we can find the
    % boundaries for this dimension
    aux = bw;
    for J = 1:ndims(bw)-1
        aux = shiftdim(aux, 1);
        aux = sum(aux, 1);
    end

    % boundaries for this dimension
    aux = find(aux);
    if length(aux) == 1
        idx(I, :) = aux(1);
    else
        idx(I, :) = aux([1 end]);
    end
    
    % extend the boundaries with a margin
    idx(I, 1) = max(1, idx(I, 1)-m);
    idx(I, 2) = min(size(bw, 1), idx(I, 2)+m);
    
    % shift image to the next dimension
    bw = shiftdim(bw, 1);
    
end
