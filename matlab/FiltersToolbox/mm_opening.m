function seg = mm_opening(seg, RAD, N)
% MM_OPENING Mathematical Morphology operator: Opening (erosion+dilation)
% with removal of small components
%
% SEG2 = mm_opening(SEG)
%
%   SEG is a 3D volume with a binary segmentation.
%
%   SEG2 is a smoothed version of SEG, after applying an opening
%   morphological operator (erosion + dilation with radius 1).
%
% ... = mm_opening(..., RAD, N)
%
%   RAD is a scalar with the radius in voxels of the erosion and dilation.
%
%   N is a scalar with the number of largest components to keep after the
%   erosion. Eroding tends to disconnect a segmented object into several
%   components, many usually noise-like. If N<Inf, right after erosion all
%   but the N largest connected components are kept. By default, N=Inf and
%   no components are removed.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
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
narginchk(1, 3);
nargoutchk(0, 1);

% defaults
if (nargin < 2 || isempty(RAD))
    RAD = 1;
end
if (nargin < 3 || isempty(N))
    N = Inf;
end

% erode the segmentation
seg = itk_imfilter('bwerode', seg, RAD, 1);

% keep only the largest or a few of the largest components
if (N < Inf)
    % connected components
    cc = bwconncomp(seg);
    
    % number of components we are going to keep
    N = min(N, length(cc));
    
    % size of each component
    n = cellfun(@(x) length(x), cc.PixelIdxList);
    [~, idx] = sort(n, 'descend');
    seg(:) = 0;
    seg(cc.PixelIdxList{idx(1:N)}) = 1;
end

% erode the segmentation
seg = itk_imfilter('bwdilate', seg, 5, 1);
