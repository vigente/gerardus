function bw = bwrmsmallcomp(bw, nobj)
% BWRMSMALLCOMP  Remove small connected components from a segmentation
%
%   This function is useful to remove segmentation noise.
%
% BW2 = bwrmsmallcomp(BW, NOBJ)
%
%   BW is a binary segmentation.
%
%   NOBJ is a scalar with the number of connected components we want to
%   keep. By default, NOBJ=1, and only the largest connected components is
%   kept.
%
%   BW2 is a binary segmentation with the same size as BW, but where all
%   connected components except for the NOBJ largest ones have been
%   removed.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
% Version: 0.1.1
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
if (nargin < 2 || isempty(nobj))
    nobj = 1;
end

% get connected components
cc = bwconncomp(bw);

% keep only the largest component to remove background noise
%
% note: it's better to clear the whole image, and then add the largest
% components, than trying to delete the smaller components. The latter
% doesn't remove all the noise, for some reason.
len = cellfun(@length, cc.PixelIdxList);
[~, idx] = sort(len, 2, 'descend');
bw = zeros(size(bw), 'uint8');
if ~isempty(idx)
    idx = cat(1, cc.PixelIdxList{idx(1:nobj)});
    bw(idx) = 1;
end
