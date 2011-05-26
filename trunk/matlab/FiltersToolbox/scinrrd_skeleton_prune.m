function nrrdsk = scinrrd_skeleton_prune(nrrdsk, cc, minlen)
% SCINRRD_SKELETON_PRUNE  Prune branches in a segmentation's skeletonization
%
% NRRDPR = SCINRRD_SKELETON_PRUNE(NRRDSK, CC, MINLEN)
%
%   NRRDSK is an SCI NRRD struct. NRRD.data contains the result of running
%   a skeletonization algorithm on a binary segmentation, e.g.
%
%     >> nrrdsk.data = itk_imfilter('skel', nrrd);
%
%   CC is a struct produced by function skeleton_label().
%
%   MINLEN is a scalar with the minimum length in voxels for a leaf. Any
%   leaf shorter than MINLEN will be pruned.
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
% Copyright © 2011 University of Oxford
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
error(nargchk(3, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (isempty(minlen))
    minlen = 0;
end

% loop each branch
for I = 1:cc.NumObjects
    
    % if leaf branch is too short...
    if (cc.IsLeaf(I) && length(cc.PixelIdxList{I}) < minlen)
        
        % chop it off
        nrrdsk.data(cc.PixelIdxList{I}) = 0;
    end
    
end
