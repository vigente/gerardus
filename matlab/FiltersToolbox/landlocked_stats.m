function stats = landlocked_stats(lab, d, dict)
% LANDLOCKED_STATS  Assess whether sections of voxels are surrounded by
% other sections (land) or by background (water)
%
% STATS = LANDLOCKED_STATS(LAB, D, DICT)
%
%   LAB is a segmentation where background voxels have value 0, and the
%   rest of voxels are partitioned in sections, each with a different
%   value. All the voxels within the same section share the same value.
%
%   D is a sparse or full distance matrix where D(i,j)>0 means that voxels
%   i and j are connected.
%
%   DICT is a sparse vector to convert from image linear indices to matrix
%   D indices.
%
%   D and DICT can be computed with function seg2dmat().
%
%   STATS is a struct with the results:
%
%     STATS.islandlocked: bool to tell whether the corresponding section is
%                         landlocked
%
%     STATS.nbound:       number of voxels in the outer boundary of the
%                         section
%
%     STATS.nwater:       number of voxels in the outer boundary that are
%                         touching the background

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.1.0
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
error(nargchk(3, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% figure out whether the data is 2D or 3D, because if it's 2D, a landlocked
% voxel has degree 8, but if it's 3D, it needs degree 26
if (size(lab, 3) > 1)
    % data is 3D
    degmax = 26;
else
    % data is 2D
    degmax = 8;
end

% number of labels
N = max(lab(:));

% compute degree of each voxel in the segmentation
deg = sum(d>0, 2);

if (max(deg) > degmax)
    error('Too many connections in D. Each voxel can only be connected to the adjacent ones')
end

% init output
stats.islandlocked = false(1, N);
stats.nbound = zeros(1, N);
stats.nwater = zeros(1, N);

% loop each label
for I = 1:N
    % distance matrix indices of voxels in the current label
    idx = dict(lab == I);
    
    % number of voxels that are on the outer boundary of the label, but not
    % touching other labels
    stats.nwater(I) = nnz(deg(idx) ~= degmax);
    
    % if all the voxels have maximum degree, then the label is landlocked
    stats.islandlocked(I) = stats.nwater(I) == 0;
    
    % create a smaller distance matrix for only the voxels in the label
    dlab = d(idx, idx);
    
    % compute degree of each voxel in the label if the label is
    % disconnected
    deglab = sum(dlab>0, 2);
    
    % total number of voxels in the outer boundary of the label, whether
    % they touch other labels or not
    stats.nbound(I) = nnz(deglab ~= degmax);
end
