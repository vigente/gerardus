function nrrdsk = scinrrd_skeleton_prune(nrrdsk, nrrd, minlen, lratio)
% SCINRRD_SKELETON_PRUNE  Prune branches in a segmentation's skeletonization
%
% This function prunes the leaves of a segmentation skeleton in two steps:
%
%   1. Remove very short leaves
%
%   2. Remove branches created by artifacts in the segmentation. Spurious
%      branches are those roughly as long as the local radius of the main
%      branch they are attached to
%
%
% NRRDPR = SCINRRD_SKELETON_PRUNE(NRRDSK, [], MINLEN)
%
%   This syntax runs step 1 only.
%
%   NRRDSK is an SCI NRRD struct. NRRDSK.data contains the result of
%   running a skeletonization algorithm on a binary segmentation NRRD,
%   e.g.
%
%     >> nrrdsk = nrrd;
%     >> nrrdsk.data = itk_imfilter('skel', nrrd);
%
%   MINLEN is a scalar with the minimum length in voxels for a leaf. Any
%   leaf shorter than MINLEN will be pruned. By default, MINLEN = 5 voxel.
%
%
% NRRDPR = SCINRRD_SKELETON_PRUNE(NRRDSK, NRRD, MINLEN, LRATIO)
%
%   This syntax runs steps 1 and 2.
%
%   NRRD is the SCI NRRD struct mentioned above.
%
%   LRATIO is a scalar. Branches with BL/R < LRATIO will be pruned, where
%   BL is the skeleton branch length and R is the estimated maximum radius
%   of the main branch the leaf is connected to.
%
%
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
% Version: 0.2.2
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
error(nargchk(3, 4, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(minlen))
    minlen = 5;
end
if (nargin < 4 || isempty(lratio))
    lratio = 1.2;
end

%% Step 1: pruning of very short leaf branches

% split skeleton in branches
[~, cc] = skeleton_label(nrrdsk.data, [], [nrrdsk.axis.spacing]);

% loop to prune branches
for I = 1:cc.NumObjects

    if ~(cc.IsLeaf(I) && length(cc.PixelIdxList{I}) < minlen)
        % skip branch
        continue
    end
    
    % get indices from current branch and its neighbours
    idx = cat(1, cc.PixelIdxList{[I cc.BranchNeighbours{I}]});
    
    % linear index => r, c, s indices
    [r, c, s] = ind2sub(size(nrrdsk.data), idx);
    
    % coordinates of a box 
    boxmin = [min(r), min(c), min(s)];
    boxmax = [max(r), max(c), max(s)];
    
    % crop box from whole skeleton segmentation
    box = nrrdsk.data(boxmin(1):boxmax(1), boxmin(2):boxmax(2), ...
        boxmin(3):boxmax(3));
    
    % compute connected components. The current branch and its neighbours
    % should form only 1, but because voxels from other branches can be
    % contained within the box, the number of connected components in the
    % box can be higher
    boxcc = bwconncomp(box);
    
    % save the number of connected components in the box for later
    boxnumcc = boxcc.NumObjects;
    
    % remove the current branch's voxels from the box
    [r, c, s] = ind2sub(size(nrrdsk.data), cc.PixelIdxList{I});
    r = r - boxmin(1) + 1;
    c = c - boxmin(2) + 1;
    s = s - boxmin(3) + 1;
    idx = sub2ind(size(box), r, c, s);
    box(idx) = 0;
    
    % recompute connected components
    boxcc = bwconncomp(box);
    
    % if the number of connected components hasn't changed, then we remove
    % the whole current branch from the skeleton segmentation
    if (boxcc.NumObjects == boxnumcc)
        nrrdsk.data(cc.PixelIdxList{I}) = 0;
    else
        % otherwise, we remove the branch except for the bifurcation voxel
        nrrdsk.data(cc.PixelIdxList{I}(cc. Degree{I} < 3)) = 0;
    end
    
end

%% Step 2: pruning of leaf branches produced by segmentation artifacts

% skip if no full segmentation is provided
if isempty(nrrd)
    return
end

% separate segmentation in branch candidates using the cleaned up skeleton
[~, cc] = skeleton_label(nrrdsk.data, nrrd.data, [nrrd.axis.spacing]);

% init vector for maximum eigenvalues
maxeigd = zeros(1, cc.NumObjects);

% loop each branch
for I = 1:cc.NumObjects
    
    % get index of voxel that connects this branch to other branches
    idx0 = cc.PixelIdxList{I}(cc.Degree{I} > 2);
    
    % we are only interested in leaves that are attached to other branches
    % as potential candidates for pruning
    if (~cc.IsLeaf(I) || isempty(idx0))
        continue
    end
    
    % real world coordinates of the voxel
    [r, c, s] = ind2sub(size(nrrd.data), idx0);
    xyz0 = scinrrd_index2world([r, c, s], nrrd.axis);

    % labels that are connected to this leaf
    idx = cc.BranchNeighbours{I};
    if (isempty(idx)) % skip if this branch is not connected
        continue
    end
    
    % coordinates of skeleton points that are connected to this leaf
    idx = unique(cat(1, cc.PixelIdxList{idx}));

    [r, c, s] = ind2sub(size(nrrd.data), idx);
    xyz = scinrrd_index2world([r, c, s], nrrd.axis);

    % compute eigenvectors of the neighbour branches' skeletons. We are
    % going to use the largest as an estimate of the main vessels direction
    eigv = pts_pca(xyz');
    
    % intersect 3D image with a plane that contains the first skeleton
    % voxel, and that is orthogonal to the main vessel
    [im, gx, gy, gz, xyz0idx] = scinrrd_intersect_plane(nrrd, xyz0, ...
        eigv(:, 1));
    
    % r, c index => linear index
    % indices in im of the rotation voxel
    xyz0idx = sub2ind(size(im), xyz0idx(1), xyz0idx(2));
    
    % convert NaN to background voxels
    im(isnan(im)) = 0;
    
    % compute connected components on the intersecting plane
    ccim = bwconncomp(im);
    
    % find the connected component that contains the bifurcation voxel
    FOUND = false;
    for CCI = 1:ccim.NumObjects
        if (any(ccim.PixelIdxList{CCI} == xyz0idx))
            FOUND = true;
            break
        end
    end
    
    if (~FOUND)
        error('Assertion fail: Root voxel not found')
    end
    
    % extract coordinates of voxels that belong to the connected component
    % from above
    xyz = [gx(ccim.PixelIdxList{CCI}), ...
        gy(ccim.PixelIdxList{CCI}), ...
        gz(ccim.PixelIdxList{CCI})];
    
    % compute PCA on the intersected voxel coordinates
    [~, eigd] = pts_pca(xyz');
    
    % we are just interested in the largest eigenvalue, that will allow us
    % to estimate the maximum thickness of the vessel at this point
    maxeigd(I) = eigd(1);
    
end

% list of branches that we want to prune
idxprune = (cc.BranchLenght ./ sqrt(4*maxeigd) < lratio) & cc.IsLeaf;

% loop to prune branches
for I = find(idxprune)
    
    % get indices from current branch and its neighbours
    idx = cat(1, cc.PixelIdxList{[I cc.BranchNeighbours{I}]});
    
    % linear index => r, c, s indices
    [r, c, s] = ind2sub(size(nrrdsk.data), idx);
    
    % coordinates of a box 
    boxmin = [min(r), min(c), min(s)];
    boxmax = [max(r), max(c), max(s)];
    
    % crop box from whole skeleton segmentation
    box = nrrdsk.data(boxmin(1):boxmax(1), boxmin(2):boxmax(2), ...
        boxmin(3):boxmax(3));
    
    % compute connected components. The current branch and its neighbours
    % should form only 1, but because voxels from other branches can be
    % contained within the box, the number of connected components in the
    % box can be higher
    boxcc = bwconncomp(box);
    
    % save the number of connected components in the box for later
    boxnumcc = boxcc.NumObjects;
    
    % remove the current branch's voxels from the box
    [r, c, s] = ind2sub(size(nrrdsk.data), cc.PixelIdxList{I});
    r = r - boxmin(1) + 1;
    c = c - boxmin(2) + 1;
    s = s - boxmin(3) + 1;
    idx = sub2ind(size(box), r, c, s);
    box(idx) = 0;
    
    % recompute connected components
    boxcc = bwconncomp(box);
    
    % if the number of connected components hasn't changed, then we remove
    % the whole current branch from the skeleton segmentation
    if (boxcc.NumObjects == boxnumcc)
        nrrdsk.data(cc.PixelIdxList{I}) = 0;
    else
        % otherwise, we remove the branch except for the bifurcation voxel
        nrrdsk.data(cc.PixelIdxList{I}(cc. Degree{I} < 3)) = 0;
    end
    
end
