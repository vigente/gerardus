function nrrdsk = scinrrd_skeleton_prune(nrrdsk, nrrd, minlen, lratio, alphamax, p)
% SCINRRD_SKELETON_PRUNE  Prune branches in a segmentation's skeletonization
%
% This function prunes the leaves of a segmentation skeleton in three
% steps:
%
%   1. Erode clumps of voxels to look like branches
%
%   2. Prune very short leaves (iteratively until no more leaves can be
%      pruned)
%
%   3. Prune leaves created by artifacts in the segmentation. Spurious
%      leaves are those roughly as long as the local radius of the main
%      branch they are attached to (iteratively until no more leaves can
%      be pruned)
%
%
% NRRDPR = SCINRRD_SKELETON_PRUNE(NRRDSK)
% NRRDPR = SCINRRD_SKELETON_PRUNE(NRRDSK, [], MINLEN)
%
%   This syntax runs steps 1 and 2 only.
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
% NRRDPR = SCINRRD_SKELETON_PRUNE(NRRDSK, NRRD, MAXCLUMP, MINLEN, LRATIO, ALPHAMAX, P)
%
%   This syntax runs steps 1, 2 and 3.
%
%   NRRD is the binary segmentation mentioned above.
%
%   LRATIO is a scalar. Leaves with L/R < LRATIO will be pruned, where
%   L is the length from the bifurcation to the tip of the leaf. By
%   default, LRATIO=1.2.
%
%   ALPHAMAX and P are merging parameters. See the help of function
%   skeleton_label for details. If merging is not enabled, then no leaves
%   will be pruned in step 3.
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
% Version: 0.4.1
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
error(nargchk(1, 6, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

if ~isstruct(nrrdsk)
    error('NRRDSK must be an SCI NRRD struct')
end

% defaults
if (nargin < 2)
    nrrd = [];
end
if (nargin < 3 || isempty(minlen))
    minlen = 5;
end
if (nargin < 4 || isempty(lratio))
    lratio = 1.2;
end
if (nargin < 5 || isempty(alphamax))
    alphamax = -Inf; % don't merge by default
end
if (nargin < 6 || isempty(p))
    p = 1.0; % don't smooth by default
end

%% Step 1: removal of big clumps of voxels

% if we have no voxels to prune, don't waste time processing
if (nnz(nrrdsk.data) == 0)
    return
end

% get sparse matrix of distances between voxels. To label the skeleton we
% don't care about the actual distances, just need to know which voxels are
% connected to others. Actual distances are needed to parameterize the
% branches, though
[dsk, dictsk, idictsk] = seg2dmat(nrrdsk.data, 'seg', ...
    [nrrdsk.axis.spacing]);

% find bifurcation voxels

% compute degree of each skeleton voxel
deg = sum(dsk > 0, 2);

% get distance matrix index of the bifurcation voxels
bifidx = deg >= 3;

% matrix index => image index
bifidx = idictsk(bifidx);

% label connected components of skeleton branches

% remove bifurcation voxels from original skeleton
nrrdsk.data(bifidx) = 0;

% get connected components in the image
cc = bwconncomp(nrrdsk.data);

% make size vector always have size(3)
if length(cc.ImageSize) == 2
    cc.ImageSize = [cc.ImageSize 1];
end

% connectivity between branches and bifurcation clumps,
% and find which branches should be merged together

% tag each branch with its label
nrrdsk.data = labelmatrix(cc);

% create empty image volume and add only bifurcation voxels
sk2 = zeros(cc.ImageSize, 'uint8');
sk2(bifidx) = 1;

% compute connected components to obtain clumps of bifurcation voxels
bifcc = bwconncomp(sk2);

% loop all the bifurcation clumps, to find which branches are neighbours of
% each other
for I = 1:bifcc.NumObjects

    %% find which branches are connected to which bifurcation clumps

    % get voxels in the bifurcation clump
    % linear index -> row, col, slice
    [r, c, s] = ind2sub(cc.ImageSize, bifcc.PixelIdxList{I});

    % get a box 1 voxel bigger than the clump
    r0 = max(1, min(r) - 1);
    c0 = max(1, min(c) - 1);
    s0 = max(1, min(s) - 1);
    rend = min(cc.ImageSize(1), max(r) + 1);
    cend = min(cc.ImageSize(2), max(c) + 1);
    send = min(cc.ImageSize(3), max(s) + 1);

    % extract that box from the volume with the branches
    boxbr = nrrdsk.data(r0:rend, c0:cend, s0:send);
    
    % create a box for the bifurcation clump
    boxbif = zeros(size(boxbr), class(sk2));
    r = r - r0 + 1;
    c = c - c0 + 1;
    s = s - s0 + 1;
    boxbif(sub2ind(size(boxbif), r, c, s)) = 1;
    
    % create a box with the same size
    box  = zeros(size(boxbif));
    
    % tag as TODO=2 branch voxels
    box(boxbr ~= 0) = 2;
    
    % add to the vox the bifurcation clump voxels
    box(boxbif == 1) = 1;
    
    % dilate the clump 1 voxel
    box = bwregiongrow(box, 2, [], 1);
    
    % get the branches that the dilated bifurcation clump overlaps
    idx = double(boxbr(box == 1));

    % because in sk bifurcation voxels were removed, we can have "0" values
    % in idx. Remove them
    idx = idx(idx ~= 0);
    
    % number of branches
    N = length(idx);
    
    % in order to merge, the bifurcation needs to have at least 2
    % branches
    if (N < 2)
        continue
    end
    
    % get all combinations of pairs of branches that share this clump
    idx = nchoosek(idx, 2);
    
    % get voxel indices for bifurcation clump
    bif = bifcc.PixelIdxList{I};
    
    % duplicate the list. bif0 is going to be the list of superfluous
    % bifurcation voxels
    bif0 = bif;
    
    % loop each pair of branches combination
    for J = 1:size(idx, 1)
        
        % get voxel indices for each branch
        br0 = cc.PixelIdxList{idx(J, 1)};
        br1 = cc.PixelIdxList{idx(J, 2)};
        
        % merge and sort both branches and the bifurcation clump
        br = sort_branch([br0(:); bif(:); br1(:)], ...
            dsk, dictsk, idictsk);
        
        % remove from the list of superfluous bifurcation voxels those
        % that are needed to connect the two branches
        bif0 = setdiff(bif0, br);
        
    end
    
    % remove from the skeleton the list of superfluous bifurcation
    % voxels
    bifcc.PixelIdxList{I} = setdiff(bif, bif0);
        
end

% convert labelled skeleton to binary mask
nrrdsk.data = uint8(nrrdsk.data ~= 0);

% add bifurcation voxels to skeleton binary mask
nrrdsk.data(cat(1, bifcc.PixelIdxList{:})) = 1;

%% Step 2: pruning of very short leaf branches
while (1)
    
    % compute skeleton labelling
    [~, cc] = skeleton_label(nrrdsk.data, [], [nrrdsk.axis.spacing]);
   
    % get number of voxels in each branch
    n = cellfun(@(x) length(x), cc.PixelIdxList);
    
    % find leaf-branches that are shorter than the minimum length
    idx1 = find(n < minlen & cc.IsLeaf);
    
    % remove short branches from the segmentation
    nrrdsk.data(cat(1, cc.PixelIdxList{idx1})) = 0;
    
    % recompute the skeleton labelling
    [~, ~, bifcc, mcon] = skeleton_label(nrrdsk.data, [], [nrrdsk.axis.spacing]);
   
    % find bifurcation clusters that are connected to 0 or 1 branches
    idx2 = find(sum(mcon, 1) < 2);
    
    % remove those bifurcation clumps, because they are not connecting
    % branches, they are either floating alone in space, or terminating a
    % branch
    nrrdsk.data(cat(1, bifcc.PixelIdxList{idx2})) = 0;
    
    % if no bifurcation clumps were found to be removed, stop the
    % algorithm, because that means that no new short leaf-braches can be
    % found either
    if (isempty(idx2))
        break;
    end
    
end


%% Step 3: pruning of leaf branches produced by segmentation artifacts

% skip if no full segmentation is provided
if isempty(nrrd)
    return
end

% repeat the process until no branches are removed
atleastonepruning = true;
while (atleastonepruning)

    % if there are no prunings in this iterations, we stop
    atleastonepruning = false;
    
    % label the segmentation using multiple merging at every bifurcation
    % clump
    [nrrd.data, cc, ~, mcon, ~, cc2] = ...
        skeleton_label(nrrdsk.data, nrrd.data, [nrrd.axis.spacing], ...
        alphamax, p, false);
    
    % measure the stats of every merged branch
    stats2 = scinrrd_seg2label_stats(nrrd, cc2, p);
    
    % loop each branch in the list of merged branches
    for I = 1:cc2.NumObjects

        % we consider the current branch a main branch, and any branch
        % coming out of it, a secondary branch
        
        % compute the major radius of the main branch
        r = sqrt(4 * stats2.Var(2, I));
        
        % loop consecutive pairs of segments in the main branch
        for J = 1:length(cc2.MergedBranches{I})-1
            
            % get the bifurcation clump between both segments
            bif = cc2.MergedBifClumps{I}(J);
            
            % get secondary branches attached to the bifurcation clump
            % (usually, there's only 1 secondary branch, but there can be more)
            bn = setdiff(find(mcon(:, bif)), cc2.MergedBranches{I}(J:J+1));
            
            % keep only secondary branches that are leaf-branches
            bn = bn(cc.IsLeaf(bn));
            
            % skip to next segment if there are no valid leaf secondary
            % branches
            if isempty(bn)
                continue
            end
            
            % compute distance between the first and last voxels in the
            % secondary branch, which gives a better estimate of whether the
            % branch protudes from the main vessel or not, than the branches
            % length (because the branch can be bent)
            len = zeros(size(bn));
            for K = 1:length(len)
                idx = cc.PixelIdxList{bn(K)}([1 end]);
                [r1, c1, s1] = ind2sub(size(nrrdsk.data), idx);
                xyz = scinrrd_index2world([r1(:), c1(:), s1(:)], ...
                    nrrdsk.axis);
                len(K) = sqrt(sum(diff(xyz).^2));
            end
            
            % get branches that have to be pruned
            toremove = (len / r) <= lratio;
            atleastonepruning = atleastonepruning || any(toremove);
            
            % if length of the secondary branch is not much larger than the
            % radius of the main branch, we assume that the secondary branch is
            % a segmentation artifact, and remove it
            nrrdsk.data(cat(1, ...
                cc.PixelIdxList{bn(toremove)})) = 0;
            
        end
        
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort a set of voxels so that they form a branch
% idx: voxel indices
% d: distance matrix
function [idx, param] = sort_branch(idx, d, dict, idict)

% number of voxels in the current branch
N = length(idx);

% degenerate case in which the branch has only one voxel
if (N == 1)
    param = 0;
    return
end

% create a "branch distance matrix" for only the voxels in the branch
idx2 = full(dict(idx));
dbr = d(idx2, idx2);

% compute Dijkstra's shortest path from an arbitrary voxel in the
% branch to every other voxel
[dp, ~] = dijkstra(dbr, 1);

% the furthest voxel should be one of the ends of the branch
[~, v0] = max(dp);

% compute shortest path to all other voxels in the branch
%
% note that even for apparently "wire" branches, sometimes we get small
% cycles of 1 voxel, and instead of a "linear" shortest-path, we have a
% tree, so in the next step, some of the voxels are going to be thrown
% away
[dp, parents] = dijkstra(dbr, v0);
[~, v1] = max(dp);

% backtrack the whole branch in order from the furthest point to
% the original extreme point
idx = [];
param = [];
J = 1;
while (v1 ~= 0)
    idx(J) = idict(idx2(v1));
    param(J) = dp(v1);
    v1 = parents(v1);
    J = J + 1;
end


% make sure that cc.PixelIdxList{I} is a column vector
idx = idx(:);

% reorder voxels so that the parameterization increases
% monotonically
idx = idx(end:-1:1);
param = param(end:-1:1);

end
