function [sk, cc2, dsk, dictsk, idictsk] = skeleton_label(sk, im, res, alphamax, p)
% SKELETON_LABEL  Give each branch of a skeleton a different label, and
% sort the voxels within each branch
%
% [LAB, CC] = SKELETON_LABEL(SK)
%
%   SK is a 3D segmentation mask. SK is assumed to be a skeleton resulting
%   from some kind of thinning algorithm, e.g. the C++ program
%   skeletonize3DSegmentation provided by Gerardus. That is, the
%   segmentation looks like a series of 1 voxel-thick branches connected
%   between them (cycles are allowed).
%
%   LAB is an image of the same dimensions of SK where the value of each
%   voxel is the label of the branch it belongs to.
%
%   Voxels in bifurcations are labelled as the nearest branch. Voxels
%   completely surrounded by other bifurcation voxels are assigned the
%   nearest label arbitrarily. Some voxels may remain unlabelled, e.g. if
%   they are not connected to any others in the segmentation.
%
%   CC is a struct like those provided by Matlab's function bwconncomp().
%   Each vector in CC.PixelIdxList{i} has the indices of all the voxels in
%   branch i. These indices are image indices, so the row, column, slice
%   coordinates of e.g. branch 5 can be obtained doing
%
%     [r, c, s] = ind2sub(size(sk), cc.PixelIdxList{5});
%
%   We have added new fields to CC:
%
%     CC.IsolatedBifurcationPixelIdx{i}: A list of bifurcation voxels in
%                       the image that are sourrounded only by other
%                       bifurcation voxels. Isolated bifurcation voxels
%                       don't belong to any branch, but they belong to the
%                       skeleton
%
%     CC.PixelParam{i}: parameterization values for the voxels in
%                       cc.PixelIdxList{i}. The parameterization is
%                       computed as the accumulated chord distance between
%                       consecutive branch voxels. This parameterization
%                       can be used, e.g. for spline interpolation.
%                       Branches that contain a loop cannot be
%                       parameterized as a 1-D spline, and thus
%                       CC.PixelParam(i)=NaN
%
%     CC.IsLeaf(i):     flags indicating whether each section is a "leaf",
%                       i.e. whether it's a branch with a free extreme
%
%     CC.IsLoop(i):     flags indicating whether each section is a "loop",
%                       i.e. whether it's a branch around a hole
%
%     CC.BranchLength(i): chord-length of each skeleton branch. This is a
%                       parameterization of the branch's skeleton as a
%                       1-D spline. Branches that contain a loop cannot
%                       be parameterized
%
%     CC.Degree{i}:     degree of each branch voxel, i.e. how many voxels
%                       it is connected to
%
%     CC.BranchNeighboursLeft{i}: index of the branches that are neighbours
%                       to branch i on its "left"
%
%     CC.BranchNeighboursRight{i}: index of the branches that are
%                       neighbours to branch i on its "right"
%
%     CC.MergedBranches{i}: List of branches in the pre-merged skeleton
%                      that were merged to create branch i
%
%   There may be a small discrepancy between cc.PixelIdxList and the voxels
%   labelled in LAB. The reason is that sometimes the same bifurcation
%   voxel can be shared by two or more branches.
%
%
% [LAB, CC, D, DICT, IDICT] = SKELETON_LABEL(SK)
%
%   D is the sparse distance matrix between skeleton points.
%
%   DICT, IDICT are vectors to convert between image and distance matrix
%   indices. i=DICT(I), I=IDICT(i), where I is an image linear index, and i
%   is a distance matrix index, i.e. you use them e.g. SK(I) or D(i, :).
%
% ... = SKELETON_LABEL(SK, IM, RES)
%
%   IM is the original 3D the skeleton SK was computed from. In this case,
%   LAB has the labelling of IM. If you want to extract the labelling for
%   SK, just run
%
%     >> LAB .* SK
%
%   RES is a 3-vector with the voxel size as [row, column, slice]. By
%   default, RES=[1 1 1].
%
% ... = SKELETON_LABEL(SK, IM, RES, ALPHAMAX, P)
%
%   ALPHAMAX is a scalar. ALPHAMAX >= 0 means that each branch end will be
%   merged with one of its neighbours if the angle between both is 
%   <= ALPHAMAX (radians). This is a way to get the same label for all the
%   segments that result from a long branch having lateral secondary
%   branches. By default ALPHAMAX = -1, so no merging is performed. As a
%   rule of thumb, ALPHAMAX = 10/180*pi seems to work well.
%
%   P is a scalar in [0, 1]. To compute the angle between branches, an
%   approximating or smoothing cubic spline is fit to the skeleton voxels
%   using csaps(..., P). P=0 is the smoothest spline (a line with the least
%   squares approximation), while P=1 is a rugged spline (the spline
%   interpolated the voxels). Adequate values of P depend on the image
%   resolution, so it's difficult to propose a formula. As a rule of thumb,
%   it seems that if resolution is in the order 1e-n, then a good value for
%   P=1-1e-n. For example, if resolution is in the order 1e-5, 
%   P=1-1e-5=0.999999. By default, P=1 and no smoothing is performed.
%
%   If you want to see how branches are being merged and smooth, uncomment
%   the DEBUG block at the end of internal function angle_btw_branches()
%   below.
%
%
% See also: skeleton_plot, scinrrd_skeleton_prune.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.11.10
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
error(nargchk(1, 5, nargin, 'struct'));
error(nargoutchk(0, 5, nargout, 'struct'));

% defaults
if (nargin < 2)
    im = [];
end
if (nargin < 3 || isempty(res))
    res = [1 1 1];
end
if (nargin < 4 || isempty(alphamax))
    alphamax = -Inf; % don't merge by default
end
if (nargin < 5 || isempty(p))
    p = 1.0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% split skeleton into branches and bifurcation voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get sparse matrix of distances between voxels. To label the skeleton we
% don't care about the actual distances, just need to know which voxels are
% connected to others. Actual distances are needed to parameterize the
% branches, though
[dsk, dictsk, idictsk] = seg2dmat(sk, 'seg', res);

% compute degree of each skeleton voxel
deg = sum(dsk > 0, 2);

% get distance matrix index of the bifurcation voxels
bifidx = deg >= 3;

% convert distance matrix index to image index
bifidx = idictsk(bifidx);

% remove bifurcation voxels from image
sk(bifidx) = 0;

% get connected components in the image
cc = bwconncomp(sk);

% make size vector always have size(3)
if length(cc.ImageSize) == 2
    cc.ImageSize = [cc.ImageSize 1];
end

% init output
cc.IsolatedBifurcationPixelIdx = [];

% loop each bifurcation voxel (working with matrix indices, not image
% indices)
for v = find(deg >= 3)'
    
    % get all its neighbours
    vn = find(dsk(v, :));
    
    % keep only neighbours that are not bifurcation voxels
    vn = vn(full(deg(vn)) < 3);
    
    % if this is a bifurcation voxel that is completely surroundered by
    % other bifurcation voxels, add it to the corresponding list
    if isempty(vn)
        
        cc.IsolatedBifurcationPixelIdx = ...
            [cc.IsolatedBifurcationPixelIdx ; idictsk(v)];
        
    else % otherwise, add it to all branches that it touches
        
        % get labels of branches that are neighbours to current bifurcation
        % voxel
        idx = find(...
            arrayfun(@(I) any(ismember(idictsk(vn), ...
            cc.PixelIdxList{I})), 1:cc.NumObjects));
        
        % add current bifurcation voxel to neighbour branches
        for I = 1:length(idx)
            cc.PixelIdxList{idx(I)} = ...
                union(cc.PixelIdxList{idx(I)}, idictsk(v));
        end
    
    end

end

% make sure that lists of voxels are column vectors
cc.PixelIdxList = cellfun(@(x) x(:), cc.PixelIdxList, ...
    'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sort the skeleton voxels in each branch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init outputs
cc.PixelParam = cell(1, cc.NumObjects);
cc.IsLeaf = false(1, cc.NumObjects);
cc.IsLoop = false(1, cc.NumObjects);
cc.BranchLength = zeros(1, cc.NumObjects);
cc.Degree = cell(1, cc.NumObjects);
cc.BranchBifurcationPixelIdx = cell(1, cc.NumObjects);

% create the part of the NRRD struct necessary to convert to real world
% coordinates
if (isempty(res))
    nrrdaxis.spacing(1) = 1;
    nrrdaxis.spacing(2) = 1;
    nrrdaxis.spacing(3) = 1;
else
    nrrdaxis.spacing(1) = res(1);
    nrrdaxis.spacing(2) = res(2);
    nrrdaxis.spacing(3) = res(3);
end
nrrdaxis.min = deal(nrrdaxis.spacing / 2);
nrrdaxis.size(1) = cc.ImageSize(1);
nrrdaxis.size(2) = cc.ImageSize(2);
nrrdaxis.size(3) = cc.ImageSize(3);

% loop each branch in the skeleton
for I = 1:cc.NumObjects

    % list of skeleton voxels in current branch
    br = cc.PixelIdxList{I};
    
    % number of voxels in the current branch
    N = length(br);
    
    % degenerate case in which the branch has only one voxel
    if (N == 1)
        cc.PixelParam{I} = 0;
        cc.IsLeaf(I) = false;
        cc.IsLoop(I) = false;
        cc.BranchLength(I) = 0;
        cc.Degree{I} = 0;
        cc.BranchBifurcationPixelIdx{I} = [];
        continue
    end

%     % DEBUG: coordinates of branch voxels
%     [r, c, s] = ind2sub(size(sk), br);
%     xyz = scinrrd_index2world([r c s], nrrdaxis);

    % create a "branch distance matrix" for only the voxels in the branch
    bri = dictsk(br);
    dbr = dsk(bri, bri);
    
    % compute Dijkstra's shortest path from an arbitrary voxel in the
    % branch to every other voxel
    [d, ~] = dijkstra(dbr, 1);
    
    % the furthest voxel should be one of the ends of the branch (unless
    % there's a loop)
    [~, v0] = max(d);
    
    % compute shortest path to all other voxels in the branch
    %
    % note that even for apparently "wire" branches, sometimes we get small
    % cycles of 1 voxel, and instead of a "linear" shortest-path, we have a
    % tree, so in the next step, some of the voxels are going to be thrown
    % away
    [d, parents] = dijkstra(dbr, v0);
    [~, v1] = max(d);
    
    % backtrack the whole branch in order from the furthest point to
    % the original extreme point
    cc.PixelIdxList{I} = [];
    cc.PixelParam{I} = [];
    J = 1;
    while (v1 ~= 0)
        cc.PixelIdxList{I}(J) = idictsk(bri(v1));
        cc.PixelParam{I}(J) = d(v1);
        v1 = parents(v1);
        J = J + 1;
    end
    
    % the branch's skeleton can have two shapes: a "wire" or a "loop". In
    % principle, we could try to break up loops, but there's not a real
    % good reason to do so, because does it really make sense to straighten
    % a loop as if it were a "wire"? A loop can have the following topology
    %
    %                         --
    %                         ||
    %                         ||
    %                       ------
    %
    % In this case, the loop cannot even be broken in a meaningful way.

    % if almost all the voxels have been put into cc.PixelIdxList{I}, it
    % means that we have a "wire" branch
    cc.IsLoop(I) = length(cc.PixelIdxList{I}) + 3 < N;
    if (cc.IsLoop(I))
        
        % we are not going to sort loops, so recover the unordered voxels
        cc.PixelIdxList{I} = idictsk(bri);
        
        % and clear up the parameterisation
        cc.PixelParam{I} = [];
        cc.BranchLength(I) = nan;
        
    else % the branch is a "wire", we keep the parameterisation and ordering
        
        % make sure that cc.PixelIdxList{I} is a column vector
        cc.PixelIdxList{I} = cc.PixelIdxList{I}(:);
        
        % reorder voxels so that the parameterization increases
        % monotonically
        cc.PixelIdxList{I} = cc.PixelIdxList{I}(end:-1:1);
        cc.PixelParam{I} = cc.PixelParam{I}(end:-1:1);
        
        % extract length of branch
        cc.BranchLength(I) = cc.PixelParam{I}(end);
        
    end
    
    % force to be a column vector
    cc.PixelParam{I} = cc.PixelParam{I}(:);
    
    % degree of each total branch voxel
    cc.Degree{I} = full(deg(dictsk(cc.PixelIdxList{I})));
    
    % each branch can have 0, 1, or 2 bifurcation voxels:
    %   0: if the branch is floating in the air without any neighbours
    %   1: if the branch is a left, so it's connected only on one end
    %   2: if the branch is connected on both ends
    vend = dictsk(cc.PixelIdxList{I}([1 end]));
    cc.BranchBifurcationPixelIdx{I} = idictsk(full(vend(deg(vend) > 2)));
    
    % a branch is a leaf is it has at least one voxel with degree==1 (tip of a
    % branch) or 0 (single voxel floating in space)
    cc.IsLeaf(I) = any(full(vend(deg(vend) < 2)));
    
end

% % Debug (2D image) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for I = 1:cc.NumObjects
%     % list of voxels in current branch
%     br = cc.PixelIdxList{I};
%     
%     % convert image linear indices to image r, c, s indices
%     [r, c, s] = ind2sub(size(sk), br);
%     
%     % plot branch in order
%     plot(c, r, 'w')
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get neighbours of each branch and angles between them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, we need to get the neighbours of the isolated bifurcation voxels.
% The reason is that two branches can be neighbours, but not touch because
% they have an isolated voxel in between. Knowning the neighbours of
% isolated voxels will be useful when we are looking for branch neighbours
cc.IsolatedBifurcationPixelNeighbours = cell(1, ...
    length(cc.IsolatedBifurcationPixelIdx));

for I = 1:length(cc.IsolatedBifurcationPixelIdx)
    
    % get current isolated bifurcation voxel
    v = cc.IsolatedBifurcationPixelIdx(I);
    
    % get neighbour voxels
    vn = idictsk(dsk(dictsk(v), :) > 0);
    
    % branches that contain any voxel in the list of neighbours
    idx = find(cellfun(@(x) any(ismember(vn, x)), ...
        cc.PixelIdxList));
    
    % add list of neighbour branches to the corresponding isolated voxel
    cc.IsolatedBifurcationPixelNeighbours{I} = idx(:)';
    
end

% init left and right neighbours
cc.BranchNeighboursLeft = cell(1, cc.NumObjects);
cc.BranchNeighboursRight = cell(1, cc.NumObjects);

% search for branch left and right neighbours
for I = 1:cc.NumObjects
    
    % get voxels at both ends of the branch. The first one is the "left"
    % end of the branch, and the last one is the "right" end
    vl = cc.PixelIdxList{I}(1);
    vr = cc.PixelIdxList{I}(end);
    
    % get neighbour voxels on each end
    vln = idictsk(dsk(dictsk(vl), :) > 0);
    vrn = idictsk(dsk(dictsk(vr), :) > 0);
    
    % branches that contain any voxel in the list of neighbours
    cc.BranchNeighboursLeft{I} = find(cellfun(@(x) ...
        any(ismember(vln, x)), cc.PixelIdxList));
    cc.BranchNeighboursRight{I} = find(cellfun(@(x) ...
        any(ismember(vrn, x)), cc.PixelIdxList));
    
    % if any of the neighbour voxels are isolated bifurcated voxels, then
    % they are not going to be in any branch, but we know which branches
    % they are touching. We add those branches as neighbours to the current
    % branch
    cc.BranchNeighboursLeft{I} = union(cc.BranchNeighboursLeft{I}, ...
        cell2mat(cc.IsolatedBifurcationPixelNeighbours(...
        ismember(cc.IsolatedBifurcationPixelIdx, vln))));
    cc.BranchNeighboursRight{I} = union(cc.BranchNeighboursRight{I}, ...
        cell2mat(cc.IsolatedBifurcationPixelNeighbours(...
        ismember(cc.IsolatedBifurcationPixelIdx, vrn))));
    
    % remove current branch, so that it isn't its own neighbour
    cc.BranchNeighboursLeft{I} = setdiff(cc.BranchNeighboursLeft{I}, I);
    cc.BranchNeighboursRight{I} = setdiff(cc.BranchNeighboursRight{I}, I);
    
end

% init outputs
cc.BranchMergeCandidateLeft = zeros(1, cc.NumObjects);
cc.BranchMergeCandidateLeftAngle = inf(1, cc.NumObjects);
cc.BranchMergeCandidateRight = zeros(1, cc.NumObjects);
cc.BranchMergeCandidateRightAngle = inf(1, cc.NumObjects);

% loop again to compute the angles between neighbour branches. This cannot
% be part of the previous loop, because beforehand we need to know the
% correspondences between all branches
for I = 1:cc.NumObjects
    
    % loops are not going to be merged
    if (cc.IsLoop(I))
        continue
    end

    % neighbours to the left of this branch
    for J = cc.BranchNeighboursLeft{I}(:)'
        
        % concatenate indices of neighbour branches
        if (cc.IsLoop(J))
            
            alpha = Inf;
            
        elseif (any(cc.BranchNeighboursLeft{J} == I)) % connected to left of neighbour
            
            % compute angle between branches
            alpha = angle_btw_branches(cc.PixelIdxList{I}(end:-1:1), ...
                cc.PixelIdxList{J}, cc.ImageSize, nrrdaxis, p);
            
        elseif (any(cc.BranchNeighboursRight{J} == I)) % connected to right of neighbour
            
            alpha = angle_btw_branches(cc.PixelIdxList{J}, ...
                cc.PixelIdxList{I}, cc.ImageSize, nrrdaxis, p);
            
        else % assertion fail, a neighbour must be connected on the left or right
            error(['Assertion fail: Branch ' num2str(J) ...
                ' not connected to neighbour ' num2str(I)])
        end
        
        % update the candidate to be merged
        if (alpha < cc.BranchMergeCandidateLeftAngle(I) ...
                && alpha < alphamax)
            cc.BranchMergeCandidateLeftAngle(I) = alpha;
            cc.BranchMergeCandidateLeft(I) = J;
        end
            
    end
    
    % neighbours to the right of this branch
    for J = cc.BranchNeighboursRight{I}(:)'
        
        % concatenate indices of neighbour branches
        if (cc.IsLoop(J))
            
            alpha = Inf;
            
        elseif (any(cc.BranchNeighboursLeft{J} == I)) % connected to left of neighbour
            
            alpha = angle_btw_branches(cc.PixelIdxList{I}, ...
                cc.PixelIdxList{J}, cc.ImageSize, nrrdaxis, p);
            
        elseif (any(cc.BranchNeighboursRight{J} == I)) % connected to right of neighbour
            
            alpha = angle_btw_branches(cc.PixelIdxList{I}, ...
                cc.PixelIdxList{J}(end:-1:1), cc.ImageSize, nrrdaxis, p);
            
        else % assertion fail, a neighbour must be connected on the left or right
            error(['Assertion fail: Branch ' num2str(J) ...
                ' not connected to neighbour ' num2str(I)])
        end
        
        % update the candidate to be merged
        if (alpha < cc.BranchMergeCandidateRightAngle(I) ...
                && alpha < alphamax)
            cc.BranchMergeCandidateRightAngle(I) = alpha;
            cc.BranchMergeCandidateRight(I) = J;
        end
            
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute merging buckets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (alphamax ~= -Inf)

    % branches that want to merge with somebody on the left or on the right
    idxabl = find(cc.BranchMergeCandidateLeft);
    idxabr = find(cc.BranchMergeCandidateRight);
    idxab = [idxabl idxabr];
    
    % candidates to be merged
    idxba = [cc.BranchMergeCandidateLeft(idxabl) ...
        cc.BranchMergeCandidateRight(idxabr)];
    
    % check that if "a" wants to merge with "b", "b" also wants to merge
    % with "a"
    ok = cc.BranchMergeCandidateLeft(idxba) == idxab ...
        | cc.BranchMergeCandidateRight(idxba) == idxab;
    
    % init list of merging buckets
    cc2.MergedBranches = cell(0);
    
    % create merging buckets. Each bucket will contain an unsorted list of all
    % branches that need to be merged together
    for I = find(ok)
        
        % look if any of these branches are already going to be merged
        % to others
        idx = find(cellfun(@(x) any(ismember([idxab(I) idxba(I)], x)), ...
            cc2.MergedBranches));

        if (isempty(idx))
            
            % create new merging bucket
            cc2.MergedBranches = [cc2.MergedBranches {[idxab(I) idxba(I)]}];
            
        elseif (length(idx) == 1)
            
            % add this branch to the merging bucket it belongs to
            cc2.MergedBranches{idx} = ...
                union(cc2.MergedBranches{idx}, [idxab(I) idxba(I)]);
            
        elseif (length(idx) == 2)
            
            % add the candidates to the first bucket
            cc2.MergedBranches{idx(1)} = ...
                union(cc2.MergedBranches{idx(1)}, [idxab(I) idxba(I)]);
            
            % empty the other buckets into the first
            cc2.MergedBranches{idx(1)} = ...
                union(cc2.MergedBranches{idx(1)}, ...
                cc2.MergedBranches{idx(2)});
            
            % delete the second bucket
            cc2.MergedBranches(idx(2)) = [];
            
        else
            
            error('Assertion fail: Branches belong to more than 2 merging buckets')
            
        end
        
    end
    
    % add branches that are not going to be merged, one branck per bucket
    cc2.MergedBranches = [ ...
        num2cell(setdiff(1:cc.NumObjects, [cc2.MergedBranches{:}])) ...
        cc2.MergedBranches];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge branches together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (alphamax ~= -Inf)
    
    % add to cc2 struct the fields that are common to all branches
    cc2.Connectivity = cc.Connectivity;
    cc2.ImageSize = cc.ImageSize;
    cc2.NumObjects = length(cc2.MergedBranches);
    
    % init rest of the fields
    cc2.PixelIdxList = cell(1, cc2.NumObjects);
    cc2.PixelParam = cell(1, cc2.NumObjects);
    cc2.IsLeaf = false(1, cc2.NumObjects);
    cc2.IsLoop = false(1, cc2.NumObjects);
    cc2.BranchLength = zeros(1, cc2.NumObjects);
    cc2.Degree = cell(1, cc2.NumObjects);
    
    % loop every merging bucket
    for I = 1:cc2.NumObjects
        
        % get union of voxels from all the branches that are going to be
        % merged
        %
        % concatenate all the isolated bifurcation voxels. We could look up
        % which isolated voxels touch the branches in the bucket, but it's
        % going to be faster to just concatenate all, and let Dijkstra's
        % algorithm prune out the free riders
        br = [unique(cat(1, cc.PixelIdxList{cc2.MergedBranches{I}})) ; ...
            cc.IsolatedBifurcationPixelIdx(:)];
        
        % create a "branch distance matrix" for only the voxels in the
        % branch
        bri = dictsk(br);
        dbr = dsk(bri, bri);
        
        % compute Dijkstra's shortest path from an arbitrary voxel in the
        % branch to every other voxel
        [d, ~] = dijkstra(dbr, 1);
        
        % remove infinite values (not connected voxels), so that we can
        % compute the maximum distance
        d(isinf(d)) = nan;
        
        % the furthest voxel should be one of the ends of the branch
        % (the voxel's index is br(v0))
        [~, v0] = max(d);
        
        % compute shortest path to all other voxels in the branch
        %
        % note that even for apparently "wire" branches, sometimes we get
        % small cycles of 1 voxel, and instead of a "linear" shortest-path,
        % we have a tree, so in the next step, some of the voxels are going
        % to be thrown away
        [d, parents] = dijkstra(dbr, v0);
        d(isinf(d)) = nan;
        [~, v1] = max(d);
        
        % backtrack the whole branch in order from the furthest point to
        % the original extreme point
        J = 1;
        while (v1 ~= 0)
            cc2.PixelIdxList{I}(J) = br(v1);
            cc2.PixelParam{I}(J) = d(v1);
            v1 = parents(v1);
            J = J + 1;
        end
        
        % invert the order of the points, so that the parameterization
        % grows
        cc2.PixelIdxList{I} = cc2.PixelIdxList{I}(end:-1:1);
        cc2.PixelParam{I} = cc2.PixelParam{I}(end:-1:1);
        
        % make list of pixels column vector
        cc2.PixelIdxList{I} = cc2.PixelIdxList{I}(:);
        cc2.PixelParam{I} = cc2.PixelParam{I}(:);
        
        % if any merged branch is a leaf, then the total branch has to be a
        % leaf too
        cc2.IsLeaf(I) = any(cc.IsLeaf(cc2.MergedBranches{I}));
        
        % loops cannot merge, so a loop can only occur for single branches
        cc2.IsLoop(I) = any(cc.IsLoop(cc2.MergedBranches{I}));
        
        % degree of each total branch voxel
        cc2.Degree{I} = full(deg(dictsk(cc2.PixelIdxList{I})));
        
        % branch length of the total branch is obtained from the
        % parameterisation
        cc2.BranchLength(I) = cc2.PixelParam{I}(end);
        
    end
    
    % find the isolated voxels that are now part of a merged branch
    idx = ismember(cc.IsolatedBifurcationPixelIdx, ...
        cat(1, cc2.PixelIdxList{:}));
    
    % new isolated voxels are those that are not part of any branch
    cc2.IsolatedBifurcationPixelIdx = cc.IsolatedBifurcationPixelIdx(~idx);
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% this code copied from "get neighbours of each branch and angles 
    %% between them" above
    
    % first, we need to get the neighbours of the isolated bifurcation
    % voxels. The reason is that two branches can be neighbours, but not
    % touch because they have an isolated voxel in between. Knowning the
    % neighbours of isolated voxels will be useful when we are looking for
    % branch neighbours
    cc2.IsolatedBifurcationPixelNeighbours = cell(1, ...
        length(cc2.IsolatedBifurcationPixelIdx));
    
    for I = 1:length(cc2.IsolatedBifurcationPixelIdx)
        
        % get current isolated bifurcation voxel
        v = cc2.IsolatedBifurcationPixelIdx(I);
        
        % get neighbour voxels
        vn = idictsk(dsk(dictsk(v), :) > 0);
        
        % branches that contain any voxel in the list of neighbours
        idx = find(cellfun(@(x) any(ismember(vn, x)), ...
            cc2.PixelIdxList));
        
        % add list of neighbour branches to the corresponding isolated
        % voxel
        cc2.IsolatedBifurcationPixelNeighbours{I} = idx(:)';
        
    end
    
    % init left and right neighbours
    cc2.BranchNeighboursLeft = cell(1, cc2.NumObjects);
    cc2.BranchNeighboursRight = cell(1, cc2.NumObjects);
    
    % search for branch left and right neighbours
    for I = 1:cc2.NumObjects
        
        % get voxels at both ends of the branch. The first one is the
        % "left" end of the branch, and the last one is the "right" end
        vl = cc2.PixelIdxList{I}(1);
        vr = cc2.PixelIdxList{I}(end);
        
        % get neighbour voxels on each end
        vln = idictsk(dsk(dictsk(vl), :) > 0);
        vrn = idictsk(dsk(dictsk(vr), :) > 0);
        
        % branches that contain any voxel in the list of neighbours
        cc2.BranchNeighboursLeft{I} = find(cellfun(@(x) ...
            any(ismember(vln, x)), cc2.PixelIdxList));
        cc2.BranchNeighboursRight{I} = find(cellfun(@(x) ...
            any(ismember(vrn, x)), cc2.PixelIdxList));
        
        % if any of the neighbour voxels are isolated bifurcated voxels,
        % then they are not going to be in any branch, but we know which
        % branches they are touching. We add those branches as neighbours
        % to the current branch
        cc2.BranchNeighboursLeft{I} = union(...
            cc2.BranchNeighboursLeft{I}, ...
            cell2mat(cc2.IsolatedBifurcationPixelNeighbours(...
            ismember(cc2.IsolatedBifurcationPixelIdx, vln))));
        cc2.BranchNeighboursRight{I} = union(...
            cc2.BranchNeighboursRight{I}, ...
            cell2mat(cc2.IsolatedBifurcationPixelNeighbours(...
            ismember(cc2.IsolatedBifurcationPixelIdx, vrn))));
        
        % remove current branch, so that it isn't its own neighbour
        cc2.BranchNeighboursLeft{I} = ...
            setdiff(cc2.BranchNeighboursLeft{I}, I);
        cc2.BranchNeighboursRight{I} = ...
            setdiff(cc2.BranchNeighboursRight{I}, I);
        
    end
    
    %% end code copied from above
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % sort the struct fields and overwrite description of the skeleton
    cc2 = orderfields(cc2, [2:5 11 6:10 13 14 1 12]);
    cc2 = rmfield(cc2, 'IsolatedBifurcationPixelNeighbours');
    
else % in case we are not merging
    
    fn = fieldnames(cc);
    cc2 = rmfield(cc, fn([11 12 15:end]));
    cc2.MergedBranches = num2cell(1:cc2.NumObjects);
    
end

clear cc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% label skeleton voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set amount of memory allowed for labels (we need label "0" for the
% background, and labels "1", "2", ..., "cc.NumObjects" for each component
%
% if a segmentation is provided to be labelled too, we need an extra "TODO"
% label
if isempty(im)
    req_bits = ceil(log2(cc2.NumObjects + 1));
else
    req_bits = ceil(log2(cc2.NumObjects + 2));
end
if (req_bits == 1)
    lab_class = 'boolean';
elseif (req_bits < 8)
    lab_class = 'uint8';
elseif (req_bits < 16)
    lab_class = 'uint16';
elseif (req_bits < 32)
    lab_class = 'uint32';
elseif (req_bits < 64)
    lab_class = 'uint64';
else
    lab_class = 'uint64';
    warning('Too many labels. There are more than 2^64 labels. Using uint64, some labels will be lost');
end
if ~strcmp(lab_class, class(sk)) % we need a different class than the input image so that we can represent all labels
    if strcmp(lab_class, 'boolean')
        sk = false(size(sk));
    else
        sk = zeros(size(sk), lab_class);
    end
else
    sk = sk * 0;
end

% give each voxel in the skeleton its label
for lab = 1:cc2.NumObjects
    sk(cc2.PixelIdxList{lab}) = lab;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% label whole segmentation voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isempty(im))
    
    % convert im data type to be the same as the seeds, so that we are
    % guaranteed that all label values can be represented in the output
    % image
    im = cast(im, class(sk));
    
    % label for voxels that need to be labelled
    TODO = cast(max(sk(:))+1, class(sk));
    
    % mark segmentation voxels with the highest label, so that we know that
    % they still need to be labelled, but belong to the segmentation as
    % opposed to the background
    im(im > 0) = TODO;
    
    % init algorithm by copying the labelled skeleton to the segmentation
    idx = sk > 0;
    im(idx) = sk(idx);
  
    % region grow algorithm to label all segmented voxels
    % overwrite skeleton image
    sk = bwregiongrow(im, TODO, res);

    % in some very particular cases, a small patch of voxels may be left
    % unlabelled. We are just going to remove them from the segmentation
    sk(sk == TODO) = 0;
    
end

% % Debug (2D image) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% hold off
% imagesc(sk)
% axis xy
% 
% % find tree leaves: convert matrix indices to image linear indices
% idx = idictsk(deg == 1);
% 
% % convert image linear indices to row, col, slice
% [r, c, s] = ind2sub(size(sk), idx);
% 
% % plot leaves
% hold on
% plot(c, r, 'ro')
% 
% % find tree bifurcations: convert matrix indices to image linear indices
% idx = idictsk(deg >= 3);
% 
% % convert image linear indices to row, col, slice
% [r, c, s] = ind2sub(size(sk), idx);
% 
% % plot leaves
% plot(c, r, 'go')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% angle_btw_branches(): compute angle between two concatenated branches
function alpha = angle_btw_branches(idx0, idx1, sz, nrrdaxis, p)

% make sure that indices are vertical vectors
idx0 = idx0(:);
idx1 = idx1(:);

% concatenate branches, avoiding to duplicate the joint point
if (idx0(end) == idx1(1))

    idx = [idx0(1:end-1) ; idx1];
    
else
    idx = [idx0 ; idx1];
    
end

% index of joint voxel
bifidx = size(idx0, 1);

% real world coordinates of voxels
[r, c, s] = ind2sub(sz, idx);
xyz = scinrrd_index2world([r, c, s], nrrdaxis);
        
% compute spline parameterization (Lee's centripetal scheme)
t = cumsum([0; (sum((xyz(2:end, :) - xyz(1:end-1, :)).^2, 2)).^.25]);

% compute cubic smoothing spline
pp = csaps(t', xyz', p);

% get curve tangent to the left and right of the joint point
% (note that dx0 doesn't necessarily correspond to br0, as opposed
% to br1)
dx0 = ppval(pp, pp.breaks(bifidx-1:bifidx));
dx1 = ppval(pp, pp.breaks(bifidx:bifidx+1));

dx0 = dx0(:, 2) - dx0(:, 1);
dx1 = dx1(:, 2) - dx1(:, 1);

% compute angle between tangent vectors
alpha = acos(dot(dx0, dx1) / norm(dx0) / norm(dx1));

% % DEBUG
% hold off
% plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'b')
% hold on
% plot3(xyz(bifidx, 1), xyz(bifidx, 2), xyz(bifidx, 3), 'b*')
% fnplt(pp, 'r')
% axis xy equal
% view(2)
% pause

end
