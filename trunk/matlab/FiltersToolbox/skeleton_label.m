function [sk, cc, dsk, dictsk, idictsk] = skeleton_label(sk, im, res, alphamax, p)
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
%     CC.BifurcationPixelIdx: A list of voxels in the image that are
%                             bifurcation voxels
%
%     CC.MergedBranches{i}: only available if branch merging is selected by
%                           setting ALPHA>=0. Gives the list of branches in
%                           the pre-merged skeleton that were merged to
%                           create branch i
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
%     CC.Degree{i}:      degree of each skeleton voxel, i.e. how many voxel
%                        it is connected to
%
%     CC.BranchLength(i): chord-length of each skeleton branch. This is a
%                         parameterization of the branch's skeleton as a
%                         1-D spline. Branches that contain a loop cannot
%                         be parameterize, and thus CC.IsLoop(i)=NaN.
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
%   resolution, so it's difficult to propose a formula. For resolution in
%   the order of 2.5e-5, P=.999999 seems to give good results (note that
%   for small resolution, P=.999999 gives a very different result to
%   P=1.0). For resolution in the order of 1, P=0.8 seems to give good
%   results. By default, P=0.8.
%
%
% See also: skeleton_plot.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.11.5
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
    p = .8;
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

% loop each bifurcation voxel (working with matrix indices, not image
% indices)
for v = find(deg >= 3)'
    
    % get all its neighbours that are not bifurcation points
    vn = find(dsk(v, :));
    vn = vn(deg(vn) < 3);
    
    % if this is a bifurcation voxel that is completely surroundered by
    % other bifurcation voxels, we are not going to label it yet
    if isempty(vn)
        continue
    end
    
    % get labels of branches that are neighbours to current bifurcation
    % voxel
    idx = find(...
        arrayfun(@(I) any(ismember(idictsk(vn), cc.PixelIdxList{I})), ...
        1:cc.NumObjects));

    % add current bifurcation voxel to neighbour branches
    for I = 1:length(idx)
        cc.PixelIdxList{idx(I)} = ...
            union(cc.PixelIdxList{idx(I)}, idictsk(v));
    end
    
end

% make sure that lists of voxels are column vectors
cc.PixelIdxList = cellfun(@(x) x(:), cc.PixelIdxList, ...
    'UniformOutput', false);

% init output
cc.BifurcationPixelIdx = [];

% add bifurcation voxels that are surrounded only by other bifurcation
% voxels to the list of Bifurcation voxels, but don't add them to any
% branch
for v = find(deg >= 3)'
    
    % get all its neighbours
    vn = find(dsk(v, :), 1);

    % process only bifurcation voxels surrounded by bifurcation voxels
    if isempty(vn)
        continue
    end
    
    cc.BifurcationPixelIdx = union(cc.BifurcationPixelIdx, idictsk(v));
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sort the skeleton voxels in each branch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init output
cc.IsLeaf = false(1, cc.NumObjects);
cc.IsLoop = false(1, cc.NumObjects);
cc.BranchLength = zeros(1, cc.NumObjects);

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
        cc.IsLoop(I) = false;
        continue
    end

    % coordinates of branch voxels
    [r, c, s] = ind2sub(size(sk), br);
    xyz = scinrrd_index2world([r c s], nrrdaxis);

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
    
end

% degree of each total branch voxel
cc.Degree = cellfun(@(x) ...
    full(deg(dictsk(x))), cc.PixelIdxList, 'UniformOutput', false);

% we only consider as bifurcation voxels those at the end of each branch
% that have degree >= 3
cc.BifurcationPixelIdx = cellfun(@(x) ...
    x([1 ; end]), cc.PixelIdxList, 'UniformOutput', false);
cc.BifurcationPixelIdx = unique(cell2mat(cc.BifurcationPixelIdx));
cc.BifurcationPixelIdx = cc.BifurcationPixelIdx(...
    full(deg(dictsk(cc.BifurcationPixelIdx))) > 2);

% a branch is a leaf is it has at least one voxel with degree==1 (tip of a
% branch) or 0 (single voxel floating in space)
cc.IsLeaf = cellfun(@(x) any(x == 1 | x == 0), cc.Degree);


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

% initialize cells to contain branch neighbours
cc.BranchNeighbours = cell(1, cc.NumObjects);

% loop bifurcation voxels
for I = 1:length(cc.BifurcationPixelIdx)

    % index of bifurcation voxel (matrix index)
    idx = dictsk(cc.BifurcationPixelIdx(I));
    
    % get voxel's neighbours (image index)
    vn = idictsk(dsk(idx, :) > 0);
    
    % include voxel in the list of its own voxel neighbours
    vn = [vn; cc.BifurcationPixelIdx(I)];
    
    % branches that contain any voxel in the list of neighbours
    idx = find(cellfun(@(x) any(ismember(vn, x)), ...
        cc.PixelIdxList));
    
    % add all branches to the list of neighbours
    for J = idx
        cc.BranchNeighbours{J} = union(cc.BranchNeighbours{J}, ...
            idx);
    end
    
end

% a branch cannot be its own neighbour
cc.BranchNeighbours = arrayfun(@(I) setdiff(cc.BranchNeighbours{I}, I), ...
    1:numel(cc.BranchNeighbours), 'UniformOutput', false);

% init left and right neighbours
cc.BranchNeighboursLeft = cell(1, cc.NumObjects);
cc.BranchNeighboursRight = cell(1, cc.NumObjects);
cc.BranchNeighboursLeftAngle = cell(1, cc.NumObjects);
cc.BranchNeighboursRightAngle = cell(1, cc.NumObjects);

% split branch's neighbours into neighbours from the "left"
% (cc.PixelParam{i}(1)) or from the "right" (cc.PixelParam{i}(end))
for I = 1:cc.NumObjects
    
    % loop branch's neighbours that are not loops. To avoid duplicating
    % operations, we only look at neighbours with a larger label
    for J = cc.BranchNeighbours{I}(cc.BranchNeighbours{I} > I)
        
        % get list of indices in both skeleton branches
        idx0 = cc.PixelIdxList{I};
        idx1 = cc.PixelIdxList{J};

        % real world coordinates of voxels
        [r0, c0, s0] = ind2sub(cc.ImageSize, idx0);
        [r1, c1, s1] = ind2sub(cc.ImageSize, idx1);
        xyz0 = scinrrd_index2world([r0, c0, s0], nrrdaxis);
        xyz1 = scinrrd_index2world([r1, c1, s1], nrrdaxis);
        
        % compute distances between the first and last voxel in the
        % current and neihgbour branches
        d = dmatrix(xyz0([1 end], :)', xyz1([1 end], :)');
        
        % get matrix index with the minimum distance
        [dmin, idmin] = min(d(:));
        
        % depending on where the branches are joined, they are neighbours
        % on the left or on the right of each other
        if (idmin == 1) % 1st voxel branch <=> 1st voxel neighbour
            
            % index of bifurcation voxel
            bifidx = size(xyz0, 1);
                
            % concatenate voxels
            if (dmin == 0)
                xyz = [xyz0(end:-1:1, :); xyz1(2:end, :)];
            else
                xyz = [xyz0(end:-1:1, :); xyz1];
            end
            
        elseif (idmin == 2) % end voxel branch <=> 1st voxel neighbour
            
            % index of bifurcation voxel
            bifidx = size(xyz0, 1);
                
            % concatenate voxels
            if (dmin == 0)
                xyz = [xyz0; xyz1(2:end, :)];
            else
                xyz = [xyz0; xyz1];
            end
            
        elseif (idmin == 3) % 1st voxel branch <=> end voxel neighbour
            
            % index of bifurcation voxel
            bifidx = size(xyz1, 1);
                
            % concatenate voxels
            if (dmin == 0)
                xyz = [xyz1; xyz0(2:end, :)];
            else
                xyz = [xyz1; xyz0];
            end
            
        elseif (idmin == 4) % end voxel branch <=> end voxel neighbour
            
            % index of bifurcation voxel
            bifidx = size(xyz1, 1);
                
            % concatenate voxels
            if (dmin == 0)
                xyz = [xyz1; xyz0(end-1:-1:1, :)];
            else
                xyz = [xyz1; xyz0(end:-1:1, :)];
            end
            
        end

        % compute spline parameterization (Lee's centripetal scheme)
        t = cumsum([0; (sum((xyz(2:end, :) - xyz(1:end-1, :)).^2, 2)).^.25]);
        
        % compute cubic smoothing spline
        pp = csaps(t', xyz', p);
            
%         % DEBUG
%         hold off
%         fnplt(pp)
%         hold on
%         plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'o', 'LineWidth', 1)
%         foo = ppval(pp, pp.breaks(bifidx));
%         plot3(foo(1), foo(2), foo(3), 'ro', 'LineWidth', 1)
%         plot3(xyz(bifidx, 1), xyz(bifidx, 2), xyz(bifidx, 3), 'go', 'LineWidth', 1)
%         axis equal
%         view(2)
        
        % get curve tangent to the left and right of the bifurcation point
        % (note that dx0 doesn't necessarily correspond to br0, as opposed
        % to br1)
        dx0 = ppval(pp, pp.breaks(bifidx-1:bifidx));
        dx1 = ppval(pp, pp.breaks(bifidx:bifidx+1));
        
        dx0 = dx0(:, 2) - dx0(:, 1);
        dx1 = dx1(:, 2) - dx1(:, 1);
        
        % compute angle between derivatives. If the current branch or the
        % neighbour are loop, then we are not going to merge them
        if any(cc.IsLoop([I J]))
            alpha = Inf;
        else
            alpha = acos(dot(dx0, dx1) / norm(dx0) / norm(dx1));
        end

        if (idmin == 1) % 1st voxel branch <=> 1st voxel neighbour
            
            % add to left and right neighbours lists
            cc.BranchNeighboursLeft{I} = [cc.BranchNeighboursLeft{I} J];
            cc.BranchNeighboursLeft{J} = [cc.BranchNeighboursLeft{J} I];
            
            cc.BranchNeighboursLeftAngle{I} = ...
                [cc.BranchNeighboursLeftAngle{I} alpha];
            cc.BranchNeighboursLeftAngle{J} = ...
                [cc.BranchNeighboursLeftAngle{J} alpha];
            
        elseif (idmin == 2) % end voxel branch <=> 1st voxel neighbour
            
            % add to left and right neighbours lists
            cc.BranchNeighboursRight{I} = [cc.BranchNeighboursRight{I} J];
            cc.BranchNeighboursLeft{J} = [cc.BranchNeighboursLeft{J} I];
            
            cc.BranchNeighboursRightAngle{I} = ...
                [cc.BranchNeighboursRightAngle{I} alpha];
            cc.BranchNeighboursLeftAngle{J} = ...
                [cc.BranchNeighboursLeftAngle{J} alpha];
            
        elseif (idmin == 3) % 1st voxel branch <=> end voxel neighbour
            
            % add to left and right neighbours lists
            cc.BranchNeighboursLeft{I} = [cc.BranchNeighboursLeft{I} J];
            cc.BranchNeighboursRight{J} = [cc.BranchNeighboursRight{J} I];
            
            cc.BranchNeighboursLeftAngle{I} = ...
                [cc.BranchNeighboursLeftAngle{I} alpha];
            cc.BranchNeighboursRightAngle{J} = ...
                [cc.BranchNeighboursRightAngle{J} alpha];
            
        elseif (idmin == 4) % end voxel branch <=> end voxel neighbour
            
            % add to left and right neighbours lists
            cc.BranchNeighboursRight{I} = [cc.BranchNeighboursRight{I} J];
            cc.BranchNeighboursRight{J} = [cc.BranchNeighboursRight{J} I];
            
            cc.BranchNeighboursRightAngle{I} = ...
                [cc.BranchNeighboursRightAngle{I} alpha];
            cc.BranchNeighboursRightAngle{J} = ...
                [cc.BranchNeighboursRightAngle{J} alpha];
            
        end
        
    end
    
end

% remove neighbours field, that is now redundant
cc = rmfield(cc, 'BranchNeighbours');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge branches together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (alphamax ~= -Inf)
    
    % init list of pairs of branches to merge
    cc2.MergedBranches = cell(0);
    
    % ignore branches that are loops
    for I = find(~cc.IsLoop)
        
        % get the left neighbour that is best aligned to this branch
        [alphamin, idx] = min(cc.BranchNeighboursLeftAngle{I});
        
        % if there is no neighbour or the alignment is larger than the
        % angle threshold, skip this branch
        if (isempty(idx) || alphamin > alphamax)
            vn = [];
        else
            vn = cc.BranchNeighboursLeft{I}(idx);
        end

        % get the right neighbour that is best aligned to this branch
        [alphamin, idx] = min(cc.BranchNeighboursRightAngle{I});
        
        % if there is no neighbour or the alignment is larger than the
        % angle threshold, skip this branch
        if ~(isempty(idx) || alphamin > alphamax)
            vn = [vn cc.BranchNeighboursRight{I}(idx)];
        end

        % check each potencial merging
        for J = 1:length(vn)
            
            % we need to check whether the inverse relationship also
            % exists, i.e. if "I" wants to merge with "vn", that "vn" also
            % wants to merge with "I"
            [~, idxvn] = min(cc.BranchNeighboursLeftAngle{vn(J)});
            if (I == cc.BranchNeighboursLeft{vn(J)}(idxvn))
                % I -> vn and vn -> I, so we are going to merge this
                % branches (nothing needs to be done in this "if" block)
            else
                [~, idxvn] = min(cc.BranchNeighboursRightAngle{vn(J)});
                if (I == cc.BranchNeighboursRight{vn(J)}(idxvn))
                    % I -> vn and vn -> I, so we are going to merge this
                    % branches (nothing needs to be done in this "if"
                    % block)
                else
                    % I -> vn, but not viceversa, so we are skipping this
                    % branch and not merging
                    continue
                end
            end
            
            % the two labels that have to be merged
            v = sort([I vn(J)]);
            
            % look if any of these branches are already going to be merged
            % to others
            idx = find(cellfun(@(x) any(ismember(v, x)), ...
                cc2.MergedBranches), 1);
            
            if (isempty(idx))
                % create new merging bucket for this branch
                cc2.MergedBranches = [cc2.MergedBranches {v}];
            else
                % add this branch to the merging bucket it belongs to
                cc2.MergedBranches{idx} = ...
                    union(cc2.MergedBranches{idx}, v);
            end
            
        end
        
    end
    
    % add branches that are not going to be merged
    cc2.MergedBranches = [...
        num2cell(setdiff(1:cc.NumObjects, [cc2.MergedBranches{:}])) ...
        cc2.MergedBranches];
    
    % add to cc2 struct the fields that are common to all branches
    cc2.Connectivity = cc.Connectivity;
    cc2.ImageSize = cc.ImageSize;
    cc2.NumObjects = length(cc2.MergedBranches);
    
    % loop every merging bucket
    for I = 1:cc2.NumObjects
        
        % get union of voxels from all the branches that are going to be merged
        br = unique(cat(1, cc.PixelIdxList{cc2.MergedBranches{I}}));
        
        % get number of voxels
        N = length(br);
        
        % concatenate all the bifurcation voxels that are only surrounded by
        % other bifurcation voxels. The reason is that these voxels connect
        % neighbour branches that are not touching. It's too complicated to
        % keep track of which voxel touches which branch, so it's going to be
        % simpler to just concatenate all of them, and let Dijkstra's algorithm
        % throw away those that are not connected
        br = [br ; cc.BifurcationPixelIdx'];
        
        % create a "branch distance matrix" for only the voxels in the branch
        bri = dictsk(br);
        dbr = dsk(bri, bri);
        
        % compute Dijkstra's shortest path from an arbitrary voxel in the
        % branch to every other voxel
        [d, ~] = dijkstra(dbr, 1);
        
        % remove voxels that are not connected to the branch
        idx = ~isinf(d);
        %     br = br(idx);
        d = d(idx);
        dbr = dbr(idx, idx);
        bri = bri(idx);
        
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
        cc2.PixelIdxList{I} = [];
        cc2.PixelParam{I} = [];
        J = 1;
        while (v1 ~= 0)
            cc2.PixelIdxList{I}(J) = idictsk(bri(v1));
            cc2.PixelParam{I}(J) = d(v1);
            v1 = parents(v1);
            J = J + 1;
        end
        
        % real world coordinates
        [r, c, s] = ind2sub(cc2.ImageSize, cc2.PixelIdxList{I}');
        xyz = scinrrd_index2world([r, c, s], nrrdaxis);
        
        % recompute parameterization for the skeleton voxels (chord length)
        cc2.PixelParam{I} = ...
            cumsum([0; (sum((xyz(2:end, :) - xyz(1:end-1, :)).^2, 2)).^.5]);
        
        % if any merged branch is a leaf, then the total branch has to be a
        % leaf too
        cc2.IsLeaf(I) = any(cc.IsLeaf(cc2.MergedBranches{I}));
        
        % loops are prevented from merging, hence if we are here, the merged
        % branch cannot contain a loop
        cc2.IsLoop(I) = false;
        
        % degree of each total branch voxel
        cc2.Degree{I} = deg(dictsk(cc2.PixelIdxList{I}));
        
        % branch length of the total branch is obtained from the
        % parameterisation
        cc2.BranchLength(I) = cc2.PixelParam{I}(end);
        
    end
    
    % replace labels in the skeleton to reflect merging
    cc = cc2;
    clear cc2

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% label skeleton voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set amount of memory allowed for labels (we need label "0" for the
% background, and labels "1", "2", ..., "cc.NumObjects" for each component
%
% if a segmentation is provided to be labelled too, we need an extra "TODO"
% label
if isempty(im)
    req_bits = ceil(log2(cc.NumObjects + 1));
else
    req_bits = ceil(log2(cc.NumObjects + 2));
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
for lab = 1:cc.NumObjects
    sk(cc.PixelIdxList{lab}) = lab;
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

