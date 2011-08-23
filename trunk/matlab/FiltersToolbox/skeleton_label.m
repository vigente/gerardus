function [sk, cc, bifcc, mcon, madj, cc2, mmerge] = skeleton_label(sk, im, res, alphamax, p, SINGLEMERGE)
% SKELETON_LABEL  Give each branch of a skeleton a different label, and
% sort the voxels within each branch
%
% [LAB, CC, BIFCC, MCON, MADJ] = SKELETON_LABEL(SK, [], RES)
%
%   SK is a 3D segmentation mask. SK is assumed to be a skeleton resulting
%   from some kind of thinning algorithm, e.g. 
%     >> sk = itk_filter('skel', im);
%   (itk_filter() is a Matlab function provided by Gerardus). That is, the
%   segmentation looks like a series of 1 voxel-thick branches connected
%   between them (cycles are allowed).
%
%   RES is a 3-vector with the voxel size as [row, column, slice]. By
%   default, RES=[1 1 1].
%
%   LAB is an image of the same dimensions of SK where the value of each
%   voxel is the label of the branch it belongs to.
%
%   Voxels in bifurcations are labelled as the nearest branch. Voxels
%   completely surrounded by other bifurcation voxels are assigned the
%   nearest label arbitrarily.
%
%   CC is a struct like the one provided by Matlab's function bwconncomp(),
%   but with some extra fields
%
%     CC.PixelIdxList{i}: list of connected voxels that form the i-th
%                       branch. These indices are image indices, so the
%                       row, column, slice coordinates of e.g. branch 5 can
%                       be obtained doing
%
%                         [r, c, s] = ind2sub(cc.ImageSize, cc.PixelIdxList{5});
%
%     CC.PixelParam{i}: parameterization values for the voxels in
%                       cc.PixelIdxList{i}. The parameterization is
%                       computed as the accumulated chord distance between
%                       consecutive branch voxels. This parameterization
%                       can be used, e.g. for spline interpolation.
%
%     CC.IsLeaf(i):     flags indicating whether each section is a "leaf",
%                       i.e. whether it's a branch with a free extreme
%
%     CC.BranchLength(i): chord-length of each skeleton branch. This is a
%                       parameterization of the branch's skeleton as a
%                       1-D spline.
%
%     CC.Degree{i}:     degree of each branch voxel, i.e. how many voxels
%                       it is connected to.
%
%   BIFCC is a struct a bit like CC, but only with field PixelIdxList
%
%     BIFCC.PixelIdxList{i}: list of connected voxels that form the i-th
%                       bifurcation clump.
%
%   MCON is a boolean sparse matrix, where rows = branches, columns =
%   bifurcation clumps. MCON(5, 2)==true means that branch 5 is connected
%   to bifurcation clump 2.
%
%   MADJ is a square sparse matrix where MADJ(7, 3)==10 means that branches
%   7 and 3 and connected through the bifurcation clump 10.
%
%
% [LAB, ...] = SKELETON_LABEL(SK, IM, RES)
%
%   IM is an array with the same size as SK, and contains the whole
%   segmentation. SK is the skeleton of IM.
%
%   If IM is provided, then LAB will contain IM labelled, instead of SK
%   labelled. IM is labelled using a region grow algorithm on the labelled
%   skeleton. If you want to extract the labelling for SK, just run
%
%     >> LAB .* SK
%
% [..., CC2] = SKELETON_LABEL(SK, IM, RES, ALPHAMAX, P, SINGLEMERGE)
%
%   With this syntax you can align branches that are well aligned with each
%   other.
%
%   ALPHAMAX is an angle in radians. ALPHAMAX >= 0 means that branch
%   merging is performed. If the angle between those two branches is <=
%   ALPHAMAX, then they are merged. By default ALPHAMAX = -Inf, so no
%   merging is performed.
%
%   P is a scalar in [0, 1]. To compute the angle between branches, an
%   approximating or smoothing cubic spline is fit to the skeleton voxels
%   using csaps(..., P). P=0 is the smoothest spline (a straight line with
%   the least squares approximation), while P=1 is a rugged spline (the
%   spline interpolated the voxels). Adequate values of P depend on the
%   image resolution, so it's difficult to propose a formula. As a rule of
%   thumb, it seems that if resolution is in the order 1e-n, then a good
%   value for P=1-1e-n. For example, if resolution is in the order 1e-5, 
%   P=1-1e-5=0.999999. By default, P=1 and no smoothing is performed.
%
%   SINGLEMERGE is a boolean flag. If SINGLEMERGE=true, then at each
%   bifurcation clump only the two branches with the smallest angle between
%   them are considered for merging. If SINGLEMERGE=false, then any pair of
%   branches with an angle smaller than ALPHAMAX will be merged. By
%   default, SINGLEMERGE=true.
%
%   Note: If you want to see how branches are being merged and smoothed,
%   uncomment the DEBUG block at the end of internal function
%   angle_btw_branches().
%
%   CC2 is a struct like CC, but with two additional fields
%
%     CC2.MergedBranches{i}: List of branches in the pre-merged skeleton
%                      that were merged to create branch i
%
%     CC2.MergedBifClumps{i}: List of bifurcation clumps that were merged
%                      to create branch i
%
%
% See also: skeleton_plot, scinrrd_skeleton_prune.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.13.0
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
error(nargoutchk(0, 8, nargout, 'struct'));

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
if (nargin < 6 || isempty(SINGLEMERGE))
    SINGLEMERGE = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% split skeleton into branches and bifurcation voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get sparse matrix of distances between voxels. To label the skeleton we
% don't care about the actual distances, just need to know which voxels are
% connected to others. Actual distances are needed to parameterize the
% branches, though
[dsk, dictsk, idictsk] = seg2dmat(sk, 'seg', res);

%% find bifurcation voxels

% compute degree of each skeleton voxel
deg = sum(dsk > 0, 2);

% get distance matrix index of the bifurcation voxels
bifidx = deg >= 3;

% matrix index => image index
bifidx = idictsk(bifidx);

%% label connected components of skeleton branches

% remove bifurcation voxels from original skeleton
sk(bifidx) = 0;

% get connected components in the image
cc = bwconncomp(sk);

% make size vector always have size(3)
if length(cc.ImageSize) == 2
    cc.ImageSize = [cc.ImageSize 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sort the skeleton voxels in each branch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init outputs
cc.PixelParam = cell(1, cc.NumObjects);
cc.IsLeaf = false(1, cc.NumObjects);
cc.BranchLength = zeros(1, cc.NumObjects);
cc.Degree = cell(1, cc.NumObjects);

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

    % sort the voxels in the branch
    [cc.PixelIdxList{I}, cc.PixelParam{I}] = ...
        sort_branch(cc.PixelIdxList{I}, dsk, dictsk, idictsk);
    
    % extract length of branch
    cc.BranchLength(I) = cc.PixelParam{I}(end);
    
    % degree of each total branch voxel
    cc.Degree{I} = full(deg(dictsk(cc.PixelIdxList{I})));
    
    % a branch is a leaf is it has at least one voxel with degree==1 (tip
    % of a branch) or 0 (single voxel floating in space)
    cc.IsLeaf(I) = any(cc.Degree{I}([1 end]) < 2);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% connectivity between branches and bifurcation clumps,
%% and find which branches should be merged together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tag each branch with its label
sk = labelmatrix(cc);

% create empty image volume and add only bifurcation voxels
sk2 = zeros(size(sk), 'uint8');
sk2(bifidx) = 1;

% compute connected components to obtain clumps of bifurcation voxels
bifcc = bwconncomp(sk2);

% init matrix to describe the connection between branches
%
%   row = branch index
%   col = bifurcation clump index
mcon = boolean(sparse(cc.NumObjects, bifcc.NumObjects));

% init adjacency matrix that says which branch to connected to which
% branches, via which bifurcation clump
madj = sparse(cc.NumObjects, cc.NumObjects);

% init matrix to keep track of which branches are merged
if (alphamax >= 0)
    mmerge = sparse(cc.NumObjects, cc.NumObjects);
else
    mmerge = [];
end

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
    boxbr = sk(r0:rend, c0:cend, s0:send);
    
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
    
    % add the neighbour connections to the connection and adjacency
    % matrices
    mcon(idx, I) = true;
    madj(idx, idx) = I;
    
    %% measure angle between pairs of branches
    
    if (alphamax >= 0)
        
        % number of branches
        N = length(idx);
        
        % in order to merge, the bifurcation needs to have at least 2
        % branches
        if (N < 2)
            continue
        end
        
        % preserve the list of branches for later
        idx0 = idx;
        
        % get all combinations of pairs of branches that share this clump
        idx = nchoosek(idx, 2);
        
        % get voxel indices for bifurcation clump
        bif = bifcc.PixelIdxList{I};
        
        % loop each pair of branches combination
        alpha = zeros(1, size(idx, 1));
        for J = 1:size(idx, 1)
            
            % get voxel indices for each branch
            br0 = cc.PixelIdxList{idx(J, 1)};
            br1 = cc.PixelIdxList{idx(J, 2)};
            
            % merge and sort both branches and the bifurcation clump
            br = sort_branch([br0(:); bif(:); br1(:)], ...
                dsk, dictsk, idictsk);
            
            % compute angle between the branches if they are merged
            alpha(J) = angle_btw_branches(br, bif, nrrdaxis, p);
            
        end
        
        if (SINGLEMERGE) % only 2 branches can be merged per bifurcation clump
            
            % get the two branches with the smallest angle
            [alphamin, J] = min(alpha);
            
            % if the smallest angle is small enough, we mark these two branches
            % to be merged via the current bifurcation clump (note that we
            % don't know yet which label the merged branches will have, we only
            % know that they need to be merged)
            if (alphamin <= alphamax)
                mmerge(idx(J, 1), idx(J, 2)) = true;
                mmerge(idx(J, 2), idx(J, 1)) = true;
            end
            
        else % any pair of suitable branches can be merged
            
            % create small matrix to write the angles between the branches
            % at this bifurcation clump
            malpha = zeros(N);
            
            % populate the matrix with the computed angles
            aux = nchoosek(1:N, 2);
            malpha(sub2ind([N N], aux(:, 1), aux(:, 2))) = alpha;
            malpha = malpha + malpha';
            
            % make the main diagonal Inf so that it doesn't count when we
            % look for the minimum angle
            malpha(1:N+1:end) = Inf;
            
            % loop all branches at this bifurcation clump
            for J = 1:N
                
                % find the best-aligned branch for current branch
                [alphamin, idxmin] = min(malpha(J, :));
                
                % if the aligment is larger than the threshold angle, then
                % we are not going to merge these branches
                if (alphamin > alphamax)
                    continue
                end
                
                % if branch A wants to merge with branch B, check that
                % branch B also wants to merge with branch A
                [~, idxmin2] = min(malpha(idxmin, :));
                
                if (idxmin2 == J)
                    
                    % if so, add both branches to the list of branches to
                    % merge
                    mmerge(idx0(J), idx0(idxmin)) = true;
                    
                end
                
            end
            
        end
        
    end

end

% mark branches to be merged with themselves, so that we can keep track of
% single branches in mmerge
if (alphamax >= 0)
    mmerge(1:size(mmerge, 1)+1:numel(mmerge)) = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge branches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list of branches
if (alphamax >= 0)
    br = 1:cc.NumObjects;
else
    br = [];
end

% init struct to contain the merged branches
cc2.Connectivity = cc.Connectivity;
cc2.ImageSize = cc.ImageSize;
cc2.NumObjects = 0;
cc2.PixelIdxList = [];
cc2.PixelParam = [];
cc2.MergedBranches = [];
cc2.MergedBifClumps = [];
cc2.IsLeaf = [];
cc2.BranchLength = [];
cc2.Degree = [];

% loop branches
I = 0; % label of the current merged branch
while (~isempty(br))
    
    % new label for merged branches
    I = I + 1;

    % get all the branches the first branch in the list is connected to,
    % and sort them forming a chain
    idx = sort_connected_branches(br(1), mmerge);
    
    % keep track of which branches have been merged
    cc2.MergedBranches{I} = idx;
    mmerge(idx, idx) = mmerge(idx, idx) * I;
    
    % remove merged branches from the list of branches
    N = length(br);
    br = setdiff(br, idx);
    if (N == length(br))
        error('Assertion fail: No branches removed, we have entered an infinite loop')
    end
    
    % find all intermediate bifurcation clumps that connect together the
    % branches
    bifidx = full(madj(sub2ind(size(madj), idx(1:end-1), idx(2:end))))';
    if any(bifidx == 0)
        error('Assertion fail: Branches are connected but have no bifurcation clump between them')
    end
    cc2.MergedBifClumps{I} = bifidx';
    
    % get voxels from all the branches and the bifurcation clumps
    v = [cat(1, cc.PixelIdxList{idx}); cat(1, bifcc.PixelIdxList{bifidx})];

    % sort the voxels so that they form a new branch, and get the
    % parameterisation
    [cc2.PixelIdxList{I}, cc2.PixelParam{I}] = ...
        sort_branch(v, dsk, dictsk, idictsk);
    cc2.PixelIdxList{I} = cc2.PixelIdxList{I}(:);
    cc2.PixelParam{I} = cc2.PixelParam{I}(:);
    
    % get some more info about the merged branch
    cc2.Degree{I} = full(deg(dictsk(cc2.PixelIdxList{I})));
    cc2.IsLeaf(I) = any(cc2.Degree{I} < 2);
    cc2.BranchLength(I) = cc2.PixelParam{I}(end);
    
end

% get number of merged branches
cc2.NumObjects = length(cc2.PixelIdxList);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% label skeleton voxels or segmentation voxels, merged or not merged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% skeleton labelling

if (alphamax >= 0) % merging
    
    % add a "TODO" label for intersection voxels
    cc2.NumObjects = cc2.NumObjects + 1;
    cc2.PixelIdxList(end+1) = {[]};

    % label all branch skeleton voxels
    sk = labelmatrix(cc2);
    TODO = sk(1)*0 + cc2.NumObjects;

    % add the bifurcation voxels that are not part of branches
    sk(setdiff(cat(1, bifcc.PixelIdxList{:}), ...
        cat(1, cc2.PixelIdxList{:}))) = TODO;
    
else % not merging

    % add a "TODO" label for intersection voxels
    cc.NumObjects = cc.NumObjects + 1;
    cc.PixelIdxList(end+1) = {[]};

    % label all branch skeleton voxels
    sk = labelmatrix(cc);
    TODO = sk(1)*0 + cc.NumObjects;

    % add the bifurcation voxels that are not part of branches
    sk(setdiff(cat(1, bifcc.PixelIdxList{:}), ...
        cat(1, cc.PixelIdxList{:}))) = TODO;

end

% if a whole segmentation is provided, then we are also going to segment it
if (~isempty(im))
    
    sk((im ~= 0) & (sk == 0)) = TODO;
    
end

% region grow algorithm to extend branch labels
sk = bwregiongrow(sk, TODO, res);

% in some very particular cases, a small patch of voxels may be left
% unlabelled. We are just going to remove them from the segmentation
sk(sk == TODO) = 0;

if (alphamax >= 0) % merging
    
    % remove the empty list of voxels we added for the TODO label
    cc2.PixelIdxList(end) = [];
    cc2.NumObjects = cc2.NumObjects - 1;
    
else % not merging

    % remove the empty list of voxels we added for the TODO label
    cc.PixelIdxList(end) = [];
    cc.NumObjects = cc.NumObjects - 1;
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort a set of branches
% idx: index of one of the branches in the branch set
% d: adjacency matrix
function idx = sort_connected_branches(idx, d)

% type double is necessary, otherwise distances will be
% wrong, because booleans can only take 0 or 1 values
if (~strcmp(class(d), 'double'))
    error('Adjacency matrix must be of type double')
end

% compute Dijkstra's shortest path from the branch to any other branch
[dp, ~] = dijkstra(d, idx);

% indices of all branches connected to the input branch, and the input
% branch itself
idx = find(~isinf(dp));

% if there's no need to sort the branch, just return all the branches
if (length(idx) <= 2)
    return
end

% three branches don't need to be sorted if they form a cycle
if (length(idx) == 3)
    if (d(idx(1), idx(2)) ...
            && d(idx(1), idx(3)) ...
            && d(idx(2), idx(3)))
        return
    end
end

% the largest distance indicates one of the branches at the end of the
% merge
dp(isinf(dp)) = 0;
[~, v0] = max(dp);

% compute shortest path from v0 to all other branches
[dp, parents] = dijkstra(d, v0);
dp(isinf(dp)) = 0;
[~, v1] = max(dp);

% get the whole string of branches from v1 to v0
idx = [];
while (v1 ~= 0)
    idx(end+1) = v1;
    v1 = parents(v1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% angle_btw_branches(): compute angle between two concatenated branches
function alpha = angle_btw_branches(br, bif, nrrdaxis, p)

% real world coordinates of voxels
[r, c, s] = ind2sub([nrrdaxis.size], br);
xyz = scinrrd_index2world([r, c, s], nrrdaxis);
        
% compute spline parameterization (Lee's centripetal scheme)
t = cumsum([0; (sum((xyz(2:end, :) - xyz(1:end-1, :)).^2, 2)).^.25]);

% compute cubic smoothing spline
pp = csaps(t', xyz', p);

% compute where are the bifurcation voxels within the sorted branch
[~, bifidx] = intersect(br, bif);

% compute the middle parameterisation value for the bifurcation voxels
tbif = mean(t(bifidx));

% compute mean step size in the parameterisation
dt = mean(diff(t));

% compute tanget vectors 1 step before and after the bifurcation junction
v0 = ppval(pp, [tbif-dt tbif]);
v1 = ppval(pp, [tbif tbif+dt]);

v0 = v0(:, 2) - v0(:, 1);
v1 = v1(:, 2) - v1(:, 1);

% compute angle between tangent vectors
alpha = acos(dot(v0, v1) / norm(v0) / norm(v1));

% % DEBUG
% aux = ppval(pp, tbif);
% hold off
% fnplt(pp, 'r')
% hold on
% plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'b')
% plot3(aux(1), aux(2), aux(3), 'bo')
% axis xy equal

end
