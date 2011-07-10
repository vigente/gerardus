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
%                         parameterization of the branch's skeleton as a
%                         1-D spline. Branches that contain a loop cannot
%                         be parameterize, and thus CC.IsLoop(i)=NaN
%
%     CC.BranchNeighbours{i}: branch labels that share a bifurcation point
%                             with i
%
%     CC.BifurcationPixelIdx{i}: indices of the bifurcation voxels for
%                                branch i. Empty if the branch is "floating
%                                in the air" and doesn't have any
%                                bifurcation points
%
%     CC.Degree{i}:      degree of each skeleton voxel, i.e. how many voxel
%                        it is connected to
%
%     CC.MergedBranches{i}: only available if branch merging is selected by
%                           setting ALPHA>=0. Gives the list of branches in
%                           the pre-merged skeleton that were merged to
%                           create branch i
%
%   There may be a small discrepancy between cc.PixelIdxList and the voxels
%   labelled in LAB. The reason is that:
%
%     - We impose each branch to have two termination voxels. Termination
%     voxels can be leaves or bifurcation points. Bifurcation points can be
%     shared by more than one branch, so in that case they need to be
%     repeated in CC, while in SK the bifurcation voxel is assigned to only
%     one branch.
%
%     - Some bifurcation voxels are completely surrounded by other
%     bifurcation voxels. The former are not included in CC, but they are
%     given a label in SK.
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
% Version: 0.10.1
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
    alphamax = -1; % don't merge by default
end
if (nargin < 5 || isempty(p))
    p = .8;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% label skeleton voxels
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

for lab = 1:cc.NumObjects
    % give each voxel in the skeleton its label
    sk(cc.PixelIdxList{lab}) = lab;
    
    % initialize cells to contain branch neighbours
    cc.BranchNeighbours{lab} = [];
end

% keep log of voxels that have been given a label
withlab = false(size(deg));
withlab(deg < 3) = true;

% loop each bifurcation voxel (working with matrix indices, not image
% indices)
for v = find(deg >= 3)'
    % get all its neighbours that are not bifurcation points
    vn = find(dsk(v, :));
    vn = vn(deg(vn) < 3);
    
    % if this is a bifurcation voxel that is completely surroundered by
    % other bifurcation voxels, we are not going to label it
    if isempty(vn)
        continue
    end
    
    % get all its neighbour's labels
    vnlab = sk(idictsk(vn));
    
    for I = 1:length(vnlab)
        % add the bifurcation point to each branch that it finishes,
        % without duplicating voxels
        cc.PixelIdxList{vnlab(I)} = ...
            union(cc.PixelIdxList{vnlab(I)}, idictsk(v));
        sk(idictsk(v)) = vnlab(I);
    end
    
    % record this in the log
    withlab(v) = true;
end

% loop each bifurcation voxel again to obtain neighbours of each branch. We
% cannot use the previous loop because first we need to attach each
% bifurcation voxel to one or more branches, until all of them have a
% label. Only when all of them have a label, can we obtain the label
% neighbours of each branch. Note that both loops are slightly different.
% The one above looks for neighbours with degree 2 to the bifurcation
% voxel. The loop below looks for any neighbour
for v = find(deg >= 3)'
    
    % skip if this voxel doesn't belong to any branch
    if ~sk(idictsk(v))
        continue
    end
    
    % get all its neighbours
    vn = find(dsk(v, :));

    % if this is a bifurcation voxel that is completely surroundered by
    % other bifurcation voxels, we skip it
    if isempty(deg(vn) < 3)
        continue
    end
    
    % get all its neighbour's labels
    vnlab = sk(idictsk(vn));
    vnlab = vnlab(vnlab ~= 0);

    for I = 1:length(vnlab)
        % add list of neighbour labels to each label
        cc.BranchNeighbours{vnlab(I)} = union( ...
            cc.BranchNeighbours{vnlab(I)}, setdiff(vnlab, vnlab(I))');
    end
    
end

for I = 1:cc.NumObjects
    
    % make sure that cc.PixelIdxList{I} is a column vector
    cc.PixelIdxList{I} = cc.PixelIdxList{I}(:);
    
end

% % loop any voxels that have been left without a label, and assign them one
% % of a neighbour's
% for v = find(~withlab)'
%     % get its neighbours that are not bifurcation voxels
%     vn = find(dsk(v, :));
%     
%     % get its first neighbour's label
%     vnlab = sk(idictsk(vn(1)));
%     
%     % if the label is 0, that means that this voxel is surrounded by voxels
%     % that have not been labelled yet, and we leave it unlabelled
%     if vnlab
%         % give it the neighbour's label (we need every voxel in the tree to
%         % have a label, otherwise we cannot extend the segmentation)
%         % but don't add it to the list of voxels in the branch (we want
%         % every branch to contain no more than one bifurcation point)
%         sk(idictsk(v)) = vnlab;
% %         cc.PixelIdxList{vnlab} = [cc.PixelIdxList{vnlab}; idictsk(v)];
%         
%         % record this in the log
%         withlab(v) = true;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sort the skeleton voxels in each branch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init output
cc.IsLeaf = false(1, cc.NumObjects);
cc.IsLoop = false(1, cc.NumObjects);
cc.BranchLength = nan(1, cc.NumObjects);

% loop each branch in the skeleton
for I = 1:cc.NumObjects

    % list of skeleton voxels in current branch
    br = cc.PixelIdxList{I};
    
    % number of voxels in the current branch
    N = length(br);
    
    % degenerate case in which the branch has only one voxel
    if (N == 1)
        cc.PixelParam{I} = 0;
        continue
    end

    % image index => skeleton distance matrix index
    idx = dictsk(br);
    
    % image indices of skeleton voxels that either terminate or connect the
    % branch to the rest of the skeleton. We call these extreme points
    idx = br(deg(idx) ~= 2);
    
    % if the branch is a loop that is not connected to any other branches,
    % then idx==[]. In this case, any voxel can be used 
    if isempty(idx)
        idx = br(1);
        
        % image index => skeleton distance matrix index
        idx = dictsk(idx);
    
        cc.IsLeaf(I) = true;
        cc.IsLoop(I) = true;
    else
        % image index => skeleton distance matrix index
        idx = dictsk(idx);
    
        % if any extreme point has degree==1, it means that this branch is a
        % "leaf", i.e. one extreme is not connected to the rest of the skeleton
        cc.IsLeaf(I) = any(deg(idx) == 1);
    end
   
    % select any extreme point
    v0 = idx(1);
    
    % create a "branch distance matrix" for only the voxels in the branch
    dbr = 0 * dsk;
    aux = dictsk(br);
    dbr(aux, aux) = dsk(aux, aux);
    
    % the branch's skeleton can have two shapes: a "wire" or a "loop". If
    % any of the extreme points is connected to more than 1 voxel in the
    % branch distance matrix, it means that we have a loop. In principle,
    % we could try to break up loops, but there's not a real good reason to
    % do so, because does it really make sense to straighten a loop as if
    % it were a "wire"? A loop can have the following topology
    %
    %                         --
    %                         ||
    %                         ||
    %                       ------
    %
    % In this case, the loop cannot even be broken in a meaningful way.
    %
    % Thus, if we find a loop, we flag the branch as such, and leave the
    % skeleton voxels unsorted
    for J = 1:length(idx)
        % branch voxels the extreme point is connected to
        idxconn = find(dbr(idx(J), :));
        
        % the branch contains a loop
        if (length(idxconn) > 1)
            cc.IsLoop(I) = true;
            break
        end
        
    end
    
    % sort the skeleton voxels in the branch only if the branch is a "wire"
    % without loops
    if (~cc.IsLoop(I))
        % compute shortest distance from the extreme voxel to every other
        % voxel in the branch, and reuse variable
        [dbr, parents] = dijkstra(dbr, v0);
        
        % convert Inf values in dbr to 0
        dbr(isinf(dbr)) = 0;
        
        % get the voxel that is furthest from the origin
        [~, v1] = max(dbr);
        
        % backtrack the whole branch in order from the furthest point to
        % the original extreme point
        cc.PixelIdxList{I}(:) = 0;
        cc.PixelParam{I} = cc.PixelIdxList{I};
        cc.PixelIdxList{I}(1) = idictsk(v1);
        cc.PixelParam{I}(1) = dbr(v1);
        J = 2;
        while (v1 ~= v0)
            v1 = parents(v1);
            cc.PixelIdxList{I}(J) = idictsk(v1);
            cc.PixelParam{I}(J) = dbr(v1);
            J = J + 1;
        end
        
        % reorder voxels so that the parameterization increases monotonically
        cc.PixelIdxList{I} = cc.PixelIdxList{I}(end:-1:1);
        cc.PixelParam{I} = cc.PixelParam{I}(end:-1:1);
        
        % extract length of each branch
        cc.BranchLength(I) = cc.PixelParam{I}(end);
    else
        % the branch contains a loop, so it doesn't make sense to
        % parameterize or sort the skeleton voxels as a line
        cc.BranchLength(I) = nan;
        cc.PixelParam{I} = nan;
    end
    
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
%% get degree of each voxel in the skeleton
%% find bifurcation voxels for each branch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop each branch in the skeleton
for I = 1:cc.NumObjects

    % get degree of each voxel in the skeleton
    cc.Degree{I} = full(deg(dictsk(cc.PixelIdxList{I})));
    
    % find the bifurcation point indices
    idx = cc.PixelIdxList{I}(deg(dictsk(cc.PixelIdxList{I})) > 2);
    
    if ((length(idx) > 1) && cc.IsLeaf(I))
        warning(['Branch ' num2str(I) ' is a leaf but has more than 1 bifurcation voxel'])
    end
    
    % if there are no bifurcation points, this branch is floating in the
    % air
    cc.BifurcationPixelIdx{I} = idx;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge branches together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (alphamax > 0)

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
    
    % now we can loop knowing exactly the list of bifurcation voxels that
    % terminate each branch
    v = full(dictsk(unique(cat(1, cc.BifurcationPixelIdx{:}))))';
    br2merge = [];
    for I = 1:length(v)
        
        % get all its neighbours
        vn = dsk(v(I), :) > 0;
        
        % get all its neighbour's labels
        vnlab = unique(sk(idictsk(vn)));
        
        % bifurcation voxels surrounded only by bifurcation voxels don't belong
        % to any branch, and have label 0, that needs to be removed
        vnlab = vnlab(vnlab ~= 0);
        
        % remove branches that are loops, we are not going to merge those in
        % any case
        vnlab = vnlab(~cc.IsLoop(vnlab));
        
        % skip if we don't have at least 2 branches to consider merging
        if (length(vnlab) < 2 )
            continue
        end
        
        % compute all unique combinations of labels
        br0br1 = nchoosek(vnlab, 2);
        alpha = zeros(size(br0br1, 1), 1);
        
        % loop each unique combination of neighbouring branches
        for J = 1:size(br0br1, 1)
            
            % get list of indices in both skeleton branches
            idx0 = cc.PixelIdxList{br0br1(J, 1)};
            idx1 = cc.PixelIdxList{br0br1(J, 2)};
            
            % real world coordinates
            [r0, c0, s0] = ind2sub(cc.ImageSize, idx0);
            [r1, c1, s1] = ind2sub(cc.ImageSize, idx1);
            xyz0 = scinrrd_index2world([r0, c0, s0], nrrdaxis);
            xyz1 = scinrrd_index2world([r1, c1, s1], nrrdaxis);
            
            % compute distances between the first and last voxel in the current
            % and neihgbour branches
            d = dmatrix(xyz0([1 end], :)', xyz1([1 end], :)');
            
            % get matrix index with the minimum distance
            [dmin, idmin] = min(d(:));
            
            % concatenate br0 and br1 so that they are joined by the
            % bifurcation voxel
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
            
%             % DEBUG
%             hold off
%             fnplt(pp)
%             hold on
%             plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'o', 'LineWidth', 1)
%             foo = ppval(pp, pp.breaks(bifidx));
%             plot3(foo(1), foo(2), foo(3), 'ro', 'LineWidth', 1)
%             plot3(xyz(bifidx, 1), xyz(bifidx, 2), xyz(bifidx, 3), 'go', 'LineWidth', 1)
%             axis equal
%             view(2)
            
            % get curve tangent to the left and right of the bifurcation point
            % (note that dx0 doesn't necessarily correspond to br0, as opposed
            % to br1)
            dx0 = ppval(pp, pp.breaks(bifidx-1:bifidx));
            dx1 = ppval(pp, pp.breaks(bifidx:bifidx+1));
            
            dx0 = dx0(:, 2) - dx0(:, 1);
            dx1 = dx1(:, 2) - dx1(:, 1);
            
            % compute angle between derivatives
            alpha(J) = acos(abs(dot(dx0, dx1) / norm(dx0) / norm(dx1)));
            
        end
        
        % find minimum angle
        [alphamin, idx] = min(alpha);
        
        % if the minimum angle is lower than the angle threshold, then we
        % are going to merge the branches
        % as branches can share bifurcation points, it is likely that we are
        % going to get repeat mergings
        if (alphamin <= alphamax)
            br2merge = [br2merge; br0br1(idx, :)];
        end
        
    end
    
    % remove repeated merging operations
    br2merge = unique(sort(br2merge, 2), 'rows');
    
    % add to cc2 struct the fields that are common to all branches
    cc2.Connectivity = cc.Connectivity;
    cc2.ImageSize = cc.ImageSize;
    
    % create vector of new labels for the branches that are going to be merged
    newlab = zeros(1, cc.NumObjects);

    % add new labels for branches that are not going to be merged
    idx = setdiff(1:cc.NumObjects, unique(br2merge(:)));
    LAB = 0;
    for I = 1:length(idx)
        
        % create new label for cc2
        LAB = LAB + 1;
        
        % add fields to cc2
        cc2.MergedBranches{LAB} = idx(I);
        cc2.BifurcationPixelIdx{LAB} = cc.BifurcationPixelIdx{idx(I)};
        cc2.PixelIdxList{LAB} = cc.PixelIdxList{idx(I)};
        cc2.PixelParam{LAB} = cc.PixelParam{idx(I)};
        cc2.BranchNeighbours{LAB} = cc.BranchNeighbours{idx(I)};
        cc2.IsLeaf(LAB) = cc.IsLeaf(idx(I));
        cc2.IsLoop(LAB) = cc.IsLoop(idx(I));
        cc2.Degree{LAB} = cc.Degree{idx(I)};
        cc2.BranchLength(LAB) = cc.BranchLength(idx(I));
    end

    for I = 1:size(br2merge, 1)
        
        % 2 branches that have to be merged
        idx = br2merge(I, :);
        
        % if one of the branches to merge already has a label assigned ...
        if (newlab(idx(1)) ~= 0)
            
            % ... give that label to the other branch too
            newlab(idx(2)) = newlab(idx(1));
            
            % get list of indices in both skeleton branches
            idx0 = cc2.PixelIdxList{newlab(idx(1))};
            idx1 = cc.PixelIdxList{br2merge(I, 2)};
            
        elseif (newlab(idx(2)) ~= 0)
            
            % ... give that label to the other branch too
            newlab(idx(1)) = newlab(idx(2));
            
            % get list of indices in both skeleton branches
            idx0 = cc.PixelIdxList{br2merge(I, 1)};
            idx1 = cc2.PixelIdxList{newlab(idx(2))};
            
        else
            
            % ... else create a new label for both of them
            LAB = LAB + 1;
            newlab(idx) = LAB;
            
            % get list of indices in both skeleton branches
            idx0 = cc.PixelIdxList{br2merge(I, 1)};
            idx1 = cc.PixelIdxList{br2merge(I, 2)};
            
        end
        
        % save the label we are giving to the merging branches
        nowlab = newlab(idx(1));
        
        if (length(cc2.MergedBranches) < LAB)
            cc2.MergedBranches{LAB} = [];
            cc2.BifurcationPixelIdx{LAB} = [];
        end
        
        % merge the branches and put them in a new cc struct
        
        % real world coordinates
        [r0, c0, s0] = ind2sub(cc.ImageSize, idx0);
        [r1, c1, s1] = ind2sub(cc.ImageSize, idx1);
        xyz0 = scinrrd_index2world([r0, c0, s0], nrrdaxis);
        xyz1 = scinrrd_index2world([r1, c1, s1], nrrdaxis);
        
        % compute distances between the first and last voxel in the current
        % and neihgbour branches
        d = dmatrix(xyz0([1 end], :)', xyz1([1 end], :)');
        
        % get matrix index with the minimum distance
        [dmin, idmin] = min(d(:));
        
        % concatenate br0 and br1 so that they are joined by the
        % bifurcation voxel
        if (idmin == 1) % 1st voxel branch <=> 1st voxel neighbour
            if (dmin == 0)
                idx = [idx0(end:-1:1); idx1(2:end)];
                xyz = [xyz0(end:-1:1, :); xyz1(2:end, :)];
            else
                idx = [idx0(end:-1:1); idx1];
                xyz = [xyz0(end:-1:1, :); xyz1];
            end
        elseif (idmin == 2) % end voxel branch <=> 1st voxel neighbour
            if (dmin == 0)
                idx = [idx0; idx1(2:end)];
                xyz = [xyz0; xyz1(2:end, :)];
            else
                idx = [idx0; idx1];
                xyz = [xyz0; xyz1];
            end
        elseif (idmin == 3) % 1st voxel branch <=> end voxel neighbour
            if (dmin == 0)
                idx = [idx1; idx0(2:end)];
                xyz = [xyz1; xyz0(2:end, :)];
            else
                idx = [idx1; idx0];
                xyz = [xyz1; xyz0];
            end
        elseif (idmin == 4) % end voxel branch <=> end voxel neighbour
            if (dmin == 0)
                idx = [idx1; idx0(end-1:-1:1)];
                xyz = [xyz1; xyz0(end-1:-1:1, :)];
            else
                idx = [idx1; idx0(end:-1:1)];
                xyz = [xyz1; xyz0(end:-1:1, :)];
            end
        end
        
        cc2.PixelIdxList{nowlab} = idx;
        
        % recompute parameterization for the skeleton voxels (chord length)
        cc2.PixelParam{nowlab} = ...
            cumsum([0; (sum((xyz(2:end, :) - xyz(1:end-1, :)).^2, 2)).^.5]);
        
        % note which branches from the original cc have contributed to this
        % total branch
        if (isempty(cc2.MergedBranches{nowlab}))
            cc2.MergedBranches{nowlab} = br2merge(I, :);
        else
            cc2.MergedBranches{nowlab} = ...
                union(cc2.MergedBranches{nowlab}, br2merge(I, :));
        end
        
        % total branch neighbours are the union of neighbours from each
        % sub-branch
        cc2.BranchNeighbours{nowlab} = unique(cat(2, ...
            cc.BranchNeighbours{br2merge(I, :)}));
        
        % if any merged branch is a leaf, then the total branch has to be a
        % leaf too
        cc2.IsLeaf(nowlab) = any(cc.IsLeaf(cc2.MergedBranches{nowlab}));
        
        % loops are prevented from merging, hence if we are here, the merged
        % branch cannot contain a loop
        cc2.IsLoop(nowlab) = false;
        
        % degree of each total branch voxel
        cc2.Degree{nowlab} = deg(dictsk(cc2.PixelIdxList{nowlab}));
        
        % branch length of the total branch is obtained from the
        % parameterisation
        cc2.BranchLength(nowlab) = cc2.PixelParam{nowlab}(end);
        
        % bifurcation voxels in the total branch are the union  of bifurcation
        % voxels in each sub-branch
        cc2.BifurcationPixelIdx{nowlab} = ...
            union(cc2.BifurcationPixelIdx{nowlab}, ...
            cat(1, cc.BifurcationPixelIdx{br2merge(I, :)}));
        
    end
    
    % number of total branches in the merged skeleton
    cc2.NumObjects = LAB;
    
    % replace labels in the skeleton to reflect merging
    cc = cc2;
    for I = 1:cc.NumObjects
        % give each voxel in the skeleton its label
        sk(cc.PixelIdxList{I}) = I;
    end

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% label all voxels in the original segmentation, not only skeleton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isempty(im))
    
    % convert im data type to be the same as the seeds, so that we are
    % guaranteed that all label values can be represented in the output
    % image
    im = cast(im, class(sk));
    
    % label for voxels that need to be labelled
    TODO = cast(N+1, class(sk));
    
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

