function [sk, cc, dsk, dictsk, idictsk] = skeleton_label(sk, im, res)
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
%                         be parameterize, and thus CC.IsLoop(i)=NaN.
%
%     CC.BifurcationPixelIdx(i): index of the bifurcation voxel for
%                                branches that are leaves. Branches that
%                                are not leaves get a 0 index.
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
% See also: skeleton_plot.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.7.0
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
error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(0, 5, nargout, 'struct'));

% defaults
if (nargin < 2)
    im = [];
end
if (nargin < 3 || isempty(res))
    res = [1 1 1];
end

%% Label skeleton voxels

% get sparse matrix of distances between voxels. To label the skeleton we
% don't care about the actual distances, just need to know which voxels are
% connected to others. Actual distances are needed to parameterize the
% branches, though
[dsk, dictsk, idictsk] = seg2dmat(sk, 'seg', res);

% compute degree of each skeleton voxel
deg = sum(dsk > 0, 2);

% get distance matrix index of the bifurcation voxels
idx = deg >= 3;

% convert distance matrix index to image index
idx = idictsk(idx);

% remove bifurcation voxels from image
sk(idx) = 0;

% get connected components in the image
cc = bwconncomp(sk);

% set amount of memory allowed for labels
req_bits = ceil(log2(cc.NumObjects));
if (req_bits == 1)
    lab_class = 'bool';
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
    sk = zeros(size(sk), lab_class);
else
    sk = sk * 0;
end

% give each voxel in the image its label
for lab = 1:cc.NumObjects
    sk(cc.PixelIdxList{lab}) = lab;
end

% keep log of voxels that have been given a label
withlab = false(size(deg));
withlab(deg < 3) = true;

% loop each bifurcation voxel (working with matrix indices, not image
% indices)
for v = find(deg >= 3)'
    % get its neighbours that are not bifurcation voxels
    vn = find(dsk(v, :));
    vn = vn(deg(vn) < 3);
    
    % if this is a bifurcation voxel that is completely surroundered by
    % other bifurcation voxels, we are not going to label it
    if ~isempty(vn)
        % get all its neighbour's labels
        vnlab = sk(idictsk(vn));
        
        % add the bifurcation point to each branch that it finishes,
        % without duplicating voxels
        for I = 1:length(vnlab)
            cc.PixelIdxList{vnlab(I)} = ...
                union(cc.PixelIdxList{vnlab(I)}, idictsk(v));
            sk(idictsk(v)) = vnlab(I);
        end
        
        % record this in the log
        withlab(v) = true;
    end
end

% if cc.PixelIdxList{I} had only 1 element, now it is a row vector, instead
% of a column vector. Make sure that they are all column vectors
for I = 1:cc.NumObjects
    cc.PixelIdxList{I} = cc.PixelIdxList{I}(:);
end

% loop any voxels that have been left without a label, and assign them one
% of a neighbour's
for v = find(~withlab)'
    % get its neighbours that are not bifurcation voxels
    vn = find(dsk(v, :));
    
    % get its first neighbour's label
    vnlab = sk(idictsk(vn(1)));
    
    % if the label is 0, that means that this voxel is surrounded by voxels
    % that have not been labelled yet, and we leave it unlabelled
    if vnlab
        % give it the neighbour's label (we need every voxel in the tree to
        % have a label, otherwise we cannot extend the segmentation)
        % but don't add it to the list of voxels in the branch (we want
        % every branch to contain no more than one bifurcation point)
        sk(idictsk(v)) = vnlab;
%         cc.PixelIdxList{vnlab} = [cc.PixelIdxList{vnlab}; idictsk(v)];
        
        % record this in the log
        withlab(v) = true;
    end
end

%% Label all voxels in the original segmentation, not only skeleton

if (~isempty(im))
    
    % check that the image data type size is large enough to add the TODO
    % label
    MAXLAB = zeros(1, class(sk));
    MAXLAB(1) = Inf; % this casts Inf to the image data type. E.g., if 
                     % image is uint8, then MAXLAB==255

    % get highest label in the skeleton
    N = max(sk(:));
    if (N >= MAXLAB)
        error('Current image data type size doesn''t allow to add a new label. Cast to a larger data type, e.g. uint16 instead of uint8')
    end

    % convert im data type to be the same as the seeds, so that we are
    % guaranteed that all label values can be represented in the output
    % image
    im = cast(im, class(sk));
    
    % label for voxels that need to be labelled
    TODO = N+1;
    
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

%% Sort the skeleton voxels in each branch

% init output
cc.IsLeaf = false(1, N);
cc.IsLoop = false(1, N);
cc.BranchLenght = nan(1, N);

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
        [dbr, p] = dijkstra(dbr, v0);
        
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
            v1 = p(v1);
            cc.PixelIdxList{I}(J) = idictsk(v1);
            cc.PixelParam{I}(J) = dbr(v1);
            J = J + 1;
        end
        
        % reorder voxels so that the parameterization increases monotonically
        cc.PixelIdxList{I} = cc.PixelIdxList{I}(end:-1:1);
        cc.PixelParam{I} = cc.PixelParam{I}(end:-1:1);
        
        % extract length of each branch
        cc.BranchLenght(I) = cc.PixelParam{I}(end);
    else
        % the branch contains a loop, so it doesn't make sense to
        % parameterize or sort the skeleton voxels as a line
        cc.BranchLenght(I) = nan;
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

%% find bifurcation voxels for each branch

cc.BifurcationPixelIdx = zeros(1, cc.NumObjects);

% loop each branch in the skeleton
for I = 1:cc.NumObjects

    % skip branches that are not leaves
    if (~cc.IsLeaf(I))
        continue
    end
    
    % find the bifurcation point index
    idx = cc.PixelIdxList{I}(deg(dictsk(cc.PixelIdxList{I})) > 2);
    
    if (length(idx) > 1)
        warning(['Branch ' num2str(I) ' has more than 1 bifurcation voxel'])
    end
    
    % if there are no bifurcation points, this branch is floating in the
    % air
    if (~isempty(idx))
        cc.BifurcationPixelIdx(I) = idx;
    end
    
end
