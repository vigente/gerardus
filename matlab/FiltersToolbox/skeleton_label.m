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
%   We have added a new field to CC:
%
%     CC.PixelParam{i} contains the parameterization values for the voxels
%     in cc.PixelIdxList{i}. The parameterization is computed as the
%     accumulated chord distance between consecutive branch voxels. This
%     parameterization can be used, e.g. for spline interpolation.
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
% Version: 0.4.0
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

% compute degree of each voxel
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
        
        % add the bifurcation point to each branch that it finishes
        for I = 1:length(vnlab)
            cc.PixelIdxList{vnlab(I)} = ...
                [cc.PixelIdxList{vnlab(I)}; idictsk(v)];
            sk(idictsk(v)) = vnlab(I);
        end
        
        % record this in the log
        withlab(v) = true;
    end
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

% loop each branch in the skeleton
for I = 1:cc.NumObjects

    % list of voxels in current branch
    br = cc.PixelIdxList{I};
    
    % number of voxels in the current branch
    N = length(br);
    
    % degenerate case in which the branch has only one voxel
    if (N == 1)
        cc.PixelParam{I} = 0;
        continue
    end
    
    % find the two extreme voxels of the current branch
    idx = dictsk(br(deg(dictsk(br)) ~= 2));
    
    % create a distance matrix for only the voxels in the branch
    dbr = 0 * dsk;
    aux = dictsk(br);
    dbr(aux, aux) = dsk(aux, aux);

    % save the initial extreme point for later
    v0 = idx(1);
    
    % compute shortest distance from the extreme voxel to every other voxel
    % in the branch, and reuse variable
    [dbr, p] = dijkstra(dbr, v0);
    
    % convert Inf values in dbr to 0
    dbr(isinf(dbr)) = 0;
    
    % get the voxel that is furthest from the origin
    [~, idx] = max(dbr);
    
    % it can happen that 4 voxels are connected forming a rhombus. This
    % will make that v1 doesn't necessarily correspond to the branches
    % extreme point. This is a degenerate case, and shouldn't affect the
    % results
    v1 = idx;
    
    % backtrack the whole branch in order from the furthest point to the
    % original extreme point
    cc.PixelIdxList{I}(:) = 0;
    cc.PixelParam{I} = cc.PixelIdxList{I};
    cc.PixelIdxList{I}(1) = idictsk(v1);
    cc.PixelParam{I}(1) = dbr(v1);
    J = 2;
    while (idx ~= v0)
        idx = p(idx);
        cc.PixelIdxList{I}(J) = idictsk(idx);
        cc.PixelParam{I}(J) = dbr(idx);
        J = J + 1;
    end
    
    % in the degenerate case explained above, it's possible to lose a voxel
    % in the branch reordering. In this case, we have an spurious 0 index,
    % that we want to remove
    cc.PixelIdxList{I} = cc.PixelIdxList{I}(cc.PixelIdxList{I} ~= 0);
    
    % reorder voxels so that the parameterization increases monotonically
    cc.PixelIdxList{I} = cc.PixelIdxList{I}(end:-1:1);
    cc.PixelParam{I} = cc.PixelParam{I}(end:-1:1);
    
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
