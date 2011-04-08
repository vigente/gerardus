function [sk, cc] = skeleton_label(sk)
% SKELETON_LABEL  Give each branch of a skeleton a different label
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
%   Each vector in CC.PixelIdxList{i} has the list of image indices of the
%   voxels with label i.

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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
error(nargchk(1, 1, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% get sparse matrix of distances between voxels. We don't care about the
% actual distances, just need to know which voxels are connected to others
[d, ~, idict] = seg2dmat(sk, 'seg');

% compute degree of each voxel
deg = sum(d > 0, 2);

% get distance matrix index of the bifurcation voxels
idx = deg >= 3;

% convert distance matrix index to image index
idx = idict(idx);

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
    vn = find(d(v, :));
    vn = vn(deg(vn) < 3);
    
    % if this is a bifurcation voxel that is completely surroundered by
    % other bifurcation voxels, we are not going to label it
    if ~isempty(vn)
        % get its first neighbour's label
        vnlab = sk(idict(vn(1)));
        
        % give it the neighbour's label
        cc.PixelIdxList{vnlab} = [cc.PixelIdxList{vnlab}; idict(vn)];
        sk(idict(v)) = vnlab;
        
        % record this in the log
        withlab(v) = true;
    end
end

% loop any voxels that have been left without a label, and assign them one
% of a neighbour's
for v = find(~withlab)'
    % get its neighbours that are not bifurcation voxels
    vn = find(d(v, :));
    
    % get its first neighbour's label
    vnlab = sk(idict(vn(1)));
    
    % if the label is 0, that means that this voxel is surrounded by voxels
    % that have not been labelled yet, and we leave it unlabelled
    if vnlab
        % give it the neighbour's label
        cc.PixelIdxList{vnlab} = [cc.PixelIdxList{vnlab}; idict(vn)];
        sk(idict(v)) = vnlab;
        
        % record this in the log
        withlab(v) = true;
    end
end

% % Debug (2D image) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% hold off
% imagesc(sk)
% 
% % find tree leaves: convert matrix indices to image linear indices
% idx = idict(deg == 1);
% 
% % convert image linear indices to row, col, slice
% [r, c, s] = ind2sub(size(sk), idx);
% 
% % plot leaves
% hold on
% plot(c, r, 'ro')
% 
% % find tree bifurcations: convert matrix indices to image linear indices
% idx = idict(deg >= 3);
% 
% % convert image linear indices to row, col, slice
% [r, c, s] = ind2sub(size(sk), idx);
% 
% % plot leaves
% plot(c, r, 'go')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
