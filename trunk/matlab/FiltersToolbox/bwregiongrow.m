function im = bwregiongrow(im, seed, res)
% BWREGIONGROW  Region grow segmentation of binary image from multiple
% seeds
%
% LAB = BWREGIONGROW(IM, SEED)
%
%   IM is a 2D matrix or 3D array with a binary segmentation.
%
%   SEED has the same size as IM, and contains the initial seed values for
%   the segmentation. For example, if SEED(45)=3, this means that voxel 45
%   in IM has label 3.
%
%   LAB has the same size as IM, and has the label values that partition
%   the original binary image IM. Each partition has a different label.
%   Partitions are computed with a region grow algorithm that expands the
%   labels from the seeds.
%
%   At each iteration, partitions grow 1 voxel until the whole IM is
%   labelled.
%
% LAB = BWREGIONGROW(..., RES)
%
%   RES is a 2-vector (in 2D) or 3-vector (in 3D) with the voxel size in
%   each dimension. By default, it is assumed that RES=[1, 1, 1]. Voxel
%   size is used to compute distances between voxels in the labelling
%   process.

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
error(nargchk(2, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(res))
    res = [1 1 1];
end

% maximum number of voxels connected to each voxel, depending on whether
% the image is 2D or 3D
if (size(im, 3) == 1)
    nnsize = 8;
else
    nnsize = 26;
end    

% compute distance matrix for the whole image segmentation so that we have
% voxel connectivity
[d, dict, idict] = seg2dmat(im, 'seg', res);

% the highest value for the image data type
MAXLAB = zeros(1, class(seed));
MAXLAB(1) = Inf; % this casts Inf to the image data type. E.g., if image is 
                 % uint8, then MAXLAB==255

% get highest label in the skeleton
N = max(seed(:));
if (N >= MAXLAB)
    error('Current image data type doesn''t allow processing. You need to use a larger data type, e.g. uint16 instead of uint8')
end

% label for voxels that need to be labelled
TODO = N+1;

% mark segmentation voxels with the highest label, so that we know that
% they still need to be labelled, but belong to the segmentation as opposed
% to the background
im(im > 0) = TODO;

% init algorithm by copying the labelled skeleton to the segmentation
idx = seed > 0;
im(idx) = seed(idx);

% % DEBUG: plot algorithm initialization
% imagesc(im)
% axis equal xy

% loop until all voxels have been labelled
imiter = 1;
while (any(imiter(:)))
    
    % get indices of all voxels already labelled
    idxlab = (im(:) ~= TODO) & (im(:) ~= 0);
    
    % image index => d matrix index
    idxlab = full(dict(idxlab));
    
    % create another distance matrix for only voxels that have been
    % labelled
    diter = d(idxlab, idxlab);
    
    % compute degree of labelled voxels
    degiter = full(sum(diter>0, 2));
    
    % find which labelled voxels are at the boundary. Boundary voxels in 2D
    % are those with degree < 8, and degree < 26 in 3D
    idxperim = degiter < nnsize;
    
    % keep only degree of boundary voxels
    degiter = degiter(idxperim);
    
    % diter index => d matrix index
    idxperim = idxlab(idxperim);
    
%     % DEBUG: plot boundary voxels
%     aux = idict(idxperim);
%     auxim = 0 * im;
%     auxim(aux) = 1;
%     imagesc(auxim)
%     axis xy equal
    
    % get degree of the boundary voxels in the whole segmentation
    deg = full(sum(d(idxperim, :)>0, 2));
    
    % voxels with deg > degiter are connected to non-labelled voxels, and
    % hence need to be propagated. We can forget about the rest of the
    % boundary voxels
    idxperim = idxperim(deg > degiter);
    
    % create mask with the voxels that need to be labelled in this
    % iteration, i.e. voxels that are next to the boundary and haven't been
    % labelled yet
    imiter = false(size(im));
    % loop boundary voxels
    for I = 1:length(idxperim)
        % get indices of voxels attached to this boundary voxel
        idxtodo = idict(d(idxperim(I), :)>0);
        
        % tag voxels that are candidates to be labelled in this iteration
        % (some of them may have been labelled already)
        imiter(idxtodo) = true;
    end
    
    % remove as candidates those that have a label already. After this, the
    % remaining candidates all need to be labelled
    idxtodo = find(imiter);
    idxtodo = idxtodo(im(idxtodo) == TODO);
    
    % loop candidates that need to be labelled
    new_labels = zeros(1, length(idxtodo));
    for I = 1:length(idxtodo)
        % get d matrix indices of voxels adjacent to this one
        idxnn = find(d(dict(idxtodo(I)), :) > 0);
        
        % keep only adjacent voxels that are labelled already
        idxnn = idxnn(im(idict(idxnn)) ~= TODO);
        
        % of all the adjacent voxels with labels, choose the closest one
        [~, idx] = min(d(dict(idxtodo(I)), idxnn));
        
        % save label for this candidate voxel
        new_labels(I) = im(idict(idxnn(idx)));
        
    end
    
    % actually label the voxels (it's important to wait until now, because
    % if we label the voxels within the loop, we get spill overs of some
    % regions into other regions
    im(idxtodo) = new_labels;
    
%     % DEBUG: plot voxels to label
%     imagesc(im)
%     axis equal xy
%     pause
    
end
