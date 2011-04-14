function stats = scinrrd_seg2label_stats(nrrd, cc, STRAIGHT, d, dict)
% SCINRRD_SEG2LABEL_STATS  Shape stats for each label in a multi-label
% segmentation
%
% STATS = SCINRRD_SEG2LABEL_STATS(NRRD, CC)
%
%   NRRD is an SCI NRRD struct with a labelled segmentation mask of
%   different structures, e.g. blobs, tubes, etc.
%
%   CC is a struct produced by function skeleton_label() with the list of
%   skeleton voxels that belongs to each label, and the parameterization
%   vector for the skeleton.
%
%   The labels can be created, e.g.
%
%     >> [nrrd.data, cc] = skeleton_label(sk.data, seg.data, [seg.axis.spacing]);
%
%     where seg.data is a binary segmentation mask of all the image,
%     sk.data is a binary segmentation mask of the skeleton.
%
%   STATS is a struct with the shape parameters computed for each label in
%   the segmentation. The provided measures are:
%
%     STATS.VAR: Variance in the three principal components of the cloud
%                of voxels that belong to each label. These are the
%                ordered eigenvalues obtained from computing Principal
%                Component Analysis on the voxel coordinates.
%
% STATS = SCINRRD_SEG2LABEL_STATS(..., STRAIGHT)
%
%   STRAIGHT is a flag. If STRAIGHT=true, then each label is warped so that
%   its skeleton becomes a straight line before computing the shape stats.
%
%   This is very useful if you have bent tubular structures in the image
%   (e.g. blood vessels) and you want to estimate the length and
%   cross-section diameters of the tube. By default, STRAIGHT=false.
%
% STATS = SCINRRD_SEG2LABEL_STATS(..., D, DICT)
%
%   D is the sparse distance matrix for adjacent voxels in the
%   segmentation. DICT is the dictionary vector to convert between image
%   linear indices and distance matrix indices. They can be computed using
%
%     >> [d, dict] = seg2dmat(seg.data, 'seg', [seg.axis.spacing]);
%
%   If they are not provided externally, they are computed internally by
%   the function.
%
% See also: scinrrd_seg2voxel_stats.

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
error(nargchk(2, 5, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% number of skeleton voxels that we are going to use to warp every
% segmentation voxel
K = 3;

% defaults
if (nargin < 3 || isempty(STRAIGHT))
    STRAIGHT=false;
end
if (nargin < 5 || isempty(d) || isempty(dict))
    % compute distance matrix for the whole image segmentation; this is a
    % very sparse matrix because only adjacent voxels get a connection
    [d, dict] = seg2dmat(nrrd.data, 'seg', [nrrd.axis.spacing]);
end

% init matrix to store the eigenvalues of each branch
eigd = zeros(3, cc.NumObjects);

% loop every branch
for I = 1:cc.NumObjects
    
    %% find the 3 closest skeleton points for each branch voxel
    
    % list of voxels in current branch
    br = find(nrrd.data == I);
    
    % list of voxels that are part of the skeleton in the branch
    sk = cc.PixelIdxList{I};
    
    % add skeleton voxels to branch voxels (bifurcation points belong to
    % the branch CC, but they may not be in the branch LAB)
    br = union(br, sk);
    
    % coordinates of branch voxels
    [r, c, s] = ind2sub(size(nrrd.data), br);
    x = scinrrd_index2world([r, c, s], nrrd.axis)';
    
    % create a distance matrix for only the voxels in the branch (this
    % includes the skeleton)
    aux = dict(br);
    dbr = d(aux, aux);
    dictbr = sparse(br, ones(length(br), 1), 1:length(br));
    
    % get the indices for skeleton voxels only in the branch
    idxsk = dictbr(sk);
    
    % compute distances from all branch voxels to skeleton points in the
    % branch
    dbrsk = dijkstra(dbr, idxsk);
    
    % sort the skeleton points so that we can select the closest skeleton
    % points
    [~, idxbrsk] = sort(dbrsk, 1, 'ascend');
    idxbrsk = idxbrsk(1:K, :);

    %% straighten all branch voxels using a local similarity transformation
    
    if (STRAIGHT)
        % create a straightened section of the skeleton of the same length
        % and with the same spacing between voxels
        xskp = [cc.PixelParam{I}' ; zeros(2, length(sk))];
        
        % init variable for warped branch points
        z = zeros(size(xskp, 1), length(br));
        
        % loop every branch voxel
        for J = 1:length(br)
            
            % convert indices of closest skeleton voxels to indices in the
            % branch, and then extract their real world coordinates
            ysk = x(:, idxsk(idxbrsk(:, J)));
            
            % get the coordinates of the corresponding stretched skeleton
            % points
            ysk2 = xskp(:, idxbrsk(:, J));
            
            % compute similarity transformation between the original and
            % straightened skeleton points
            [~, ~, t] = procrustes(ysk2', ysk');
            
            % extract real world coordinates of the voxel we are trying to
            % warp
            x2 = x(:, J);
            
            % interpolate point
            z(:, J) = (t.b * x2' * t.T + t.c(1, :))';
        end
        
    else % STRAIGHT=false
        
        z = x;
        
    end
    
    %% compute statistics about the shape of the branch
    
    % compute eigenvalues of branch
    [~, eigd(:, I)] = pts_pca(z);
    
end

% create output struct
stats.var = eigd;
