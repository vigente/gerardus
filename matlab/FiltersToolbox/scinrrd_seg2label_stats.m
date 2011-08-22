function stats = scinrrd_seg2label_stats(nrrd, cc, p)
% SCINRRD_SEG2LABEL_STATS  Shape stats for each object in a multi-label
% segmentation; objects can be straightened with an skeleton or medial line
% before computing some measures
%
%   This function was developed to measure the dimensions of different
%   objects found in a segmentation. 
%
%   Each voxel in the input segmentation has to be labelled as belonging to
%   the different objects or the background.
%
%   Boundary voxels are counted to assess whether each label is next to the
%   background ("water"), or interior ("landlocked").
%
%   This function then computes Principal Component Analysis (PCA) on the
%   voxels of each object, to estimate the variance in the 3 principal
%   directions, var(1), var(2) and var(3).
%
%   If the object is e.g. a cylindrical vessel with elliptical
%   cross-section, the length and the 2 main diamters of the cylinder can
%   be estimated as
%
%      L = sqrt(12 * var(1))
%      d1 = sqrt(16 * var(2))
%      d2 = sqrt(16 * var(3))
%
%   Because vessels are quite often curved or bent, and this gives
%   misleading values of variance, the function can also straighten the
%   objects using a skeleton or medial line prior to computing PCA.
%
%   If the branches are straightened, note that the returned variance
%   values are not necessarily in largest to smallest order, but:
%
%     var(1): eigenvector that is closest to the straightened skeleton
%     var(2): largest of the remaining eigenvalues
%     var(3): smallest of the remaining eigenvalues
%
%   Boundary voxel counting is performed without straightening the labels.
%
% STATS = SCINRRD_SEG2LABEL_STATS(NRRD)
%
%   NRRD is an SCI NRRD struct with a labelled segmentation mask.
%
%   All voxels in nrrd.data with value 0 belong to the background. All
%   voxels with value 1 belong to object 1, value 2 corresponds to object
%   2, and so on.
%
%   STATS is a struct with the shape parameters computed for each object in
%   the segmentation. The measures provided are
%
%     STATS.var: variance in the three principal components of the cloud
%                of voxels that belong to each object. These are the
%                ordered eigenvalues obtained from computing Principal
%                Component Analysis on the voxel coordinates.
%
%     STATS.islandlocked: bool to tell whether the corresponding section is
%                landlocked
%
%     STATS.nbound: number of voxels in the outer boundary of the section
%
%     STATS.nwater: number of voxels in the outer boundary that are
%                touching the background
%
%     STATS.nvox: number of voxels in the branch
%
%     STATS.vol: volume of the branch (in m^3) units
%
%     STATS.dice: Dice's coefficient (relative volume overlap) between the
%                branch and the cylinder that corresponds to the STATS.var
%                values. Dice's coefficient = 2 * |A & B| / (|A| + |B|)
%
%     STATS.jaccard: Jaccard index (relative volume overlap) between the
%                branch and the cylinder that corresponds to the STATS.var
%                values. Jaccard's coefficient = 2 * |A & B| / |A | B|
%
%   If an object has no voxels, the corresponding STATS values are NaN.
%
% STATS = SCINRRD_SEG2LABEL_STATS(NRRD, CC)
%
%   If CC is provided, then the function straightens each object before
%   computing PCA.
%
%   CC is a struct produced by function skeleton_label() with the list of
%   skeleton voxels that belongs to each object, and the parameterization
%   vector for the skeleton.
%
%   The labels can be created, e.g.
%
%     >> nrrd = seg;
%     >> [nrrd.data, cc] = skeleton_label(sk, seg.data, [seg.axis.spacing]);
%
%     where seg is an NRRD struct with a binary segmentation, and sk is the
%     corresponding skeleton, that can be computed using
%
%     >> sk = itk_imfilter('skel', seg.data);
%
% STATS = SCINRRD_SEG2LABEL_STATS(..., P)
%
%   P is a scalar in [0, 1]. To straighten branches, an approximating or
%   smoothing cubic spline is fit to the skeleton voxels using csaps(...,
%   P). P=0 is the smoothest spline (a line with the least squares
%   approximation), while P=1 is a rugged spline (the spline interpolated
%   the voxels). Adequate values of P depend on the image resolution, so
%   it's difficult to propose a formula. For resolution in the order of
%   2.5e-5, P=.999999 seems to give good results (note that for small
%   resolution, P=.999999 gives a very different result to P=1.0). For
%   resolution in the order of 1, P=0.8 seems to give good results. By
%   default, P=1 and no smotthing is performed.
%
% See also: skeleton_label, seg2dmat, scinrrd_seg2voxel_stats.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.7.0
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
error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 2)
    cc = [];
end
if (nargin < 3 || isempty(p))
    p = 1.0;
end

% figure out whether the data is 2D or 3D, because if it's 2D, a landlocked
% voxel has degree 8, but if it's 3D, it needs degree 26
if (size(nrrd.data, 3) > 1)
    % data is 3D
    degmax = 26;
    degbox = ones(3, 3, 3, 'uint8');
    degbox(2, 2, 2) = 0;
else
    % data is 2D
    degmax = 8;
    degbox = ones(3, 'uint8');
    degbox(2, 2) = 0;
end

if (p < 0 || p > 1)
    error('P must be a scalar in [0, 1]')
end

% number of objects
N = max(nrrd.data(:));
if (~isempty(cc) && (N ~= cc.NumObjects))
    error('If CC is provided, then it must have one element per object in NRRD')
end

% compute degree of each voxel in the segmentation
deg = nrrd.data ~= 0;
deg = uint8(convn(deg, degbox, 'same')) .* uint8(deg);

% init output
stats.var = zeros(3, N);
stats.islandlocked = true(1, N);
stats.nbound = nan(1, N);
stats.nwater = nan(1, N);
stats.nvox = nan(1, N);
stats.vol = nan(1, N);
stats.dice = zeros(1, N);
stats.jaccard = zeros(1, N);

% loop every branch
for I = 1:N

    % list of voxels in current branch
    br = find(nrrd.data == I);
    
    % count number of voxels
    stats.nvox(I) = length(br);
    
    % if there are no voxels, skip this branch
    if (isempty(br))
        continue
    end
    
    %% compute boundary stats
    
    % number of voxels that are touching the background
    stats.nwater(I) = nnz(deg(br) ~= degmax);
    
    % if all the voxels have maximum degree, then the label is landlocked
    stats.islandlocked(I) = stats.nwater(I) == 0;
    
    % crop the part of the segmentation that contains the branch
    [r, c, s] = ind2sub(cc.ImageSize, br);
    
    from = min([r c s]);
    to = max([r c s]);
    
    deglab = nrrd.data(from(1):to(1), from(2):to(2), from(3):to(3));
    
    % keep only voxels of the current branch
    deglab = uint8(deglab) == I;
    
    % compute degree of each voxel in the label if the label had been
    % disconnected from all other labels
    deglab = convn(deglab, degbox, 'same') .* deglab;
    
    % total number of voxels in the outer boundary of the label, whether
    % they touch other labels or not
    stats.nbound(I) = nnz(deglab ~= degmax & deglab ~= 0);
    
    %% compute eigenvalues using PCA
        
    if (~isempty(cc))
        % list of voxels that are part of the skeleton in the branch
        sk = cc.PixelIdxList{I};
        
        % add skeleton voxels to the branch, in case they are not already
        br = union(sk, br);
    end
    
    % coordinates of branch voxels
    [r, c, s] = ind2sub(size(nrrd.data), br);
    xi = scinrrd_index2world([r, c, s], nrrd.axis)';
    
    % straighten all branch voxels
    if (~isempty(cc) && all(~isnan(cc.PixelParam{I})) ...
            && (length(sk) > 2) && (length(br) > 2))
        
        % coordinates of skeleton voxels
        [r, c, s] = ind2sub(size(nrrd.data), sk);
        x = scinrrd_index2world([r, c, s], nrrd.axis)';
        
        % smooth skeleton
        if (p < 1)
            
            % compute spline parameterization for interpolation (Lee's
            % centripetal scheme)
            t = cumsum([0 (sum((x(:, 2:end) - x(:, 1:end-1)).^2, 1)).^.25]);
            
            % compute cubic smoothing spline
            pp = csaps(t, x, p);
            
            % sample spline
            x = ppval(pp, t);
            
            % recompute skeleton parameterisation (chord length)
            cc.PixelParam{I} = ...
                cumsum([0 sqrt(sum((x(:, 2:end) - x(:, 1:end-1)).^2, 1))])';
            
        end
    
        % create a straightened section of the skeleton of the same length
        % and with the same spacing between voxels
        y0 = [cc.PixelParam{I}' ; zeros(2, length(sk))];

        % middle point in the parameterisation
        y0m = y0(:, end) / 2;
        
        % compute rigid transformation to align straight line with skeleton
        [~, y, t] = procrustes(x', y0', 'Scaling', false);
        y = y';
        y0m = (y0m' * t.T + t.c(1, :))';

        % straighten vessel using B-spline transform
        yi = itk_pstransform('bspline', x', y', xi', [], 6)';

        % compute eigenvalues of branch (most of the time we are going to
        % get 3 eigenvalues, but not always, e.g. if we have only two
        % voxels in the branch)
        [eigv, stats.var(:, I)] = pts_pca(yi);
        
        % find the eigenvector that is aligned with the straightened
        % skeleton, that's going to be our "eigenvalue 1", whether it's the
        % largest one or not. The reason is that we are going to always
        % assume that "eigenvalue 1" can be used to estimate the length of
        % the cylinder.
        yv = y(:, end) - y(:, 1);
        yv = yv / norm(yv);
        [~, idx] = max(abs(dot(eigv, ...
            repmat(yv, 1, size(eigv, 2)), 1)));
        
        % create index vector to reorder the eigenvalues and eigenvectors
        idx = [idx 1:idx-1 idx+1:3];
        stats.var(:, I) = stats.var(idx, I);
        eigv = eigv(:, idx);
        
    else % don't straighten objects
        
        yi = xi;
        
        % compute middle point
        y0m = median(yi, 2);
        
        % compute eigenvalues of branch (most of the time we are going to
        % get 3 eigenvalues, but not always, e.g. if we have only two
        % voxels in the branch)
        [eigv, stats.var(:, I)] = pts_pca(yi);
        
    end

    %% convert voxel coordinates to segmentation mask and create cylinder 
    %% segmentation mask
    
    % translate and rotate segmentation voxels so that they are on the
    % X-axis centered around 0
    yi = eigv' * (yi - repmat(y0m, 1, size(yi, 2)));
    
    % compute dimensions of the cylinder
    L = sqrt(12 * stats.var(1, I));
    r1 = sqrt(4 * stats.var(2, I));
    r2 = sqrt(4 * stats.var(3, I));
    
    % vertices of the box that contains the cylinder
    cylbox = [-L/2 -r1 -r2; L/2 r1 r2]';
  
    % real world coordinates => voxel indices
    idx = scinrrd_world2index(yi', nrrd.axis, false)';
    idx0m = scinrrd_world2index([0 0 0], nrrd.axis, false)';
    cylbox = scinrrd_world2index(cylbox', nrrd.axis, false)';
    
    % box that contains both segmentation and cylinder
    box = [min([idx cylbox], [], 2) max([idx cylbox], [], 2)];
    
    % length of boxes
    % Note: boxsz is a length in voxel units and number of voxels
    boxlen = box(:, 2) - box(:, 1) + 1;
    cylboxlen = cylbox(:, 2) - cylbox(:, 1) + 1;
    
    % translate index coordinates so that they are centered around (0,0,0)
    idx = idx - repmat(idx0m, 1, size(idx, 2));
    cylbox = cylbox - repmat(idx0m, 1, size(cylbox, 2));
    box = box - repmat(idx0m, 1, size(box, 2));
    
%     % DEBUG: plot branch and boxes' corners
%     hold off
%     plot3(idx(2, :), idx(1, :), idx(3, :), '.')
%     hold on
%     plot3(cylbox(2, :), cylbox(1, :), cylbox(3, :), 'ok')
%     plot3(box(2, :), box(1, :), box(3, :), 'xg')
%     axis equal xy
    
    % compute regular grid at integer positions over vessel cross section;
    % this is where the ellipse will go
    [gr, ~, gs] = ndgrid(floor(box(1, 1)):ceil(box(1, 2)), 0, ...
        floor(box(3, 1)):ceil(box(3, 2)));

    % find which points in the grid belong inside the ellipse
    
    % polar coordinates of the grid points
    theta = atan2(gs, gr);
    r = sqrt(gs.^2 + gr.^2);
    
    % distance from the origin to the ellipse along the line that connects
    % the origin with each grid point
    rel = cylboxlen(1) * cylboxlen(3) / 4 ...
        ./ sqrt((cylboxlen(3) / 2 * cos(theta)).^2 ...
        + (cylboxlen(1) / 2 * sin(theta)).^2);
    
    % points inside the ellipse
    % Note: to see the ellipse the right way
    %   >> imagesc(squeeze(imin)');
    %   >> axis equal ij
    % because in imin, rows => y-axis, columns => z-axis
    imin = uint8(r <= rel);
    
%     % DEBUG: plot cylinder voxels
%     for J = round(cylbox(2, 1)):round(cylbox(2, 2))
%         plot3(J*ones(numel(find(r<=rel)), 1),  ...
%             gr(r <= rel), gs(r <= rel), 'r.')
%     end
    
    % number of voxels for the cylinder
    cylboxn = round(cylboxlen) + 1;
    
    % number of voxels for the long axis of whole box
    boxn = ceil(boxlen) + 1;
    
    % if number of left over voxels is not even, then we make the box 1
    % voxel larger
    boxn = boxn + mod(boxn - cylboxn, 2);
    
    % cylinder mask is going to be a padding of zeros, a repetition of the
    % ellipse, and another padding of zeros
    mask = [ ...
        zeros(size(imin, 1), (boxn(2) - cylboxn(2))/2, size(imin, 3), 'uint8') ...
        repmat(imin, [1, cylboxn(2), 1]) ...
        zeros(size(imin, 1), (boxn(2) - cylboxn(2))/2, size(imin, 3), 'uint8') ...
        ];
    
    % translate segmentation to (1, 1, 1), so that we can use the
    % coordinates as voxel indices of the _whole_ box (not the cylinder
    % box)
    idx = idx - repmat(round(min(idx, [], 2)) - 1, 1, size(idx, 2));

    % convert the coordinates of vessel voxels to a segmentation mask
    im = zeros(size(mask), 'uint8');
    im(sub2ind(size(im), round(idx(1, :)), ...
        round(idx(2, :)), round(idx(3, :)))) = 1;

    % fill holes
    se = strel('ball', 3, 3);
    im = imerode(imdilate(im, se), se);
    
%     % DEBUG: save the segmentation masks as NRRD files so that they can be
%     % inspected with Seg3D
%     foonrrd = scinrrd_im2nrrd(mask);
%     scinrrd_save('/tmp/foomask.mat', foonrrd);
%     foonrrd = scinrrd_im2nrrd(im);
%     scinrrd_save('/tmp/fooim.mat', foonrrd);
    
    % compute overlap between vessel and cylinder
    stats.dice(I) = 2 * nnz(mask & im) / (nnz(mask) + nnz(im));
    stats.jaccard(I) = nnz(mask & im) / nnz(mask | im);
    
end

% compute volume of the branch
stats.vol = stats.nvox * prod([nrrd.axis.spacing]);
