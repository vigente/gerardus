function stats = scinrrd_seg2label_stats(nrrd, cc, p, STRAIGHT)
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
% STATS = SCINRRD_SEG2LABEL_STATS(NRRD, CC)
%
%   NRRD is an SCI NRRD struct with a labelled segmentation mask. All
%   voxels in nrrd.data with value 0 belong to the background. All voxels
%   with value 1 belong to object 1, value 2 corresponds to object 2, and
%   so on.
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
%   STATS is a struct with the shape parameters computed for each object in
%   the segmentation. The measures provided are
%
%     STATS.Var: variance in the three principal components of the cloud
%                of voxels that belong to each object. These are the
%                ordered eigenvalues obtained from computing Principal
%                Component Analysis on the voxel coordinates.
%
%     STATS.IsLandlocked: bool to tell whether the corresponding section is
%                landlocked
%
%     STATS.NBound: number of voxels in the outer boundary of the section
%
%     STATS.NWater: number of voxels in the outer boundary that are
%                touching the background
%
%     STATS.NVox: number of voxels in the branch
%
%     STATS.Vol: volume of the branch (in m^3) units
%
%     STATS.CylDivergence: 75%-quantile of distances between actual voxels
%                on the branch's surface, and the corresponding voxel on
%                the predicted cylinder. This value sometimes cannot be
%                estimated in branches with very few voxels, and a NaN is
%                returned
%
%
%
% STATS = SCINRRD_SEG2LABEL_STATS(..., P, STRAIGHT)
%
%   P is a scalar in [0, 1]. To straighten branches, an approximating or
%   smoothing cubic spline is fitted to the skeleton voxels using
%   csaps(..., P). P=0 is the smoothest spline (a line with the least
%   squares approximation), while P=1 is a rugged spline (the spline
%   interpolated the voxels). Adequate values of P depend on the image
%   resolution, so it's difficult to propose a formula. For resolution in
%   the order of 2.5e-5, P=.999999 seems to give good results (note that
%   for small resolution, P=.999999 gives a very different result to
%   P=1.0). For resolution in the order of 1, P=0.8 seems to give good
%   results. By default, P=1 and no smoothing is performed.
%
%   STRAIGHT is a boolean flag. If STRAIGHT==true, then branches are
%   straightened using the skeleton before computing PCA. By default,
%   STRAIGHT=true.
%
%
% See also: skeleton_label, seg2dmat, scinrrd_seg2voxel_stats.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.8.2
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
error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 3 || isempty(p))
    p = 1.0;
end
if (nargin < 4 || isempty(STRAIGHT))
    STRAIGHT = true;
end

% figure out whether the data is 2D or 3D, because if it's 2D, a landlocked
% voxel has degree 8, but if it's 3D, it needs degree 26
if (size(nrrd.data, 3) > 1)
    % data is 3D
    degmax = 26;
else
    % data is 2D
    degmax = 8;
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
%
% to compute the degree, we want to convolve with
%
% 1 1 1
% 1 0 1
% 1 1 1
%
% Because this array is not linearly separable (rank ~= 1), we use instead
%
% 1 1 1
% 1 1 1 = convn([1 1 1], [1 1 1]')
% 1 1 1
%
% and then have to substract "1" from the degree of each segmented voxel
%
% Also, we have to mask out everything outside the segmentation, because
% the convolution has a "blurring" effect at the edges of the object
deg = convn(single(nrrd.data ~= 0), ones(3, 1, 1, 'single'), 'same');
deg = convn(deg, ones(1, 3, 1, 'single'), 'same');
if (size(nrrd.data, 3) > 1)
    deg = convn(deg, ones(1, 1, 3, 'single'), 'same');
end
idx = deg > 0;
deg(idx) = deg(idx) - 1;
deg = uint8(deg);
deg = deg .* uint8(nrrd.data ~= 0);

% init output
stats.Var = zeros(3, N);
stats.IsLandlocked = true(1, N);
stats.NBound = nan(1, N);
stats.NWater = nan(1, N);
stats.NVox = nan(1, N);
stats.Vol = nan(1, N);
stats.CylDivergence = nan(1, N);

% get coordinates and labels of the segmented voxels
idxlab = find(nrrd.data);
lab = nonzeros(nrrd.data);

% sort the label values
[lab, idx] = sort(lab);
idxlab = idxlab(idx);

% find where each label begins, and add a last index for an inexistent
% label; that index will be used to know where the last label ends
idxlab0 = [0 ; find(diff(lab)) ; length(lab)] + 1;

% free some memory
clear lab

% loop every branch
for I = 1:N

    %% compute number of voxels

    % list of voxels in current branch. The reason why we are not doing a
    % simple br = find(nrrd.data == I); is because for a large volume,
    % that's a comparatively very slow operation
    br = idxlab(idxlab0(I):idxlab0(I+1)-1);
    
    % count number of voxels
    stats.NVox(I) = length(br);
    
    % if there are no voxels, skip this branch
    if (isempty(br))
        continue
    end

    %% compute boundary stats

    % number of voxels that are touching the background
    stats.NWater(I) = nnz(deg(br) ~= degmax);
    
    % if all the voxels have maximum degree, then the label is landlocked
    stats.IsLandlocked(I) = stats.NWater(I) == 0;
    
    % crop the part of the segmentation that contains the branch
    [r, c, s] = ind2sub(size(nrrd.data), br);
    
    from = min([r c s], [], 1);
    to = max([r c s], [], 1);
    
    % keep only voxels of the current branch
    deglab0 = (nrrd.data(from(1):to(1), from(2):to(2), from(3):to(3)) == I);
    
    % compute degree of each voxel in the label if the label had been
    % disconnected from all other labels
    deglab = convn(single(deglab0), ones(3, 1, 1, 'single'), 'same');
    deglab = convn(deglab, ones(1, 3, 1, 'single'), 'same');
    if (size(nrrd.data, 3) > 1)
        deglab = convn(deglab, ones(1, 1, 3, 'single'), 'same');
    end
    idx = deglab0 > 0;
    deglab(idx) = deglab(idx) - 1;
    deglab = uint8(deglab .* deglab0);
    
    % total number of voxels in the outer boundary of the label, whether
    % they touch other labels or not
    stats.NBound(I) = nnz(deglab ~= degmax & deglab ~= 0);

    %% compute eigenvalues using PCA

    if (STRAIGHT)
        % list of voxels that are part of the skeleton in the branch
        sk = cc.PixelIdxList{I};
        
        % add skeleton voxels to the branch, in case they are not already
        br = union(sk, br);
    end
    
    % coordinates of branch voxels
    [r, c, s] = ind2sub(size(nrrd.data), br);
    xi = scinrrd_index2world([r, c, s], nrrd.axis)';
    
    % straighten all branch voxels
    if (STRAIGHT && all(~isnan(cc.PixelParam{I})) ...
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
        yi = itk_pstransform('bspline', x', y', xi', [], 5)';

        % compute eigenvalues of branch (most of the time we are going to
        % get 3 eigenvalues, but not always, e.g. if we have only two
        % voxels in the branch)
        [eigv, stats.Var(:, I)] = pts_pca(yi);
        
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
        stats.Var(:, I) = stats.Var(idx, I);
        eigv = eigv(:, idx);
        
    else % don't straighten objects
        
        yi = xi;
        
        % compute middle point
        y0m = median(yi, 2);
        
        % compute eigenvalues of branch (most of the time we are going to
        % get 3 eigenvalues, but not always, e.g. if we have only two
        % voxels in the branch)
        [eigv, aux] = pts_pca(yi);
        
        % if we don't have 3 distinct eigenvectors/eigenvalues, we can skip
        % this branch
        if (length(aux) < 3)
            continue
        end
        stats.Var(1:length(aux), I) = aux;
        
    end
    
    % numeric errors can cause the appearance of small negative
    % eigenvalues. We make those values 0 to avoid errors below
    stats.Var(stats.Var(:, I) < 0, I) = 0;

    %% convert voxel coordinates to segmentation mask and create cylinder
    %% segmentation mask

    % translate and rotate segmentation voxels so that they are on the
    % X-axis centered around 0
    yi = eigv' * (yi - repmat(y0m, 1, size(yi, 2)));
    
    % compute dimensions of the cylinder
    L = sqrt(12 * stats.Var(1, I));
    r1 = sqrt(4 * stats.Var(2, I));
    r2 = sqrt(4 * stats.Var(3, I));
    
    % polar coordinates of the segmentation voxels
    theta = atan2(yi(3, :), yi(2, :));
    r = sqrt(yi(2, :).^2 + yi(3, :).^2);
    
    % distance from the origin to the elliptical perimeter of the cylinder
    % for each value of the angle
    rel = r1 * r2 ./ sqrt((r2 * cos(theta)).^2 + (r1 * sin(theta)).^2);
    
    % find segmentation voxels that are within the cylinder
    isin = ...
        (yi(1, :) >= -L/2) & (yi(1, :) <= L/2) ...
        & (r <= rel);
    
    % compute overlap between vessel and cylinder
    stats.CylOverlap(I) = nnz(isin) / length(isin);
    
    % compute Delaunay triangulation of the segmentation points
    tri = DelaunayTri(yi');
    
    % sometimes for small branches the points can be colinear, in which
    % case it's not possible to compute a 3D triangulation. In that case,
    % we skip to the next branch
    warning('off', 'MATLAB:TriRep:EmptyTri3DWarnId')
    if (isempty(tri.Triangulation))
        
        warning('on', 'MATLAB:TriRep:EmptyTri3DWarnId')
        stats.CylDivergence(I) = nan;
        continue
        
    end
    warning('on', 'MATLAB:TriRep:EmptyTri3DWarnId')
    
    % get all edges
    e = tri.edges;
    
    % compute length of each edge
    len = sqrt(sum((yi(:, e(:, 1)) - yi(:, e(:, 2))).^2, 1));
    
    % get length of voxel diagonal
    len0 = sqrt(sum([nrrd.axis.spacing].^2));
    
    % find tetrahedra where at least an edge is very long
    badtetra = tri.edgeAttachments(e(len > 1.1 * len0, :));
    badtetra = unique([badtetra{:}]);
    
    % create new triangulation, removing the tetrahedra with long edges
    aux = tri.Triangulation;
    aux(badtetra, :) = [];
    if isempty(aux)
        
        stats.CylDivergence(I) = nan;
        
    else
        
        warning('off', 'MATLAB:TriRep:PtsNotInTriWarnId')
        tri = TriRep(aux, yi(1, :)', yi(2, :)', yi(3, :)');
        warning('on', 'MATLAB:TriRep:PtsNotInTriWarnId')
        
        % find voxels that are on the surface of the segmentation
        triboundary = freeBoundary(tri);
        idx = unique(triboundary(:));
        
        % get voxel actual distance to the central line normalized by the
        % expected distance if the segmentation is a cylinder
        stats.CylDivergence(I) = quantile(abs(r - rel), .75);
        
%         % DEBUG: plot the mesh
%         hold off
%         trisurf(triboundary, yi(1,:), yi(2,:), yi(3,:))
    end
    
end

% compute volume of the branch
stats.Vol = stats.NVox * prod([nrrd.axis.spacing]);
