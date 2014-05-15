function [stats, idx] = scinrrd_seg2voxel_stats(nrrd, RAD, idx)
% SCINRRD_SEG2VOXEL_STATS  Shape stats for each voxel in a segmentation
% based on a windowed neighbourhood.
%
% [STATS, IDX] = scinrrd_seg2voxel_stats(NRRD, RAD, IDX)
%
%   NRRD is an SCI NRRD struct with a binary segmentation mask of different
%   structures, e.g. blobs, tubes, etc.
%
%   This function computes certain parameters (STATS) that can be used to
%   decide whether a voxel belongs to a blob, a tube, etc.
%
%   RAD is a scalar with the window size (in real world coordinates) used
%   to define local neighbourhoods for computations. The local
%   neighbourhood is a cube of side RAD.
%
%   STATS is a struct with the computed parameters for voxel IDX. STATS is
%   not computed for all the voxels in the neighbourhood, but for the
%   connected component that includes the target voxel in that
%   neighbourhood.
%
%   IDX is a vector with linear index values of the segmented voxels. If
%   you want to read the voxels, you can do
%
%   >> nrrd.data(idx)
%
%   The row, column and slice coordinates of the segmented voxels are given
%   by
%
%   >> [r, c, s] = ind2sub(size(nrrd.data), idx);
%
%   By default, IDX is computed internally as all the voxels different from
%   0 in the segmentation mask. But a much faster way to segment vessels is
%   to first skeletonize the segmentation mask (using the C++ program
%   "skeletonize3DSegmentation" in Gerardus), and only compute STATS on the
%   skeleton.
%
%   >> nrrd = scimat_load('im.mat');
%   >> nrrdsk = scimat_load('imsk.mat');
%   >> idxsk = find(nrrdsk.data);
%   >> stats = scinrrd_seg2voxel_stats(nrrd, 25e-4, idxsk);
%
%   The measures provided by STATS are:
%
%     STATS.VAR1: Variance of the connected component in the direction of
%                 maximum variability (i.e. largest eigenvalue). If this
%                 corresponds to the length L of a cylinder, then
%                 L = sqrt(12 * VAR1).
%
%     STATS.VAR2: Second eigenvalue. If this corresponds to the the
%                 elliptical cross-section of a cylinder with major radius
%                 R, then
%                 R = 2 * sqrt(VAR2).
%
%     STATS.VAR3: Third eigenvalue. Likewise for minor radius r, then
%                 r = 2 * sqrt(VAR3).
%
%     STATS.VOL:  Volume of the connected component.
%
%     STATS.DC:   Distance between the centroid of the connected component
%                 and the target voxel.
%
% See also: scinrrd_seg2label_stats.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011-2014 University of Oxford
% Version: 0.1.2
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
narginchk(2, 3);
nargoutchk(0, 2);

% get radius size in voxels in every dimension
x0 = scinrrd_index2world([1 1 1], nrrd.axis);
irad = round(scinrrd_world2index(x0+RAD*ones(1,3), nrrd.axis) - [1 1 1]);

% find all the voxels that are part of the segmentation, if they are not
% provided by the user
if (nargin < 3 || isempty(idx))
    idx = find(nrrd.data);
end

%initialize output matrices for the eigenvalues
stats.var1 = zeros(size(idx));
stats.var2 = stats.var1;
stats.var3 = stats.var1;
stats.vol = stats.var1;
stats.dc = stats.var1;

% compute voxel volume
vol0 = prod([nrrd.axis.spacing]);

%loop every voxel that is part of the segmentation (skip the rest)
for I = 1:length(idx)
    [R, C, S] = ind2sub(size(nrrd.data), idx(I));

    % extract a neighbourhood around the target voxel
    Rmin = max(R-irad(1), 1);
    Rmax = min(R+irad(1), size(nrrd.data, 1));
    Cmin = max(C-irad(2), 1);
    Cmax = min(C+irad(2), size(nrrd.data, 2));
    Smin = max(S-irad(3), 1);
    Smax = min(S+irad(3), size(nrrd.data, 3));
    im = nrrd.data(...
        Rmin:Rmax, ...
        Cmin:Cmax, ...
        Smin:Smax ...
        );
    
    % compute connected components
    cc = bwconncomp(im);
    
    % find which connected component the target voxel belongs to
    idxtarg = sub2ind(size(im), R - Rmin + 1, C - Cmin + 1, S - Smin + 1);
    
    Jtarg = [];
    for J = 1:cc.NumObjects
        if (any(find(cc.PixelIdxList{J} == idxtarg)))
            Jtarg = J;
            break
        end
    end
    if isempty(Jtarg)
        error(['Assertion error: Target point [' num2str(R) ', ' num2str(C) ', ' ...
            num2str(S) ...
            '] is not in any connected component']);
    end
    
    % convert the index values, that are indices of im (neighbourhood), to
    % indices in nrrd.data (whole image)
    [Rcc, Ccc, Scc] = ind2sub(size(im), cc.PixelIdxList{Jtarg});
    Rcc = Rcc+Rmin-1;
    Ccc = Ccc+Cmin-1;
    Scc = Scc+Smin-1;
    
    % convert indices to real world coordinates
    x = scinrrd_index2world([Rcc, Ccc, Scc], nrrd.axis);
    
    % compute centroid of the points
    xmean = mean(x);
    
    % compute coordinates of the target point
    xtarg = scinrrd_index2world([R, C, S], nrrd.axis);
    
    % compute distance between target point and centroid
    stats.dc(I) = sqrt(sum((xmean - xtarg).^2));
    
    % compute PCA on the points, to find the directions of maximum
    % variability
    [~, d] = pts_pca(x');
    
    % the 2nd and 3rd eigenvectors give as the variance values in the two
    % main axis of the orthogonal plane. It can be shown that for an
    % ellipse, each radius of the ellipse r = sqrt(var)*2
    %
    % the 1st eigenvector gives the variance along the vessel. Again, it
    % can be shown that the length L = sqrt(12*var)
    if (length(d) < 1)
        stats.var1(I) = 0;
    else
        stats.var1(I) = d(1);
    end
    if (length(d) < 2)
        stats.var2(I) = 0;
    else
        stats.var2(I) = d(2);
    end
    if (length(d) < 3)
        stats.var3(I) = 0;
    else
        stats.var3(I) = d(3);
    end

    % compute actual volume of connected component
    stats.vol(I) = length(cc.PixelIdxList{Jtarg}) * vol0;
    
end
