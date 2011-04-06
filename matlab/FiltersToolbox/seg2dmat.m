function d = seg2dmat(im)
% SEG2DMAT  Local neighbourhood distance matrix between segmentation voxels
%
% D = SEG2DMAT(IM)
%
%   D is a sparse matrix where D(i,j) gives the distance between the i-th
%   and j-th voxels in the binary segmentation IM.
%
%   The index values i, j are computed with SUB2IND(), in the usual Matlab
%   way if you reshape IM into a vector, IM(:).
%
%   A 26-neighbourhood is assumed. That is, a voxel is only connected to
%   the 26 voxels that form a cube around it. That's why the sparse matrix
%   representation is convenient.
%
%   This function is fully vectorized and is thus very fast.
%
%   This function can be used with 2D images instead of 3D volumes too,
%   although some of the intermediate steps are not as efficient
%   memory-wise as they could be.
%
% See also: im2imat.

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
error(nargoutchk(0, 1, nargout, 'struct'));

% total number of voxels in the image
N = numel(im);

% image size
sz = size(im);
if (length(sz) == 2)
    sz(3) = 1;
end

% get linear indices of segmented voxels
idx = find(im);

% convert to r, c, s indices
[r, c, s] = ind2sub(sz, idx);

% compute neighbourhood with connectivity 26 around origin
[gr, gc, gs] = ndgrid(-1:1, -1:1, -1:1);

% compute distances from each point to the origin
dlocal = sqrt(gr.^2 + gc.^2 + gs.^2);

% convert volume of distances and coordinates into vectors
dlocal = dlocal(:);
gr = gr(:);
gc = gc(:);
gs = gs(:);

% each row has the coordinates of a neighbourhood around a segmented voxel
r = repmat(gr', length(idx), 1) + repmat(r, 1, 27);
c = repmat(gc', length(idx), 1) + repmat(c, 1, 27);
s = repmat(gs', length(idx), 1) + repmat(s, 1, 27);

% convert out of range subscripts into NaNs, so that they are not
% considered later when computing distances
r(r < 1 | r > sz(1)) = nan;
c(c < 1 | c > sz(2)) = nan;
s(s < 1 | s > sz(3)) = nan;

% reduce the r, c, s matrices to a single matrix with linear indices for
% the neighbours. Now, we have vector idx and matrix nn. For example:
% idx(5) = 34864
% nn(5,:) = [NaN NaN NaN NaN NaN NaN NaN NaN NaN 34395 34396 34397 ...
%            34863 34864 34865 35331 35332 35333 108807 108808 ...
%            108809 109275 109276 109277 109743 109744 109745]
% means that voxel 34864 has only 18 neighbours (voxel 34864 is in the
% first slice, so the neighbourhood is constrained that way)
nn = sub2ind(sz, r, c, s);

% we are going to remove the central column, because we don't need to
% connect a voxel to itself with distance 0
nn = nn(:, [1:13 15:end]);
dlocal = dlocal([1:13 15:end]);

% replicate local distance and index vectors so that it corresponds with
% the nn matrix
dlocal = repmat(dlocal', length(idx), 1);
idx = repmat(idx, 1, 26);

% find the not out of range connections
ok = ~isnan(nn);

% remove connections to out of range voxels
idx = idx(ok);
nn = nn(ok);
dlocal = dlocal(ok);

% create sparse matrix for distances between all voxels in the image
d = sparse(idx, nn, dlocal, N, N);

% remove the connections to voxels that are not segmented
notidx = find(~im);
d(notidx, :) = 0;
d(:, notidx) = 0;
