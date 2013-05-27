function [d2, idx, p, d] = scimat_dmatrix_imgeodesic(scimat, idx)
% SCIMAT_DMATRIX_IMGEODESIC  Compute a distance/adjacency matrix for a set
% of scattered points using the image intensity to approximate the
% geodesics
%
% D = scimat_dmatrix_imgeodesic(SCIMAT, IDX)
%
%   SCIMAT is SCI MAT struct that contains a grayscale image. The live wire
%   will try to avoid bright (white) areas, and instead select darker
%   voxels.
%
%   IDX is a 3-column matrix with the index-coordinates of the points in
%   the scattered set.
%
%   In particular, the algorithm looks for the shortest paths between pairs
%   of points in IDX, where the Euclidean distance between two adjacent
%   voxels is weighted by their mean intensity. Note that if the intensity
%   difference between background and object is small, the algorithm may
%   prefer a shorter path even if it crosses the object. Also, backgrounds
%   with intensity = 0 can produce unexpected results. Take into account
%   that in that case, any path on the background will have the same
%   cost=0.
%
%   D is a square matrix. D(i,j) is the estimated geodesic distance between
%   points i and j. This distance is purely geometric, not weighted by
%   intensity values.
%
% [D, IDX2, P, DINT] = scimat_dmatrix_imgeodesic(...)
%
%   IDX2 is a vector with the linear indices that correspond to the
%   multiple subscript indices in IDX.
%
%   P is the predecessors matrix. Each row corresponds to a point in IDX.
%   The path from target node IDX2(j) to source node IDX2(i) can be
%   backtracked using
%
%     pth = graphpred2path(p(I, :), idx2(J));
%
%  Note that, for efficiency reasons, the predecessors path is only
%  guaranteed for the nodes in IDX. For other nodes, the path may not have
%  been computed or may be suboptimal.
%
%  DINT is a distance matrix with the intensity-weighted distances from the
%  nodes in IDX to all the voxels in the image. As with P, distances are
%  only guaranteed to be correct and minimal for the nodes in IDX.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.1.0
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
narginchk(2, 2);
nargoutchk(0, 4);

% voxel size
res = [scimat.axis.spacing];

% convert index coordinates to linear index coordinates
idx = sub2ind(size(scimat.data), idx(:, 1), idx(:, 2), idx(:, 3));

% extract a sparse matrix with the distances between adjacent voxels in the
% MRI volume weighted by the mean voxel intensity
ds = im2dmatrix(double(scimat.data), res(1:ndims(scimat.data)));

% create another sparse matrix that contains only the geodesic distances,
% without any intensity weighting
ds2 = im2dmatrix(ones(size(scimat.data)), res(1:ndims(scimat.data)));

% shortest paths between landmarks and every other voxel in the image
% (actually, the algorithm stops when all the target landmarks have been
% found)
%
% Note that Euclidean distances between voxels are weighted by voxel
% intensities, so they are not actual pure geodesic shortest path distances
%
% Piggyback the Euclidean-distance-only sparse matrix so that we have the
% geodesic distances too
[d, p, d2] = dijkstra(ds, idx, idx', ds2);

% we only need to keep the distances between landmarks
d2 = d2(:, idx);
