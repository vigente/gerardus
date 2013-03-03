function [d, points] = scimat_dmatrix_thickslice(scimat, K)
% SCIMAT_DMATRIX_THICKSLICE  Compute a distance/adjacency matrix for a
% segmentation that consists of scattered points in slices wide apart
%
% [D, POINTS] = scimat_dmatrix_thickslice(SCIMAT, K)
%
%   SCIMAT is a structure with the segmentation. The segmentation is
%   assumed to represent scattered points, in slices that are wide apart
%   (because the voxels are very thick in the z-axis, or because only a few
%   slices of the whole volume have been segmented).
%
%   The adjacency matrix is built in two steps:
%
%   First, slices are done independently. Each point in the slice is
%   connected only to its K nearest neighbours.
%
%   Second, each slice is connected to the slices immediately above and
%   below it. Each point is again connected to its K nearest neighbours
%   above, and K nearest neighbours below.
%
%   K is a scalar with the number of nearest neighbours that are
%   considered. By default, K=3.
%
%   D is a matrix PxP, assuming there are P points in SCIMAT. D(i,j) is the
%   Euclidean distance between points i and j only if they are connected.
%   Otherwise, D(i,j)=0.
%
%   POINTS is a Px3 matrix with the real world (x, y, z)-coordinates of the
%   segmented points.
%
% Example:
%
%    scimat = scinrrd_load('test/thick-slice-points.mat');
%    [d, points] = scimat_dmatrix_thickslice(scimat, 2);
%    gplot3d(d, points, 'o')
%    hold on
%    gplot3d(d, points)
%    axis equal
%
% See also: thickslice_collate_sax_la

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
narginchk(1, 2);
nargoutchk(0, 2);

% defaults
if (nargin < 2 || isempty(K))
    K = 3;
end

% number of slices
N = size(scimat.data, 3);

% indices of user selected points
idx = find(scimat.data);

% in real world coordinates
[r, c, s] = ind2sub(size(scimat.data), idx);
points = scinrrd_index2world([r c s], scimat.axis);

% total number of points
Ntot = length(r);

% distance matrix between all pairs of points
d = zeros(Ntot);

% split points by slice
idxs = cell(1, N);
for I = 1:N
    idxs{I} = find(s==I);
end

% remove slices without any points. This way, we can have points in only a
% few slices of a high resolution image
idxs = idxs(~cellfun(@isempty, idxs));
N = length(idxs);

% find nearest neighbours within each slice
for I = 1:N
    
    % points in this slice
    Ns = length(idxs{I});
    
    % skip slices without points
    if (Ns == 0)
        continue
    end
    
    % distance matrix for points within this slice
    ds = dmatrix(points(idxs{I}, :)');
    
    % compute local neighbourhood for the slice and update the big distance
    % matrix with it
    d = update_dmatrix(d, ds, idxs{I}, idxs{I}, min([K, Ns-1]));
    
end

% find neighbours between slices
for I = 1:N-1
    
    % points in each slice
    Ns1 = length(idxs{I});
    Ns2 = length(idxs{I+1});
    
    % skip if either slice has no points
    if (Ns1 == 0 || Ns2 == 0)
        continue
    end
    
    % distance matrix for points from one slice to the other
    ds = dmatrix(points(idxs{I}, :)', points(idxs{I+1}, :)');
    
    % compute local neighbourhood for slice I to slice I+1 and update the
    % big distance matrix with it
    d = update_dmatrix(d, ds, idxs{I}, idxs{I+1}, min([K, Ns1-1, Ns2-1]));
    
    % compute local neighbourhood for slice I+1 to slice I and update the
    % big distance matrix with it
    d = update_dmatrix(d, ds', idxs{I+1}, idxs{I}, min([K, Ns1-1, Ns2-1]));
    
end

end

% this function also in thickslice_collate_sax_la()
function d = update_dmatrix(d, ds, idx1, idx2, K)

if (size(ds, 1) ~= length(idx1) || size(ds, 2) ~= length(idx2))
    error('Internal error: Slice distance matrix has incongruent dimensions with the index vectors')
end

% sort distances from smaller to larger in each row
[ds, idxto] = sort(ds, 2, 'ascend');

% the first column will be distance 0 of each point to itself
ds(:, 1) = [];
idxto(:, 1) = [];

% transfer only the distances and connections of the K-nearest
% neighbours to the global distance matrix (making sure that the output
% distance matrix is symmetric)
idxto = idxto(:, 1:K);
idxto = idx2(idxto);
idxfrom = repmat(idx1, 1, K);
d(sub2ind(size(d), idxfrom, idxto)) = ds(:, 1:K);
d(sub2ind(size(d), idxto, idxfrom)) = ds(:, 1:K);

end
