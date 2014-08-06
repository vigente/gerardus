function scimat = scimat_surface2seg(scimat, tri, x)
% scimat_surface2seg  Convert triangular mesh into segmentation.
%
%   scimat_surface2seg sets to "1" those voxels in a segmentation volume
%   that contain or are intersected by the triangles in a triangular mesh.
%   The method works by sampling each triangle in the mesh uniformly, using
%   barycentric coordinates. Samples are close enough to each other so that
%   each voxel will have at least two samples. Then sample coordinates are
%   rounded to the closest voxel centres, and the voxels are set.
%
% SCIMAT2 = scimat_surface2seg(SCIMAT, TRI, X)
%
%   SCIMAT is a struct to hold the segmentation (see "help scimat" for
%   details). SCIMAT.axis must be provided. SCIMAT.axis is important
%   because it tells us where the segmentation volume begings and ends, and
%   the voxel size. SCIMAT.data is not necessary, and if it is not
%   provided, a new one will be created. If it is provided, previously set
%   voxels will not be cleared.
%
%   (TRI, X) describe the triangular mesh.
% 
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh.
%
%   X is a 2-column matrix. X(i, :) contains the xy-coordinates of the
%   i-th node in the mesh.
%
%   SCIMAT2 is the output segmentation.
%
% See also: surface_tridomain, surface_param, surface_interp.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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
narginchk(3, 3);
nargoutchk(0, 1);

% smallest dimension of the voxel size
lmin = min([scimat.axis.spacing]);

% if the scimat struct has no data field, create one
if (~isfield(scimat, 'data'))
    scimat.data = zeros(size([scimat.axis.size]), 'uint8');
end

% loop triangles
Ntri = size(tri, 1);
for I = 1:Ntri
    
    % vertices of current triangle
    v = x(tri(I, :), :);
    
    % sample triangle
    xi = sample_triangle(v, lmin);
    
    % convert real world coordinates to indices
    idx = round(scimat_world2index(xi, scimat));
    
    % linearize indices
    idx = sub2ind(size(scimat.data), idx(:, 1), idx(:, 2), idx(:, 3));
    
    % set voxels that belong to the triangle to 1
    scimat.data(idx) = 1;
    
end

end

% given a triangle, sample it uniformly using barycentric coordinates
function xi = sample_triangle(x, lmin)

% nomenclature
v1 = x(1, :);
v2 = x(2, :);
v3 = x(3, :);

% compute number of points between triangle vertices. We want to have at
% least 2 samples within each voxel
dmax = max([norm(v1-v2) norm(v2-v3) norm(v1-v3)]);
N = ceil(1 + 2 * dmax / lmin);

% sample the triangle
a1 = linspace(0, 1, N);
a2 = linspace(0, 1, N);

% barycentric coordinate variables
[a1, a2] = meshgrid(a1, a2);
a1 = a1(:);
a2 = a2(:);
idx = a1 + a2 <= 1;
a1 = a1(idx);
a2 = a2(idx);
a3 = 1 - a1 - a2;

% point coordinates
xi = [a1 a1 a1] .* repmat(v1, length(a1), 1) ...
    + [a2 a2 a2] .* repmat(v2, length(a2), 1) ...
    + [a3 a3 a3] .* repmat(v3, length(a3), 1);

% % plot sampled triangle
% hold off
% plot3(xi(:, 1), xi(:, 2), xi(:, 3), '.')
% hold on
% plot3(x([1:end 1], 1), x([1:end 1], 2), x([1:end 1], 3))

end
