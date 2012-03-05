function [tri, triboundary] = bwmesh(im, res)
% BWMESH  Tetrahedral volumetric and triangular surface mesh of a binary
% segmentation
%
% [TRI, TRIBOUNDARY] = BWMESH(IM, RES)
%
%   IM is a 3D array with a binary segmentation.
%
%   RES is a 3-vector with the size of the voxel in each coordinate. By
%   default, RES=[1 1 1].
%
%   TRI is a 4-column array where each row contains the indices of 4 nodes
%   forming a tetrahedron. Indices correspond to linear indices of the IM
%   array voxels.
%
%   Example: If you have a SCI NRRD volume
%
%     >> [tri, trisurface] = bwmesh(nrrd.data, [nrrd.axis.spacing]);
%
%   To plot the volumetric mesh you can use (but note that this is very
%   slow even for small meshes)
%
%     >> [r, c, s] = ind2sub(size(nrrd.data), find(nrrd.data));
%     >> x = scinrrd_index2world([r c s], nrrd.axis);
%     >> tetramesh(tri, x(:, 1), x(:, 2), x(:, 3))
%
%   TRIBOUNDARY is a 3-column array. Each row has the indices of 3 points
%   that form a triangle on the surface of the mesh. To plot the surface
%   mesh (usually very fast) you can use
%
%     >> trisurf(triboundary, x(:, 1), x(:, 2), x(:, 3), 'EdgeColor', 'none')
%     >> axis xy equal
%     >> camlight('headlight')
%     >> lighting gouraud
%
%   To obtain a list of indices of the points on the surface,
%
%     >> idx = unique(triboundary(:));
%
%
% See also: pts_mesh, TriRep, DelaunayTri

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.1.3
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.


% check arguments
error(nargchk(1, 2, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% defaults
if (nargin < 2 || isempty(res))
    res = [1 1 1];
end

% auxiliary struct so that we can use a function already in the toolbox
for I = 3:-1:1
    nrrdaxis(I).spacing = res(I);
    nrrdaxis(I).min = 0;
    nrrdaxis(I).size = size(im, I);
end

% real world coordinates of the voxel centres
[r, c, s] = ind2sub(size(im), find(im));
x = scinrrd_index2world([r c s], nrrdaxis);

% length of voxel diagonal
len = sqrt(sum(res.^2));

% compute mesh
if (nargout > 1)
    [tri, triboundary] = pts_mesh(x, 1.75 * len);
else
    tri = pts_mesh(x, 1.75 * len);
end
