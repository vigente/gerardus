function [tri, triboundary] = bwmesh(im, res)
% BWMESH  Tetrahedral volumetric and triangular surface mesh of a binary
% segmentation
%
% This is a very basic function for toy examples. For a much more
% sophisticated surface mesher, see cgal_meshseg, also in Gerardus.
%
% [TRI, TRIBOUNDARY] = BWMESH(IM, RES)
%
%   IM is a 3D array with a binary segmentation.
%
%   RES is a 3-vector with the size of the voxel in each coordinate. By
%   default, RES=[1 1 1].
%
%   TRI is a 4-column array where each row contains the indices of 4 nodes
%   forming a tetrahedron. E.g. TRI(2, :)=[23 178 3 90] means that the
%   first tetrahedron contains nodes 23, 178, 3 and 90. Nodes are indexed
%   in the order given by find(IM). For instance, the segmentation
%
%     0 1 1 0
%     1 0 0 0
%     0 1 0 1
%
%   has 5 nodes, with corresponding indices 1, 2, 3, 4 and 5.
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
%   TRIBOUNDARY is a 3-column array. Each row has the indices of 3 nodes
%   that form a triangle on the surface of the mesh. To plot the surface
%   mesh (usually very fast) you can use
%
%     >> trisurf(triboundary, x(:, 1), x(:, 2), x(:, 3), 'EdgeColor', 'none')
%     >> axis xy equal
%     >> camlight('headlight')
%     >> lighting gouraud
%
%   To obtain a list of all nodes on the surface,
%
%     >> idx = unique(triboundary(:));
%
%
% See also: cgal_meshseg, pts_mesh, TriRep, delaunay.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011-2013 University of Oxford
% Version: 0.1.6
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
narginchk(1, 2);
nargoutchk(0, 2);

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
