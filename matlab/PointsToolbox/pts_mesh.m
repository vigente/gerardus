function [tri, triboundary] = pts_mesh(x, maxlen)
% PTS_MESH  Tetrahedral volumetric and triangular surface mesh of a cloud
% of points
%
% [TRI, TRIBOUNDARY] = PTS_MESH(X, MAXLEN)
%
%   X is a 3-column array where each row has the Cartesian coordinates of a
%   point from the cloud of points. If you have a segmentation in SCI NRRD
%   format, you can compute the coordinates of the segmented voxels as
%
%     >> [r, c, s] = ind2sub(size(nrrd.data), find(nrrd.data));
%     >> x = scinrrd_index2world([r c s], nrrd.axis);
%
%   MAXLEN is a scalar with the maximum length allowed for an edge
%   connecting two nodes in the mesh. By default, MAXLEN=Inf, so the mesh
%   will be the convex hull of the cloud of points. This is usually *not*
%   the wanted result. A sensible value of MAXLEN is e.g. 1.75 * L, where L
%   is the length of the voxel diagonal.
%
%   TRI is a 'TriRep' struct with the description of the mesh. To plot the
%   volumetric mesh you can use (but note that this is very slow even for
%   small meshes)
%
%     >> tetramesh(tri)
%
%   TRIBOUNDARY is a 3-row array. Each column has the indices of 3 points
%   that form a triangle on the surface of the mesh. To plot the surface
%   mesh (usually very fast) you can use
%
%     >> trisurf(triboundary, x(:, 1), x(:, 2), x(:, 3))
%
%   To obtain a list of indices of the points on the surface,
%
%     >> idx = unique(triboundary(:));
%
%
% See also: bwmesh, TriRep, DelaunayTri

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
error(nargchk(1, 2, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

% defaults
if (nargin < 2 || isempty(maxlen))
    maxlen = Inf;
end

% compute Delaunay triangulation of the segmentation points
tri = DelaunayTri(x);

% sometimes the points can be colinear, in which case it's not possible to
% compute a 3D triangulation. In that case, we return an empty
% triangulation
warning('off', 'MATLAB:TriRep:EmptyTri3DWarnId')
if (isempty(tri.Triangulation))
    
    warning('on', 'MATLAB:TriRep:EmptyTri3DWarnId')
    tri = [];
    triboundary = [];
    return
    
end
warning('on', 'MATLAB:TriRep:EmptyTri3DWarnId')

% get all edges
e = tri.edges;

% compute length of each edge
len = sqrt(sum((x(e(:, 1), :) - x(e(:, 2), :)).^2, 2));

% find tetrahedra where at least an edge is very long
badtetra = tri.edgeAttachments(e(len > maxlen, :));
badtetra = unique([badtetra{:}]);

% create new triangulation, removing the tetrahedra with long edges
aux = tri.Triangulation;
aux(badtetra, :) = [];
if isempty(aux)
    
    tri = [];
    triboundary = [];
    return
    
else
    
    warning('off', 'MATLAB:TriRep:PtsNotInTriWarnId')
    tri = TriRep(aux, x(:, 1), x(:, 2), x(:, 3));
    warning('on', 'MATLAB:TriRep:PtsNotInTriWarnId')
    
    % find voxels that are on the surface of the segmentation
    triboundary = freeBoundary(tri);
    
%     % DEBUG: plot the mesh
%     hold off
%     trisurf(triboundary, x(:, 1), x(:, 2), x(:, 3))
end
