function [tri, x] = addpoint2tri(tri, x, xi)
% ADDPOINT2TRI  Add one or more nodes to a triangular mesh
%
% This is a not very sophisticated function that allows to add new points
% to an exisiting triangular mesh. Note that it doesn't deal with the case
% when an inserted point is identical to an existing node, or with the case
% when an inserted point is equally closest to more than one facet.
%
% [TRI2, X2] = ADDPOINT2TRI(TRI, X, XI)
%
%   TRI is a 3-column matrix describing a triangular mesh. TRI(i,:) are the
%   indices of the 3 nodes that form the i-th triangle of the mesh.
%
%   X is a 3-column matrix with the xyz-coordinates of the nodes.
%
%   XI is a 3-column matrix with the xyz-coordinates of the points we want
%   to insert into the mesh. Each point becomes a node. Iteratively, the
%   closest point to the mesh is chosen. The closest facet to that node is
%   removed and replaced by 3 new triangles, all originating from the new
%   node, and ending on the boundary of the previous facet.
%
%   TRI2, X2 describe the new triangular mesh, after all the points in XI
%   have been added.

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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% TODO: deal with the case when an inserted point is identical to an
% existing node

% TODO: deal with the case when an inserted point is equally closest to
% more than one facet

% check arguments
narginchk(3, 3);
nargoutchk(0, 2);

if (size(x, 2) ~= size(xi, 2) || size(x, 2) ~= 3)
    error('X and XI must have 3 columns')
end

% loop until we have added all the points
while(~isempty(xi))
    
    % compute distances from all the points to the mesh
    [f, d] = closest_trifacet(tri, x, xi);
    
    % find closest point to the mesh. That's the point we are going to
    % insert
    [~, xiidx] = min(d);
    
    % facet that will be removed
    fidx = f(xiidx);
    
    % move the new point to the list of mesh nodes
    x(end+1, :) = xi(xiidx, :);
    xi(xiidx, :) = [];

    % for convenience, index of the new node
    newidx = size(x, 1);
    
    % nodes of the triangle that is closest to the point to be inserted
    % tri(f(xiidx), :)
    
    % add 3 new triangles to the mesh. If the new node P is closest to
    % triangle (a,b,c), then we add triangles (P,b,c), (a,P,c), (a,b,P) and
    % delete the original triangle
    tri(end+1, :) = [newidx tri(fidx, 2:3)];
    tri(end+1, :) = [tri(fidx, 1) newidx tri(fidx, 3)];
    tri(end+1, :) = [tri(fidx, 1:2) newidx];
    tri(fidx, :) = [];
    
end
