function [tri, x] = remove_vertex_from_tri(tri, x, idx)
% REMOVE_VERTEX_FROM_TRI  Remove vertices from a triangular mesh
%
% This function removes vertices from a mesh and all triangles they are
% part of. It also relabels vertex indices accordingly.
%
% [TRI2, X2] = remove_vertex_from_tri(TRI, X, IDX)
%
%   TRI is a 3-column matrix. Each row represents the indices of the three
%   vertices that form a triangle. TRI as a whole represents the closed
%   surface.
%
%   X is a 3-column matrix. Each row represents the Cartesian coordinates
%   of a vertex on the surface, indexed by TRI values.
%
%   IDX is a vector of vertex indices.
%
%   TRI2 and X2 are the matrices that describe the mesh after the vertices
%   have been removed.
%
% See also: addpoint2tri.

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

% check arguments
narginchk(3, 3);
nargoutchk(0, 2);

% map between node indices map(10) = 6 means that vertex 10 will be renamed
% as 6
map = 1:size(x, 1);

% remove duplicates
idx = unique(idx);

% number of indices
N = length(idx);

% loop every index to remove
for I = 1:N
    
    % if we remove one vertex, all vertices after it decrease their index
    map(idx(I):end) = map(idx(I):end) - 1;
    
end

% find triangles connected to the vertex that are going to be removed
[I, J] = find(ismember(tri, idx));

% remove triangles
tri(I, :) = [];

% remove vertices
x(idx, :) = [];

% rename vertices in triangles
tri = map(tri);
