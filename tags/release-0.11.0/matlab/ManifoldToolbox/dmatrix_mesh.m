function [d, dtot] = dmatrix_mesh(tri, x)
% DMATRIX_MESH  Sparse distance and shortest-path distance matrices between
% the nodes of a mesh
%
% [D, DTOT] = dmatrix_mesh(TRI)
% [D, DTOT] = dmatrix_mesh(TRI, X)
%
%   TRI is a matrix where each row contains the indices of the nodes that
%   form an element in the mesh. Thus, for a triangulation, TRI has 3
%   columns, for a tetrahedral mesh, TRI has 4 columns and so on.
%
%   X is a matrix where each row contains the real world coordinates of a
%   mesh node. Thus, for 2D points, X has 2 columns, for 3D points, 3
%   columns and so on. If X is not provided, then distances between
%   connected vertices are assumed to be 1.
%
%   D is a sparse matrix. D(i,j) is the length of the edge between nodes i
%   and j. If the nodes are not connected, D(i,j) is a "hole" in the
%   matrix. Note that although under Matlab sparse matrix's definition this
%   technically means D(i,j)=0, by convention in function dijkstra(), this
%   means for us D(i,j)=Inf.
%
%   DTOT is a full matrix. DTOT(i,j) is the shortest-path length between
%   nodes i and j, using Dijkstra's algorithm. Note that if node j cannot
%   be reached from node i at all, then DTOT(i,j)=Inf. For nodes connected
%   by edges, D(i,j)=DTOT(i,j).
%
% See also: dmatrix_con, dmatrix_sphmesh, dmatrix.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.2.0
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

% count the number of edges in the mesh
Nedge = size(tri, 1) * (size(tri, 2) - 1);

% rearrange tri so that we have a 2-column matrix with a list of all the
% edges
e = zeros(Nedge, 2);
for I = 1:(size(tri, 2) - 1)
    e(size(tri, 1)*(I-1)+1:size(tri, 1)*I, :) = tri(:, I:I+1);
end
I=I+1;
e(size(tri, 1)*(I-1)+1:size(tri, 1)*I, :) = tri(:, [I 1]);

% sort the vertex indices, so that the small one goes before the larger one
e = sort(e, 2, 'ascend');

% remove duplicated edges
e = unique(e, 'rows');

% compute length of each edge
if (nargin > 1)
    d = x(e(:, 1), :) - x(e(:, 2), :);
    d = sqrt(sum(d.*d, 2));
else
    d = ones(size(e, 1), 1);
end

% create a sparse matrix with the distances (and make it symmetric)
d = sparse([e(:, 1); e(:, 2)], [e(:, 2); e(:, 1)], [d; d]);

% if requested by the user, compute the full distance matrix using
% Dijkstra's shortest path algorithm
if (nargout > 1)
    dtot = dijkstra(d, 1:length(d));
end
