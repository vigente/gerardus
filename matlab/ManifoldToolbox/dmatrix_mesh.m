function [d, dtot] = dmatrix_mesh(tri, x, totmethod, options)
% DMATRIX_MESH  Sparse distance and shortest-path distance matrices between
% the nodes of a mesh
%
% D = dmatrix_mesh(TRI)
% D = dmatrix_mesh(TRI, X)
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
% [D, DTOT] = dmatrix_mesh(TRI, X, 'dijkstra')
% [D, DTOT] = dmatrix_mesh(TRI, X, 'fastmarching', OPTIONS)
%
%   DTOT is a full matrix. DTOT(i,j) is the shortest-path length between
%   nodes i and j. Note that if node j cannot be reached from node i at
%   all, then DTOT(i,j)=Inf. For nodes directly connected by edges,
%   D(i,j)=DTOT(i,j).
%
%   'dijkstra': Dijkstra's shortest-path algorithm is used. 
%
%   'fastmarching': (only valid for triangular meshes), distances are
%   computed as the solution to an Eikonal equation by propagating a front
%   with uniform speed using the Fast Marching method. OPTIONS is a struct
%   with parameters for the Fast Marching method. From the help of
%   perform_fast_marching_mesh():
%
%       OPTIONS.end_points:  stop when these points are reached
%       OPTIONS.nb_iter_max: stop when a given number of iterations is
%       reached.
%       OPTIONS.constraint_map: A column vector with size(X,1) elements to
%       reduce the set of explored points. Only points with current
%       distance smaller than L will be expanded. Set some entries of L to
%       -Inf to avoid any exploration of these points.
%
%   Both Dijkstra and Fast Marching use similar Fibonacci heap
%   implementations, and have the same computational complexity. That said,
%   the Fast Marching implementation we have needs to loop for each vertex,
%   and is 1.5 times slower (7.2 min vs. 4.7 min) in a mesh with 5381 nodes
%   and 10758 elements. However, Dijkstra's method suffers greatly from
%   metrication errors in regular meshes.
%
%   For Dijkstra's algorithm, we use our own modification of Joshua
%   Tenenbaum’s implementation (http://isomap.stanford.edu/), that uses
%   Fibonacci heaps implemented in C by John Boyer. For the Fast Marching
%   method, we use Gabriel Peyre’s implementation of Kimmel and Sethian
%   (1998), that uses Fibonacci heaps implemented by John-Mark Gurney.
%
% R. Kimmel and J. A. Sethian, “Computing geodesic paths on manifolds,”
% Proceedings of the National Academy of Sciences, vol. 95, no. 15, pp.
% 8431–8435, 1998.
%
% See also: dmatrix_con, dmatrix_sphmesh, dmatrix, dijkstra,
% perform_fast_marching_mesh.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.4.3
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
narginchk(1, 4);
nargoutchk(0, 2);

% defaults
if (nargin < 3 || isempty(totmethod))
    totmethod = 'dijkstra';
end
if (nargin < 4)
    options = [];
end

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

% if requested by the user, compute the full distance matrix
if (nargout > 1)
    switch totmethod
        case 'dijkstra'
            
            narginchk(1, 3); % this syntax does not accept options input variable
            dtot = dijkstra(d, 1:length(d));
            
        case 'fastmarching'
            
            if (~isfield(options, 'constraint_map') ...
                    || isempty(options.constraint_map))
                % if the user doesn't limit the neighbourhood, we assume a
                % full distance matrix
                dtot = zeros(size(d));
                options.constraint_map = [];
            else
                % if the user limits the size of the neighbourhood, we are
                % going to assume that it's more efficient to store the
                % distance matrix in sparse form. This of course will not
                % be true if the neighbourhood is quite close to the size
                % of the mesh
                dtot = sparse(d);
            end
            
            % number of vertices
            N = length(d);
            
            % loop the columns of the distance matrix
            for I = 1:N
                aux = perform_fast_marching_mesh(x, tri, I, options);
                auxidx = ~isinf(aux) & (aux ~= 0);
                dtot(auxidx, I) = aux(auxidx);
            end
            
            % symmetrize the distance matrix by replacing the bottom
            % triangular half with the top half
            dtot = triu(dtot);
            dtot = dtot + dtot';
            
        otherwise
            error('Total distance method not implemented')
    end
end
