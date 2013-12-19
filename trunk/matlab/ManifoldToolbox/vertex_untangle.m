function [x0, exitflag] = vertex_untangle(tri, x, idx0)
% VERTEX_UNTANGLE  Find 2D coordinates for a free vertex surrounded
% by a polygon of fixed counter-clockwise sorted vertices such that the
% free vertex is untangled.
%
% This function is an enhanced implementation of the 2D case (i.e. each
% simplex is a triangle) of the "Optimization-based mesh untangling"
% method proposed by Freitag and Plassmann (2000).
%
% Note that the original paper contains some errors in the formulation of
% the linear programming problem. The correct problem formulation is
%
%     max b^T \pi (equivalenty: min -b^T \pi)
%     such that -A^T \pi <= c
%
% where b = [0 0 1], A^T = [ax ay -2*ones(N, 1)], N=number of
% triangles, and ax, ay, c as given in the paper.
%
% We check whether the free vertex is already untangled. In that case, we
% do not relocate it (this is an untangling algorithm, not a mesh
% improvement or mesh smoothing algorithm).
%
% X0 = vertex_untangle(TRI, X);
%
%   Let TRI, X describe a closed fan triangular mesh. That is, we have a
%   central or free vertex connected to N neighbours. The neighbours are
%   connected between them to form a closed perimeter. The perimeter is
%   assumed to be oriented in counter-clockwise orientation (that is, if a
%   solution exists, the areas of the triangles will all be positive).
%
%   We want to find the coordinates X0 for the free vertex, such that the
%   edges that connect vertices don't overlap. The N neighbours remain
%   fixed. If the free vertex was already producing any overlaps (i.e. the
%   mesh was entangled), this can be seen as an untangling algorithm.
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh. The orientation provided by the triangles
%   will be assumed by the algorithm to the counter-clockwise, i.e. the
%   algorithm will produce a solution where all triangles have positive
%   signed areas.
%
%   X is a 2-column matrix. X(i, :) contains the xy-coordinates of the
%   i-th node in the mesh. Note that Freitag and Plassmann (2000) also
%   provided the solution for the 3D case (i.e. each simplex is a
%   tetrahedron), but we haven't implemented it here yet.
%
%   X0 is a 2-vector with the optimal coordinates of the central (a.k.a.
%   free) vertex, in a linear programming sense, 
%
% [X0, EXITFLAG] = vertex_untangle(...)
%
%   EXITFLAG is the exit condition of the algorithm.
%
%     2: Free vertex was not tangled, no relocation performed.
%     1: Linear programming algorithm converged to a relocation solution.
%     0: Linear programming algorithm stopped because of maximum number of
%        iterations reached.
%    -1: Linear programming algorithm found no valid solution exists.
%
% L. A. Freitag and P. Plassmann, "Local optimization-based simplicial mesh
% untangling and improvement", International Journal for Numerical Methods
% in Engineering, 49(1):109-125, 2000.
%
% See also: sphtri_untangle, linprog, surfreorient.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.2.3
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
nargoutchk(0, 3);

if (size(tri, 2) ~= 3)
    error('TRI must have 3 columns')
end
if (size(x, 2) ~= 2)
    error('X must have 2 columns')
end
if (~isscalar(idx0))
    error('IDX must be a scalar')
end

% number of vertices in the neighbourhood
N = size(x, 1);

% adjacency matrix for the mesh. Edges are directed, thus the matrix is not
% symmetric. We do it this way to preserve the counter-clockwise
% orientation of the neighbourhood, i.e. the positive sign of the areas
d = sparse(N, N);
d(sub2ind([N, N], tri(:, 1), tri(:, 2))) = 1;
d(sub2ind([N, N], tri(:, 2), tri(:, 3))) = 1;
d(sub2ind([N, N], tri(:, 3), tri(:, 1))) = 1;

% number of neighbours for each vertex
nn = full(sum((d | d') ~= 0));

% valid neighbourhood check: the central vertex is connected to every other
% vertex
if (nn(idx0) ~= size(x, 1)-1)
    error(['Invalid neighbourhood: central vertex connected to ' ...
        num2str(nn(idx0)) ' vertices instead of ' num2str(size(x, 1) - 1)])
end
% valid neighbourhood check: non-central vertices must be connected to the
% central vertex, and to another two vertices each
if (any(nn([1:idx0-1 idx0+1:end]) ~= 3))
    error('Invalid neighbourhood: at least one boundary vertex not connected to another 3 vertices')
end

% get sorted list of vertices that form the perimeter
v = sort_perim_vertices(d, idx0);

% % DEBUG: plot boundary edges, marking the first segment
% subplot(2, 2, 1)
% hold off
% gplot(d, x)
% hold on
% plot(x(v([1:end 1]), 1), x(v([1:end 1]), 2), 'r')
% plot(x(v(1), 1), x(v(1), 2), 'ro')
% plot(x(v(2), 1), x(v(2), 2), 'rp')

% compute the optimal coordinates
% [x1, triarea, exitflag, output] ...
[x1, ~, exitflag] ...
    = linear_programming_coordinates(x, idx0, v);

% % DEBUG: plot solution of current orientation
% subplot(2, 2, 2)
% aux = x;
% aux(idx0, :) = x1;
% hold off
% gplot(d, aux)
% hold on
% plot(aux(v([1:end 1]), 1), aux(v([1:end 1]), 2), 'r')
% plot(aux(v(1), 1), aux(v(1), 2), 'ro')
% plot(aux(v(2), 1), aux(v(2), 2), 'r*')

% % DEBUG: plot solution of opposite orientation
% subplot(2, 2, 4)
% aux = x;
% aux(idx0, :) = x2;
% hold off
% gplot(d, aux)
% hold on
% plot(aux(v([1:end 1]), 1), aux(v([1:end 1]), 2), 'r')
% plot(aux(v(1), 1), aux(v(1), 2), 'ro')
% plot(aux(v(2), 1), aux(v(2), 2), 'r*')

%% choose the best solution

% if unsuccessful, we return NaN as a solution
if (exitflag < 0)
    
    x0 = nan(1, 2);
    exitflag = -1;

% else, we return the new location
else
    
    x0 = x1;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary functions

% d: the sparse adjacency matrix with the directed connections between the
% perimeter vertices
%
% idx0: index of the central vertex
%
% v: ordered list of perimeter vertices in counter-clockwise orientation
% (so that the signed area is positive)
function v = sort_perim_vertices(d, idx0)

% clear up the connections from the perimeter to the central vertex
d(idx0, :) = 0;
d(:, idx0) = 0;

% Dijkstra's shortest paths from one arbitrary vertex in the perimeter to
% all the others
if (idx0 == 1) 
    va = 2; 
else
    va = 1; 
end
[l, p] = dijkstra(d, va);

% furtherst vertex from origin vertex
l(isinf(l)) = 0;
[~, vb] = max(l);

% path from va to vb
v = graphpred2path(p, vb);

% invert path
v = v(end:-1:1);

end

% set up the Freitag and Plassmann (2000) linear programming problem and
% compute the solution. Note that the formulation in the original paper is
% wrong, as has been corrected here.
%
% x: 2D vertex coordinates
%
% v0: index of the free vertex
%
% v: indices of the perimeter vertices, sorted in counter-clockwise
% orientation (so that the signed area is positive)
%
% lb, ub: lower and upper bounds for the linear programming algorithm
function [x0, triarea, exitflag, output] ...
    = linear_programming_coordinates(x, v0, v)

% for convenience, express the boundary as a list of edges
e = [v(1:end)' v([2:end 1])'];

% number of triangles
Ntri = size(e, 1);

% terms for the linear programming problem
xi = x(e(:, 1), 1);
xj = x(e(:, 2), 1);
yi = x(e(:, 1), 2);
yj = x(e(:, 2), 2);
ax = yi - yj;
ay = xj - xi;
c = xi .* yj - xj .* yi;
A = [ax'; ay'; -2*ones(1, Ntri)];

% compute triangle areas (this formula comes from expanding the determinant
% form of a triangle's area)
triarea = 0.5 * (ax * x(v0, 1) + ay * x(v0, 2) + c);

% if all the areas are positive, then there's no tangling, and we won't
% relocate the free vertex (this is not a smoothing algorithm, thus the
% free vertex will only be moved if it's tangled)
if all(triarea > 0)
    x0 = x(v0, :);
    exitflag = 2;
    output = [];
    return
end

% the solution cannot be outside the convex hull of the boundary vertices.
% Thus, it cannot be outside of a rectangular box containing the boundary
% vertices. We can use this as the lower and upper bounds for the linear
% programming solution
lb = min(x(v, :));
ub = max(x(v, :));

% the third component of the linear programming solution vector = minimum
% area. We bound the area value between 0 and infinite. We could use as the
% upper bound the convex hull's area, but it's not worth wasting time
% computing it.
lb = [lb 0]';
ub = [ub Inf]';

% otherwise, we solve linear programming problem to try to untangle the
% free vertex
options = optimset('Display', 'off', 'Simplex', 'on', 'LargeScale', 'off');
[x0, ~, exitflag, output] = linprog([0 0 -1]', -A', c, [], [], ...
    lb, ub, [], options);

% the two first elements are the target coordinates of the central vertex
x0 = x0(1:2)';

end
