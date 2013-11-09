function [x0, exitflag] = vertex_untangle(tri, x)
% VERTEX_UNTANGLE  Location of a 2D vertex at the centre of a closed fan
% triangular mesh so that we avoid edge overlapping.
%
% This function is an enhanced implementation of the 2D case (i.e. each
% simplex is a triangle) of the "Optimization-based mesh untangling"
% solution proposed by Freitag and Plassmann (2000). The differences are:
%
%   - The original paper contains an error in the sign of "c" in the linear
%     programming inequalities. The correct problem formulation is
%
%       max b^T \pi (equivalenty: min -b^T \pi)
%       such that A^T \pi <= -c
%
%   - Instead of the simplex method, we use Matlab's linprog default,
%   'interior-point'
%   (http://www.mathworks.co.uk/help/optim/ug/linprog.html).
%
%   - The paper assumes that the triangles are correctly oriented such that
%     the correct position of the central voxel makes all areas positive.
%     In practice, we don't know that. Thus, this function computes 
%     solutions for one orientation of the perimeter vertices and the
%     opposite, and chooses whichever solution is correct. If both are
%     correct, the one with the largest minimum area is chosen.
%
% Let TRI, X describe a closed fan triangular mesh. That is, we have a
% central or free vertex connected to N neighbours. The neighbours are
% connected between them to form a closed perimeter.
%
% We want to find the coordinates X0 for the free vertex, such that the
% edges that connect vertices don't overlap. The N neighbours remain fixed.
% If the free vertex was already producing any overlaps (i.e. the mesh was
% entangled), this can be seen as an untangling algorithm.
%
% X0 = vertex_untangle(TRI, X);
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh.
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
%     1: The algorithm has converged.
%     0: Maximum number of iterations reached.
%    -1: No valid solution.
%
% L. A. Freitag and P. Plassmann, "Local optimization-based simplicial mesh
% untangling and improvement", International Journal for Numerical Methods
% in Engineering, 49(1):109-125, 2000.
%
% See also: sphtri_untangle, linprog.

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
narginchk(2, 2);
nargoutchk(0, 3);

if (size(tri, 2) ~= 3)
    error('TRI must have 3 columns')
end
if (size(x, 2) ~= 2)
    error('X must have 2 columns')
end

% adjacency/distance matrix for the mesh
d = dmatrix_mesh(tri);

% number of neighbours for each vertex
nn = full(sum(d ~= 0));

% the central vertex must be connected to every other vertex. Find which
% one it is
idx0 = find(nn == size(d, 1)-1);
if (length(idx0) ~= 1)
    error('This neighbourhood does not have a unique central vertex that is adjacent to all other vertices')
end

% the solution cannot be outside the convex hull of the boundary vertices.
% Thus, it cannot be outside of a rectangular box containing the boundary
% vertices. We can use this as the lower and upper bounds for the linear
% programming solution
idxn = [1:idx0-1 idx0+1:size(d, 1)];
lb = min(x(idxn, :));
ub = max(x(idxn, :));

% we don't give bounds to the third component of the vector, that is the
% minimum area. We need to turn the bound vectors to column format
lb = [lb 0]';
ub = [ub Inf]';

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

% compute the optimal coordinates  for perimeter with current orientation
[x1, triarea1, exitflag1, output1] ...
    = linear_programming_coordinates(x, v, lb, ub);

% compute the optimal coordinates  for perimeter with opposite orientation
[x2, triarea2, exitflag2, output2] ...
    = linear_programming_coordinates(x, v(end:-1:1), lb, ub);

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

% if both orientations were unsuccessful, we return NaN as a solution
if ((exitflag1 < 0) && (exitflag2 < 0))
    
    x0 = nan(1, 2);
    exitflag = -1;

% if both were successful, we choose the solution that has the largest min
% triangle area
elseif ((exitflag1 >= 0) && (exitflag2 >= 0))

    if (min(triarea1) > min(triarea2))
        x0 = x1;
        exitflag = exitflag1;
    else
        x0 = x2;
        exitflag = exitflag2;
    end

% if only one orientation was successful
elseif (exitflag1 >= 0)
    
    x0 = x1;
    exitflag = exitflag1;
    
elseif (exitflag2 >= 0)
    
    x0 = x2;
    exitflag = exitflag2;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary functions

% d is the sparse adjacency matrix with the connections between the
% perimeter vertices
%
% v is a list of the 
function v = sort_perim_vertices(d, idx0)

% clear up the connections from the perimeter to the central vertex
d(idx0, :) = 0;
d(:, idx0) = 0;

% Dijkstra's shortest paths from one arbitrary vertex in the perimeter to
% all the others
if (idx0 == 1) va = 2; else va = 1; end
[l, p] = dijkstra(d, va);

% arbitrarily choose one vertex at the antipodes in the perimeter 
l(isinf(l)) = 0;
[~, vb] = max(l);

% path from va to vb
branch1 = graphpred2path(p, vb);

% remove the connections in branch1, so that we can find the other branch
V = size(d, 1);
d(sub2ind([V V], branch1(1:end-1)', branch1(2:end)')) = 0;
d(sub2ind([V V], branch1(2:end)', branch1(1:end-1)')) = 0;

% recompute Dijkstra's shortest path
[~, p] = dijkstra(d, va);
if (isinf(d(vb)))
    error('Gerardus:Assertion', 'Assertion fail: Antipode vertex only reachable through one branch')
end

% 2nd path from va to vb
branch2 = graphpred2path(p, vb);

% combine both branches into final solution
v = [branch1 branch2(end-1:-1:2)];
    
end

% set up the Freitag and Plassmann (2000) linear programming problem and
% compute the solution
function [x0, triarea, exitflag, output] ...
    = linear_programming_coordinates(x, v, lb, ub)

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
A = [ax'; ay'; ones(1, Ntri)];

% solve linear programming problem
options = optimset('Display', 'off');
[x0, ~, exitflag, output] = linprog(-[0 0 1]', A', -c, [], [], lb, ub, [], ...
    options);

% the two first elements are the target coordinates of the central vertex
x0 = x0(1:2)';

% compute triangle areas
triarea = zeros(1, Ntri);
for I = 1:Ntri
    m = [[x(e(I, :), :); x0] ones(3, 1)];
    triarea(I) = 0.5 * det(m);
end

end
