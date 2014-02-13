function [u, fval, exitflag, info] = vertices_untangle(tri, x, isFree)
% VERTICES_UNTANGLE Find coordinates for a set of free vertices surrounded
% by a polygon of fixed counter-clockwise sorted vertices such that the
% free vertices are untangled.
%
% This function is an extension of the work of Freitag and Plassmann
% (2000), who deal only with the case of a single free vertex.
%
% We pose the untangling problem as a Quadratic Constrained Quadratic
% Program (QCQP) with multiple quadratic constraints. To solve this QCQP,
% we use SCIP [1] by Achterber (2009) with the Matlab interface provided by
% Jonathan Currie in OptiToolbox [2].
%
% U = vertices_untangle(TRI, X, ISFREE)
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh. The orientation provided by the triangles
%   will be assumed by the algorithm to the counter-clockwise, i.e. the
%   algorithm will produce a solution where all triangles have positive
%   signed areas.
%
%   X is a 2-column matrix. X(i, :) contains the xy-coordinates of the
%   i-th node in the mesh.
%
%   ISFREE is a logical vector with one element per point in X.
%   ISFREE(i)=true means that X(i, :) is a free vertex. In that case, the
%   coordinates of X(i, :) are irrelevant, as they will be ignored and
%   recomputed by the untangling algorithm.
%
%   U is a two-colum matrix with the optimal coordinates found by the
%   algorithm for X(ISFREE, :). Thus, to obtain an untangled mesh, run
%   X(ISFREE, :) = U.
%
% [U, FVAL, EXITFLAG, INFO] = vertices_untangle(...)
%
%   FVAL is a scalar with the objective value at the solution.
%
%   EXITFLAG is a scalar with the exit status of the SCIP algorithm.
%      1: optimum found
%     -1: solution does not exist
%
%   INFO is a struct with information about the particulars of the SCIP
%   algorithm.
%
%
% L. A. Freitag and P. Plassmann, "Local optimization-based simplicial mesh
% untangling and improvement", International Journal for Numerical Methods
% in Engineering, 49(1):109-125, 2000.
%
% T. Achterberg, "SCIP: Solving constraint integer programs", Mathematical
% Programming Computation 1(1), pp. 1-41, 2009.
%
% [1] SCIP Solving Constraint Integer Programs, Konrad-Zuse-Zentrum für
% Informationstechnik Berlin, Division Scientific Computing, Department
% Optimization, http://scip.zib.de/.
%
% [2] OptiToolbox. A free Matlab Toolbox for Optimization, by Jonathan
% Currie, http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/Solvers/SCIP.
%
% See also: vertex_untangle, opti_scip.

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
nargoutchk(0, 4);

% number of vertices and triangles
N = size(x, 1);
Ntri = size(tri, 1);

% check inputs
if (size(tri, 2) ~= 3)
    error('TRI must have 3 columns')
end
if (size(x, 2) ~= 2)
    error('X must have 2 columns (i.e. points must be 2D)')
end
if (length(isFree) ~= N)
    error('There must be an element in ISFREE per each point in X')
end
if (~islogical(isFree))
    error('ISFREE must be of type logical (boolean)')
end

%% upper bound of f = max min triangle_area

% directed edges on the external polygon
e = [tri(:, 1:2); tri(:, 2:3); tri(:, [3 1])];
e = e(sum(uint8(~isFree(e)), 2) == 2, :);

% % DEBUG
% for I = 1:size(e, 1)
%     plot(x(e(I, :), 1), x(e(I, :), 2), 'r')
%     pause
% end

% directed adjacency matrix for the boundary polygon
amat = sparse(N, N);
amat(sub2ind([N N], e(:, 1), e(:, 2))) = 1;

% the first directed edge is e.g. [18 15]. That is, we know if that we ask
% for the shorted path from 15 to 18, the path will have to go all the way
% around the polygon
[~, polyidx] = graphshortestpath(amat, e(1, 2), e(1, 1));

% % DEBUG
% for I = 1:length(polyidx)-1
%     plot(x(polyidx([I I+1]), 1), x(polyidx([I I+1]), 2), 'g')
%     pause
% end

% polygon area
apoly = polyarea(x(polyidx, 1), x(polyidx, 2));

%% upper and lower bounds for the quadratic problem variables

% box that contains all the mesh points
xmin = min(x);
xmax = max(x);

% decision variable bounds, lb<=nu<=ub
lb = zeros(2*nnz(isFree)+1, 1);
ub = zeros(2*nnz(isFree)+1, 1);

lb(1:2:end) = xmin(1);
lb(2:2:end) = xmin(2);
lb(end) = 0;

ub(1:2:end) = xmax(1);
ub(2:2:end) = xmax(2);
ub(end) = apoly;

%% quadratic constraints for the quadratic problem

% quadratic constrained quadratic program: initialise constraints
qc.Q = cell(1, Ntri);
qc.l = zeros(2*N + 1, Ntri);
qc.qrl = zeros(Ntri, 1); % lower bounds for quadratic constraints
qc.qru = zeros(Ntri, 1); % upper bounds for quadratic constraints

% each triangle contributes one quadratic constraint of the form 
% - c <= nu'*Q*nu + l'*x < 2*Apoly - c
for I = 1:Ntri
    
    % coordinates of the three vertices in the triangle
    triloc = tri(I, :);
    xloc = x(triloc, :);
    isFreeloc = isFree(triloc);
    
    % depending on the number of free vertices in the triangle, we create
    % different constraints
    switch (nnz(isFreeloc))
        
        case 0 % 0 free vertices, 3 fixed vertices
            
            % this case doesn't contribute any constraints to the quadratic
            % program, as the fixed vertices cannot be moved
        
        case 1 % 1 free vertex, 2 fixed vertices
    
            % position of the free vertex
            idx = find(isFreeloc);
            
            % shift the vertices until the free vertex is the last one,
            % without changing the sign of the area
            xloc = circshift(xloc, 3 - idx);
            isFreeloc = circshift(isFreeloc, 3 - idx);
            if any(isFreeloc ~= [0 0 1]')
                error(['Assertion: Triangle ' num2str(I) ...
                    ' vertices cannot be shifted to canonic configuration'])
            end
            triloc = circshift(triloc(:), 3 - idx);
            
            % auxiliary variables to make the code more readable
            x1 = xloc(1, 1); % x-coordinate of 1st fixed vertex
            y1 = xloc(1, 2); % y-coordinate of 1st fixed vertex
            x2 = xloc(2, 1); % x-coordinate of 2nd fixed vertex
            y2 = xloc(2, 2); % y-coordinate of 2nd fixed vertex
            idx = triloc(3); % free vertex index in the mesh
            
            % quadratic constraint
            qc.Q{I} = sparse(2*N + 1, 2*N + 1); % Q = 0;
            qc.l(2*(idx-1)+[1 2], I) = [y1 - y2; x2 - x1];
            qc.l(end, I) = -2;
            qc.qrl(I) = -(x1 * y2 - x2 * y1);
            qc.qru(I) = 2*apoly + qc.qrl(I);
            
        case 2 % 2 free vertices, 1 fixed vertex
            
            % position of the fixed vertex
            idx = find(~isFreeloc);
            
            % shift the vertices until the fixed vertex is the first one,
            % without changing the sign of the area
            xloc = circshift(xloc, 1 - idx);
            isFreeloc = circshift(isFreeloc, 1 - idx);
            if any(isFreeloc ~= [0 1 1]')
                error(['Assertion: Triangle ' num2str(I) ...
                    ' vertices cannot be shifted to canonic configuration'])
            end
            triloc = circshift(triloc(:), 1 - idx);
            
            % auxiliary variables to make the code more readable
            x1 = xloc(1, 1); % x-coordinate of fixed vertex
            y1 = xloc(1, 2); % y-coordinate of fixed vertex
            
            % indices of the elements in the vector of unknown variables
            % that correspond to the coordinates of the free vertices
            idx_u1 = 2*(triloc(2)-1) + 1;
            idx_v1 = idx_u1 + 1;
            idx_u2 = 2*(triloc(3)-1) + 1;
            idx_v2 = idx_u2 + 1;
            
            % quadratic constraint
            qc.Q{I} = sparse(2*N + 1, 2*N + 1);
            qc.Q{I}([idx_u1, idx_v1, idx_u2, idx_v2], ...
                [idx_u1, idx_v1, idx_u2, idx_v2]) = [
                0   0   0   .5
                0   0   -.5 0
                0   -.5 0   0
                .5  0   0   0
                ];
            qc.l([idx_u1, idx_v1, idx_u2, idx_v2], I) = [-y1 x1 y1 -x1];
            qc.l(end, I) = -2;
            qc.qrl(I) = 0;
            qc.qru(I) = 2*apoly;
            
        case 3 % three free vertices
            
            % auxiliary variables to make the code more readable
            idx_u1 = 2*(triloc(1)-1) + 1;
            idx_v1 = idx_u1 + 1;
            idx_u2 = 2*(triloc(2)-1) + 1;
            idx_v2 = idx_u2 + 1;
            idx_u3 = 2*(triloc(3)-1) + 1;
            idx_v3 = idx_u3 + 1;
            
            % quadratic constraints
            qc.Q{I} = sparse(2*N + 1, 2*N + 1);
            qc.Q{I}([idx_u1, idx_v1, idx_u2, idx_v2, idx_u3, idx_v3], ...
                [idx_u1, idx_v1, idx_u2, idx_v2, idx_u3, idx_v3]) = [
                0 0 0 .5 0 -.5
                0 0 -.5 0 .5 0
                0 -.5 0 0 0 .5
                .5 0 0 0 -.5 0
                0 .5 0 -.5 0 0
                -.5 0 .5 0 0 0
                ];
            qc.l(end, I) = -2;
            qc.qrl(I) = 0;
            qc.qru(I) = 2*apoly;
            
        otherwise
            
            error('Assertion error: Triangle has more than 3 vertices')
            
    end
    
end


% we only need to solve the problem for the coordinates of the free
% vertices, the fixed vertices remain constant. Thus, we are going to
% remove the rows and columns that correspond to the fixed vertices in the
% constraints
idx = (find(isFree)*2-1);
idx = sort([idx; idx+1]);
for I = 1:Ntri
    
    qc.Q{I} = qc.Q{I}([idx; end], [idx; end]);
    
end
qc.l = qc.l([idx; end], :);

%% create OPTI object

% quadratic program objective function
H = [];
b = zeros(2*nnz(isFree) + 1, 1);
b(end) = 1;

% linear constraints
A = sparse([]);
rl = [];
ru = [];

% coordinates and area are continuous variables
xtype = repmat('c', 1, 2*nnz(isFree)+1);

% we don't have any special ordered set (SOS) constraints
sos = [];

% options
opts.display = 'iter'; % DEBUG
% opts.warnings = 'off'; % DEBUG
% opts.display = 'off';
opts.warnings = 'off';

% default options
% opts.maxiter = 1500;
opts.maxiter = 1;
% opts.maxnodes = 10000;
 opts.maxnodes = 500000;
% opts.maxtime = 1000;
opts.maxtime = 5000;
% opts.tolrfun = 1e-6;
% opts.objbias = 0.0;
% opts.display = 0;
opts.freq = 1;
opts.limits.solutions = 1;

opts.gamsfile = 'rcasero-problem.gams';

% solve the quadratic problem
[nu, fval, exitflag, info] ...
    = opti_scip(H, -b, A, rl, ru, lb, ub, xtype, sos, qc, opts);

% format the output
u = reshape(nu(1:end-1, :), 2, nnz(isFree))';
