function [y, stopCondition, sigma, t] ...
    = qcqp_smacof_pip(dx, y, bnd, w, con, smacof_opts, scip_opts)
% QCQP_SMACOF  Scaling by MAjorizing a COnvex Function (SMACOF) algorithm
% posed as a Quadratic Program with Quadratic Constraints (QPQC).
%
% The classic implementation of SMACOF uses an iterative algorithm that
% relies on the Guttman transform update. This function implements a
% different approach, where the Guttman transform update is replaced by
% solving a Quadratic Program with Quadratic Constraints (QPQC) using the
% SCIP solver in the Opti Toolbox.
%
% The syntax of the QCQP problem in the Opti Toolbox can be seen here
%
%   http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/Probs/QCQP
%
% The objective function is
%
%   min_y 1/2 y' * H * y + f' * y
%
% subject to:
%
%   rl <= A * y <= ul
%   lb <= y <= ub
%
% and a series of quadratic constraints, each of the form
%
%   qrl <= y' * Q * y + l' * x <= qru
%
%
% Y = qcqp_smacof(D, Y0)
%
%   D is an (N, N)-distance matrix, with distances between the points in an
%   N-point configuration. D can be full or sparse.
%
%   Y0 is an initial guess of the solution, given as an (N, P)-matrix,
%   where P is the dimensionality of the output points. Y0 can be generated
%   randomly.
%
%   Y is the solution computed by SMACOF. Y is a point configuration with
%   the same size as Y0.
%
% Y = qcqp_smacof(..., LB, UB, W, A, RL, RU, QC, SMACOF_OPTS, SCIP_OPTS)
%
%   LB, UB, W, A, RL, RU, QC, SMACOF_OPTS, SCIP_OPTS are optional inputs
%   and can be left empty.
%
%   LB, UB are matrices with the same size as Y. They provide lower and
%   upper bounds, respectively, for the corresponding coordinates. While
%   SCIP accepts unbound problems, the solution may diverge or converge
%   slowlier.
%
%   W is a weight matrix the same size as D. W_ij = 0 means that the
%   distance between points i and j does not affect the stress measure. By
%   default, all weights between connected nodes are 1.
%
%   A, RL, RU are a matrix and two vectors, respectively, that provide the
%   linear constraints for the SCIP quadratic program objective function, 
%   RL <= A * Y <= RU. By default, the algorithm runs with no linear
%   constraints.
%
%   QC is a struct with the quadratic constraints for the SCIP quadratic
%   program. Its fields are (compare with constraint formula above):
%
%     'Q': cell vector. Each cell contains the Q matrix for a quadratic
%           constraint.
%     'l': matrix. Each column contains the l vector for a quadratic
%           constraint.
%     'qrl': vector. Each element is the lower bound of a quadratic
%          constraint.
%     'qul': vector. Each element is the upper bound of a quadratic
%          constraint.
%
%   SMACOF_OPTS is a struct with tweaking parameters for the SMACOF
%   algorithm.
%
%     'MaxIter': (default = 0) Maximum number of majorization iterations we
%                allow the optimisation algorithm.
%
%     'Epsilon': (default = Inf) The algorithm will stop if
%                (SIGMA(I+1)-SIGMA(I))/SIGMA(I) < OPTS.Epsilon.
%
%     'Display': (default = 'off') Do not display any internal information.
%                'iter': display internal information at every iteration.
%
%     'TolFun':  (default = 1e-12) Termination tolerance of the stress
%                value.
%
%   SCIP_OPTS is a struct with tweaking parameters for the SCIP algorithm.
%
%     'display': (default 'off') Similar to 'Display' above.
%
%     'warnings': (default 'off')
%
%     'maxiter': (default 1500) Maximum number of iterations in some
%                internal part of the algorithm. This is _not_ the maximum
%                number of overall iterations.
%
%     'maxnodes': (default 10000) Maximum number of branching nodes. If the
%                algorithm runs out of nodes, it is unlikely that
%                increasing this number will produce better solutions.
%
%     'maxtime': (default 1000) Maximum total time allowed for the
%                algorithm to run. As far as I know, this is the closest
%                way to stop the algorithm compared to 'MaxIter' above.
%
%     'tolrfun': (default 1e-6)
%
%     'objbias': (default 0.0)
%
%
% [1] T. Dwyer, Y. Koren, and K. Marriott, "Drawing directed graphs using
% quadratic programming," IEEE Transactions on Visualization and Computer
% Graphics, vol. 12, no. 4, pp. 536�548, 2006.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.0.1
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

%% Input arguments

% check arguments
narginchk(2, 7);
nargoutchk(0, 4);

% number of points
N = size(dx, 1);

% dimensionality of the points
D = size(y, 2);
if (D ~= 2)
    error('Only implemented for 2D output')
end
        

% check inputs
if (N ~= size(dx, 2))
    error('D must be a square matrix')
end
if (N ~= size(y, 1))
    error('Y0 must have the same number of rows as D')
end

% defaults

if (nargin < 3)
    bnd = [];
end
if (nargin < 4 || isempty(w))
    % if the user doesn't provide a weight matrix, we simply assign 1 if
    % two vertices are connected, and 0 if not
    w = double(dx ~= 0);
    if (any(diag(w)) ~= 0)
        error('Assertion error: W matrix has diagonal elements that are non-zero')
    end
end
if (any(size(w) ~= [N N]))
    error('W must be a square matrix with the same size as D')
end
if (nargin < 5)
    con = [];
end

% SMACOF_OPTS defaults
if (nargin < 6 || isempty(smacof_opts) || ~isfield(smacof_opts, 'MaxIter'))
    smacof_opts.MaxIter = 100;
end
if (nargin < 6 || isempty(smacof_opts) || ~isfield(smacof_opts, 'Epsilon'))
    smacof_opts.Epsilon = 0;
end
if (nargin < 6 || isempty(smacof_opts) || ~isfield(smacof_opts, 'Display'))
    smacof_opts.Display = 'off';
end
if (nargin < 6 || isempty(smacof_opts) || ~isfield(smacof_opts, 'TolFun'))
    smacof_opts.TolFun = 1e-12;
end

% SCIP_OPTS defaults

if (nargin < 7 || isempty(scip_opts) || ~isfield(scip_opts, 'display'))
    scip_opts.display = 'off';
end
if (nargin < 7 || isempty(scip_opts) || ~isfield(scip_opts, 'warnings'))
    scip_opts.warnings = 'off';
end
if (nargin < 7 || isempty(scip_opts) || ~isfield(scip_opts, 'maxiter'))
    scip_opts.maxiter = 1500;
end
if (nargin < 7 || isempty(scip_opts) || ~isfield(scip_opts, 'maxnodes'))
    scip_opts.maxnodes = 10000;
end
if (nargin < 7 || isempty(scip_opts) || ~isfield(scip_opts, 'maxtime'))
    scip_opts.maxtime = 1000;
end
if (nargin < 7 || isempty(scip_opts) || ~isfield(scip_opts, 'tolrfun'))
    scip_opts.tolrfun = 1e-6;
end
if (nargin < 7 || isempty(scip_opts) || ~isfield(scip_opts, 'objbias'))
    scip_opts.objbias = 0.0;
end

%% Objective function: 1/2 nu' * H * nu + f' * nu

% pre-compute the weighted Laplacian matrix
V = -w;
V(1:N+1:end) = sum(w, 2);

% quadratic terms of the objective function
objfunq = cell(1, N);

% main diagonal terms
for I = 1:N
    
    objfunq{I} = sprintf(...
        '+%.6g x%d x%d + %.6g y%d y%d', ...
        0.5*full(V(I, I)), I, I, 0.5*full(V(I, I)), I, I);
        
end
objfunq{1} = [' obj: ' objfunq{1}];

% upper triangular matrix terms
for I = 1:N
    for J = I+1:N
        
        % V is sparse, so we don't add terms that are going to be
        % multiplied by 0 anyway
        if (V(I, J))
            
            % main diagonal terms
            objfunq{end+1} = sprintf(...
                '+%.6g x%d x%d + %.6g y%d y%d', ...
                full(V(I, J)), I, J, full(V(I, J)), I, J);
            
        end
        
    end
end

% the linear term of the objective function (f) has to be computed at each
% iteration of the QPQC-SMACOF algorithm. Thus, it is not computed here

%% SMACOF algorithm

% init stopCondition
stopCondition = [];

% Euclidean distances between vertices in the current solution
dy = dmatrix_con(dx, y);

% initial stress
sigma = zeros(1, smacof_opts.MaxIter+1);
sigma(1) = sum(sum(w .* (dx - dy).^2));

% display algorithm's evolution
t = zeros(1, smacof_opts.MaxIter+1); % time past from 0th iteration
tic
if (strcmp(smacof_opts.Display, 'iter'))
    fprintf('Iter\tSigma\t\t\tTime (sec)\n')
    fprintf('===================================================\n')
    fprintf('%d\t\t%.4e\t\t%.4e\n', 0, sigma(1), 0.0)
end

% auxiliary intermediate result
mwdx = -w .* dx;

% majorization loop
for I = 1:smacof_opts.MaxIter

    % auxiliary matrix B: non-main-diagonal elements
    B = mwdx ./ dy;
    B(isnan(B)) = 0;

    % auxiliary matrix B: main diagonal elements
    B(1:N+1:end) = -sum(B, 2);
    
    % the linear term (f) of the quadratic objective function 
    % 1/2 nu' * H * nu + f' * nu
    % has to be recomputed at every iteration
    f = -2 * B * y;
    
    % convert linear term to PIP format
    objfunl = cell(1, N);
    for J = 1:N
        objfunl{J} = sprintf(...
            '+%.6g x%d +%.6g y%d', ...
            f(J, 1), J, f(J, 2), J);
    end
    
    % create PIP file to describe problem
    fid = fopen('/tmp/model.pip', 'w');
    if (fid == -1)
        error('Cannot open file /tmp/model.pip to save PIP model')
    end
    fprintf(fid, '%s\n%s\n%s\n%s\n%s\n%s\n', ...
        'Maximize', objfunq{:}, objfunl{:}, bnd{:}, con{:}, 'End');
    if (fclose(fid) == -1)
        error('Cannot close file /tmp/model.pip to save PIP model')
    end
    
    % solve the quadratic problem
    % Note: fval, exitflag, info are not required by the algorithm, but
    % they are very useful for debugging
%     [y, fval, exitflag, info] ...
    sout = system('~/Documents/gerardus/programs/scip-3.0.2.linux.x86_64.gnu.opt.spx -c "read /tmp/model.pip" -c "optimize" -c "write solution /tmp/model-sol.txt" -c "quit"');
    
    % read solution
    y = read_solution('/tmp/model-sol.txt', size(y));
    
    % recompute distances between vertices in the current solution
    dy = dmatrix_con(dx, y);

    % compute stress with the current solution
    sigma(I+1) = sum(sum(w .* (dx - dy).^2));
    
    % display algorithm's evolution
    t(I+1) = toc;
    if (strcmp(smacof_opts.Display, 'iter'))
        fprintf('%d\t\t%.4e\t\t%.4e\n', I, sigma(I+1), t(I+1))
    end
    
    % check whether the stress is under the tolerance level requested by
    % the user
    if (sigma(I+1) < smacof_opts.TolFun)
        stopCondition{end+1} = 'TolFun';
    end
    
    % check whether the improvement in stress is below the user's request
    if ((sigma(I)-sigma(I+1))/sigma(I) < smacof_opts.Epsilon)
        stopCondition{end+1} = 'Epsilon';
    end

    % stop if any stop condition has been met
    if (~isempty(stopCondition))
        break;
    end
    
end

% check whether the "maximum number of iterations" stop condition has been
% met
if (I == smacof_opts.MaxIter)
    stopCondition{end+1} = 'MaxIter';
end

% prune stress and time vectors if convergence was reached before the
% maximum number of iterations
sigma(I+2:end) = [];
t(I+2:end) = [];

end

function y = read_solution(file, sz)

fid = fopen(file, 'r');
if (fid == -1)
    error(['Cannot open file ' file ' to read solution'])
end

% read contents of the file. Example result:
%
% solution status: optimal solution found
% objective value:                     468.345678663118
% x1                                                 -4 	(obj:0)
% y1                                                  2 	(obj:0)
% x2                                  0.613516669331233 	(obj:0)
% y2                                                 -4 	(obj:0)
% x3                                   1.24777861008035 	(obj:0)
% y3                                  -3.43327058768579 	(obj:0)
% x4                                                 -4 	(obj:0)
% y4                                 -0.804343353251697 	(obj:0)
% x5                                                  2 	(obj:0)
% y5                                                 -4 	(obj:0)
% x6                                  -3.31069797707353 	(obj:0)
% y6                                   1.21192102496956 	(obj:0)
% x7                                     1.643412276563 	(obj:0)
% y7                                  -3.72674560989634 	(obj:0)
% quadobjvar                           468.345678663118 	(obj:1)
c = textscan(fid, '%s%f%s', 'Headerlines', 2, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
fclose(fid);

y = nan(sz);

for I = 1:length(c{1})-1
    
    % variable name
    varname = c{1}{I};
    
    % vertex index
    idx = str2double(varname(2:end));
    
    if (c{1}{I}(1) == 'x') % this is an x-coordinate
        y(idx, 1) = c{2}(I);
    elseif (c{1}{I}(1) == 'y') % this is a y-coordinate
        y(idx, 2) = c{2}(I);
    end
    
end

end
