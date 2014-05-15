function [y, stopCondition, sigma, t] = smacof(dx, y, w, opts)
% SMACOF  Classic implementation of SMACOF (Scaling by MAjorizing a COnvex
% Function) algorithm for MDS (Multidimensional Scaling)
%
% SMACOF (Scaling by MAjorizing a COnvex Function) is an iterative method
% to solve Multidimensional Scaling problems, proposed by de Leeuw and
% Heiser [1]-[3] from the late 1970s. A modern account can be found in [4].
%
% An advantage of SMACOF over classical MDS (cmdscale) is that the former
% accepts sparse distance matrices. A disadvantage is that the solution is
% found iteratively and can get trapped in local minima.
%
% Y = smacof(D, Y0)
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
% [..., STOPCONDITION, SIGMA] = smacof(..., W, OPTS)
%
%   STOPCONDITION is a cell array with the condition that were met
%   simultaneously to make the algorithm stop.
%
%   SIGMA is a vector with the weighted stress value at each iteration:
%   SIGMA(1) corresponds to the stress of Y0, SIGMA(2) corresponds to the
%   first majorization iteration, etc. Weighted stress is given as
%
%     SIGMA = \sum_{i<j} W_ij (D_ij - DY_ij)^2
%
%   where W is a weight matrix the same size as D. W_ij = 0 means that the
%   distance between points i and j does not affect the stress measure.
%
%   OPTS is a struct with parameters for the algorithm:
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
%
% [1] J. De Leeuw, "Applications of convex analysis to multidimensional
% scaling," Recent Developments in Statistics, pp. 133�146, 1977.
%
% [2] J. De Leeuw and W. J. Heiser, "Convergence of correction matrix
% algorithms for multidimensional scaling," ser. Geometric representations
% of relational data, J. C. Lingoes, Ed. Mathesis Press, 1977, pp. 735�753.
%
% [3] ��, "Multidimensional scaling with restrictions on the
% configuration," ser. Multivariate analysis, P. R. Krishnaiah, Ed., vol.
% 5. North Holland Publishing Company, 1980, pp. 501�522.
%
% [4] J. De Leeuw and P. Mair, "Multidimensional scaling using
% majorization: SMACOF in R," Journal of Statistical Software, vol. 31, no.
% 3, 2009.
%
% See also: cmds.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.0.3
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
narginchk(2, 4);
nargoutchk(0, 4);

% number of points
N = size(dx, 1);

% check inputs
if (N ~= size(dx, 2))
    error('D must be a square matrix')
end
if (N ~= size(y, 1))
    error('Y0 must have the same number of rows as D')
end

% defaults

% if the user doesn't provide a weight matrix, we simply assign 1 if two
% vertices are connected, and 0 if not
if (nargin < 3 || isempty(w))
    w = double(dx ~= 0);
    if (any(diag(w)) ~= 0)
        error('Assertion error: W matrix has diagonal elements that are non-zero')
    end
end
if (any(size(w) ~= [N N]))
    error('W must be a square matrix with the same size as D')
end

if (nargin < 4 || isempty(opts) || ~isfield(opts, 'MaxIter'))
    opts.MaxIter = 100;
end
if (nargin < 4 || isempty(opts) || ~isfield(opts, 'Epsilon'))
    opts.Epsilon = 0;
end
if (nargin < 4 || isempty(opts) || ~isfield(opts, 'Display'))
    opts.Display = 'off';
end
if (nargin < 4 || isempty(opts) || ~isfield(opts, 'TolFun'))
    opts.TolFun = 1e-12;
end

% pre-compute the weighted Laplacian matrix
V = -w;
V(1:N+1:end) = sum(w, 2);

% pre-compute peudoinverse of V (V is square but not full rank, so we
% cannot compute the inverse)
Vinv = pinv(full(V));

% we don't need V anymore
clear V

% init stopCondition
stopCondition = [];

% Euclidean distances between vertices in the current solution
dy = dmatrix_con(dx, y);

% initial stress
sigma = zeros(1, opts.MaxIter+1);
sigma(1) = sum(sum(w .* (dx - dy).^2));

% display algorithm's evolution
t = zeros(1, opts.MaxIter+1); % time past from 0th iteration
tic
if (strcmp(opts.Display, 'iter'))
    fprintf('Iter\tSigma\t\t\tTime (sec)\n')
    fprintf('===================================================\n')
    fprintf('%d\t\t%.4e\t\t%.4e\n', 0, sigma(1), 0.0)
end

% auxiliary intermediate result
mwdx = -w .* dx;

% majorization loop
for I = 1:opts.MaxIter

    % auxiliary matrix B: non-main-diagonal elements
    B = mwdx ./ dy;
    B(isnan(B)) = 0;

    % auxiliary matrix B: main diagonal elements
    B(1:N+1:end) = -sum(B, 2);
    
    % Guttman transform update
    y = Vinv * B * y;

    % recompute distances between vertices in the current solution
    dy = dmatrix_con(dx, y);

    % compute stress with the current solution
    sigma(I+1) = sum(sum(w .* (dx - dy).^2));
    
    % display algorithm's evolution
    t(I+1) = toc;
    if (strcmp(opts.Display, 'iter'))
        fprintf('%d\t\t%.4e\t\t%.4e\n', I, sigma(I+1), t(I+1))
    end
    
    % check whether the stress is under the tolerance level requested by
    % the user
    if (sigma(I+1) < opts.TolFun)
        stopCondition{end+1} = 'TolFun';
    end
    
    % check whether the improvement in stress is below the user's request
    if ((sigma(I)-sigma(I+1))/sigma(I) < opts.Epsilon)
        stopCondition{end+1} = 'Epsilon';
    end
    
    % stop if any stop condition has been met
    if (~isempty(stopCondition))
        break;
    end
    
end

% check whether the "maximum number of iterations" stop condition has been
% met
if (I == opts.MaxIter)
    stopCondition{end+1} = 'MaxIter';
end

% prune stress and time vectors if convergence was reached before the
% maximum number of iterations
sigma(I+2:end) = [];
t(I+2:end) = [];
