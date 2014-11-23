function [y, stopCondition, sigma, sigma0, t] ...
    = cons_smacof_pip(dx, y, isFree, bnd, w, con, smacof_opts, scip_opts)
% CONS_SMACOF_PIP  SMACOF algorithm with polynomial constraints (PIP file
% format).
%
% Scaling by MAjorizing a COnvex Function (SMACOF) is an iterative solution
% to the Multidimensional Scaling (MDS) problem (see de Leeuw and Mair [1]
% for an up to date review).
%
% The classic implementation of SMACOF uses an iterative algorithm that
% relies on the Guttman transform update. Our function here implements a
% different approach, where the Guttman transform update is replaced by
% solving a Quadratic Program (QP) (as proposed by Dwyer et al [2]) with
% polynomial constraints.
%
% In this implementation, we use the SCIP binary to solve the constrained
% QP. This binary can be downloaded for Linux, MacOS X and Windows from the
% Zuse Institute Berlin (ZIB) website
%
%   http://scip.zib.de/#download
%
% The objective function is
%
%   min_y 1/2 y' * H * y + f' * y
%
% subject to the constraints and bounds provided by the user.
%
% We use the PIP file format to formulate the constrained QP and pass it to
% SCIP.
%
%   http://polip.zib.de/pipformat.php
%
%
% [Y, STOPCONDITION, SIGMA, SIGMA0, T] = cons_smacof_pip(D, Y0, ISFREE, BND, [], CON)
%
%   D is an (N, N)-distance matrix, with distances between the points in an
%   N-point configuration. D can be full or sparse. D(i,j)=0 means that
%   vertices i and j are not directly connected.
%
%   Y0 is an initial guess of the solution, given as an (N, P)-matrix,
%   where P is the dimensionality of the output points. Currently, P must
%   be either 2 or 3.
%
%   BND is a cell array with the variable bounds in PIP format, and cannot
%   be empty (otherwise SCIP doesn't return a solution). E.g.
%
%      BND = {'Bounds', ' -1 <= x1 <= 4', ' -1 <= y1 <= 2.5', ...
%             ' -1 <= x2 <= 4', ' -1 <= y2 <= 2.5'};
%
%   CON is a cell array with the problem constraints in PIP format, and
%   cannot be empty (otherwise SCIP doesn't return a solution). E.g.
%
%      CON = {'Subject to', ...
%            ' c1: -0.5 x6 y7 +0.5 x3 y7 +0.5 x7 y6 -0.5 x3 y6 -0.5 x7 y3 +0.5 x6 y3 >= 0.1'};
%
%   Y is the best solution computed by the algorithm within all iterations.
%   Y is a point configuration with the same size as Y0. If the algorithm
%   cannot find any valid solution (e.g. because we set a too short time
%   limit, or because a valid solution doesn't exist, Y will be an array of
%   NaNs).
%
%   STOPCONDITION is a cell array with a string for each stop condition
%   that made the algorithm stop at the last iteration.
%
%   SIGMA is a vector with the stress value at each iteration. Weighted
%   stress is given as
%
%     SIGMA = \sum_{i<j} W_ij (D_ij - DY_ij)^2
%
%   where W is a weight matrix the same size as D. W_ij = 0 means that the
%   distance between points i and j does not affect the stress measure.
%
%   If the algorithm could not find any solution, SIGMA is a single scalar
%   with value = Inf. Note that even if the last SIGMA value is Inf, this
%   doesn't mean that the algorithm hasn't found a valid solution in a
%   previous iteration.
%
%   SIGMA0 is a scalar with the stress value of the initial guess Y0. Note
%   that we have no guarantee that Y0 is a valid solution, so be careful
%   when comparing SIGMA0 to SIGMA. In particular, SIGMA0 may be smaller
%   than any SIGMA, but may not be a valid solution according to the
%   constraints.
%
%   T is a vector with the time between the beginning of the algorithm and
%   each iteration. Units in seconds.
%
% Y = cons_smacof_pip(..., SMACOF_OPTS, SCIP_OPTS)
%
%   SMACOF_OPTS is a struct with parameters to tweak the SMACOF algorithm.
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
%   SCIP_OPTS is a struct with parameters to tweak the SCIP algorithm.
%
%     'scipbin': (defaults = 'scip-3.1.0.linux.x86_64.gnu.opt.spx' (Linux),
%                            'scip-3.1.0.darwin.x86_64.gnu.opt.spx' (Mac),
%                            'scip-3.1.0.mingw.x86_64.intel.opt.spx.exe' (Win))
%                Name of the SCIP binary/executable. This binary should be
%                available in the system path. E.g. place it in
%                gerardus/programs. The binaries/executable can be
%                downloaded from http://scip.zib.de/#download.
%
%     'display_verblevel': (default 4) verbosity level of output (0: SCIP
%                quiet mode).
%
%     'display_freq': (default 100) frequency for displaying node
%                information lines.
%
%     'limits_absgap': (default 0.0) solving stops, if the absolute 
%                gap = |primalbound - dualbound| is below the given value.
%
%     'limits_gap': (default 0.0) solving stops, if the relative 
%                gap = |primal - dual|/MIN(|dual|,|primal|) is below the
%                given value.
%
%     'limits_time': (default 1e+20) maximal time in seconds to run.
%
%     'limits_solutions': (default -1) solving stops, if the given number
%                of solutions were found (-1: no limit).
%
%     'lp_threads': (default 0: automatic) number of threads used for
%                solving the LP.
%
%     'numerics_feastol': (default 1e-6) feasibility tolerance for
%                constraints in SCIP. If this value is changed from the
%                default, it needs to be passed to the function that
%                creates the constraints too, e.g.
%                tri_ccqp_smacof_nofold_sph_pip().
%
%
% [1] J de Leeuw, P Mair, "Multidimensional scaling using majorization:
% SMACOF in R", Journal of Statistical Software, 31(3), 2009.
%
% [2] T. Dwyer, Y. Koren, and K. Marriott, "Drawing directed graphs using
% quadratic programming," IEEE Transactions on Visualization and Computer
% Graphics, vol. 12, no. 4, pp. 536-548, 2006.
%
% See also: cmdscale, qcqp_smacof, tri_sphparam.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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

%% Input arguments

% check arguments
narginchk(6, 8);
nargoutchk(0, 5);

% start clock
tic

% number of points
N = size(dx, 1);

% dimensionality of the points
D = size(y, 2);
if ((D ~= 2) && (D ~= 3))
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
if (nargin < 3 || isempty(isFree))
    isFree = true(N, 1);
end
Nfree = nnz(isFree);

if (isempty(w))
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

% if the user doesn't allow any vertices to be optimized, then we can exit
% the function
if (Nfree == 0)
    stopCondition = 'No free vertices to optimise';
    sigma = [];
    dy = dmatrix_con(dx, y);
    sigma0 = 0.5 * sum(sum(w .* (dx - dy).^2));
    t = [];
    return;
end

% SMACOF_OPTS defaults
if (nargin < 7 || isempty(smacof_opts) || ~isfield(smacof_opts, 'MaxIter'))
    smacof_opts.MaxIter = 100;
end
if (nargin < 7 || isempty(smacof_opts) || ~isfield(smacof_opts, 'Epsilon'))
    smacof_opts.Epsilon = 0;
end
if (nargin < 7 || isempty(smacof_opts) || ~isfield(smacof_opts, 'Display'))
    smacof_opts.Display = 'off';
end
if (nargin < 7 || isempty(smacof_opts) || ~isfield(smacof_opts, 'TolFun'))
    smacof_opts.TolFun = 1e-12;
end

% SCIP_OPTS
if (nargin < 8 || isempty(scip_opts) || ~isfield(scip_opts, 'scipbin'))
    
    % default name of the SCIP binary depending on the architecture
    if (isunix)
        SCIPBIN = 'scip-3.1.0.linux.x86_64.gnu.opt.spx';
    elseif (ismac)
        SCIPBIN = 'scip-3.1.0.darwin.x86_64.gnu.opt.spx';
    elseif (ispc)
        SCIPBIN = 'scip-3.1.0.mingw.x86_64.intel.opt.spx.exe';
    else
        error('Operating system not recognized: I do not know what is the default name for the SCIP binary')
    end
end

scip_opts_comm = {};
if (nargin >= 8 && ~isempty(scip_opts))
    
    % SCIP binary
    if (isfield(scip_opts, 'scipbin'))
        SCIPBIN = scip_opts.scipbin;
    end
    
    % display options
    QUIETFLAG = [];
    OUTPUTREDIR = []; % output redirecton, e.g. "> /dev/null"
    if isfield(scip_opts, 'display_verblevel')
        if (scip_opts.display_verblevel == 0)
            QUIETFLAG = ' -q ';
            
            % bug workaround: there's a bug in
            % scip-3.1.0.linux.x86_64.gnu.opt.spx, which causes the
            % solution to be written as an empty file when the quiet flag
            % is used. As a workaround, we disable the quiet flag in that
            % case, and instead send the output to /dev/null, as suggested
            % by Stefan Vigerske
            if strcmp(SCIPBIN, 'scip-3.1.0.linux.x86_64.gnu.opt.spx')
                QUIETFLAG = [];
                OUTPUTREDIR = ' > /dev/null';
            end
            
        end
        scip_opts_comm{end+1} = [' -c "set display verblevel ' num2str(scip_opts.display_verblevel) '"'];
    end
    
    % frequency for displaying node information lines [100]
    if (isfield(scip_opts, 'display_freq'))
        scip_opts_comm{end+1} = [' -c "set display freq ' num2str(scip_opts.display_freq) '"'];
    end
    
    % limits options
    
    % solving stops, if the absolute gap = |primalbound - dualbound| is below the given value [0.0]
    if (isfield(scip_opts, 'limits_absgap'))
        scip_opts_comm{end+1} = [' -c "set limits absgap ' num2str(scip_opts.limits_absgap) '"'];
    end
    
    % solving stops, if the relative gap = |primal - dual|/MIN(|dual|,|primal|) is below the given value [0.0]
    if (isfield(scip_opts, 'limits_gap'))
        scip_opts_comm{end+1} = [' -c "set limits gap ' num2str(scip_opts.limits_gap) '"'];
    end
    
    % maximal time in seconds to run [1e+20]
    if (isfield(scip_opts, 'limits_time'))
        scip_opts_comm{end+1} = [' -c "set limits time ' num2str(scip_opts.limits_time) '"'];
    end
    
    % solving stops, if the given number of solutions were found (-1: no limit) [-1]
    if (isfield(scip_opts, 'limits_solutions'))
        scip_opts_comm{end+1} = [' -c "set limits solutions ' num2str(scip_opts.limits_solutions) '"'];
    end
    
    % lp options
    
    % number of threads used for solving the LP (0: automatic)
    if (isfield(scip_opts, 'lp_threads'))
        scip_opts_comm{end+1} = [' -c "set lp advanced threads ' num2str(scip_opts.lp_threads) '"'];
    end
    
    % numerics options
    
    % feasibility tolerance for constraints in SCIP [1e-6]
    if (isfield(scip_opts, 'numerics_feastol'))
        scip_opts_comm{end+1} = [' -c "set numerics feastol ' num2str(scip_opts.numerics_feastol) '"'];
    end
    
end

%% Objective function: 1/2 nu' * H * nu + f' * nu

% pre-compute the weighted Laplacian matrix
V = -w;
V(1:N+1:end) = sum(w, 2);

% quadratic terms of the objective function

% upper triangular matrix terms
[I, J] = find(dx);
idx = I > J;
I(idx) = [];
J(idx) = [];

% remove pairs where both vertices are fixed, as those only contribute a
% constant to the objective function and can be ignored in the optimization
idx = ~isFree(I) & ~isFree(J);
I(idx) = [];
J(idx) = [];

Nterms = length(I);
objfunq = cell(1, Nterms + Nfree);
for idx = 1:Nterms
    % upper triangular terms
    
    % 2D outputs
    if (D == 2)
        
        % xi (free), xj (free)
        if (isFree(I(idx)) && isFree(J(idx)))
            objfunq{idx} = sprintf(...
                '+%.16g x%d x%d + %.16g y%d y%d', ...
                2*full(V(I(idx), J(idx))), I(idx), J(idx), ...
                2*full(V(I(idx), J(idx))), I(idx), J(idx));
        % xi (free), xj (fixed)
        elseif (isFree(I(idx)) && ~isFree(J(idx)))
            objfunq{idx} = sprintf(...
                '+%.16g x%d + %.16g y%d', ...
                2*full(V(I(idx), J(idx))) * y(J(idx), 1), I(idx), ...
                2*full(V(I(idx), J(idx))) * y(J(idx), 2), I(idx));
        % xi (fixed), xj (free)
        elseif (~isFree(I(idx)) && isFree(J(idx)))
            objfunq{idx} = sprintf(...
                '+%.16g x%d + %.16g y%d', ...
                2*full(V(I(idx), J(idx))) * y(I(idx), 1), J(idx), ...
                2*full(V(I(idx), J(idx))) * y(I(idx), 2), J(idx));
        end
        
    % 3D outputs
    elseif (D == 3)
        
        % xi (free), xj (free)
        if (isFree(I(idx)) && isFree(J(idx)))
            objfunq{idx} = sprintf(...
                '+%.16g x%d x%d + %.16g y%d y%d + %.16g z%d z%d', ...
                2*full(V(I(idx), J(idx))), I(idx), J(idx), ...
                2*full(V(I(idx), J(idx))), I(idx), J(idx), ...
                2*full(V(I(idx), J(idx))), I(idx), J(idx));
        % xi (free), xj (fixed)
        elseif (isFree(I(idx)) && ~isFree(J(idx)))
            objfunq{idx} = sprintf(...
                '+%.16g x%d + %.16g y%d + %.16g z%d', ...
                2*full(V(I(idx), J(idx))) * y(J(idx), 1), I(idx), ...
                2*full(V(I(idx), J(idx))) * y(J(idx), 2), I(idx), ...
                2*full(V(I(idx), J(idx))) * y(J(idx), 3), I(idx));
        % xi (fixed), xj (free)
        elseif (~isFree(I(idx)) && isFree(J(idx)))
            objfunq{idx} = sprintf(...
                '+%.16g x%d + %.16g y%d + %.16g z%d', ...
                2*full(V(I(idx), J(idx))) * y(I(idx), 1), J(idx), ...
                2*full(V(I(idx), J(idx))) * y(I(idx), 2), J(idx), ...
                2*full(V(I(idx), J(idx))) * y(I(idx), 3), J(idx));
        end
        
    else
        error('Assertion fail: Output dimension D is not 2 or 3')
    end
end

% main diagonal terms
I = find(isFree);
for idx = 1:length(I)
    
    % 2D outputs
    if (D == 2)

        % only free vertices contribute to the objective function
        objfunq{Nterms + idx} = sprintf(...
            '+%.16g x%d^2 + %.16g y%d^2', ...
            full(V(I(idx), I(idx))), I(idx), ...
            full(V(I(idx), I(idx))), I(idx));
        
    % 3D outputs
    elseif (D == 3)
        
        % only free vertices contribute to the objective function
        objfunq{Nterms + idx} = sprintf(...
            '+%.16g x%d^2 + %.16g y%d^2 + %.16g z%d^2', ...
            full(V(I(idx), I(idx))), I(idx), ...
            full(V(I(idx), I(idx))), I(idx), ...
            full(V(I(idx), I(idx))), I(idx));
        
    else
        error('Assertion fail: Output dimension D is not 2 or 3')
    end
        
end
objfunq{1} = [' obj: ' objfunq{1}];

% the linear term of the objective function (f) has to be computed at each
% iteration of the QPQC-SMACOF algorithm. Thus, it is not computed here

%% SMACOF algorithm

% file name and path to save PIP model, computed solution and initial guess
% (generate unique names so that it is possible to run several instances of
% this function in parallel)
[~, aux] = fileparts(tempname);
pipfilename = [tempdir 'model-' aux '.pip']; % PIP model
solfilename = [tempdir 'model-' aux '.sol']; % output solution
sol0filename = [tempdir 'model-' aux '.sol']; % initial guess

% init stopCondition
stopCondition = [];

% Euclidean distances between vertices in the current solution
dy = dmatrix_con(dx, y);

% initial stress
sigma0 = 0.5 * sum(sum(w .* (dx - dy).^2));

% if dx, dy are sparse matrices, sigma0 will be a sparse scalar, and this
% gives an error with fprintf below
sigma0 = full(sigma0);

% vector of stress computed by the algorithm
sigma = zeros(1, smacof_opts.MaxIter);

% display algorithm's evolution
t = zeros(1, smacof_opts.MaxIter); % time past from 0th iteration
if (strcmp(smacof_opts.Display, 'iter'))
    fprintf('Iter\t\tSigma\t\t\tTime (sec)\n')
    fprintf('===================================================\n')
    fprintf('%d\t\t%.4e\t\t%.4e\n', 0, sigma0, 0.0)
end

% auxiliary intermediate result
mwdx = -w .* dx;

% initialize the storage of the best solution found by the algorithm. We
% start with Inf so that the initial guess will always be replaced by any
% valid solution. The reason is that we cannot be sure that the initial
% guess fulfills the constraints, but any valid solution from a later
% iteration does
sigmabest = Inf;
ybest = nan(size(y));

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
    objfunl = cell(1, Nfree);
    J = find(isFree);
    for idx = 1:length(J)

        % 2D output
        if (D == 2)
            
            objfunl{idx} = sprintf(...
                '+%.16g x%d +%.16g y%d', ...
                f(J(idx), 1), J(idx), f(J(idx), 2), J(idx));
            
        % 3D output
        elseif (D == 3)
            
            objfunl{idx} = sprintf(...
                '+%.16g x%d +%.16g y%d +%.16g z%d', ...
                f(J(idx), 1), J(idx), f(J(idx), 2), J(idx), ...
                f(J(idx), 3), J(idx));
            
        else
            error('Assertion fail: Output dimension D is not 2 or 3')
        end
        
    end
    
    % create PIP file to describe problem
    fid = fopen(pipfilename, 'w');
    if (fid == -1)
        error(['Cannot open file ' pipfilename ' to save PIP model'])
    end
    fprintf(fid, '%s\n%s\n%s\n%s\n%s\n%s\n', ...
        'Minimize', objfunq{:}, objfunl{:}, bnd{:}, con{:}, 'End');
    if (fclose(fid) == -1)
        error(['Cannot close file ' pipfilename ' to save PIP model'])
    end
    
    % create file for the initial guess
    write_solution(sol0filename, y);
    
    % solve the quadratic problem
    system([...
        SCIPBIN...
        QUIETFLAG ...
        ' -c "read ' pipfilename '"'...
        strcat(scip_opts_comm{:}) ...
        ' -c "optimize"'...
        ' -c "write solution ' solfilename '"'...
        ' -c "quit"' ...
        OUTPUTREDIR]);
    
    % read solution
    [aux, status] = read_solution(solfilename, size(y));
    
    % if SCIP cannot find a valid solution (e.g. we set a too short time
    % limit or it doesn't exist), then we cannot update the current
    % solution. This also means that we need to stop the optimization,
    % because the problem matrices and vectors won't change, and basically
    % we would try to solve the same problem again, not reaching a solution
    % either. Note however that it's possible that we have been finding
    % valid solutions before, and it's only at a later iteration that we
    % get to this "not valid solution" situation. Thus, we cannot assume
    % that the stop condition below these lines means "error".
    if (isempty(aux))
        stopCondition{end+1} = ['SCIP: ' status];
        sigma(I) = NaN;
        t(I) = toc;
        break;
    end

    % we only have to update the positions of the free vertices. The values
    % for fixed vertices in aux are all 0, so they need to be ignored
    y(isFree, :) = aux(isFree, :);
    
    % recompute distances between vertices in the current solution
    dy = dmatrix_con(dx, y);

    % compute stress with the current solution
    sigma(I) = 0.5 * sum(sum(w .* (dx - dy).^2));
    
    % update best solution
    if (sigma(I) < sigmabest)
        sigmabest = sigma(I);
        ybest = y;
    end
    
    % display algorithm's evolution
    t(I) = toc;
    if (strcmp(smacof_opts.Display, 'iter'))
        fprintf('%d\t\t%.4e\t\t%.4e\n', I, sigma(I), t(I))
    end
    
    % check whether the stress is under the tolerance level requested by
    % the user
    if (sigma(I) < smacof_opts.TolFun)
        stopCondition{end+1} = 'TolFun';
    end
    
    % check whether the improvement in stress is below the user's request,
    % and it's positive (don't stop if the stress gets worse, because
    % stress can go up and down in the optimization)
    if (I > 1)
        if ((sigma(I-1)-sigma(I))/sigma(I-1) < smacof_opts.Epsilon ...
                && (sigma(I-1)-sigma(I))/sigma(I-1) >= 0)
            stopCondition{end+1} = 'Epsilon';
        end
    end

    % stop if any stop condition has been met
    if (~isempty(stopCondition))
        break;
    end
    
end

% return the best solution the algorithm has found in all iterations
y = ybest;

% check whether the "maximum number of iterations" stop condition has been
% met
if (I == smacof_opts.MaxIter)
    stopCondition{end+1} = 'MaxIter';
end

% prune stress and time vectors if convergence was reached before the
% maximum number of iterations
sigma(I+1:end) = [];
t(I+1:end) = [];

end

% read SCIP solution from a text file created by SCIP
%
% file: path and file name of the solution file.
%
% sz: size of the output matrix with the solution. This parameter makes the
%     code easier to write, and it also allows to detect missing variables
%     in the solution
%
% y: matrix with the solution as a point configuration (each row is a
%    point)
function [y, status] = read_solution(file, sz)

if ((sz(2) ~= 2) && (sz(2) ~= 3))
    error('We only know how to read solutions that are sets of 2D or 3D points')
end

fid = fopen(file, 'r');
if (fid == -1)
    error(['Cannot open file ' file ' to read solution'])
end

% read status of the solution
status = fgetl(fid);
if ((isnumeric(status) && status == -1) ...
        || ~strcmp(status(1:16), 'solution status:'))
    error(['Assertion fail: File with SCIP solution does not start with string ''solution status:''. File ' file])
end
status = status(18:end);

% unless we obtained a valid solution, we don't bother with the rest of the
% file (which should be empty anyway), we exit, and the main function
% should detect the empty solution y, and create a stopCondition with the
% status returned by SCIP
if (~strcmp(status, 'optimal solution found') ...
        && ~strcmp(status, 'solution limit reached'))
    fclose(fid);
    y = [];
    return;
end

% the second line in the file should be the value of the objective function
aux = fgetl(fid);
if (~strcmp(aux(1:16), 'objective value:'))
    error(['Assertion fail: Second line in file with SCIP solution does not start with string ''objective value:''. File ' file])
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
c = textscan(fid, '%s%f%s', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
fclose(fid);

if isempty(c{1})
    error('SCIP did not return any solution')
end

% if SCIP has found a solution, read it from the file into the output
% matrix. Note that if a variable = 0, SCIP will not put it in the output
% file
y = zeros(sz);
for I = 1:length(c{1})-1
    
    % variable name
    varname = c{1}{I};
    
    % if list of output variables has finished, we can exit the loop
    if (strcmp(varname, 'quadobjvar'))
        break;
    end
    
    % vertex index
    idx = str2double(varname(2:end));
    
    if (c{1}{I}(1) == 'x') % this is an x-coordinate
        y(idx, 1) = c{2}(I);
    elseif (c{1}{I}(1) == 'y') % this is a y-coordinate
        y(idx, 2) = c{2}(I);
    elseif (c{1}{I}(1) == 'z') % this is a z-coordinate
        y(idx, 3) = c{2}(I);
    end
    
end

end

% write SCIP initial solution to a text file that SCIP can understand
%
% y: matrix with the solution as a point configuration (each row is a
%    point), including the points that are not free points
%
% file: path and file name of the solution file.
function write_solution(file, y, isFree)

% size of full solution configuration
sz = size(y);

if ((sz(2) ~= 2) && (sz(2) ~= 3))
    error('We only know how to read solutions that are sets of 2D or 3D points')
end

fid = fopen(file, 'w');
if (fid == -1)
    error(['Cannot open file ' file ' to write solution'])
end

% write solution to file. Example result:
%
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
for I = 1:sz(1)
    
    % write x-coordinate
    fprintf(fid, '%s%d\t\f.16\t%s\n', 'x', I, y(I, 1), '(obj:0)');

    % write y-coordinate
    fprintf(fid, '%s%d\t\f.16\t%s\n', 'y', I, y(I, 2), '(obj:0)');
    
    % write z-coordinate
    if (sz(2) == 3)
        fprintf(fid, '%s%d\t\f.16\t%s\n', 'z', I, y(I, 3), '(obj:0)');
    end
    
end

end
