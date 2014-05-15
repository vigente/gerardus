function [x,fval,exitflag,info] = opti_lipsol(f,A,b,Aeq,beq,lb,ub,opts)
%OPTI_LIPSOL Solve a LP using LIPSOL
%
%   min f'*x                   subject to:     A*x <= b
%    x                                         Aeq*x == beq
%                                              lb <= x <= ub
%
%   x = opti_lipsol(f,A,b,Aeq,beq,lb,ub) solves a LP where f is the objective 
%   vector, A,b,Aeq,beq are the linear constraints and lb,ub are the bounds.
%
%   x = opti_lipsol(f,...,ub,opts) uses opts to pass optiset options to the
%   solver. 
%
%   [x,fval,exitflag,info] = opti_lipsol(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR LIPSOL USING THE MEX INTERFACE
%   See supplied GNU GPL

%   Copyright (C) 2013 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 8, opts = optiset; end 
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, beq = []; end
if nargin < 4, Aeq = []; end
if nargin < 3, error('You must supply at least 3 arguments to opti_lipsol'); end

warn = strcmpi(opts.warnings,'all');

%Add in lipsolset settings
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    popts = lipsolset(opts.solverOpts);    
else
    popts = lipsolset;
end
%Add in options from optiset   
popts.maxiter = opts.maxiter;
popts.maxtime = opts.maxtime;
popts.tol = opts.tolrfun;
%Setup display level
popts.verb = dispLevel(opts.display,1,1);

%Check sparsity
if(~isempty(A) && ~issparse(A))
    if(warn), optiwarn('OPTI:NotSparseA','The A matrix should be sparse, correcting: [sparse(A)]'); end
    A = sparse(A);
end
if(~isempty(Aeq) && ~issparse(Aeq))
    if(warn), optiwarn('OPTI:NotSparseAeq','The Aeq matrix should be sparse, correcting: [sparse(Aeq)]'); end
    Aeq = sparse(Aeq);
end

%Check for bounds
if(isempty(lb)), lb = -Inf(size(f)); end
if(isempty(ub)), ub = Inf(size(f)); end
%Check for empty Aeq
if(isempty(Aeq)), Aeq = zeros(0,numel(f)); end

%Augment Inequalities as Equalities (via slack vars)
neq = size(Aeq,1);
nineq = size(A,1);
if(nineq)
    Aeq = [Aeq spalloc(neq,nineq,0);A speye(nineq)];
    beq = [beq;b];
    fn = [f;zeros(nineq,1)];
    lb = [lb;zeros(nineq,1)];
    ub = [ub;Inf(nineq,1)];
else
    fn = f;
end

%LIPSOL contains error checking
[x,fval,inf,msg] = lipsol(Aeq, beq, fn, lb, ub, popts);

%Just Keep Existing x (not augmented system)
x = x(1:numel(f));

%Add objbias if present
if(~isempty(fval) && isfield(opts,'objbias') && ~isempty(opts.objbias))
    fval = fval + opts.objbias;
end

%Assign Outputs
info.Iterations = inf(2);
info.Time = toc(t);
info.Algorithm = 'LIPSOL: Interior Point LP Solver'; 
info.Status = msg;

switch(inf(1))
    case 1
        exitflag = 1;
    case 2
        exitflag = -2;
    case -1
        exitflag = -1;
    case -3
        info.Status = 'Exceeded Maximum Iterations';
        exitflag = 0;
    case -3.5
        info.Status = 'Exceeded Maximum Time';
        exitflag = 0;
    otherwise
        info.Status = 'Unknown Exit';
        exitflag = -4;
end

%Assign Lambda
% info.Lambda = struct('ineqlin',lam.dual_row(~eq),'eqlin',lam.dual_row(eq),'bounds',lam.dual_col);