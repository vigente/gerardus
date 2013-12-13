function [x,fval,exitflag,info] = opti_qsopt(f,A,b,Aeq,beq,lb,ub,opts)
%OPTI_QSOPT Solve a LP using QSOPT
%
%   [x,fval,exitflag,info] = opti_qsopt(f,A,b,Aeq,beq,lb,ub) solves the linear
%   program min f'x where A,b are the inequality constraints, Aeq,beq are
%   the equality constraints and lb,ub are the bounds.
%
%   THIS IS A WRAPPER FOR QSOPT USING THE MEX INTERFACE
%   See supplied License

%   Copyright (C) 2011 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 8, opts = optiset; end 
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, beq = []; end
if nargin < 4, Aeq = []; end
if nargin < 3, error('You must supply at least 3 arguments to opti_qsopt'); end

warn = strcmpi(opts.warnings,'al');

%Merge linear constraints
opts.nin = size(A,1); %save no ineq
A = [A;Aeq];
b = [b;beq];

if(~issparse(A))
    if(warn)
        optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)]');
    end
    A = sparse(A);
end

%Setup Printing
opts.display = dispLevel(opts.display);
 
%MEX contains error checking
[x,fval,exitflag,dualrow,dualcol] = qsopt(f, A, b, lb, ub, opts);

%Assign Outputs
info.Iterations = [];
info.Time = toc(t);
info.Algorithm = 'QSOPT: Dual Simplex';

switch(exitflag)
    case 1
        info.Status = 'Optimal';
    case 0
        info.Status = 'Exceeded Iterations';
    case -1
        info.Status = 'Infeasible';
    case -2
        info.Status = 'Unbounded or Infeasible';
    otherwise
        info.Status = [];
end

%Add constant objective term
if(~isempty(fval) && isfield(opts,'objbias') && ~isempty(opts.objbias))
    fval = fval + opts.objbias;
end

%Assign Lambda
if(opts.nin==0)
    ineq = 0;
    dualin = [];
else
    ineq = 1:opts.nin;
    dualin = dualrow(ineq);
end
eq = (1:length(beq)) + ineq(end);
info.Lambda = struct('ineqlin',dualin,'eqlin',dualrow(eq),'bounds',dualcol);

function  print_level = dispLevel(lev)
%Return QSOPT compatible display level
switch(lower(lev))
    case'off'
        print_level = 0;
    case 'iter'
        print_level = 1;
    case 'final'
        print_level = 0; %not used
end
