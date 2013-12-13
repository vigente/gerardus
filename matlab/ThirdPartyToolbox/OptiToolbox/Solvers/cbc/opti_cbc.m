function [x,fval,exitflag,info] = opti_cbc(f,A,rl,ru,lb,ub,xtype,sos,opts)
%OPTI_CBC Solve a MILP using CBC
%
%   min f'*x      subject to:     rl <= A*x <= ru
%    x                            lb <= x <= ub
%                                 for i = 1..n: xi in Z
%                                 for j = 1..m: xj in {0,1} 
%
%   x = opti_cbc(f,A,rl,ru,lb,ub,xint) solves a MILP where f is the 
%   objective vector, A,rl,ru are the linear constraints, lb,ub are the
%   bounds and xint is a string of integer variables ('C', 'I', 'B')
%
%   x = opti_cbc(f,...,xint,sos) sos is a structure with fields type, 
%   index, weight for SOS.
%
%   x = opti_cbc(f,...,sos,opts) uses opts to pass optiset options to the
%   solver.
%
%   [x,fval,exitflag,info] = opti_cbc(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR CBC USING THE MEX INTERFACE
%   See supplied Eclipse Public License

%   Copyright (C) 2012 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 9, opts = optiset; end 
if nargin < 8, sos = []; end
if nargin < 7, xtype = repmat('C',size(f)); end
if nargin < 6, ub = []; end
if nargin < 5, lb = []; end
if nargin < 4, error('You must supply at least 4 arguments to opti_cbc'); end

warn = strcmpi(opts.warnings,'all');

%Check sparsity
if(~issparse(A))
    if(warn)
        optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)]');
    end
    A = sparse(A);
end

%Setup Printing
opts.display = dispLevel(opts.display);
    
%Check Integer Vars
if(~ischar(xtype) || length(xtype) ~= length(f))
    error('The integer string must be a char array %d x 1!',length(f));
else
    xtype = upper(xtype);
end

%Setup SOS
if(~isempty(sos))
    if(~isstruct(sos) || ~isfield(sos,'type') || ~isfield(sos,'index') || ~isfield(sos,'weight'))
        error('SOS constraints must be a structure with fields ''type'', ''index'' and ''weight''');
    end
    if(length(sos.type) == 1)
        sos.index = {sos.index};
        sos.weight = {sos.weight};
    end
end

%Get Objective Bias (cbc does not seem to process it)
if(isfield(opts,'objbias') && ~isempty(opts.objbias))
    objbias = opts.objbias;
    opts = rmfield(opts,'objbias'); %incase a later version does!
else
    objbias = 0;
end

%MEX contains error checking
[x,fval,exitflag,iter,cobj] = cbc(f, A, rl, ru, lb, ub, xtype, sos, opts);

%Add Objective Bias
fval = fval + objbias;

%Assign Outputs
info.Nodes = iter;
info.AbsGap = abs(cobj-fval);
info.RelGap = abs(cobj-fval)/(1e-1 + abs(fval));
info.Time = toc(t);
info.Algorithm = 'CBC: Branch and Cut using CLP';

switch(exitflag)
    case 1
        info.Status = 'Integer Optimal';
    case 0
        info.Status = 'Exceeded Iterations';
    case -1
        info.Status = 'Infeasible';
    case -2
        info.Status = 'Unbounded or Infeasible';
    case -5
        info.Status = 'User Exited';
    otherwise
        info.Status = [];
end


function  print_level = dispLevel(lev)
%Return CLP compatible display level
switch(lower(lev))
    case'off'
        print_level = 0;
    case 'iter'
        print_level = 1;
    case 'final'
        print_level = 1;
end