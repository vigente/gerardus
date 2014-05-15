function [x,fval,exitflag,info] = opti_lbfgsb(fun,grad,lb,ub,x0,opts)
%OPTI_LBFGSB Solve a bounded NLP using L-BFGS-B (Nocedal & Zhu)
%
%   min f(x)       subject to:   lb <= x <= ub
%    x
%
%   x = opti_lbfgsb(fun,grad,lb,ub,x0) solves a NLP where fun is the 
%   objective function and grad is the gradient of the objective. lb and ub
%   are the decision variable bounds (required) and x0 is a starting guess.
%   Avoid Infinite bounds.
%
%   x = opti_lbfgsb(fun,...,x0,opts) uses opts to pass optiset options to 
%   the solver. 
%
%   [x,fval,exitflag,info] = opti_lbfgsb(...) returns the objective value 
%   at the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR L-BFGS-B
%   See referenced BSD License

%   Copyright (C) 2012 Jonathan Currie (I2C2)

if(nargin < 6), opts = optiset; end
if(nargin < 5), error('LBFGSB requires at least 5 arguments'); end

%Setup display level
opts.display = dispLevel(opts.display);

%Ensure we have bounds
if(isempty(lb) || isempty(ub))
    error('L-BFGS-B only solves bounded NLPs - You must supply lb and ub to this function');
end
%Ensure we have a gradient
if(isempty(grad))
    error('L-BFGS-B requires a gradient function');
end

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('L-BFGS-B requires an initial guess, x0!');
end

t = tic;
% Run L-BFGS-B
[x, fval, exitflag, iter, feval] = lbfgsb(fun,grad,lb,ub,x0,opts);

%Collect Results
info.Iterations = iter;
info.FuncEvals = feval;
info.Time = toc(t);
info.Algorithm = 'L-BFGS-B: Limited Memory BFGS Bounded Optimization';

switch(exitflag)
    case 1
        info.Status = 'Optimal';
    case 0
        info.Status = 'Exceeded Iterations';
    case -1
        info.Status = 'Infeasible / Could not Converge';
    case -2
        info.Status = 'L-BFGS-B Error';
    case -5
        info.Status = 'User Exited';
    otherwise        
        info.Status = 'L-BFGS-B Error';
end
