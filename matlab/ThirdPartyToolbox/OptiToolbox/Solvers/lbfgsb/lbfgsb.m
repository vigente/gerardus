% LBFGSB  Solve a Bounded NLP using LBFGSB
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_lbfgsb() INSTEAD!
%
% lbfgsb uses the Limited Memory BFGS Bounded Optimization library.
%
%   [x,fval,exitflag,iter] = lbfgsb(fun,grad,x0,lb,ub,opts)
%
%   Input arguments:
%       fun - nonlinear function handle
%       grad - gradient of nonlinear function handle
%       x0 - initial solution guess
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%
%   Option Fields (all optional):
%       tolrfun - relative function tolerance
%       maxiter - maximum iterations
%       display - solver display level [0,1,2]
%
%   Return Status:
%       1 - optimal
%       0 - maximum iterations exceeded
%      -1 - abnormal termination
%      -2 - error on input / other
%      -5 - user exit
%
%   
%   Based in parts on the original MEX interface by Dr. Peter Carbonetto.
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)